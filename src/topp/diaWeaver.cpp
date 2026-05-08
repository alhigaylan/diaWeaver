// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Mohammed Alhigaylan $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/APPLICATIONS/diaWeaver.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerIM.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include <OpenMS/FEATUREFINDER/ElutionPeakDetection.h>
#include <OpenMS/FEATUREFINDER/FeatureFindingPeptide.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ClusterMassTracesByPrecursor.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/FORMAT/MSNumpressCoder.h>
#include <OpenMS/METADATA/SourceFile.h>

#ifdef WITH_OPENTIMS
#include <OpenMS/FORMAT/BrukerTimsFile.h>
#endif

#include <cmath>
#include <set>
#include <chrono>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace OpenMS;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_diaWeaver diaWeaver

@brief Generates pseudo spectra from DIA data by correlating MS1 precursor features with MS2 fragment traces.

This tool processes DIA (Data Independent Acquisition) data to generate pseudo MS/MS spectra
by correlating precursor elution profiles with fragment ion traces. It outputs one mzML file
per DIA window containing the reconstructed spectra.

The processing pipeline for each DIA window:
1. Peak picking on MS1 and MS2 spectra (PeakPickerIM for ion mobility data, PeakPickerHiRes otherwise)
2. FeatureFinderPeptide on MS1 to detect monoisotopic peptide precursor features
3. MassTraceExtractor on MS2 to extract fragment mass traces
4. Correlation of MS1 feature elution profiles with MS2 fragment traces using Pearson correlation
5. Assembly of correlated fragments into pseudo MS/MS spectra

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_diaWeaver.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_diaWeaver.html
*/

class TOPPDiaWeaver :
  public TOPPBase
{
public:
  TOPPDiaWeaver() :
    TOPPBase(
      "diaWeaver",
      "Split a DIA mzML into per-window MS1 and MS2 mzML files",
      false)
  {
  }

protected:

  /// Aggregate a single spectrum with its RT neighbors using Gaussian weighting.
  /// Reads from the immutable source experiment; safe to call from a parallel loop
  /// that writes picked results into a separate output experiment.
  static void aggregateSpectrum_(
    const MSExperiment& exp,
    Size center_idx,
    const PeakPickerIM& picker,
    MSSpectrum& out)
  {
    if (center_idx >= exp.size()) return;

    Param params = picker.getParameters();
    double fwhm = (double)params.getValue("aggregation:rt_FWHM");
    double cutoff = (double)params.getValue("aggregation:cutoff");
    double factor = -4.0 * std::log(2.0) / (fwhm * fwhm);

    double center_rt = exp[center_idx].getRT();

    std::vector<MSSpectrum> spectra_to_aggregate;
    std::vector<double> weights;

    // Search forward (including center)
    for (Size j = center_idx; j < exp.size(); ++j)
    {
      double rt_diff = exp[j].getRT() - center_rt;
      double weight = std::exp(factor * rt_diff * rt_diff);
      if (weight < cutoff && j != center_idx) break;
      spectra_to_aggregate.push_back(exp[j]);
      weights.push_back(weight);
    }

    // Search backward
    for (SignedSize j = static_cast<SignedSize>(center_idx) - 1; j >= 0; --j)
    {
      double rt_diff = exp[j].getRT() - center_rt;
      double weight = std::exp(factor * rt_diff * rt_diff);
      if (weight < cutoff) break;
      spectra_to_aggregate.push_back(exp[j]);
      weights.push_back(weight);
    }

    // Normalize weights
    double sum_w = 0.0;
    for (double w : weights) sum_w += w;
    for (double& w : weights) w /= sum_w;

    picker.aggregateScans(spectra_to_aggregate, weights, out);
  }

#ifdef WITH_OPENTIMS
  BrukerTimsFile::Config getBrukerConfig_()
  {
    BrukerTimsFile::Config c;
    c.calibration_tolerance = getDoubleOption_("bruker:calibration_tolerance");
    c.calibrate = (getStringOption_("bruker:calibrate") == "true");
    String mode = getStringOption_("bruker:export_mode");
    c.export_mode = (mode == "frame") ? BrukerTimsFile::Config::FRAME : BrukerTimsFile::Config::AUTO;
    c.ms1_centroid_mz_ppm = static_cast<float>(getDoubleOption_("bruker:ms1_centroid_mz_ppm"));
    c.ms1_centroid_im_pct = static_cast<float>(getDoubleOption_("bruker:ms1_centroid_im_pct"));
    c.ms1_n_neighbors         = getIntOption_("bruker:ms1_n_neighbors");
    c.ms1_min_support         = getIntOption_("bruker:ms1_min_support");
    c.ms1_max_rt_distance_sec = getDoubleOption_("bruker:ms1_max_rt_distance_sec");
    c.ms1_centroid_max_peaks  = getIntOption_("bruker:ms1_centroid_max_peaks");
    c.dia_ms2_n_neighbors = getIntOption_("bruker:dia_ms2_n_neighbors");
    c.dia_ms2_min_support = getIntOption_("bruker:dia_ms2_min_support");
    c.dia_ms2_centroid = (getStringOption_("bruker:dia_ms2_centroid") == "true");
    return c;
  }
#endif

  void registerOptionsAndFlags_() override
  {
    registerInputFile_(
      "in",
      "<file>",
      "",
      "Input DIA file (mzML or Bruker .d)",
      true);
    setValidFormats_("in", {"mzML"
#ifdef WITH_OPENTIMS
      , "d"
#endif
    });

    registerOutputPrefix_(
      "out",
      "<prefix>",
      "",
      "Output directory prefix (default: <input>_diaWindows)");

    registerFlag_(
      "save_unfragmented_precursors",
      "If set, save peaks within precursor isolation window in MS2 and apply precursor detection algorithm");

    registerFlag_(
      "keep_ms1",
      "If set, include peak-picked MS1 spectra in the output file alongside pseudo spectra");

    registerFlag_(
      "aggregate_across_scans",
      "If set, aggregate signal across neighboring scans using Gaussian weighting before peak picking. "
      "This can improve signal-to-noise for low-intensity peaks (requires IM data).", false);

    registerSubsection_("PeakPickerIM", "Parameters for ion mobility peak picking (used when input has IM data)");

    registerSubsection_("PeakPickerHiRes", "Parameters for high-resolution peak picking (used when input has no IM data)");

    registerSubsection_("FeatureFinderPeptide", "Parameters for FeatureFinderPeptide algorithm (precursor detection on MS1 and unfragmented precursors)");

    registerSubsection_("MassTraceExtractor", "Parameters for MassTraceExtractor algorithm (mass trace detection on MS2 data)");

    registerSubsection_("ClusterMassTraces", "Parameters for clustering mass traces to create pseudo spectra");

#ifdef WITH_OPENTIMS
    registerTOPPSubsection_("bruker", "Options for reading Bruker TimsTOF .d files (requires WITH_OPENTIMS)");
    registerStringOption_("bruker:export_mode", "<mode>", "auto", "Export mode: 'auto' detects DDA/DIA acquisition type, "
      "'frame' returns raw 4D frames without signal processing.", false, true);
    setValidStrings_("bruker:export_mode", {"auto", "frame"});
    registerDoubleOption_("bruker:calibration_tolerance", "<float>", 0.0, "m/z recalibration tolerance (0 = library default)", false, true);
    setMinFloat_("bruker:calibration_tolerance", 0.0);
    registerStringOption_("bruker:calibrate", "<toggle>", "false", "Enable m/z recalibration (may fail on some datasets)", false, true);
    setValidStrings_("bruker:calibrate", {"true", "false"});
    registerDoubleOption_("bruker:ms1_centroid_mz_ppm", "<float>", 5.0,
      "MS1 frame IM-centroiding m/z tolerance in ppm. Collapses the ion mobility dimension "
      "by aggregating neighboring peaks directly on the raw gridded data (Sage algorithm, Lazear 2023). "
      "Both this and bruker:ms1_centroid_im_pct must be > 0 to enable. Set to 0 to disable.", false, true);
    setMinFloat_("bruker:ms1_centroid_mz_ppm", 0.0);
    registerDoubleOption_("bruker:ms1_centroid_im_pct", "<float>", 3.0,
      "MS1 frame IM-centroiding ion mobility tolerance in percent. "
      "Both this and bruker:ms1_centroid_mz_ppm must be > 0 to enable. Set to 0 to disable.", false, true);
    setMinFloat_("bruker:ms1_centroid_im_pct", 0.0);
    registerIntOption_("bruker:ms1_n_neighbors", "<int>", 0,
      "MS1 frame aggregation: number of adjacent MS1 frames on each side to sum. "
      "0 = disabled (raw export), 1 = 3-frame sum, 2 = 5-frame sum. "
      "Applies to both DIA and DDA; ignored in FRAME export mode.", false, true);
    setMinInt_("bruker:ms1_n_neighbors", 0);
    setMaxInt_("bruker:ms1_n_neighbors", 50);
    registerIntOption_("bruker:ms1_min_support", "<int>", 0,
      "MS1 denoising: minimum occupied neighbor cells in a 3x3 (m/z x IM) grid to keep a point. "
      "Applied after aggregation. 0 = disabled, 8 = all 8 neighbors required (strictest). "
      "Only effective when ms1_n_neighbors > 0. Appropriate for dense survey runs; disable for "
      "rare-species discovery.", false, true);
    setMinInt_("bruker:ms1_min_support", 0);
    setMaxInt_("bruker:ms1_min_support", 8);
    registerDoubleOption_("bruker:ms1_max_rt_distance_sec", "<float>", 0.0,
      "Cap the RT distance (seconds) between a neighbor MS1 frame and the center frame during "
      "aggregation. 0.0 = no cap. Recommended for DDA (e.g. 5.0) where MS1 frame cadence is "
      "irregular. The center frame is always included regardless of this cap.", false, true);
    setMinFloat_("bruker:ms1_max_rt_distance_sec", 0.0);
    registerIntOption_("bruker:ms1_centroid_max_peaks", "<int>", 100000,
      "Cap on the number of centroided peaks retained per MS1 spectrum. Top-intensity peaks "
      "are kept; low-intensity tail is dropped if the limit is hit (a warning is logged in that "
      "case). Only effective when MS1 centroiding is enabled via ms1_centroid_mz_ppm/pct. Raise "
      "for aggregated MS1 (ms1_n_neighbors > 0) on dense surveys; lower to trim long-tail noise.", false, true);
    setMinInt_("bruker:ms1_centroid_max_peaks", 1);
    registerIntOption_("bruker:dia_ms2_n_neighbors", "<int>", 0,
      "DIA MS2 frame aggregation: number of adjacent frames on each side to sum per SWATH window. "
      "0 = disabled (raw export), 1 = 3-frame sum, 2 = 5-frame sum.", false, true);
    setMinInt_("bruker:dia_ms2_n_neighbors", 0);
    registerIntOption_("bruker:dia_ms2_min_support", "<int>", 1,
      "DIA MS2 denoising: minimum occupied neighbor cells in a 3x3 (m/z x IM) grid to keep a point. "
      "Only effective when bruker:dia_ms2_n_neighbors > 0.", false, true);
    setMinInt_("bruker:dia_ms2_min_support", 1);
    registerStringOption_("bruker:dia_ms2_centroid", "<toggle>", "true",
      "Apply 2D Gaussian smoothing + local maxima peak picking to the denoised DIA MS2 grid. "
      "Produces IM_CENTROIDED spectra. Only effective when bruker:dia_ms2_n_neighbors > 0.", false, true);
    setValidStrings_("bruker:dia_ms2_centroid", {"true", "false"});
#endif

    registerIntOption_(
      "threads",
      "<n>",
      1,
      "Total number of threads to use for processing",
      false);
    setMinInt_("threads", 1);

    registerIntOption_(
      "threads_outer_loop",
      "<n>",
      -1,
      "Number of threads for the outer loop (over DIA windows). Remaining threads are used for "
      "inner loop (peak picking within each window). Set to -1 to use all threads in outer loop only (no nested parallelism). "
      "Example: with 24 total threads and 4 outer threads, each window gets 6 threads for peak picking.",
      false);
  }

  Param getSubsectionDefaults_(const String& name) const override
  {
    if (name == "PeakPickerIM")
    {
      PeakPickerIM ppim;
      return ppim.getDefaults();
    }
    if (name == "PeakPickerHiRes")
    {
      PeakPickerHiRes pphr;
      return pphr.getDefaults();
    }
    if (name == "FeatureFinderPeptide")
    {
      Param combined;

      // Common parameters for all FFP sub-algorithms
      Param p_com;
      p_com.setValue("noise_threshold_int", 60.0, "Intensity threshold below which peaks are regarded as noise. Ignored if auto_noise_threshold is true.");
      p_com.setValue("chrom_peak_snr", 1.0, "Minimum signal-to-noise a mass trace should have.");
      p_com.setValue("chrom_fwhm", 5.0, "Expected chromatographic peak width (in seconds).");
      p_com.setValue("auto_noise_threshold", "true", "If true, automatically estimates the noise threshold from the input map using random scan sampling. Overrides noise_threshold_int.");
      p_com.setValidStrings("auto_noise_threshold", {"true", "false"});
      p_com.setValue("noise_estimation_n_scans", 50, "Number of scans randomly sampled to estimate the noise level when auto_noise_threshold is true.");
      p_com.setValue("noise_estimation_percentile", 80.0, "Intensity percentile used to define the noise level from sampled scans when auto_noise_threshold is true.");
      combined.insert("common:", p_com);
      combined.setSectionDescription("common", "Common parameters for all other subsections");

      // MassTraceDetection parameters
      Param p_mtd = MassTraceDetection().getDefaults();
      p_mtd.setValue("mass_error_ppm", 7.0, "Allowed mass deviation (in ppm).");
      p_mtd.setValue("min_trace_length", 5.0, "Minimum expected length of a mass trace (in seconds).");
      p_mtd.setValue("ion_mobility_tolerance", 0.01, "Allowed ion mobility deviation (in 1/k0).");
      p_mtd.setValue("reestimate_mt_sd", "false", "Enables dynamic re-estimation of m/z variance during mass trace collection stage.");
      p_mtd.setValue("quant_method", "max_height", "Method of quantification for mass traces. For LC data 'area' is recommended, 'median' for direct injection data. 'max_height' simply uses the most intense peak in the trace.");
      p_mtd.setValue("trace_termination_outliers", 2, "Mass trace extension in one direction cancels if this number of consecutive spectra with no detectable peaks is reached.");
      p_mtd.setValue("impute_zeros_missing_scans", "true", "Insert a zero-intensity data point at every scan position where no acceptable peak was found during trace extension, instead of silently skipping it. This preserves the full RT grid of each trace.");

      p_mtd.remove("noise_threshold_int");
      p_mtd.remove("chrom_peak_snr");
      p_mtd.remove("auto_noise_threshold");
      p_mtd.remove("noise_estimation_n_scans");
      p_mtd.remove("noise_estimation_percentile");
      combined.insert("mtd:", p_mtd);
      combined.setSectionDescription("mtd", "Mass Trace Detection parameters");

      // ElutionPeakDetection parameters
      Param p_epd;
      p_epd.setValue("enabled", "true", "Enable splitting of isobaric mass traces by chromatographic peak detection. Disable for direct injection.");
      p_epd.setValue("width_filtering", "off", "Enable filtering of unlikely peak widths. The fixed setting filters out mass traces outside the [min_fwhm, max_fwhm] interval (set parameters accordingly!). The auto setting filters with the 5 and 95% quantiles of the peak width distribution.");
      p_epd.setValidStrings("enabled", {"true", "false"});
      p_epd.insert("", ElutionPeakDetection().getDefaults());

      p_epd.remove("chrom_peak_snr");
      p_epd.remove("chrom_fwhm");
      combined.insert("epd:", p_epd);
      combined.setSectionDescription("epd", "Elution Profile Detection (to separate isobaric Mass Traces by elution time).");

      // FeatureFindingPeptide parameters
      Param p_ffp = FeatureFindingPeptide().getDefaults();
      p_ffp.setValue("local_rt_range", 5.0, "RT range where to look for coeluting mass traces");
      p_ffp.setValue("local_mz_range", 3.0, "MZ range where to look for isotopic mass traces");
      p_ffp.setValue("local_im_range", 0.02, "IM range where to look for isotopic mass traces");
      p_ffp.setValue("charge_lower_bound", 2, "Lowest charge state to consider");
      p_ffp.setValue("charge_upper_bound", 4, "Highest charge state to consider");
      p_ffp.setValue("remove_single_traces", "true", "Remove unassembled traces (single traces).");
      p_ffp.setValue("use_smoothed_intensities", "true", "Use Savitzky-Golay smoothed intensities (produced by ElutionPeakDetection) instead of raw intensities.");
      p_ffp.setValue("mass_defect_filtering", "true", "Filter feature hypotheses by peptide mass defect boundaries.");
      p_ffp.setValue("mass_defect_offset", 0.1, "Mass defect tolerance offset for the peptide mass defect filter.");
      p_ffp.setValue("overlapping_features", "false", "Allow low-confidence hypotheses to reuse mass traces already claimed by higher-scoring features, provided they propose a different charge state.");
      p_ffp.setValue("hypothesis_score_quantile", 0.5, "Score quantile threshold used to classify hypotheses as low-confidence when overlapping_features is true. Hypotheses scoring below this quantile may reuse traces from higher-scoring features at a different charge state.");
      p_ffp.setValue("rt_max_lag", 5, "Maximum lag (in scans) for normalised cross-correlation between isotope elution profiles.");
      p_ffp.setValue("rt_min_pearson_correlation", 0.3, "Minimum Pearson correlation between two mass trace elution profiles.");
      p_ffp.setValue("rt_peak_overlap_threshold", 0.3, "Minimum FWHM overlap proportion required between two co-eluting mass traces.");

      p_ffp.remove("chrom_fwhm");
      p_ffp.remove("report_chromatograms");
      combined.insert("ffp:", p_ffp);
      combined.setSectionDescription("ffp", "FeatureFindingPeptide parameters (assembling mass traces to charged features)");

      return combined;
    }
    if (name == "MassTraceExtractor")
    {
      Param combined;

      // Common parameters
      Param p_com;
      p_com.setValue("noise_threshold_int", 30.0, "Intensity threshold below which peaks are regarded as noise. Ignored if auto_noise_threshold is true.");
      p_com.setValue("chrom_peak_snr", 1.0, "Minimum signal-to-noise a mass trace should have.");
      p_com.setValue("chrom_fwhm", 3.0, "Expected chromatographic peak width (in seconds).");
      p_com.setValue("auto_noise_threshold", "true", "If true, automatically estimates the noise threshold from the input map using random scan sampling. Overrides noise_threshold_int.");
      p_com.setValidStrings("auto_noise_threshold", {"true", "false"});
      p_com.setValue("noise_estimation_n_scans", 50, "Number of scans randomly sampled to estimate the noise level when auto_noise_threshold is true.");
      p_com.setValue("noise_estimation_percentile", 80.0, "Intensity percentile used to define the noise level from sampled scans when auto_noise_threshold is true.");

      combined.insert("common:", p_com);
      combined.setSectionDescription("common", "Common parameters for all other subsections");

      // MassTraceDetection parameters
      Param p_mtd = MassTraceDetection().getDefaults();
      p_mtd.setValue("mass_error_ppm", 7.0, "Allowed mass deviation (in ppm).");
      p_mtd.setValue("min_trace_length", 2.0, "Minimum expected length of a mass trace (in seconds).");
      p_mtd.setValue("ion_mobility_tolerance", 0.01, "Allowed ion mobility deviation (in 1/k0).");
      p_mtd.setValue("reestimate_mt_sd", "false", "Enables dynamic re-estimation of m/z variance during mass trace collection stage.");
      p_mtd.setValue("quant_method", "max_height", "Method of quantification for mass traces. For LC data 'area' is recommended, 'median' for direct injection data. 'max_height' simply uses the most intense peak in the trace.");
      p_mtd.setValue("trace_termination_outliers", 2, "Mass trace extension in one direction cancels if this number of consecutive spectra with no detectable peaks is reached.");
      p_mtd.setValue("impute_zeros_missing_scans", "true", "Insert a zero-intensity data point at every scan position where no acceptable peak was found during trace extension, instead of silently skipping it. This preserves the full RT grid of each trace.");

      p_mtd.remove("noise_threshold_int");
      p_mtd.remove("chrom_peak_snr");
      p_mtd.remove("auto_noise_threshold");
      p_mtd.remove("noise_estimation_n_scans");
      p_mtd.remove("noise_estimation_percentile");
      combined.insert("mtd:", p_mtd);
      combined.setSectionDescription("mtd", "Mass Trace Detection parameters");

      // ElutionPeakDetection parameters
      Param p_epd;
      p_epd.setValue("enabled", "true", "Enable splitting of isobaric mass traces by chromatographic peak detection.");
      p_epd.setValue("width_filtering", "off", "Enable filtering of unlikely peak widths. The fixed setting filters out mass traces outside the [min_fwhm, max_fwhm] interval (set parameters accordingly!). The auto setting filters with the 5 and 95% quantiles of the peak width distribution.");
      p_epd.setValidStrings("enabled", {"true", "false"});
      p_epd.insert("", ElutionPeakDetection().getDefaults());

      p_epd.remove("chrom_peak_snr");
      p_epd.remove("chrom_fwhm");
      combined.insert("epd:", p_epd);
      combined.setSectionDescription("epd", "Elution Profile Detection (to separate isobaric Mass Traces by elution time).");

      return combined;
    }
    if (name == "ClusterMassTraces")
    {
      Param p;
      p.setValue("min_pearson_correlation", 0.3, "Minimal pearson correlation score to match elution profiles to each other.");
      p.setValue("max_lag", 1, "Maximal lag (e.g. by how many spectra the peak may be shifted at most).");
      p.setValue("min_nr_ions", 30, "Minimal number of ions to report a spectrum.");
      p.setValue("max_rt_apex_difference", 5.0, "Maximal difference of the apex in retention time (in seconds).");
      p.setValue("im_tolerance", 0.02, "Ion mobility tolerance for matching precursors to fragments.");
      p.setValue("nr_precursors_per_fragment", 50, "Maximum number of precursors a fragment can be assigned to.");
      p.setValue("rt_tolerance", 2.0, "RT tolerance (in seconds) for matching up mass trace points during correlation.");
      p.setValue("pearson_weight", 1.0, "Weight for the Pearson correlation component in the combined score used to rank precursor-fragment assignments.");
      p.setValue("delta_rt_weight", 1.0, "Weight for the delta RT component in the combined score used to rank precursor-fragment assignments.");
      p.setValue("delta_im_weight", 1.0, "Weight for the delta ion mobility component in the combined score used to rank precursor-fragment assignments.");
      p.setValue("max_nr_ions", 500, "Maximum number of fragment ions written per output spectrum. Assignments are ranked and truncated to this limit. Set to 0 to disable the limit.");
      p.setValue("assign_unassigned_to_all", "false", "Assign unassigned MS2 fragments to all precursors within RT range.");
      p.setValidStrings("assign_unassigned_to_all", {"true", "false"});
      p.setValue("use_combined_scores", "true", "If true, rank precursor-fragment assignments by a weighted combination of Pearson correlation, delta RT and delta IM. If false, use Pearson correlation only.");
      p.setValidStrings("use_combined_scores", {"false", "true"});
      p.setValue("output_fragment_scores", "false", "If true, output per-fragment scores as FloatDataArrays in the output spectra.");
      p.setValidStrings("output_fragment_scores", {"false", "true"});
      p.setValue("smooth_ms1", "true", "Use EPD-smoothed intensities for MS1 elution profiles in correlation scoring. Falls back to raw intensities if smoothed data is absent.");
      p.setValidStrings("smooth_ms1", {"false", "true"});
      p.setValue("smooth_ms2", "false", "Use EPD-smoothed intensities for MS2 elution profiles in correlation scoring. Falls back to raw intensities if smoothed data is absent.");
      p.setValidStrings("smooth_ms2", {"false", "true"});
      return p;
    }
    return Param();
  }

  /**
   * @brief Run FeatureFinderPeptide pipeline on a centroided MSExperiment
   * @param[in,out] ms_peakmap Input centroided peak map (will be sorted)
   * @param[in] common_param Common parameters for FFP algorithms
   * @param[in] mtd_param MassTraceDetection parameters
   * @param[in] epd_param ElutionPeakDetection parameters
   * @param[in] ffp_param FeatureFindingPeptide parameters
   * @param[out] feat_map Output feature map
   * @param[out] traces_out Output mass traces (for accessing raw intensity data)
   * @return True on success, false on error
   */
  bool runFeatureFinderPeptide_(MSExperiment& ms_peakmap,
                                const Param& common_param,
                                Param mtd_param,
                                Param epd_param,
                                Param ffp_param,
                                FeatureMap& feat_map,
                                std::vector<MassTrace>& traces_out)
  {
    if (ms_peakmap.empty())
    {
      return true; // Nothing to process
    }

    // Ensure spectra are sorted by m/z
    ms_peakmap.sortSpectra(true);

    std::vector<MassTrace> m_traces;

    // Configure and run mass trace detection
    MassTraceDetection mtdet;
    mtd_param.insert("", common_param);
    mtd_param.remove("chrom_fwhm");
    mtdet.setParameters(mtd_param);
    mtdet.run(ms_peakmap, m_traces);

    if (m_traces.empty())
    {
      OPENMS_LOG_INFO << "No mass traces detected." << std::endl;
      return true;
    }

    // Configure and run elution peak detection
    std::vector<MassTrace> m_traces_final;
    if (epd_param.getValue("enabled").toBool())
    {
      std::vector<MassTrace> split_mtraces;
      epd_param.remove("enabled");
      epd_param.insert("", common_param);
      epd_param.remove("noise_threshold_int");
      ElutionPeakDetection epdet;
      epdet.setParameters(epd_param);
      epdet.detectPeaks(m_traces, split_mtraces);
      if (epdet.getParameters().getValue("width_filtering") == "auto")
      {
        m_traces_final.clear();
        epdet.filterByPeakWidth(split_mtraces, m_traces_final);
      }
      else
      {
        m_traces_final = split_mtraces;
      }
    }
    else
    {
      m_traces_final = m_traces;
      for (Size i = 0; i < m_traces_final.size(); ++i)
      {
        m_traces_final[i].estimateFWHM(false);
      }
      if (ffp_param.getValue("use_smoothed_intensities").toBool())
      {
        OPENMS_LOG_WARN << "Without EPD, smoothing is not supported. Setting 'use_smoothed_intensities' to false!" << std::endl;
        ffp_param.setValue("use_smoothed_intensities", "false");
      }
    }

    // Configure and run feature finding
    ffp_param.insert("", common_param);
    ffp_param.remove("noise_threshold_int");
    ffp_param.remove("chrom_peak_snr");
    ffp_param.setValue("report_chromatograms", "false");

    std::vector<std::vector<MSChromatogram>> feat_chromatograms;
    FeatureFindingPeptide ffpep;
    ffpep.setParameters(ffp_param);
    ffpep.run(m_traces_final, feat_map, feat_chromatograms);

    // Filter features with zero intensity
    auto intensity_zero = [](Feature& f) { return f.getIntensity() == 0; };
    feat_map.erase(std::remove_if(feat_map.begin(), feat_map.end(), intensity_zero), feat_map.end());

    // Output the mass traces for use in clustering (contains raw intensity data)
    traces_out = m_traces_final;

    OPENMS_LOG_INFO << "FFMetabo: " << m_traces_final.size() << " traces -> "
                    << feat_map.size() << " features" << std::endl;

    return true;
  }

  /**
   * @brief Run MassTraceExtractor pipeline on a centroided MSExperiment
   * @param[in,out] ms_peakmap Input centroided peak map (will be sorted)
   * @param[in] common_param Common parameters for MTE algorithms
   * @param[in] mtd_param MassTraceDetection parameters
   * @param[in] epd_param ElutionPeakDetection parameters
   * @param[out] traces_out Output mass traces
   * @return True on success, false on error
   */
  bool runMassTraceExtractor_(MSExperiment& ms_peakmap,
                              const Param& common_param,
                              Param mtd_param,
                              Param epd_param,
                              std::vector<MassTrace>& traces_out)
  {
    if (ms_peakmap.empty())
    {
      return true; // Nothing to process
    }

    // Ensure spectra are sorted by m/z
    ms_peakmap.sortSpectra(true);

    std::vector<MassTrace> m_traces;

    // Configure and run mass trace detection
    MassTraceDetection mtdet;
    mtd_param.insert("", common_param);
    mtd_param.remove("chrom_fwhm");
    mtdet.setParameters(mtd_param);
    mtdet.run(ms_peakmap, m_traces);

    if (m_traces.empty())
    {
      OPENMS_LOG_INFO << "No mass traces detected." << std::endl;
      return true;
    }

    // Configure and run elution peak detection if enabled
    bool use_epd = epd_param.getValue("enabled").toBool();

    if (use_epd)
    {
      std::vector<MassTrace> split_mtraces;
      epd_param.remove("enabled");
      epd_param.insert("", common_param);
      epd_param.remove("noise_threshold_int");
      ElutionPeakDetection epdet;
      epdet.setParameters(epd_param);
      epdet.detectPeaks(m_traces, split_mtraces);
      if (epdet.getParameters().getValue("width_filtering") == "auto")
      {
        traces_out.clear();
        epdet.filterByPeakWidth(split_mtraces, traces_out);
      }
      else
      {
        traces_out = std::move(split_mtraces);
      }
    }
    else
    {
      traces_out = std::move(m_traces);
    }

    // Remove empty traces
    traces_out.erase(
      std::remove_if(traces_out.begin(), traces_out.end(),
                     [](const MassTrace& t) { return t.getSize() == 0; }),
      traces_out.end());

    OPENMS_LOG_INFO << "MassTraceExtractor: " << m_traces.size() << " traces -> "
                    << traces_out.size() << " final traces" << std::endl;

    return true;
  }

  ExitCodes main_(int, const char**) override
  {
    const String in = getStringOption_("in");
    const bool save_precursors = getFlag_("save_unfragmented_precursors");
    const bool keep_ms1 = getFlag_("keep_ms1");
    const bool aggregate_scans = getFlag_("aggregate_across_scans");
    const Param ppim_params = getParam_().copy("PeakPickerIM:", true);
    const Param pphr_params = getParam_().copy("PeakPickerHiRes:", true);

    // FeatureFinderPeptide parameters (for MS1 and precursor data)
    const Param ffm_common_param = getParam_().copy("FeatureFinderPeptide:common:", true);
    Param ffm_mtd_param = getParam_().copy("FeatureFinderPeptide:mtd:", true);
    Param ffm_epd_param = getParam_().copy("FeatureFinderPeptide:epd:", true);
    Param ffm_ffp_param = getParam_().copy("FeatureFinderPeptide:ffp:", true);

    // MassTraceExtractor parameters (for MS1 and MS2 data)
    const Param mte_common_param = getParam_().copy("MassTraceExtractor:common:", true);
    Param mte_mtd_param = getParam_().copy("MassTraceExtractor:mtd:", true);
    Param mte_epd_param = getParam_().copy("MassTraceExtractor:epd:", true);

    // ClusterMassTraces parameters (for pseudo spectra generation)
    const Param cluster_param = getParam_().copy("ClusterMassTraces:", true);

#ifdef _OPENMP
    const int num_threads = getIntOption_("threads");
    const int threads_outer_loop = getIntOption_("threads_outer_loop");
#endif

    String out = getStringOption_("out");
    if (out.empty())
    {
      out = File::path(in) + "/" +
            File::basename(in) + "_diaWindows";
    }
    File::makeDir(out);

    // Start timing
    auto start_time = std::chrono::high_resolution_clock::now();

    // Detect input file type
    FileTypes::Type in_type = FileHandler::getTypeByFileName(in);

    // Bruker pre-extracted window data (populated only for .d files)
    bool is_bruker = false;
    bool bruker_im_centroiding = false; // true when Bruker built-in IM centroiding was applied during loading
    DiaWeaver::WindowedExperiments bruker_ms2_windows, bruker_ms1_windows, bruker_precursor_windows;

    // ------------------------------
    // Step 1: Determine DIA windows and IM info
    // ------------------------------
    DiaWeaver::WindowMap windows;
    DiaWeaver::IMInfo im_info;

#ifdef WITH_OPENTIMS
    if (in_type == FileTypes::BRUKER_TDF)
    {
      is_bruker = true;
      OPENMS_LOG_INFO << "Bruker .d file detected. Loading via BrukerTimsFile..." << std::endl;

      BrukerTimsFile tims_file;
      tims_file.setLogType(log_type_);
      auto bruker_config = getBrukerConfig_();

      bruker_im_centroiding = (bruker_config.ms1_centroid_mz_ppm > 0.0f && bruker_config.ms1_centroid_im_pct > 0.0f);

      PeakMap bruker_exp;
      tims_file.load(in, bruker_exp, bruker_config);

      DiaWeaver::determineWindows(bruker_exp, windows);
      im_info = DiaWeaver::determineIMInfo(bruker_exp, windows);

      // Pre-extract all windows into memory (already in-memory anyway)
      DiaWeaver::extractMS2Windows(bruker_exp, windows, bruker_ms2_windows,
                                    save_precursors ? &bruker_precursor_windows : nullptr);
      DiaWeaver::extractMS1Windows(bruker_exp, windows, bruker_ms1_windows);

      OPENMS_LOG_INFO << "Loaded " << windows.size() << " DIA windows from Bruker .d file." << std::endl;
    }
#endif

    OnDiscMSExperiment on_disc;
    if (!is_bruker)
    {
      // ------------------------------
      // mzML path: use OnDiscMSExperiment for memory-efficient metadata access
      // This determines DIA windows and IM info without loading all peak data
      // ------------------------------
      OPENMS_LOG_INFO << "Opening file for metadata access..." << std::endl;

      if (!on_disc.openFile(in))
      {
        OPENMS_LOG_ERROR << "Failed to open file as indexed mzML." << std::endl;
        return INPUT_FILE_NOT_FOUND;
      }

      DiaWeaver::determineWindows(on_disc, windows);
      im_info = DiaWeaver::determineIMInfo(on_disc, windows);
    }

    // Convert map to vector for OpenMP indexed access
    std::vector<std::pair<DiaWeaver::DIAWindow, std::vector<Size>>> window_vec(
      windows.begin(), windows.end());
    const Size total_windows = window_vec.size();

    // Note: OnDiscMSExperiment provides memory-efficient on-demand spectrum loading
    // Each thread will get its own copy via firstprivate (creates separate file handles)

    OPENMS_LOG_INFO << "Processing " << total_windows << " DIA windows";
#ifdef _OPENMP
    OPENMS_LOG_INFO << " using " << num_threads << " thread(s)";
#endif
    OPENMS_LOG_INFO << " with on-disc random access..." << std::endl;

    if (im_info.available)
    {
      if (bruker_im_centroiding)
      {
        OPENMS_LOG_INFO << "Bruker built-in IM centroiding was applied during .d loading. "
                        << "Skipping PeakPickerIM mobilogram method — spectra are already centroided." << std::endl;
      }
      else
      {
        OPENMS_LOG_INFO << "Ion mobility data detected. Using PeakPickerIM (mobilogram method)." << std::endl;
        if (aggregate_scans)
        {
          OPENMS_LOG_INFO << "Aggregation across scans enabled. Signal will be boosted before peak picking." << std::endl;
        }
      }
    }
    else
    {
      OPENMS_LOG_INFO << "No ion mobility data detected. Using PeakPickerHiRes." << std::endl;
      if (aggregate_scans)
      {
        OPENMS_LOG_WARN << "aggregate_across_scans is set but no IM data detected. Aggregation will be skipped." << std::endl;
      }
    }

    // Create output file path
    String output_filepath = out + "/pseudo_spectra.mzML";
    OPENMS_LOG_INFO << "Output: " << output_filepath << std::endl;

    Size spectra_written = 0;
    { // Scope for consumer - file is finalized when consumer is destroyed
    PlainMSDataWritingConsumer consumer(output_filepath);
    consumer.setExpectedSize(0, 0);  // Unknown count, will be determined during processing

    // Set source file information
    SourceFile source_file;
    source_file.setNameOfFile(File::basename(in));
    source_file.setPathToFile(File::path(in));
    ExperimentalSettings exp_settings;
    exp_settings.setSourceFiles({source_file});
    consumer.setExperimentalSettings(exp_settings);

    // Shared counter for unique spectrum native IDs (protected by critical section)
    Size spectrum_index = 0;

    // ------------------------------
    // Step 3: Process windows in parallel with nested parallelism
    // Outer loop: over DIA windows
    // Inner loop: over spectra within each window (for peak picking)
    // ------------------------------
    Size processed = 0;

#ifdef _OPENMP
    // Store total number of threads available
    const int total_threads = num_threads;

    // Calculate outer and inner thread counts
    int outer_threads = total_threads;
    int inner_threads = 1;

    if (threads_outer_loop > 0)
    {
      // User specified nested parallelism
      outer_threads = std::min(threads_outer_loop, total_threads);
      inner_threads = std::max(1, total_threads / outer_threads);
      omp_set_nested(1);
      omp_set_dynamic(0);
      OPENMS_LOG_INFO << "Using nested parallelism: " << outer_threads << " outer threads x "
                      << inner_threads << " inner threads for peak picking." << std::endl;
    }
    else
    {
      OPENMS_LOG_INFO << "Using " << outer_threads << " threads for window processing (no nested parallelism)." << std::endl;
    }

    omp_set_num_threads(outer_threads);
#pragma omp parallel for schedule(dynamic, 1) firstprivate(on_disc)
#endif
    for (SignedSize idx = 0; idx < static_cast<SignedSize>(total_windows); ++idx)
    {
      const DiaWeaver::DIAWindow& w = window_vec[idx].first;
      const std::vector<Size>& indices = window_vec[idx].second;

      // Thread-local buffers
      MSExperiment ms2_exp;
      MSExperiment ms1_exp;
      MSExperiment precursor_exp;
      std::vector<MassTrace> ms2_traces;  // MS2 mass traces for clustering

      // Extract spectra for this window
      if (is_bruker)
      {
#ifdef WITH_OPENTIMS
#pragma omp critical (bruker_window_access)
        {
          auto it_ms2 = bruker_ms2_windows.find(w);
          if (it_ms2 != bruker_ms2_windows.end()) ms2_exp = std::move(it_ms2->second);
          auto it_ms1 = bruker_ms1_windows.find(w);
          if (it_ms1 != bruker_ms1_windows.end()) ms1_exp = std::move(it_ms1->second);
          if (save_precursors)
          {
            auto it_prec = bruker_precursor_windows.find(w);
            if (it_prec != bruker_precursor_windows.end()) precursor_exp = std::move(it_prec->second);
          }
        }
#endif
      }
      else
      {
        // mzML path: on-demand from disk (each thread has its own file handle)
        DiaWeaver::extractSingleMS2Window(on_disc, w, indices, im_info, ms2_exp,
                                           save_precursors ? &precursor_exp : nullptr);
      }

      // Apply peak picking to MS2 spectra
      // Skip mobilogram method when Bruker built-in IM centroiding was applied — spectra are already centroided
      if (!bruker_im_centroiding && aggregate_scans && im_info.available)
      {
        // Parallel aggregation + peak picking into a separate output experiment.
        // ms2_exp stays immutable (raw data) so aggregateSpectrum always reads
        // original neighbors. Each thread writes to a unique index in ms2_picked.
        MSExperiment ms2_picked;
        ms2_picked.resize(ms2_exp.size());

#pragma omp parallel num_threads(inner_threads)
        {
          PeakPickerIM picker_im;
          picker_im.setParameters(ppim_params);

#pragma omp for schedule(dynamic, 1)
          for (SignedSize s = 0; s < static_cast<SignedSize>(ms2_exp.size()); ++s)
          {
            if (ms2_exp[s].getIMPeakType() != IMPeakType::IM_CENTROIDED)
            {
              MSSpectrum aggregated;
              aggregateSpectrum_(ms2_exp, static_cast<Size>(s), picker_im, aggregated);
              picker_im.pickIMTraces(aggregated);
              ms2_picked[s] = std::move(aggregated);
            }
            else
            {
              ms2_picked[s] = ms2_exp[s];
            }
          }
        }

        ms2_exp = std::move(ms2_picked);
      }
      else if (!bruker_im_centroiding)
      {
        // Standard parallel peak picking (no aggregation)
#pragma omp parallel num_threads(inner_threads)
        {
          PeakPickerIM picker_im;
          PeakPickerHiRes picker_hr;
          if (im_info.available)
          {
            picker_im.setParameters(ppim_params);
          }
          else
          {
            picker_hr.setParameters(pphr_params);
          }

#pragma omp for schedule(dynamic, 1)
          for (SignedSize s = 0; s < static_cast<SignedSize>(ms2_exp.size()); ++s)
          {
            if (im_info.available)
            {
              if (ms2_exp[s].getIMPeakType() != IMPeakType::IM_CENTROIDED)
              {
                picker_im.pickIMTraces(ms2_exp[s]);
              }
            }
            else
            {
              if (ms2_exp[s].getType(true) != SpectrumSettings::SpectrumType::CENTROID)
              {
                MSSpectrum picked;
                picker_hr.pick(ms2_exp[s], picked);
                ms2_exp[s] = std::move(picked);
              }
            }
          }
        }
      }

      // Run MassTraceExtractor on MS2 data to get fragment traces for clustering
      if (!ms2_exp.empty())
      {
        Param mte_mtd_copy = mte_mtd_param;
        Param mte_epd_copy = mte_epd_param;
        runMassTraceExtractor_(ms2_exp, mte_common_param, mte_mtd_copy, mte_epd_copy, ms2_traces);
      }

      // Apply peak picking to precursors
      if (save_precursors && !precursor_exp.empty())
      {
        if (aggregate_scans && im_info.available)
        {
          // Parallel aggregation + peak picking (same pattern as MS2)
          MSExperiment prec_picked;
          prec_picked.resize(precursor_exp.size());

#pragma omp parallel num_threads(inner_threads)
          {
            PeakPickerIM picker_im;
            picker_im.setParameters(ppim_params);

#pragma omp for schedule(dynamic, 1)
            for (SignedSize s = 0; s < static_cast<SignedSize>(precursor_exp.size()); ++s)
            {
              if (precursor_exp[s].getIMPeakType() != IMPeakType::IM_CENTROIDED)
              {
                MSSpectrum aggregated;
                aggregateSpectrum_(precursor_exp, static_cast<Size>(s), picker_im, aggregated);
                picker_im.pickIMTraces(aggregated);
                prec_picked[s] = std::move(aggregated);
              }
              else
              {
                prec_picked[s] = precursor_exp[s];
              }
            }
          }

          precursor_exp = std::move(prec_picked);
        }
        else
        {
          // Standard parallel peak picking (no aggregation)
#pragma omp parallel num_threads(inner_threads)
          {
            PeakPickerIM picker_im;
            PeakPickerHiRes picker_hr;
            if (im_info.available)
            {
              picker_im.setParameters(ppim_params);
            }
            else
            {
              picker_hr.setParameters(pphr_params);
            }

#pragma omp for schedule(dynamic, 1)
            for (SignedSize s = 0; s < static_cast<SignedSize>(precursor_exp.size()); ++s)
            {
              if (im_info.available)
              {
                if (precursor_exp[s].getIMPeakType() != IMPeakType::IM_CENTROIDED)
                {
                  picker_im.pickIMTraces(precursor_exp[s]);
                }
              }
              else
              {
                if (precursor_exp[s].getType(true) != SpectrumSettings::SpectrumType::CENTROID)
                {
                  MSSpectrum picked;
                  picker_hr.pick(precursor_exp[s], picked);
                  precursor_exp[s] = std::move(picked);
                }
              }
            }
          }
        }
        // Run FeatureFinderPeptide on precursor data, then cluster with MS2 traces
        FeatureMap precursor_features;
        std::vector<MassTrace> precursor_traces;
        Param mtd_copy = ffm_mtd_param;
        Param epd_copy = ffm_epd_param;
        Param ffp_copy = ffm_ffp_param;

        if (runFeatureFinderPeptide_(precursor_exp, ffm_common_param, mtd_copy, epd_copy, ffp_copy, precursor_features, precursor_traces)
            && !precursor_features.empty() && !ms2_traces.empty())
        {
          MSExperiment pseudo_spectra;

          ClusterMassTracesByPrecursor clusterFragments;
          clusterFragments.setParameters(cluster_param);
          clusterFragments.run(precursor_features, precursor_traces, ms2_traces, w.lower_mz, w.upper_mz, pseudo_spectra);

          if (!pseudo_spectra.empty())
          {
#pragma omp critical (write_spectra)
            {
              for (auto& spectrum : pseudo_spectra)
              {
                spectrum.setNativeID("scan=" + String(++spectrum_index));
                consumer.consumeSpectrum(spectrum);
              }
            }
          }
        }
      }

      // MS1 already loaded for Bruker path; for mzML, extract on-demand from disk
      if (!is_bruker)
      {
        DiaWeaver::extractSingleMS1Window(on_disc, w, im_info, ms1_exp);
      }

      // Apply peak picking to MS1 spectra
      // Skip mobilogram method when Bruker built-in IM centroiding was applied — spectra are already centroided
      if (!bruker_im_centroiding && aggregate_scans && im_info.available)
      {
        // Parallel aggregation + peak picking (same pattern as MS2)
        MSExperiment ms1_picked;
        ms1_picked.resize(ms1_exp.size());

#pragma omp parallel num_threads(inner_threads)
        {
          PeakPickerIM picker_im;
          picker_im.setParameters(ppim_params);

#pragma omp for schedule(dynamic, 1)
          for (SignedSize s = 0; s < static_cast<SignedSize>(ms1_exp.size()); ++s)
          {
            if (ms1_exp[s].getIMPeakType() != IMPeakType::IM_CENTROIDED)
            {
              MSSpectrum aggregated;
              aggregateSpectrum_(ms1_exp, static_cast<Size>(s), picker_im, aggregated);
              picker_im.pickIMTraces(aggregated);
              ms1_picked[s] = std::move(aggregated);
            }
            else
            {
              ms1_picked[s] = ms1_exp[s];
            }
          }
        }

        ms1_exp = std::move(ms1_picked);
      }
      else if (!bruker_im_centroiding)
      {
        // Standard parallel peak picking (no aggregation)
#pragma omp parallel num_threads(inner_threads)
        {
          PeakPickerIM picker_im;
          PeakPickerHiRes picker_hr;
          if (im_info.available)
          {
            picker_im.setParameters(ppim_params);
          }
          else
          {
            picker_hr.setParameters(pphr_params);
          }

#pragma omp for schedule(dynamic, 1)
          for (SignedSize s = 0; s < static_cast<SignedSize>(ms1_exp.size()); ++s)
          {
            if (im_info.available)
            {
              if (ms1_exp[s].getIMPeakType() != IMPeakType::IM_CENTROIDED)
              {
                picker_im.pickIMTraces(ms1_exp[s]);
              }
            }
            else
            {
              if (ms1_exp[s].getType(true) != SpectrumSettings::SpectrumType::CENTROID)
              {
                MSSpectrum picked;
                picker_hr.pick(ms1_exp[s], picked);
                ms1_exp[s] = std::move(picked);
              }
            }
          }
        }
      }

      // Write peak-picked MS1 spectra to output file if requested
      if (keep_ms1 && !ms1_exp.empty())
      {
#pragma omp critical (write_spectra)
        {
          for (auto& spectrum : ms1_exp)
          {
            spectrum.setNativeID("scan=" + String(++spectrum_index));
            spectrum.setType(SpectrumSettings::SpectrumType::CENTROID);
            consumer.consumeSpectrum(spectrum);
          }
        }
      }

      // Run FeatureFinderPeptide on MS1 data and cluster with MS2 traces to create pseudo spectra
      if (!ms1_exp.empty() && !ms2_traces.empty())
      {
        FeatureMap ms1_features;
        std::vector<MassTrace> ms1_traces;
        Param mtd_copy = ffm_mtd_param;
        Param epd_copy = ffm_epd_param;
        Param ffp_copy = ffm_ffp_param;

        if (runFeatureFinderPeptide_(ms1_exp, ffm_common_param, mtd_copy, epd_copy, ffp_copy, ms1_features, ms1_traces)
            && !ms1_features.empty())
        {
          MSExperiment pseudo_spectra;

          // cluster MS2 fragments by precursors
          ClusterMassTracesByPrecursor clusterFragments;
          clusterFragments.setParameters(cluster_param);
          clusterFragments.run(ms1_features, ms1_traces, ms2_traces, w.lower_mz, w.upper_mz, pseudo_spectra);

          if (!pseudo_spectra.empty())
          {
            // Write pseudo spectra to output file (thread-safe via critical section)
#pragma omp critical (write_spectra)
            {
              for (auto& spectrum : pseudo_spectra)
              {
                spectrum.setNativeID("scan=" + String(++spectrum_index));
                consumer.consumeSpectrum(spectrum);
              }
            }
          }
        }
      }

      // Progress logging (thread-safe)
#ifdef _OPENMP
#pragma omp critical (progress_log)
#endif
      {
        ++processed;
        OPENMS_LOG_INFO << "Processed window " << processed << "/" << total_windows
                        << " (m/z: " << w.lower_mz << "-" << w.upper_mz << ")" << std::endl;
      }
    }

#ifdef _OPENMP
    // Restore total thread count if nested parallelism was used
    if (threads_outer_loop > 0)
    {
      omp_set_num_threads(total_threads);
    }
#endif

    spectra_written = consumer.getNrSpectraWritten();
    OPENMS_LOG_INFO << "Finished processing all windows. Wrote " << spectra_written << " pseudo spectra." << std::endl;
    } // End of consumer scope - file is finalized here

    if (spectra_written == 0)
    {
      OPENMS_LOG_WARN << "No pseudo spectra were generated." << std::endl;
      return EXECUTION_OK;
    }

    // Load the file, sort spectra by RT, and rewrite
    OPENMS_LOG_INFO << "Sorting pseudo spectra by retention time..." << std::endl;

    MSExperiment pseudo_exp;
    MzMLFile mzml;
    mzml.load(output_filepath, pseudo_exp);

    pseudo_exp.sortSpectra(false);  // Sort by RT

    // Re-assign native IDs after sorting; also set MSn spectrum type and positive polarity
    // so downstream tools (e.g. FragPipe/MSFragger) recognise these as valid MS2 spectra.
    for (Size i = 0; i < pseudo_exp.size(); ++i)
    {
      pseudo_exp[i].setNativeID("scan=" + String(i + 1));
      pseudo_exp[i].getInstrumentSettings().setScanMode(InstrumentSettings::ScanMode::MSNSPECTRUM);
      pseudo_exp[i].getInstrumentSettings().setPolarity(IonSource::Polarity::POSITIVE);
    }

    // Configure numpress compression to reduce file size
    MSNumpressCoder::NumpressConfig npconfig_mz;
    npconfig_mz.setCompression("linear");
    npconfig_mz.numpressErrorTolerance = 0.0001;  // 0.01% error tolerance for m/z

    MSNumpressCoder::NumpressConfig npconfig_int;
    npconfig_int.setCompression("slof");  // Short logged float for intensities

    MSNumpressCoder::NumpressConfig npconfig_fda;
    npconfig_fda.setCompression("slof");  // For float data arrays (e.g., ion mobility)

    mzml.getOptions().setNumpressConfigurationMassTime(npconfig_mz);
    mzml.getOptions().setNumpressConfigurationIntensity(npconfig_int);
    mzml.getOptions().setNumpressConfigurationFloatDataArray(npconfig_fda);
    mzml.getOptions().setCompression(true);  // Also apply zlib compression

    mzml.store(output_filepath, pseudo_exp);

    // Calculate and print elapsed time
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
    auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration % std::chrono::hours(1));
    auto seconds = duration % std::chrono::minutes(1);

    OPENMS_LOG_INFO << "Done. Output (numpress compressed): " << output_filepath << std::endl;
    OPENMS_LOG_INFO << "Total processing time: "
                    << hours.count() << "h "
                    << minutes.count() << "m "
                    << seconds.count() << "s" << std::endl;

    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPDiaWeaver tool;
  return tool.main(argc, argv);
}
