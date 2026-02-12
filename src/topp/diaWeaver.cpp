// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Mohammed Alhigaylan $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/APPLICATIONS/diaWeaver.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerIM.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include <OpenMS/FEATUREFINDER/ElutionPeakDetection.h>
#include <OpenMS/FEATUREFINDER/FeatureFindingMetabo.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ClusterMassTracesByPrecursor.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/FORMAT/MSNumpressCoder.h>
#include <OpenMS/METADATA/SourceFile.h>

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
2. FeatureFinderMetabo on MS1 to detect monoisotopic peptide precursor features
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

  void registerOptionsAndFlags_() override
  {
    registerInputFile_(
      "in",
      "<file>",
      "",
      "Input DIA mzML file",
      true);
    setValidFormats_("in", {"mzML"});

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

    registerSubsection_("FeatureFinderMetabo", "Parameters for FeatureFinderMetabo algorithm (precursor detection on MS1 and unfragmented precursors)");

    registerSubsection_("MassTraceExtractor", "Parameters for MassTraceExtractor algorithm (mass trace detection on MS2 data)");

    registerSubsection_("ClusterMassTraces", "Parameters for clustering mass traces to create pseudo spectra");

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
    if (name == "FeatureFinderMetabo")
    {
      Param combined;

      // Common parameters for all FFM sub-algorithms
      Param p_com;
      p_com.setValue("noise_threshold_int", 60.0, "Intensity threshold below which peaks are regarded as noise.");
      p_com.setValue("chrom_peak_snr", 1.0, "Minimum signal-to-noise a mass trace should have.");
      p_com.setValue("chrom_fwhm", 5.0, "Expected chromatographic peak width (in seconds).");
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

      p_mtd.remove("noise_threshold_int");
      p_mtd.remove("chrom_peak_snr");
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

      // FeatureFindingMetabo parameters
      Param p_ffm = FeatureFindingMetabo().getDefaults();
      p_ffm.setValue("isotope_filtering_model", "peptides", "Use peptide isotope model for filtering");
      p_ffm.setValue("local_rt_range", 5.0, "RT range where to look for coeluting mass traces");
      p_ffm.setValue("local_mz_range", 3.0, "MZ range where to look for isotopic mass traces");
      p_ffm.setValue("local_im_range", 0.02, "IM range where to look for isotopic mass traces");
      p_ffm.setValue("charge_lower_bound", 2, "Lowest charge state to consider");
      p_ffm.setValue("charge_upper_bound", 4, "Highest charge state to consider");
      p_ffm.setValue("remove_single_traces", "true", "Remove unassembled traces (single traces).");
      p_ffm.setValue("mz_scoring_13C", "true", "Use the 13C isotope peak position (~1.003355 Da) as the expected shift in m/z for isotope mass traces (highly recommended for lipidomics!). Disable for general metabolites (as described in Kenar et al. 2014, MCP.)");
      p_ffm.setValue("use_smoothed_intensities", "false", "Use LOWESS intensities instead of raw intensities.");

      p_ffm.remove("chrom_fwhm");
      p_ffm.remove("report_chromatograms");
      combined.insert("ffm:", p_ffm);
      combined.setSectionDescription("ffm", "FeatureFinder parameters (assembling mass traces to charged features)");

      return combined;
    }
    if (name == "MassTraceExtractor")
    {
      Param combined;

      // Common parameters
      Param p_com;
      p_com.setValue("noise_threshold_int", 30.0, "Intensity threshold below which peaks are regarded as noise.");
      p_com.setValue("chrom_peak_snr", 1.0, "Minimum signal-to-noise a mass trace should have.");
      p_com.setValue("chrom_fwhm", 3.0, "Expected chromatographic peak width (in seconds).");

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

      p_mtd.remove("noise_threshold_int");
      p_mtd.remove("chrom_peak_snr");
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
      p.setValue("rt_tolerance", 5.0, "RT tolerance (in seconds) for matching up mass trace points during correlation.");
      return p;
    }
    return Param();
  }

  /**
   * @brief Run FeatureFinderMetabo pipeline on a centroided MSExperiment
   * @param[in,out] ms_peakmap Input centroided peak map (will be sorted)
   * @param[in] common_param Common parameters for FFM algorithms
   * @param[in] mtd_param MassTraceDetection parameters
   * @param[in] epd_param ElutionPeakDetection parameters
   * @param[in] ffm_param FeatureFindingMetabo parameters
   * @param[out] feat_map Output feature map
   * @param[out] traces_out Output mass traces (for accessing raw intensity data)
   * @return True on success, false on error
   */
  bool runFeatureFinderMetabo_(MSExperiment& ms_peakmap,
                               const Param& common_param,
                               Param mtd_param,
                               Param epd_param,
                               Param ffm_param,
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
      if (ffm_param.getValue("use_smoothed_intensities").toBool())
      {
        OPENMS_LOG_WARN << "Without EPD, smoothing is not supported. Setting 'use_smoothed_intensities' to false!" << std::endl;
        ffm_param.setValue("use_smoothed_intensities", "false");
      }
    }

    // Configure and run feature finding
    ffm_param.insert("", common_param);
    ffm_param.remove("noise_threshold_int");
    ffm_param.remove("chrom_peak_snr");
    ffm_param.setValue("report_chromatograms", "false");

    std::vector<std::vector<MSChromatogram>> feat_chromatograms;
    FeatureFindingMetabo ffmet;
    ffmet.setParameters(ffm_param);
    ffmet.run(m_traces_final, feat_map, feat_chromatograms);

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

    // FeatureFinderMetabo parameters (for MS1 and precursor data)
    const Param ffm_common_param = getParam_().copy("FeatureFinderMetabo:common:", true);
    Param ffm_mtd_param = getParam_().copy("FeatureFinderMetabo:mtd:", true);
    Param ffm_epd_param = getParam_().copy("FeatureFinderMetabo:epd:", true);
    Param ffm_ffm_param = getParam_().copy("FeatureFinderMetabo:ffm:", true);

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

    // ------------------------------
    // Step 1: Use OnDiscMSExperiment for memory-efficient metadata access
    // This determines DIA windows and IM info without loading all peak data
    // ------------------------------
    OPENMS_LOG_INFO << "Opening file for metadata access..." << std::endl;

    // Start timing
    auto start_time = std::chrono::high_resolution_clock::now();

    OnDiscMSExperiment on_disc;
    if (!on_disc.openFile(in))
    {
      OPENMS_LOG_ERROR << "Failed to open file as indexed mzML." << std::endl;
      return INPUT_FILE_NOT_FOUND;
    }

    // Determine DIA windows from MS2 metadata (efficient - no peak data loaded)
    DiaWeaver::WindowMap windows;
    DiaWeaver::determineWindows(on_disc, windows);

    // Determine IM info (loads only representative spectra)
    DiaWeaver::IMInfo im_info = DiaWeaver::determineIMInfo(on_disc, windows);

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
      OPENMS_LOG_INFO << "Ion mobility data detected. Using PeakPickerIM (mobilogram method)." << std::endl;
      if (aggregate_scans)
      {
        OPENMS_LOG_INFO << "Aggregation across scans enabled. Signal will be boosted before peak picking." << std::endl;
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

      // Extract MS2 (on-demand from disk - each thread has its own file handle)
      DiaWeaver::extractSingleMS2Window(on_disc, w, indices, im_info, ms2_exp,
                                         save_precursors ? &precursor_exp : nullptr);

      // Apply peak picking to MS2 spectra
      if (aggregate_scans && im_info.available)
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
            MSSpectrum aggregated;
            aggregateSpectrum_(ms2_exp, static_cast<Size>(s), picker_im, aggregated);
            picker_im.pickIMTraces(aggregated);
            ms2_picked[s] = std::move(aggregated);
          }
        }

        ms2_exp = std::move(ms2_picked);
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
          for (SignedSize s = 0; s < static_cast<SignedSize>(ms2_exp.size()); ++s)
          {
            if (im_info.available)
            {
              picker_im.pickIMTraces(ms2_exp[s]);
            }
            else
            {
              MSSpectrum picked;
              picker_hr.pick(ms2_exp[s], picked);
              ms2_exp[s] = std::move(picked);
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
              MSSpectrum aggregated;
              aggregateSpectrum_(precursor_exp, static_cast<Size>(s), picker_im, aggregated);
              picker_im.pickIMTraces(aggregated);
              prec_picked[s] = std::move(aggregated);
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
                picker_im.pickIMTraces(precursor_exp[s]);
              }
              else
              {
                MSSpectrum picked;
                picker_hr.pick(precursor_exp[s], picked);
                precursor_exp[s] = std::move(picked);
              }
            }
          }
        }
        // Run FeatureFinderMetabo on precursor data (results used internally, not saved)
        FeatureMap precursor_features;
        std::vector<MassTrace> precursor_traces;
        Param mtd_copy = ffm_mtd_param;
        Param epd_copy = ffm_epd_param;
        Param ffm_copy = ffm_ffm_param;
        runFeatureFinderMetabo_(precursor_exp, ffm_common_param, mtd_copy, epd_copy, ffm_copy, precursor_features, precursor_traces);
      }

      // Extract MS1 (on-demand from disk - each thread has its own file handle)
      DiaWeaver::extractSingleMS1Window(on_disc, w, im_info, ms1_exp);

      // Apply peak picking to MS1 spectra
      if (aggregate_scans && im_info.available)
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
            MSSpectrum aggregated;
            aggregateSpectrum_(ms1_exp, static_cast<Size>(s), picker_im, aggregated);
            picker_im.pickIMTraces(aggregated);
            ms1_picked[s] = std::move(aggregated);
          }
        }

        ms1_exp = std::move(ms1_picked);
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
          for (SignedSize s = 0; s < static_cast<SignedSize>(ms1_exp.size()); ++s)
          {
            if (im_info.available)
            {
              picker_im.pickIMTraces(ms1_exp[s]);
            }
            else
            {
              MSSpectrum picked;
              picker_hr.pick(ms1_exp[s], picked);
              ms1_exp[s] = std::move(picked);
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

      // Run FeatureFinderMetabo on MS1 data and cluster with MS2 traces to create pseudo spectra
      if (!ms1_exp.empty() && !ms2_traces.empty())
      {
        FeatureMap ms1_features;
        std::vector<MassTrace> ms1_traces;
        Param mtd_copy = ffm_mtd_param;
        Param epd_copy = ffm_epd_param;
        Param ffm_copy = ffm_ffm_param;

        if (runFeatureFinderMetabo_(ms1_exp, ffm_common_param, mtd_copy, epd_copy, ffm_copy, ms1_features, ms1_traces)
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

    // Re-assign native IDs after sorting
    for (Size i = 0; i < pseudo_exp.size(); ++i)
    {
      pseudo_exp[i].setNativeID("scan=" + String(i + 1));
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
