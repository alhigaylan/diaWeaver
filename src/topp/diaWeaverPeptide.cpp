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
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/METADATA/SourceFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FASTAFile.h>
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
#include <OpenMS/ANALYSIS/ID/ProSEAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/FragmentClaimRegistry.h>
#include <OpenMS/ANALYSIS/ID/BasicProteinInferenceAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/ANALYSIS/ID/IDMergerAlgorithm.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#ifdef WITH_OPENTIMS
#include <OpenMS/FORMAT/BrukerTimsFile.h>
#endif

#include <OpenMS/CONCEPT/Constants.h>
#include <cmath>
#include <fstream>
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
@page TOPP_diaWeaverPeptide diaWeaverPeptide

@brief DIA peptide identification: pseudo-spectrum generation with integrated ProSE database search.

This tool processes DIA (Data Independent Acquisition) data by running the full diaWeaver
preprocessing pipeline and then performing peptide identification using the ProSE search engine.
A single protein database index is built once and shared across all DIA windows.

The processing pipeline:
1. Determine DIA isolation windows from the input file
2. Build a ProSE fragment index from the FASTA database (once, shared across all windows)
3. For each DIA window (sequential; ProSE parallelises internally):
   a. Peak picking on MS2 spectra (PeakPickerIM for IM data, PeakPickerHiRes otherwise)
   b. MassTraceExtractor on MS2 to extract fragment mass traces
   c. Peak picking + FeatureFinderPeptide on MS1 to detect precursor features
   d. ClusterMassTracesByPrecursor to assemble pseudo MS/MS spectra
   e. ProSE database search against the shared fragment index
4. Merge per-window PSMs, run cross-window protein inference, apply FDR
5. Write merged idXML output

@note The outer window loop is intentionally sequential: ProSE already exploits all
available threads internally via OpenMP. Nested parallelism would cause thread contention.

@note Fragment claiming across shared MassTraces (the same MS2 fragment trace appearing
in multiple pseudo spectra) is not yet implemented. This is the intended next step:
ClusterMassTracesByPrecursor will emit fragment_trace_id IntDataArrays, an
ExperimentalFragmentIndex will map trace_id to spectrum indices within each window,
and a FragmentClaimRegistry will orchestrate iterative re-scoring.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_diaWeaverPeptide.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_diaWeaverPeptide.html
*/

/// @cond TOPPCLASSES

class TOPPDiaWeaverPeptide :
  public TOPPBase
{
public:
  TOPPDiaWeaverPeptide() :
    TOPPBase(
      "diaWeaverPeptide",
      "DIA peptide identification: pseudo-spectrum generation with integrated ProSE database search.",
      false)
  {
  }

protected:

  // -------------------------------------------------------------------------
  // Helper: aggregate a single spectrum with its RT neighbors (Gaussian weighting).
  // Identical to diaWeaver; kept here to avoid a shared helper dependency.
  // -------------------------------------------------------------------------
  static void aggregateSpectrum_(
    const MSExperiment& exp,
    Size center_idx,
    const PeakPickerIM& picker,
    MSSpectrum& out)
  {
    if (center_idx >= exp.size()) return;

    Param params = picker.getParameters();
    double fwhm    = (double)params.getValue("aggregation:rt_FWHM");
    double cutoff  = (double)params.getValue("aggregation:cutoff");
    double factor  = -4.0 * std::log(2.0) / (fwhm * fwhm);
    double center_rt = exp[center_idx].getRT();

    std::vector<MSSpectrum> spectra_to_aggregate;
    std::vector<double> weights;

    for (Size j = center_idx; j < exp.size(); ++j)
    {
      double w = std::exp(factor * (exp[j].getRT() - center_rt) * (exp[j].getRT() - center_rt));
      if (w < cutoff && j != center_idx) break;
      spectra_to_aggregate.push_back(exp[j]);
      weights.push_back(w);
    }
    for (SignedSize j = static_cast<SignedSize>(center_idx) - 1; j >= 0; --j)
    {
      double w = std::exp(factor * (exp[j].getRT() - center_rt) * (exp[j].getRT() - center_rt));
      if (w < cutoff) break;
      spectra_to_aggregate.push_back(exp[j]);
      weights.push_back(w);
    }

    double sum_w = 0.0;
    for (double w : weights) sum_w += w;
    for (double& w : weights) w /= sum_w;

    picker.aggregateScans(spectra_to_aggregate, weights, out);
  }

#ifdef WITH_OPENTIMS
  BrukerTimsFile::Config getBrukerConfig_()
  {
    BrukerTimsFile::Config c;
    c.calibration_tolerance      = getDoubleOption_("bruker:calibration_tolerance");
    c.calibrate                  = (getStringOption_("bruker:calibrate") == "true");
    String mode = getStringOption_("bruker:export_mode");
    c.export_mode = (mode == "frame") ? BrukerTimsFile::Config::FRAME : BrukerTimsFile::Config::AUTO;
    c.ms1_centroid_mz_ppm        = static_cast<float>(getDoubleOption_("bruker:ms1_centroid_mz_ppm"));
    c.ms1_centroid_im_pct        = static_cast<float>(getDoubleOption_("bruker:ms1_centroid_im_pct"));
    c.ms1_n_neighbors            = getIntOption_("bruker:ms1_n_neighbors");
    c.ms1_min_support            = getIntOption_("bruker:ms1_min_support");
    c.ms1_max_rt_distance_sec    = getDoubleOption_("bruker:ms1_max_rt_distance_sec");
    c.ms1_centroid_max_peaks     = getIntOption_("bruker:ms1_centroid_max_peaks");
    c.dia_ms2_n_neighbors        = getIntOption_("bruker:dia_ms2_n_neighbors");
    c.dia_ms2_min_support        = getIntOption_("bruker:dia_ms2_min_support");
    c.dia_ms2_centroid           = (getStringOption_("bruker:dia_ms2_centroid") == "true");
    return c;
  }
#endif

  // -------------------------------------------------------------------------
  // Run the FeatureFinderPeptide pipeline on a centroided MSExperiment.
  // Identical to diaWeaver.
  // -------------------------------------------------------------------------
  bool runFeatureFinderPeptide_(MSExperiment& ms_peakmap,
                                const Param& common_param,
                                Param mtd_param,
                                Param epd_param,
                                Param ffp_param,
                                FeatureMap& feat_map,
                                std::vector<MassTrace>& traces_out)
  {
    if (ms_peakmap.empty()) return true;

    ms_peakmap.sortSpectra(true);

    std::vector<MassTrace> m_traces;
    MassTraceDetection mtdet;
    mtd_param.insert("", common_param);
    mtd_param.remove("chrom_fwhm");
    mtdet.setParameters(mtd_param);
    mtdet.run(ms_peakmap, m_traces);

    if (m_traces.empty())
    {
      OPENMS_LOG_INFO << "No MS1 mass traces detected." << std::endl;
      return true;
    }

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
      for (auto& t : m_traces_final) t.estimateFWHM(false);
      if (ffp_param.getValue("use_smoothed_intensities").toBool())
      {
        OPENMS_LOG_WARN << "Without EPD, smoothing is not supported. Setting 'use_smoothed_intensities' to false." << std::endl;
        ffp_param.setValue("use_smoothed_intensities", "false");
      }
    }

    ffp_param.insert("", common_param);
    ffp_param.remove("noise_threshold_int");
    ffp_param.remove("chrom_peak_snr");
    ffp_param.setValue("report_chromatograms", "false");

    std::vector<std::vector<MSChromatogram>> feat_chromatograms;
    FeatureFindingPeptide ffpep;
    ffpep.setParameters(ffp_param);
    ffpep.run(m_traces_final, feat_map, feat_chromatograms);

    feat_map.erase(
      std::remove_if(feat_map.begin(), feat_map.end(),
                     [](Feature& f) { return f.getIntensity() == 0; }),
      feat_map.end());

    traces_out = m_traces_final;

    OPENMS_LOG_INFO << "FeatureFinderPeptide: " << m_traces_final.size() << " traces -> "
                    << feat_map.size() << " features" << std::endl;
    return true;
  }

  // -------------------------------------------------------------------------
  // Run the MassTraceExtractor pipeline on a centroided MSExperiment.
  // Identical to diaWeaver.
  // -------------------------------------------------------------------------
  bool runMassTraceExtractor_(MSExperiment& ms_peakmap,
                              const Param& common_param,
                              Param mtd_param,
                              Param epd_param,
                              std::vector<MassTrace>& traces_out)
  {
    if (ms_peakmap.empty()) return true;

    ms_peakmap.sortSpectra(true);

    std::vector<MassTrace> m_traces;
    MassTraceDetection mtdet;
    mtd_param.insert("", common_param);
    mtd_param.remove("chrom_fwhm");
    mtdet.setParameters(mtd_param);
    mtdet.run(ms_peakmap, m_traces);

    if (m_traces.empty())
    {
      OPENMS_LOG_INFO << "No MS2 mass traces detected." << std::endl;
      return true;
    }

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

    traces_out.erase(
      std::remove_if(traces_out.begin(), traces_out.end(),
                     [](const MassTrace& t) { return t.getSize() == 0; }),
      traces_out.end());

    OPENMS_LOG_INFO << "MassTraceExtractor: " << traces_out.size() << " final MS2 traces" << std::endl;
    return true;
  }

  // -------------------------------------------------------------------------
  // Peak-pick an MSExperiment in-place (IM or HiRes, with optional aggregation).
  // -------------------------------------------------------------------------
  void peakPickInPlace_(MSExperiment& exp,
                        const PeakPickerIM& picker_im,
                        const PeakPickerHiRes& picker_hr,
                        bool im_available,
                        bool aggregate_scans,
                        bool bruker_im_centroiding,
                        int inner_threads)
  {
    if (exp.empty()) return;

    if (!bruker_im_centroiding && aggregate_scans && im_available)
    {
      MSExperiment picked;
      picked.resize(exp.size());
#pragma omp parallel for num_threads(inner_threads) schedule(dynamic, 1)
      for (SignedSize s = 0; s < static_cast<SignedSize>(exp.size()); ++s)
      {
        if (exp[s].getIMPeakType() != IMPeakType::IM_CENTROIDED)
        {
          MSSpectrum agg;
          aggregateSpectrum_(exp, static_cast<Size>(s), picker_im, agg);
          PeakPickerIM local_picker = picker_im;
          local_picker.pickIMTraces(agg);
          picked[s] = std::move(agg);
        }
        else
        {
          picked[s] = exp[s];
        }
      }
      exp = std::move(picked);
    }
    else if (!bruker_im_centroiding)
    {
#pragma omp parallel for num_threads(inner_threads) schedule(dynamic, 1)
      for (SignedSize s = 0; s < static_cast<SignedSize>(exp.size()); ++s)
      {
        if (im_available)
        {
          if (exp[s].getIMPeakType() != IMPeakType::IM_CENTROIDED)
          {
            PeakPickerIM local_picker = picker_im;
            local_picker.pickIMTraces(exp[s]);
          }
        }
        else
        {
          if (exp[s].getType(true) != SpectrumSettings::SpectrumType::CENTROID)
          {
            MSSpectrum picked;
            picker_hr.pick(exp[s], picked);
            exp[s] = std::move(picked);
          }
        }
      }
    }
  }

  // -------------------------------------------------------------------------
  // Parameter registration
  // -------------------------------------------------------------------------
  void registerOptionsAndFlags_() override
  {
    // Input
    registerInputFile_("in", "<file>", "", "Input DIA file (mzML or Bruker .d)", true);
    setValidFormats_("in", {"mzML"
#ifdef WITH_OPENTIMS
      , "d"
#endif
    });

    registerInputFile_("database", "<file>", "",
      "Input protein sequence database in FASTA format. A single fragment index is built "
      "from this database and reused across all DIA windows.");
    setValidFormats_("database", ListUtils::create<String>("fasta"));

    // Output
    registerOutputFile_("out_idxml", "<file>", "",
      "Merged output idXML with cross-window protein inference and FDR filtering. "
      "All per-window PSMs are pooled before inference.");
    setValidFormats_("out_idxml", ListUtils::create<String>("idXML"));

    registerOutputFile_("out_mzml", "<file>", "",
      "Plain pseudo spectra mzML: all diaWeaver pseudo spectra written as-is, "
      "streamed to disk window-by-window without holding the full experiment in memory.", false);
    setValidFormats_("out_mzml", ListUtils::create<String>("mzML"));

    registerOutputFile_("out_annotated_mzml", "<file>", "",
      "Annotated pseudo spectra mzML: one spectrum per identified pseudo spectrum, "
      "containing only peaks matched to b/y ions of FDR-passing PSMs. Each peak carries "
      "'fragment_annotation' (StringDataArray), 'fragment_trace_id', and "
      "'fragment_window_id' (IntegerDataArrays) for ion accounting. Written after FDR filtering.", false);
    setValidFormats_("out_annotated_mzml", ListUtils::create<String>("mzML"));

    registerOutputFile_("out_orphan_mzml", "<file>", "",
      "Orphan peaks mzML: one spectrum per pseudo spectrum, containing only MS2 peaks "
      "NOT matched to any FDR-passing PSM. Each peak carries 'fragment_trace_id' and "
      "'fragment_window_id' IntegerDataArrays for downstream ion accounting. "
      "Written after FDR filtering.", false);
    setValidFormats_("out_orphan_mzml", ListUtils::create<String>("mzML"));

    registerOutputFile_("out_debug_tsv", "<file>", "",
      "Debug TSV: one row per PSM hit, written after FDR scoring but before score "
      "threshold filtering. Columns: spectrum_native_id, sequence, hyperscore, q_value. "
      "If FDR is not applied, q_value is reported as NA and hyperscore is the raw "
      "ln(hyperscore). Useful for inspecting the full ranked hit list.", false);
    setValidFormats_("out_debug_tsv", ListUtils::create<String>("tsv"));

    // Preprocessing flags (identical to diaWeaver)
    registerFlag_("save_unfragmented_precursors",
      "Also run FeatureFinderPeptide on peaks within the precursor isolation window.");

    registerFlag_("aggregate_across_scans",
      "Aggregate signal across neighboring scans using Gaussian weighting before peak picking "
      "(requires IM data).", false);

    // Subsections shared with diaWeaver
    registerSubsection_("PeakPickerIM",
      "Parameters for ion mobility peak picking (used when input has IM data)");
    registerSubsection_("PeakPickerHiRes",
      "Parameters for high-resolution peak picking (used when input has no IM data)");
    registerSubsection_("FeatureFinderPeptide",
      "Parameters for FeatureFinderPeptide algorithm (precursor detection on MS1)");
    registerSubsection_("MassTraceExtractor",
      "Parameters for MassTraceExtractor algorithm (fragment trace detection on MS2)");
    registerSubsection_("ClusterMassTraces",
      "Parameters for clustering mass traces into pseudo spectra");

    // ProSE search parameters
    Param search_algo_params_with_subsection;
    search_algo_params_with_subsection.insert("Search:", ProSEAlgorithm().getDefaults());
    registerFullParam_(search_algo_params_with_subsection);

    registerIntOption_("min_unique_fragments", "<int>", 3,
      "Minimum number of unclaimed fragment trace IDs required to retain a PSM after "
      "the iterative fragment claiming pass.", false);
    setMinInt_("min_unique_fragments", 1);

#ifdef WITH_OPENTIMS
    registerTOPPSubsection_("bruker", "Options for reading Bruker TimsTOF .d files (requires WITH_OPENTIMS)");
    registerStringOption_("bruker:export_mode", "<mode>", "auto",
      "Export mode: 'auto' detects DDA/DIA, 'frame' returns raw 4D frames.", false, true);
    setValidStrings_("bruker:export_mode", {"auto", "frame"});
    registerDoubleOption_("bruker:calibration_tolerance", "<float>", 0.0,
      "m/z recalibration tolerance (0 = library default)", false, true);
    setMinFloat_("bruker:calibration_tolerance", 0.0);
    registerStringOption_("bruker:calibrate", "<toggle>", "false",
      "Enable m/z recalibration", false, true);
    setValidStrings_("bruker:calibrate", {"true", "false"});
    registerDoubleOption_("bruker:ms1_centroid_mz_ppm", "<float>", 5.0,
      "MS1 IM-centroiding m/z tolerance in ppm.", false, true);
    setMinFloat_("bruker:ms1_centroid_mz_ppm", 0.0);
    registerDoubleOption_("bruker:ms1_centroid_im_pct", "<float>", 3.0,
      "MS1 IM-centroiding ion mobility tolerance in percent.", false, true);
    setMinFloat_("bruker:ms1_centroid_im_pct", 0.0);
    registerIntOption_("bruker:ms1_n_neighbors", "<int>", 0,
      "MS1 frame aggregation: adjacent frames on each side to sum.", false, true);
    setMinInt_("bruker:ms1_n_neighbors", 0);
    setMaxInt_("bruker:ms1_n_neighbors", 50);
    registerIntOption_("bruker:ms1_min_support", "<int>", 0,
      "MS1 denoising: minimum occupied neighbor cells to keep a point.", false, true);
    setMinInt_("bruker:ms1_min_support", 0);
    setMaxInt_("bruker:ms1_min_support", 8);
    registerDoubleOption_("bruker:ms1_max_rt_distance_sec", "<float>", 0.0,
      "RT cap for MS1 neighbor aggregation (0.0 = no cap).", false, true);
    setMinFloat_("bruker:ms1_max_rt_distance_sec", 0.0);
    registerIntOption_("bruker:ms1_centroid_max_peaks", "<int>", 100000,
      "Cap on centroided peaks retained per MS1 spectrum.", false, true);
    setMinInt_("bruker:ms1_centroid_max_peaks", 1);
    registerIntOption_("bruker:dia_ms2_n_neighbors", "<int>", 0,
      "DIA MS2 frame aggregation: adjacent frames on each side to sum.", false, true);
    setMinInt_("bruker:dia_ms2_n_neighbors", 0);
    registerIntOption_("bruker:dia_ms2_min_support", "<int>", 1,
      "DIA MS2 denoising: minimum occupied neighbor cells.", false, true);
    setMinInt_("bruker:dia_ms2_min_support", 1);
    registerStringOption_("bruker:dia_ms2_centroid", "<toggle>", "true",
      "Apply 2D Gaussian smoothing + peak picking to the DIA MS2 grid.", false, true);
    setValidStrings_("bruker:dia_ms2_centroid", {"true", "false"});
#endif

    registerIntOption_("threads", "<n>", 1,
      "Number of threads for peak picking and ProSE search.", false);
    setMinInt_("threads", 1);

    registerIntOption_("threads_outer_loop", "<n>", -1,
      "Number of threads for the outer loop (over DIA windows). Remaining threads are used for "
      "inner loop (peak picking within each window). Set to -1 to use all threads in the outer "
      "loop only (no nested parallelism). Example: with 24 total threads and 4 outer threads, "
      "each window gets 6 threads for peak picking.", false);
  }

  // -------------------------------------------------------------------------
  // Subsection defaults (identical to diaWeaver for shared subsections)
  // -------------------------------------------------------------------------
  Param getSubsectionDefaults_(const String& name) const override
  {
    if (name == "PeakPickerIM")
    {
      return PeakPickerIM().getDefaults();
    }
    if (name == "PeakPickerHiRes")
    {
      return PeakPickerHiRes().getDefaults();
    }
    if (name == "FeatureFinderPeptide")
    {
      Param combined;
      Param p_com;
      p_com.setValue("noise_threshold_int", 60.0, "Intensity threshold below which peaks are regarded as noise.");
      p_com.setValue("chrom_peak_snr", 1.0, "Minimum signal-to-noise a mass trace should have.");
      p_com.setValue("chrom_fwhm", 5.0, "Expected chromatographic peak width (in seconds).");
      combined.insert("common:", p_com);
      combined.setSectionDescription("common", "Common parameters for all other subsections");

      Param p_mtd = MassTraceDetection().getDefaults();
      p_mtd.setValue("mass_error_ppm", 7.0, "Allowed mass deviation (in ppm).");
      p_mtd.setValue("min_trace_length", 5.0, "Minimum expected length of a mass trace (in seconds).");
      p_mtd.setValue("ion_mobility_tolerance", 0.01, "Allowed ion mobility deviation (in 1/k0).");
      p_mtd.setValue("reestimate_mt_sd", "false", "Enables dynamic re-estimation of m/z variance.");
      p_mtd.setValue("quant_method", "max_height", "Quantification method for mass traces.");
      p_mtd.setValue("trace_termination_outliers", 2, "Cancel trace extension after this many consecutive empty spectra.");
      p_mtd.setValue("impute_zeros_missing_scans", "true", "Insert zero-intensity points at empty scan positions.");
      p_mtd.remove("noise_threshold_int");
      p_mtd.remove("chrom_peak_snr");
      p_mtd.remove("auto_noise_threshold");
      p_mtd.remove("noise_estimation_n_scans");
      p_mtd.remove("noise_estimation_percentile");
      combined.insert("mtd:", p_mtd);
      combined.setSectionDescription("mtd", "Mass Trace Detection parameters");

      Param p_epd;
      p_epd.setValue("enabled", "true", "Enable elution peak detection.");
      p_epd.setValue("width_filtering", "off", "Filter unlikely peak widths.");
      p_epd.setValidStrings("enabled", {"true", "false"});
      p_epd.insert("", ElutionPeakDetection().getDefaults());
      p_epd.remove("chrom_peak_snr");
      p_epd.remove("chrom_fwhm");
      combined.insert("epd:", p_epd);
      combined.setSectionDescription("epd", "Elution Profile Detection");

      Param p_ffp = FeatureFindingPeptide().getDefaults();
      p_ffp.setValue("local_rt_range", 5.0, "RT range where to look for coeluting mass traces");
      p_ffp.setValue("local_mz_range", 3.0, "MZ range where to look for isotopic mass traces");
      p_ffp.setValue("local_im_range", 0.02, "IM range where to look for isotopic mass traces");
      p_ffp.setValue("charge_lower_bound", 2, "Lowest charge state to consider");
      p_ffp.setValue("charge_upper_bound", 4, "Highest charge state to consider");
      p_ffp.setValue("remove_single_traces", "true", "Remove unassembled traces.");
      p_ffp.setValue("use_smoothed_intensities", "true", "Use Savitzky-Golay smoothed intensities.");
      p_ffp.setValue("mass_defect_filtering", "true", "Filter by peptide mass defect boundaries.");
      p_ffp.setValue("mass_defect_offset", 0.1, "Mass defect tolerance offset.");
      p_ffp.setValue("overlapping_features", "false", "Allow low-confidence hypotheses to reuse traces.");
      p_ffp.setValue("hypothesis_score_quantile", 0.5, "Score quantile threshold for low-confidence hypotheses.");
      p_ffp.setValue("rt_max_lag", 5, "Maximum lag for cross-correlation.");
      p_ffp.setValue("rt_min_pearson_correlation", 0.3, "Minimum Pearson correlation.");
      p_ffp.setValue("rt_peak_overlap_threshold", 0.3, "Minimum FWHM overlap proportion.");
      p_ffp.remove("chrom_fwhm");
      p_ffp.remove("report_chromatograms");
      combined.insert("ffp:", p_ffp);
      combined.setSectionDescription("ffp", "FeatureFindingPeptide parameters");

      return combined;
    }
    if (name == "MassTraceExtractor")
    {
      Param combined;
      Param p_com;
      p_com.setValue("noise_threshold_int", 30.0, "Intensity threshold below which peaks are regarded as noise.");
      p_com.setValue("chrom_peak_snr", 1.0, "Minimum signal-to-noise a mass trace should have.");
      p_com.setValue("chrom_fwhm", 3.0, "Expected chromatographic peak width (in seconds).");
      combined.insert("common:", p_com);
      combined.setSectionDescription("common", "Common parameters");

      Param p_mtd = MassTraceDetection().getDefaults();
      p_mtd.setValue("mass_error_ppm", 7.0, "Allowed mass deviation (in ppm).");
      p_mtd.setValue("min_trace_length", 2.0, "Minimum expected length of a mass trace (in seconds).");
      p_mtd.setValue("ion_mobility_tolerance", 0.01, "Allowed ion mobility deviation (in 1/k0).");
      p_mtd.setValue("reestimate_mt_sd", "false", "Dynamic re-estimation of m/z variance.");
      p_mtd.setValue("quant_method", "max_height", "Quantification method.");
      p_mtd.setValue("trace_termination_outliers", 2, "Cancel trace extension after this many consecutive empty spectra.");
      p_mtd.setValue("impute_zeros_missing_scans", "true", "Insert zero-intensity points at empty scan positions.");
      p_mtd.remove("noise_threshold_int");
      p_mtd.remove("chrom_peak_snr");
      p_mtd.remove("auto_noise_threshold");
      p_mtd.remove("noise_estimation_n_scans");
      p_mtd.remove("noise_estimation_percentile");
      combined.insert("mtd:", p_mtd);
      combined.setSectionDescription("mtd", "Mass Trace Detection parameters");

      Param p_epd;
      p_epd.setValue("enabled", "true", "Enable elution peak detection.");
      p_epd.setValue("width_filtering", "off", "Filter unlikely peak widths.");
      p_epd.setValidStrings("enabled", {"true", "false"});
      p_epd.insert("", ElutionPeakDetection().getDefaults());
      p_epd.remove("chrom_peak_snr");
      p_epd.remove("chrom_fwhm");
      combined.insert("epd:", p_epd);
      combined.setSectionDescription("epd", "Elution Profile Detection");

      return combined;
    }
    if (name == "ClusterMassTraces")
    {
      Param p;
      p.setValue("min_pearson_correlation", 0.3, "Minimal Pearson correlation to match elution profiles.");
      p.setValue("max_lag", 1, "Maximal lag for cross-correlation.");
      p.setValue("min_nr_ions", 30, "Minimal number of ions to report a spectrum.");
      p.setValue("max_rt_apex_difference", 5.0, "Maximal retention time difference at apex (seconds).");
      p.setValue("im_tolerance", 0.02, "Ion mobility tolerance for precursor-fragment matching.");
      p.setValue("nr_precursors_per_fragment", 50, "Maximum number of precursors a fragment can be assigned to.");
      p.setValue("rt_tolerance", 2.0, "RT tolerance for mass trace point correlation (seconds).");
      p.setValue("pearson_weight", 1.0, "Weight for Pearson correlation in combined score.");
      p.setValue("delta_rt_weight", 1.0, "Weight for delta RT in combined score.");
      p.setValue("delta_im_weight", 1.0, "Weight for delta IM in combined score.");
      p.setValue("max_nr_ions", 500, "Maximum fragment ions per output spectrum (0 = no limit).");
      p.setValue("assign_unassigned_to_all", "false", "Assign unmatched MS2 fragments to all nearby precursors.");
      p.setValidStrings("assign_unassigned_to_all", {"false", "true"});
      p.setValue("use_combined_scores", "true", "Rank assignments by combined score (Pearson + RT + IM).");
      p.setValidStrings("use_combined_scores", {"false", "true"});
      p.setValue("output_fragment_scores", "false", "Write per-fragment scores as FloatDataArrays.");
      p.setValidStrings("output_fragment_scores", {"false", "true"});
      p.setValue("smooth_ms1", "true", "Use EPD-smoothed intensities for MS1 elution profiles.");
      p.setValidStrings("smooth_ms1", {"false", "true"});
      p.setValue("smooth_ms2", "false", "Use EPD-smoothed intensities for MS2 elution profiles.");
      p.setValidStrings("smooth_ms2", {"false", "true"});
      return p;
    }
    return Param();
  }

  // -------------------------------------------------------------------------
  // Main
  // -------------------------------------------------------------------------
  ExitCodes main_(int, const char**) override
  {
    const String in         = getStringOption_("in");
    const String database   = getStringOption_("database");
    const String out_idxml          = getStringOption_("out_idxml");
    const String out_mzml           = getStringOption_("out_mzml");
    const String out_annotated_mzml = getStringOption_("out_annotated_mzml");
    const String out_orphan_mzml    = getStringOption_("out_orphan_mzml");
    const String out_debug_tsv      = getStringOption_("out_debug_tsv");

    if (out_idxml.empty())
    {
      OPENMS_LOG_ERROR << "No output specified. Provide -out_idxml." << std::endl;
      return ILLEGAL_PARAMETERS;
    }

    const bool save_precursors  = getFlag_("save_unfragmented_precursors");
    const bool aggregate_scans  = getFlag_("aggregate_across_scans");

    const Param ppim_params     = getParam_().copy("PeakPickerIM:", true);
    const Param pphr_params     = getParam_().copy("PeakPickerHiRes:", true);

    const Param ffm_common_param = getParam_().copy("FeatureFinderPeptide:common:", true);
    Param ffm_mtd_param          = getParam_().copy("FeatureFinderPeptide:mtd:", true);
    Param ffm_epd_param          = getParam_().copy("FeatureFinderPeptide:epd:", true);
    Param ffm_ffp_param          = getParam_().copy("FeatureFinderPeptide:ffp:", true);

    const Param mte_common_param = getParam_().copy("MassTraceExtractor:common:", true);
    Param mte_mtd_param          = getParam_().copy("MassTraceExtractor:mtd:", true);
    Param mte_epd_param          = getParam_().copy("MassTraceExtractor:epd:", true);

    const Param cluster_param    = getParam_().copy("ClusterMassTraces:", true);

    // ProSE search parameters: defer FDR to post-merge (always, since we have
    // multiple windows; per-window FDR would be statistically meaningless).
    Param search_params = getParam_().copy("Search:", true);
    const double user_psm_fdr     = static_cast<double>(search_params.getValue("FDR:PSM"));
    const double user_protein_fdr = static_cast<double>(search_params.getValue("FDR:protein"));
    const String decoy_prefix     = search_params.getValue("decoy_prefix").toString();
    search_params.setValue("FDR:PSM",     0.0);
    search_params.setValue("FDR:protein", 0.0);

    const Size min_unique_fragments_ = static_cast<Size>(getIntOption_("min_unique_fragments"));

#ifdef _OPENMP
    const int num_threads = getIntOption_("threads");
    const int threads_outer_loop = getIntOption_("threads_outer_loop");
    omp_set_num_threads(num_threads);
#endif

    auto start_time = std::chrono::high_resolution_clock::now();

    // ------------------------------------------------------------------
    // Step 1: Load FASTA and build the shared ProSE SearchContext once.
    // ------------------------------------------------------------------
    OPENMS_LOG_INFO << "[diaWeaverPeptide] Loading FASTA database: " << database << std::endl;
    std::vector<FASTAFile::FASTAEntry> fasta_db;
    FASTAFile().load(database, fasta_db);

    if (fasta_db.empty())
    {
      OPENMS_LOG_ERROR << "FASTA database is empty: " << database << std::endl;
      return INPUT_FILE_EMPTY;
    }
    OPENMS_LOG_INFO << "[diaWeaverPeptide] Loaded " << fasta_db.size() << " protein sequences." << std::endl;

    ProSEAlgorithm prose;
    prose.setLogType(log_type_);
    prose.setParameters(search_params);

    OPENMS_LOG_INFO << "[diaWeaverPeptide] Building fragment index (shared across all windows)..." << std::endl;
    ProSEAlgorithm::SearchContext ctx = prose.prepareContext(fasta_db);
    OPENMS_LOG_INFO << "[diaWeaverPeptide] Fragment index built. " << ctx.fragment_index.getPeptides().size()
                    << " peptide entries indexed." << std::endl;

    // ------------------------------------------------------------------
    // Step 2: Determine DIA windows.
    // ------------------------------------------------------------------
    FileTypes::Type in_type = FileHandler::getTypeByFileName(in);

    bool is_bruker = false;
    bool bruker_im_centroiding = false;
    DiaWeaver::WindowedExperiments bruker_ms2_windows, bruker_ms1_windows, bruker_precursor_windows;

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

      DiaWeaver::WindowMap windows_tmp;
      DiaWeaver::determineWindows(bruker_exp, windows_tmp);
      DiaWeaver::extractMS2Windows(bruker_exp, windows_tmp, bruker_ms2_windows,
                                    save_precursors ? &bruker_precursor_windows : nullptr);
      DiaWeaver::extractMS1Windows(bruker_exp, windows_tmp, bruker_ms1_windows);
      OPENMS_LOG_INFO << "Loaded " << windows_tmp.size() << " DIA windows from Bruker .d file." << std::endl;
    }
#endif

    DiaWeaver::WindowMap windows;
    DiaWeaver::IMInfo im_info;
    OnDiscMSExperiment on_disc;

    if (!is_bruker)
    {
      OPENMS_LOG_INFO << "Opening mzML for metadata access..." << std::endl;
      if (!on_disc.openFile(in))
      {
        OPENMS_LOG_ERROR << "Failed to open file as indexed mzML." << std::endl;
        return INPUT_FILE_NOT_FOUND;
      }
      DiaWeaver::determineWindows(on_disc, windows);
      im_info = DiaWeaver::determineIMInfo(on_disc, windows);
    }
#ifdef WITH_OPENTIMS
    else
    {
      // Re-determine windows from the Bruker experiment for indexed access
      // (bruker_ms2_windows already contains the per-window experiments;
      // build the window map from its keys so the loop below can iterate it).
      for (const auto& kv : bruker_ms2_windows)
      {
        windows[kv.first] = {};  // indices unused for Bruker path
      }
      // im_info is derived from the pre-extracted windows
      im_info.available = !bruker_ms2_windows.empty() &&
                          bruker_ms2_windows.begin()->first.hasIonMobility();
    }
#endif

    std::vector<std::pair<DiaWeaver::DIAWindow, std::vector<Size>>> window_vec(
      windows.begin(), windows.end());
    const Size total_windows = window_vec.size();

    OPENMS_LOG_INFO << "[diaWeaverPeptide] Processing " << total_windows << " DIA windows sequentially "
                    << "(ProSE parallelises internally)." << std::endl;

    if (im_info.available)
    {
      OPENMS_LOG_INFO << (bruker_im_centroiding
        ? "Bruker IM centroiding applied during loading; skipping PeakPickerIM mobilogram."
        : "Ion mobility data detected. Using PeakPickerIM.") << std::endl;
    }
    else
    {
      OPENMS_LOG_INFO << "No ion mobility data. Using PeakPickerHiRes." << std::endl;
    }

    // ------------------------------------------------------------------
    // Step 3: Per-window preprocessing + pseudo-spectrum generation + search.
    //
    // The outer loop is sequential: ProSE's scoreSpectraAgainstIndex_ already
    // uses #pragma omp parallel for internally, consuming all available threads.
    // ------------------------------------------------------------------

    // Accumulate per-window identification results for cross-window merging.
    std::vector<std::vector<ProteinIdentification>> all_prot_ids;
    std::vector<PeptideIdentificationList> all_pep_ids;
    // Keyed by native_id: preprocessed pseudo spectra (post-searchWithClaiming) held for
    // post-FDR orphan and annotated mzML writing.
    std::unordered_map<String, MSSpectrum> all_pseudo_spectra;
    // Accumulated claimed trace keys from all per-window FragmentClaimRegistries.
    // Used in orphan mzML building to exclude exactly those traces that were claimed.
    FragmentClaimRegistry all_claimed;

    // Pickers passed as const-ref to peakPickInPlace_(), which creates per-inner-thread
    // copies internally — safe to share across outer threads.
    PeakPickerIM    picker_im;  picker_im.setParameters(ppim_params);
    PeakPickerHiRes picker_hr;  picker_hr.setParameters(pphr_params);

#ifdef _OPENMP
    const int total_threads = num_threads;
    int outer_threads = total_threads;
    int inner_threads = 1;

    if (threads_outer_loop > 0)
    {
      outer_threads = std::min(threads_outer_loop, total_threads);
      inner_threads = std::max(1, total_threads / outer_threads);
      omp_set_nested(1);
      omp_set_dynamic(0);
      OPENMS_LOG_INFO << "[diaWeaverPeptide] Nested parallelism: "
                      << outer_threads << " outer x " << inner_threads << " inner threads." << std::endl;
    }
    else
    {
      OPENMS_LOG_INFO << "[diaWeaverPeptide] " << outer_threads
                      << " threads for window processing (no nested parallelism)." << std::endl;
    }
    omp_set_num_threads(outer_threads);
#else
    const int inner_threads = 1;
#endif

    Size processed = 0;
    Size total_pseudo_spectra = 0;

    // Set up streaming consumer for plain pseudo spectra (optional).
    // PlainMSDataWritingConsumer finalizes the mzML file when it goes out of scope,
    // so it must enclose the entire window loop.
    std::unique_ptr<PlainMSDataWritingConsumer> plain_consumer;
    if (!out_mzml.empty())
    {
      plain_consumer = std::make_unique<PlainMSDataWritingConsumer>(out_mzml);
      plain_consumer->setExpectedSize(0, 0);
      SourceFile sf;
      sf.setNameOfFile(File::basename(in));
      sf.setPathToFile(File::path(in));
      ExperimentalSettings es;
      es.setSourceFiles({sf});
      plain_consumer->setExperimentalSettings(es);
    }

    // prose is copied per-thread via firstprivate: each outer thread operates on its
    // own ProSEAlgorithm instance (identical parameters, independent mutable state).
    // ctx (SearchContext / FragmentIndex) is shared read-only — concurrent reads safe.
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) firstprivate(on_disc, prose)
#endif
    for (SignedSize idx = 0; idx < static_cast<SignedSize>(total_windows); ++idx)
    {
      const DiaWeaver::DIAWindow& w = window_vec[idx].first;
      const std::vector<Size>& indices = window_vec[idx].second;

      auto fmt_mz = [](double mz) -> String {
        char buf[32];
        std::snprintf(buf, sizeof(buf), "%.4f", mz);
        String s(buf);
        auto last = s.find_last_not_of('0');
        if (last != String::npos) s = s.substr(0, last + 1);
        if (!s.empty() && s.back() == '.') s += '0';
        return s;
      };
      const String window_label = "window_" + fmt_mz(w.lower_mz) + "_" + fmt_mz(w.upper_mz);

      MSExperiment ms2_exp, ms1_exp, precursor_exp;
      std::vector<MassTrace> ms2_traces;

      // --- 3a. Extract raw spectra for this window ---
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
        DiaWeaver::extractSingleMS2Window(on_disc, w, indices, im_info, ms2_exp,
                                           save_precursors ? &precursor_exp : nullptr);
      }

      // --- 3b. Peak pick MS2 ---
      peakPickInPlace_(ms2_exp, picker_im, picker_hr, im_info.available,
                       aggregate_scans, bruker_im_centroiding, inner_threads);

      // --- 3c. MassTraceExtractor on MS2 ---
      {
        Param mte_mtd_copy = mte_mtd_param;
        Param mte_epd_copy = mte_epd_param;
        runMassTraceExtractor_(ms2_exp, mte_common_param, mte_mtd_copy, mte_epd_copy, ms2_traces);
      }

      if (ms2_traces.empty())
      {
#pragma omp critical (progress_log)
        {
          OPENMS_LOG_INFO << "[diaWeaverPeptide] No MS2 traces for window " << window_label
                          << ". Skipping." << std::endl;
          ++processed;
        }
        continue;
      }

      // --- 3d. Extract MS1 ---
      if (!is_bruker)
      {
        DiaWeaver::extractSingleMS1Window(on_disc, w, im_info, ms1_exp);
      }

      // --- 3e. Peak pick MS1 ---
      peakPickInPlace_(ms1_exp, picker_im, picker_hr, im_info.available,
                       aggregate_scans, bruker_im_centroiding, inner_threads);

      // --- 3f. Optional precursor window processing ---
      MSExperiment pseudo_precursor;
      if (save_precursors && !precursor_exp.empty())
      {
        peakPickInPlace_(precursor_exp, picker_im, picker_hr, im_info.available,
                         aggregate_scans, bruker_im_centroiding, inner_threads);

        FeatureMap precursor_features;
        std::vector<MassTrace> precursor_traces;
        Param mtd_copy = ffm_mtd_param, epd_copy = ffm_epd_param, ffp_copy = ffm_ffp_param;
        if (runFeatureFinderPeptide_(precursor_exp, ffm_common_param, mtd_copy, epd_copy, ffp_copy,
                                     precursor_features, precursor_traces)
            && !precursor_features.empty())
        {
          ClusterMassTracesByPrecursor clusterer;
          clusterer.setParameters(cluster_param);
          clusterer.run(precursor_features, precursor_traces, ms2_traces,
                        w.lower_mz, w.upper_mz, pseudo_precursor, static_cast<Int>(idx));
        }
      }

      // --- 3g. FeatureFinderPeptide on MS1 ---
      FeatureMap ms1_features;
      std::vector<MassTrace> ms1_traces;
      MSExperiment pseudo_ms1;
      {
        Param mtd_copy = ffm_mtd_param, epd_copy = ffm_epd_param, ffp_copy = ffm_ffp_param;
        if (runFeatureFinderPeptide_(ms1_exp, ffm_common_param, mtd_copy, epd_copy, ffp_copy,
                                     ms1_features, ms1_traces)
            && !ms1_features.empty())
        {
          ClusterMassTracesByPrecursor clusterer;
          clusterer.setParameters(cluster_param);
          clusterer.run(ms1_features, ms1_traces, ms2_traces,
                        w.lower_mz, w.upper_mz, pseudo_ms1, static_cast<Int>(idx));
        }
      }

      // Merge pseudo spectra from MS1 and optional precursor paths into one PeakMap.
      MSExperiment pseudo_spectra;
      for (auto& s : pseudo_ms1)      pseudo_spectra.addSpectrum(std::move(s));
      for (auto& s : pseudo_precursor) pseudo_spectra.addSpectrum(std::move(s));

      if (pseudo_spectra.empty())
      {
#pragma omp critical (progress_log)
        {
          OPENMS_LOG_INFO << "[diaWeaverPeptide] No pseudo spectra generated for window "
                          << window_label << ". Skipping." << std::endl;
          ++processed;
        }
        continue;
      }

      // Assign globally unique native IDs. idx is loop-private (OMP parallel for),
      // so "window=<idx>_scan=<si>" is collision-free across threads.
      for (Size si = 0; si < pseudo_spectra.size(); ++si)
        pseudo_spectra[si].setNativeID("window=" + String(idx) + "_scan=" + String(si));

      // Stream plain pseudo spectra to disk (thread-safe via critical section).
      if (plain_consumer)
      {
#pragma omp critical (write_spectra)
        {
          for (auto& spec : pseudo_spectra)
            plain_consumer->consumeSpectrum(spec);
        }
      }

      // --- 3h. ProSE search with iterative fragment claiming ---
      // prose is a per-thread firstprivate copy; ctx (FragmentIndex) is shared read-only.
      std::vector<ProteinIdentification> window_prot_ids;
      PeptideIdentificationList window_pep_ids;
      FragmentClaimRegistry window_registry;

      const ProSEAlgorithm::ExitCodes ec =
        prose.searchWithClaiming(pseudo_spectra, ctx, window_prot_ids, window_pep_ids,
                                 window_registry, min_unique_fragments_);

      // Always collect preprocessed pseudo spectra for post-FDR orphan/annotated output.
      // Merge this window's claim registry into the global one for orphan building.
#pragma omp critical (collect_pseudo_spectra)
      {
        for (auto& spec : pseudo_spectra)
          all_pseudo_spectra[spec.getNativeID()] = spec;
        all_claimed.merge(window_registry);
      }

      if (ec != ProSEAlgorithm::ExitCodes::EXECUTION_OK || window_prot_ids.empty() || window_pep_ids.empty())
      {
#pragma omp critical (progress_log)
        {
          if (ec != ProSEAlgorithm::ExitCodes::EXECUTION_OK)
            OPENMS_LOG_WARN << "[diaWeaverPeptide] ProSE non-OK exit code ("
                            << static_cast<int>(ec) << ") for window " << window_label << std::endl;
          else
            OPENMS_LOG_INFO << "[diaWeaverPeptide] No PSMs for window " << window_label << std::endl;
          ++processed;
        }
        continue;
      }

      window_prot_ids[0].setPrimaryMSRunPath({window_label});
      window_prot_ids[0].getSearchParameters().db = database;

#pragma omp critical (results_collect)
      {
        total_pseudo_spectra += pseudo_spectra.size();
        OPENMS_LOG_INFO << "[diaWeaverPeptide] Window " << window_label << ": "
                        << window_pep_ids.size() << " PSMs, "
                        << window_prot_ids[0].getHits().size() << " protein hits." << std::endl;
        all_prot_ids.push_back(std::move(window_prot_ids));
        all_pep_ids.push_back(std::move(window_pep_ids));
        ++processed;
      }
    } // end window loop

#ifdef _OPENMP
    if (threads_outer_loop > 0)
      omp_set_num_threads(total_threads);
#endif

    OPENMS_LOG_INFO << "[diaWeaverPeptide] All " << total_windows << " windows processed. "
                    << total_pseudo_spectra << " pseudo spectra searched across all windows." << std::endl;

    if (all_pep_ids.empty())
    {
      OPENMS_LOG_WARN << "[diaWeaverPeptide] No PSMs identified across any window. "
                      << "Writing empty output." << std::endl;
      std::vector<ProteinIdentification> empty_prot;
      PeptideIdentificationList empty_pep;
      FileHandler().storeIdentifications(out_idxml, empty_prot, empty_pep, {FileTypes::IDXML});
      return EXECUTION_OK;
    }

    // ------------------------------------------------------------------
    // Step 4: Cross-window merge, protein inference, FDR.
    // ------------------------------------------------------------------
    OPENMS_LOG_INFO << "[diaWeaverPeptide] Merging " << all_pep_ids.size()
                    << " per-window result sets..." << std::endl;

    IDMergerAlgorithm merger;
    for (Size i = 0; i < all_prot_ids.size(); ++i)
    {
      merger.insertRuns(all_prot_ids[i], all_pep_ids[i]);
    }
    ProteinIdentification merged_proteins;
    PeptideIdentificationList merged_peptides;
    merger.returnResultsAndClear(merged_proteins, merged_peptides);

    std::vector<ProteinIdentification> merged_prot_ids = {std::move(merged_proteins)};

    // Record all window labels as the merged run's MS run paths.
    StringList all_window_labels;
    for (const auto& pv : all_prot_ids)
    {
      if (!pv.empty())
      {
        StringList run_paths;
        pv[0].getPrimaryMSRunPath(run_paths);
        all_window_labels.insert(all_window_labels.end(), run_paths.begin(), run_paths.end());
      }
    }
    merged_prot_ids[0].setPrimaryMSRunPath(all_window_labels);
    merged_prot_ids[0].getSearchParameters().db = database;

    // Check whether decoys are present (needed for FDR).
    bool has_decoys = false;
    for (const auto& ph : merged_prot_ids[0].getHits())
    {
      if (ph.metaValueExists("target_decoy") &&
          ph.getMetaValue("target_decoy").toString() == "decoy")
      {
        has_decoys = true;
        break;
      }
    }

    // Helper: write debug TSV of all scored PSM hits.
    // Written after FDR assigns q-values (before score-threshold filtering) so
    // the full ranked list is visible. If FDR was not applied, q_value is "NA".
    auto write_debug_tsv = [&](bool has_qvalues)
    {
      if (out_debug_tsv.empty()) return;
      std::ofstream tsv(out_debug_tsv);
      if (!tsv)
      {
        OPENMS_LOG_WARN << "[diaWeaverPeptide] Cannot open debug TSV for writing: "
                        << out_debug_tsv << std::endl;
        return;
      }
      tsv << "spectrum_native_id\tRT\tIM\tprecursor_mz\tsequence\ttarget_decoy\thyperscore\tq_value\n";
      const String orig_score_key = "ln(hyperscore)_score";
      for (const auto& pi : merged_peptides)
      {
        const String& native_id = pi.getSpectrumReference();
        const double  rt        = pi.getRT();
        const double  prec_mz   = pi.getMZ();
        const String  im_str    = pi.metaValueExists(Constants::UserParam::IM)
                                    ? String(static_cast<double>(pi.getMetaValue(Constants::UserParam::IM)))
                                    : "NA";
        for (const auto& hit : pi.getHits())
        {
          double hyperscore, qval;
          if (has_qvalues && hit.metaValueExists(orig_score_key))
          {
            hyperscore = static_cast<double>(hit.getMetaValue(orig_score_key));
            qval       = hit.getScore();
          }
          else
          {
            hyperscore = hit.getScore();
            qval       = -1.0;
          }
          tsv << native_id << "\t"
              << rt        << "\t"
              << im_str    << "\t"
              << prec_mz   << "\t"
              << hit.getSequence().toString() << "\t"
              << (hit.isDecoy() ? "decoy" : "target") << "\t"
              << hyperscore << "\t";
          if (qval >= 0.0) tsv << qval;
          else             tsv << "NA";
          tsv << "\n";
        }
      }
      OPENMS_LOG_INFO << "[diaWeaverPeptide] Debug PSM TSV written to " << out_debug_tsv
                      << std::endl;
    };

    // PSM-level FDR — applied before protein inference so BPIA sees only
    // confident PSMs (mirrors ProSE's per-file FDR-before-merge approach).
    if (user_psm_fdr > 0.0)
    {
      if (!has_decoys)
      {
        OPENMS_LOG_WARN << "[diaWeaverPeptide] FDR:PSM requested but no decoy PSMs found. "
                        << "Enable Search:decoys or provide a FASTA with decoy proteins. "
                        << "Skipping PSM FDR filtering." << std::endl;
        write_debug_tsv(false);
      }
      else
      {
        OPENMS_LOG_INFO << "[diaWeaverPeptide] Applying PSM FDR at " << user_psm_fdr * 100
                        << "% threshold..." << std::endl;
        FalseDiscoveryRate fdr_tool;
        Param fdr_params = fdr_tool.getParameters();
        fdr_params.setValue("use_all_hits", "true");
        fdr_params.setValue("add_decoy_peptides", "true");
        fdr_tool.setParameters(fdr_params);
        fdr_tool.apply(merged_peptides);
        write_debug_tsv(true);
        IDFilter::filterHitsByScore(merged_peptides, user_psm_fdr);
        IDFilter::removeEmptyIdentifications(merged_peptides);
        OPENMS_LOG_INFO << "[diaWeaverPeptide] " << merged_peptides.size()
                        << " PSMs retained after PSM FDR." << std::endl;
      }
    }
    else
    {
      write_debug_tsv(false);
    }

    // If protein FDR is not requested, strip decoys before BPIA so protein
    // inference runs on target-only PSMs (same as ProSE merged flow).
    // If protein FDR is requested, decoys must survive so picked-protein FDR
    // has target/decoy protein pairs to work with.
    if (user_protein_fdr == 0.0)
    {
      IDFilter::removeDecoyHits(merged_prot_ids);
      IDFilter::removeDecoyHits(merged_peptides);
      IDFilter::removeEmptyIdentifications(merged_peptides);
      IDFilter::removeUnreferencedProteins(merged_prot_ids, merged_peptides);
    }

    // Cross-window protein inference via BasicProteinInferenceAlgorithm.
    OPENMS_LOG_INFO << "[diaWeaverPeptide] Running cross-window protein inference on "
                    << merged_peptides.size() << " PSMs..." << std::endl;
    BasicProteinInferenceAlgorithm bpia;
    bpia.run(merged_peptides, merged_prot_ids);

    // Protein-level picked-protein FDR
    if (user_protein_fdr > 0.0)
    {
      if (!has_decoys)
      {
        OPENMS_LOG_WARN << "[diaWeaverPeptide] FDR:protein requested but no decoy proteins found. "
                        << "Skipping protein FDR filtering." << std::endl;
      }
      else
      {
        OPENMS_LOG_INFO << "[diaWeaverPeptide] Applying picked-protein FDR at "
                        << user_protein_fdr * 100 << "% threshold..." << std::endl;
        FalseDiscoveryRate fdr_tool;
        fdr_tool.applyPickedProteinFDR(merged_prot_ids[0], decoy_prefix, true);
        IDFilter::filterHitsByScore(merged_prot_ids, user_protein_fdr);
        OPENMS_LOG_INFO << "[diaWeaverPeptide] " << merged_prot_ids[0].getHits().size()
                        << " proteins retained after protein FDR." << std::endl;
      }
    }

    // Remove decoys and dangling references from final output.
    IDFilter::removeDecoyHits(merged_prot_ids);
    IDFilter::removeDecoyHits(merged_peptides);
    IDFilter::removeEmptyIdentifications(merged_peptides);
    IDFilter::removeUnreferencedProteins(merged_prot_ids, merged_peptides);
    IDFilter::removeDanglingProteinReferences(merged_peptides, merged_prot_ids);
    // Protein groups built by BPIA may reference decoy accessions that were removed above.
    IDFilter::updateProteinGroups(merged_prot_ids[0].getProteinGroups(), merged_prot_ids[0].getHits());
    IDFilter::updateProteinGroups(merged_prot_ids[0].getIndistinguishableProteins(), merged_prot_ids[0].getHits());

    // ------------------------------------------------------------------
    // Step 5: Write output.
    // ------------------------------------------------------------------
    OPENMS_LOG_INFO << "[diaWeaverPeptide] Final result: " << merged_peptides.size() << " PSMs, "
                    << merged_prot_ids[0].getHits().size() << " proteins." << std::endl;

    // --- 5a. Annotated pseudo spectra mzML ---
    // Group surviving PeptideIdentifications by spectrum native ID (getSpectrumReference()).
    // All peptide hits sharing the same native ID originated from the same pseudo spectrum
    // and are written as one output spectrum containing only b/y-matched peaks.

    // Helper: return the fragment_trace_id and fragment_window_id IntegerDataArrays.
    // Throws if either is absent — both are unconditionally written by
    // ClusterMassTracesByPrecursor and must always be present on every pseudo spectrum.
    struct TraceArrayPair
    {
      const MSSpectrum::IntegerDataArray* trace_id_arr;
      const MSSpectrum::IntegerDataArray* window_id_arr;
    };
    auto getTraceArrays = [](const MSSpectrum& s) -> TraceArrayPair
    {
      TraceArrayPair result{nullptr, nullptr};
      for (const auto& arr : s.getIntegerDataArrays())
      {
        if      (arr.getName() == "fragment_trace_id")  result.trace_id_arr = &arr;
        else if (arr.getName() == "fragment_window_id") result.window_id_arr = &arr;
      }
      if (result.trace_id_arr == nullptr || result.window_id_arr == nullptr)
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Pseudo spectrum '" + s.getNativeID() +
            "' is missing 'fragment_trace_id' or 'fragment_window_id' IntegerDataArray. "
            "These arrays must be present on all pseudo spectra produced by ClusterMassTracesByPrecursor.");
      return result;
    };

    // Build native_id → PSM index map needed for annotated mzML label construction.
    std::map<String, std::vector<Size>> native_id_to_pids;
    if (!out_annotated_mzml.empty())
    {
      for (Size pid = 0; pid < merged_peptides.size(); ++pid)
      {
        const String& ref = merged_peptides[pid].getSpectrumReference();
        if (!ref.empty()) native_id_to_pids[ref].push_back(pid);
      }
    }

    if (!out_annotated_mzml.empty())
    {
      MSExperiment annotated_exp;

      for (const auto& [native_id, pid_indices] : native_id_to_pids)
      {
        auto spec_it = all_pseudo_spectra.find(native_id);
        if (spec_it == all_pseudo_spectra.end()) continue;
        const MSSpectrum& orig = spec_it->second;

        // Build mz → annotation label map from all FDR-passing hits for this spectrum.
        std::map<double, String> mz_to_label;
        double spectrum_rt = orig.getRT();

        for (Size pid : pid_indices)
        {
          for (const PeptideHit& hit : merged_peptides[pid].getHits())
          {
            const String seq = hit.getSequence().toString();
            for (const PeptideHit::PeakAnnotation& pa : hit.getPeakAnnotations())
            {
              String label = seq + "-" + pa.annotation;
              if (pa.charge > 1) label += "+" + String(pa.charge);
              auto it = mz_to_label.find(pa.mz);
              if (it == mz_to_label.end()) mz_to_label[pa.mz] = label;
              else                         it->second += ";" + label;
            }
          }
        }
        if (mz_to_label.empty()) continue;

        // Build sorted mz vector from original spectrum for binary search.
        std::vector<double> orig_mzs;
        orig_mzs.reserve(orig.size());
        for (const auto& pk : orig) orig_mzs.push_back(pk.getMZ());

        const auto [trace_id_arr, window_id_arr] = getTraceArrays(orig);

        MSSpectrum out_spec;
        out_spec.setNativeID(native_id);
        out_spec.setRT(spectrum_rt);
        out_spec.setMSLevel(2);
        if (!orig.getPrecursors().empty())
          out_spec.getPrecursors().push_back(orig.getPrecursors()[0]);

        MSSpectrum::StringDataArray annot_arr;  annot_arr.setName("fragment_annotation");
        MSSpectrum::IntegerDataArray out_tid;   out_tid.setName("fragment_trace_id");
        MSSpectrum::IntegerDataArray out_wid;   out_wid.setName("fragment_window_id");

        for (const auto& [mz, label] : mz_to_label)
        {
          // Find the peak in the original spectrum by exact mz.
          auto lb = std::lower_bound(orig_mzs.begin(), orig_mzs.end(), mz);
          if (lb == orig_mzs.end() || *lb != mz) continue;
          const Size pk_idx = static_cast<Size>(lb - orig_mzs.begin());

          Peak1D pk; pk.setMZ(mz); pk.setIntensity(orig[pk_idx].getIntensity());
          out_spec.push_back(pk);
          annot_arr.push_back(label);
          out_tid.push_back((*trace_id_arr)[pk_idx]);
          out_wid.push_back((*window_id_arr)[pk_idx]);
        }

        if (out_spec.empty()) continue;
        out_spec.getStringDataArrays().push_back(std::move(annot_arr));
        if (!out_tid.empty()) out_spec.getIntegerDataArrays().push_back(std::move(out_tid));
        if (!out_wid.empty()) out_spec.getIntegerDataArrays().push_back(std::move(out_wid));
        annotated_exp.addSpectrum(std::move(out_spec));
      }

      annotated_exp.sortSpectra(true);
      OPENMS_LOG_INFO << "[diaWeaverPeptide] Writing " << annotated_exp.size()
                      << " annotated pseudo spectra to: " << out_annotated_mzml << std::endl;
      MzMLFile().store(out_annotated_mzml, annotated_exp);
    }

    if (!out_orphan_mzml.empty())
    {
      MSExperiment orphan_exp;

      for (const auto& [native_id, orig] : all_pseudo_spectra)
      {
        const auto [trace_id_arr, window_id_arr] = getTraceArrays(orig);

        MSSpectrum orphan_spec;
        orphan_spec.setNativeID(native_id);
        orphan_spec.setRT(orig.getRT());
        orphan_spec.setMSLevel(2);
        if (!orig.getPrecursors().empty())
          orphan_spec.getPrecursors() = orig.getPrecursors();

        MSSpectrum::IntegerDataArray out_tid; out_tid.setName("fragment_trace_id");
        MSSpectrum::IntegerDataArray out_wid; out_wid.setName("fragment_window_id");

        for (Size i = 0; i < orig.size(); ++i)
        {
          const auto key = FragmentClaimRegistry::makeKey((*window_id_arr)[i], (*trace_id_arr)[i]);
          if (all_claimed.isClaimed(key)) continue;
          orphan_spec.push_back(orig[i]);
          out_tid.push_back((*trace_id_arr)[i]);
          out_wid.push_back((*window_id_arr)[i]);
        }

        if (orphan_spec.empty()) continue;
        if (!out_tid.empty()) orphan_spec.getIntegerDataArrays().push_back(std::move(out_tid));
        if (!out_wid.empty()) orphan_spec.getIntegerDataArrays().push_back(std::move(out_wid));
        orphan_exp.addSpectrum(std::move(orphan_spec));
      }

      orphan_exp.sortSpectra(true);
      OPENMS_LOG_INFO << "[diaWeaverPeptide] Writing " << orphan_exp.size()
                      << " orphan pseudo spectra to: " << out_orphan_mzml << std::endl;
      MzMLFile().store(out_orphan_mzml, orphan_exp);
    }

    // --- 5b. Peptide/protein identifications ---
    OPENMS_LOG_INFO << "[diaWeaverPeptide] Writing output: " << out_idxml << std::endl;
    FileHandler().storeIdentifications(out_idxml, merged_prot_ids, merged_peptides,
                                       {FileTypes::IDXML});

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    auto hours   = std::chrono::duration_cast<std::chrono::hours>(duration);
    auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration % std::chrono::hours(1));
    auto seconds = duration % std::chrono::minutes(1);
    OPENMS_LOG_INFO << "[diaWeaverPeptide] Total processing time: "
                    << hours.count() << "h "
                    << minutes.count() << "m "
                    << seconds.count() << "s" << std::endl;

    return EXECUTION_OK;
  }
};

/// @endcond

int main(int argc, const char** argv)
{
  TOPPDiaWeaverPeptide tool;
  return tool.main(argc, argv);
}
