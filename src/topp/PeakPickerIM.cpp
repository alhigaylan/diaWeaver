// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Author: Mohammed Alhigaylan $
// $Maintainer: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerIM.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>

#include <cmath>

#ifdef WITH_OPENTIMS
#include <OpenMS/FORMAT/BrukerTimsFile.h>
#endif

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_PeakPickerIM PeakPickerIM

@brief A tool for peak detection in the ion mobility dimension for mzML and Bruker .d files.

<center>
<table>
<tr>
<th ALIGN = "center"> pot. predecessor tools </td>
<td VALIGN="middle" ROWSPAN=2> &rarr; PeakPickerIM &rarr;</td>
<th ALIGN = "center"> pot. successor tools </td>
</tr>
<tr>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FileConverter </td>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool operating on MS peak data @n (in mzML format)</td>
</tr>
</table>
</center>

This tool applies peak picking in the ion mobility dimension to raw LC-IMS-MS data.
The input file can be an mzML file containing ion mobility data in concatenated format
(where each spectrum contains an ion mobility float data array) or a Bruker TimsTOF .d
directory (requires OpenMS built with WITH_OPENTIMS).

Three peak picking methods are available:
- @b mobilogram: Picks peaks along the ion mobility dimension using a peak picker.
- @b cluster: Clusters peaks in the ion mobility dimension.
- @b traces: Picks peaks using ion mobility elution profiles.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_PeakPickerIM.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_PeakPickerIM.html

For the parameters of the algorithm section see the algorithm documentation: @ref OpenMS::PeakPickerIM "PeakPickerIM"

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPPeakPickerIM : public TOPPBase
{
public:
  TOPPPeakPickerIM() :
      TOPPBase("PeakPickerIM", "Applies PeakPickerIM to an mzML or Bruker .d file", false)
  {}

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
    int center_ms_level = exp[center_idx].getMSLevel();

    std::vector<MSSpectrum> spectra_to_aggregate;
    std::vector<double> weights;

    // Search forward (including center), same MS level only
    for (Size j = center_idx; j < exp.size(); ++j)
    {
      if (exp[j].getMSLevel() != center_ms_level) continue;
      double rt_diff = exp[j].getRT() - center_rt;
      double weight = std::exp(factor * rt_diff * rt_diff);
      if (weight < cutoff && j != center_idx) break;
      spectra_to_aggregate.push_back(exp[j]);
      weights.push_back(weight);
    }

    // Search backward, same MS level only
    for (SignedSize j = static_cast<SignedSize>(center_idx) - 1; j >= 0; --j)
    {
      if (exp[j].getMSLevel() != center_ms_level) continue;
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
    registerInputFile_("in", "<file>", "", "Input file (mzML or Bruker .d)");
    setValidFormats_("in", { "mzML",
#ifdef WITH_OPENTIMS
      "d",
#endif
    });

    registerOutputFile_("out", "<file>", "", "Output mzML file");
    setValidFormats_("out", { "mzML" });

    registerStringOption_("processOption", "<name>", "inmemory",
                          "Whether to load all data and process them in-memory or process on-the-fly (lowmemory) without loading the whole file into memory first",
                          false, true);
    setValidStrings_("processOption", { "inmemory", "lowmemory" } );

    registerStringOption_("method", "<name>", "mobilogram",
                          "Method to pick peaks in IM dimension", false, true);
    setValidStrings_("method", { "mobilogram", "cluster", "traces" } );

    registerFlag_("aggregate_across_scans",
                  "If set, aggregate signal across neighboring scans using Gaussian weighting before peak picking. "
                  "This can improve signal-to-noise for low-intensity peaks.", false);

    addEmptyLine_();
    registerSubsection_("algorithm", "Algorithm parameters for PeakPickerIM (organized into pickIMTraces, pickIMCluster, pickIMElutionProfiles, aggregation).");

#ifdef WITH_OPENTIMS
    registerTOPPSubsection_("bruker", "Options for reading Bruker TimsTOF .d files (requires WITH_OPENTIMS)");
    registerStringOption_("bruker:export_mode", "<mode>", "frame", "Export mode: 'auto' detects DDA/DIA acquisition type, "
      "'frame' returns raw 4D frames without signal processing.", false, true);
    setValidStrings_("bruker:export_mode", {"auto", "frame"});
    registerDoubleOption_("bruker:calibration_tolerance", "<float>", 0.0, "m/z recalibration tolerance (0 = library default)", false, true);
    setMinFloat_("bruker:calibration_tolerance", 0.0);
    registerStringOption_("bruker:calibrate", "<toggle>", "false", "Enable m/z recalibration (may fail on some datasets)", false, true);
    setValidStrings_("bruker:calibrate", {"true", "false"});
    registerDoubleOption_("bruker:ms1_centroid_mz_ppm", "<float>", 0.0,
      "MS1 frame IM-centroiding m/z tolerance in ppm. Collapses the ion mobility dimension "
      "by aggregating neighboring peaks directly on the raw gridded data (Sage algorithm, Lazear 2023). "
      "Both this and ms1_centroid_mz_ppm must be > 0 to enable. Suggested value: 5.0. "
      "When enabled, this replaces the PeakPickerIM algorithm for MS1 frames.", false, true);
    setMinFloat_("bruker:ms1_centroid_mz_ppm", 0.0);
    registerDoubleOption_("bruker:ms1_centroid_im_pct", "<float>", 0.0,
      "MS1 frame IM-centroiding ion mobility tolerance in percent. Both this and ms1_centroid_mz_ppm "
      "must be > 0 to enable. Suggested value: 3.0.", false, true);
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
      "DIA MS2 frame aggregation: 0 = raw per-frame export, 1 = 3-frame sum, 2 = 5-frame sum. "
      "Switches the entire DIA-MS2 export pipeline regardless of ms2_centroid_algo.",
      false, true);
    setMinInt_("bruker:dia_ms2_n_neighbors", 0);
    registerIntOption_("bruker:dia_ms2_min_support", "<int>", 1,
      "DIA MS2 denoising: minimum occupied neighbor cells in a 3x3 (m/z x IM) grid to keep a point "
      "(center cell excluded from count). Applied after frame aggregation. Only effective when "
      "dia_ms2_n_neighbors > 0. Set to 0 to disable denoising (useful for pure centroiding "
      "without noise filtering).", false, true);
    setMinInt_("bruker:dia_ms2_min_support", 0);
    registerStringOption_("bruker:dia_ms2_centroid", "<toggle>", "false",
      "Apply 2D Gaussian smoothing + local maxima peak picking to the denoised DIA MS2 grid. "
      "Produces IM_CENTROIDED spectra with sub-bin (m/z, IM) precision. Only effective when dia_ms2_n_neighbors > 0.", false, true);
    setValidStrings_("bruker:dia_ms2_centroid", {"true", "false"});

    // Hill-based centroiding (IM-axis trace linking + valley splitting).
    registerStringOption_("bruker:ms1_centroid_algo", "<algo>", "off",
      "MS1 centroiding algorithm. 'off' = no IM-axis centroiding. 'greedy2d' = legacy 2D box "
      "clustering using ms1_centroid_mz_ppm/pct. 'hillbased' = IM-axis hill detection.",
      false, true);
    setValidStrings_("bruker:ms1_centroid_algo", {"off", "greedy2d", "hillbased"});
    registerStringOption_("bruker:ms2_centroid_algo", "<algo>", "off",
      "MS2 centroiding algorithm. Takes precedence over dia_ms2_centroid.",
      false, true);
    setValidStrings_("bruker:ms2_centroid_algo", {"off", "greedy2d", "hillbased"});
    registerDoubleOption_("bruker:ms2_centroid_mz_ppm", "<float>", 20.0,
      "HillBased MS2 m/z linking tolerance in ppm. Default 20.0 is DIA-PASEF-tuned. "
      "Set to 0 to refuse to run HillBased MS2 (the algo helper falls back to Off).",
      false, true);
    setMinFloat_("bruker:ms2_centroid_mz_ppm", 0.0);
    registerDoubleOption_("bruker:centroid_valley_factor", "<float>", 1.3,
      "HillBased: hill valley factor (hvf). Smaller = more aggressive splitting.",
      false, true);
    setMinFloat_("bruker:centroid_valley_factor", 1.0);
    registerIntOption_("bruker:ms1_centroid_min_hill_length", "<int>", 1,
      "HillBased MS1: minimum number of IM scans a hill must span. Default 1 keeps single-IM-scan ions.",
      false, true);
    setMinInt_("bruker:ms1_centroid_min_hill_length", 1);
    registerIntOption_("bruker:ms2_centroid_min_hill_length", "<int>", 2,
      "HillBased MS2: minimum number of IM scans a hill must span. Default 2 is "
      "DIA-PASEF-tuned. DDA-PASEF users should override to 1 (most DDA fragments "
      "are seen in only one IM scan; min=2 drops ~93% of DDA peaks).",
      false, true);
    setMinInt_("bruker:ms2_centroid_min_hill_length", 1);
    registerIntOption_("bruker:centroid_max_scan_gap", "<int>", 0,
      "HillBased: max consecutive empty IM scans a hill may bridge (0 = strict).",
      false, true);
    setMinInt_("bruker:centroid_max_scan_gap", 0);
    registerStringOption_("bruker:isotopic_prefilter", "<toggle>", "false",
      "MS1 + DIA-MS2 isotopic-partner prefilter applied after aggregation (or after raw "
      "extraction), before the centroider. Drops peaks lacking an isotopic partner at "
      "m/z ± C13C12_MASSDIFF / q (q in {1..5}) within ± isotopic_prefilter_tol_ppm AND "
      "|Δscan_id| <= 1. Not applied to DDA-MS2.",
      false, true);
    setValidStrings_("bruker:isotopic_prefilter", {"true", "false"});
    registerDoubleOption_("bruker:isotopic_prefilter_tol_ppm", "<float>", 50.0,
      "ppm tolerance for isotopic-partner matching by the prefilter.",
      false, true);
    setMinFloat_("bruker:isotopic_prefilter_tol_ppm", 0.0);

    registerStringOption_("bruker:expose_hill_bounds", "<toggle>", "false",
      "HillBased: attach hill bounding-box arrays per centroided spectrum for visual QC.",
      false, true);
    setValidStrings_("bruker:expose_hill_bounds", {"true", "false"});
#endif
  }

  Param getSubsectionDefaults_(const String& section) const override
  {
    if (section == "algorithm")
    {
      OpenMS::PeakPickerIM picker_defaults;
      Param p = picker_defaults.getDefaults();
      Param combined;
      combined.insert("pickIMTraces:",         p.copy("pickIMTraces:", true));
      combined.insert("pickIMCluster:",        p.copy("pickIMCluster:", true));
      combined.insert("pickIMElutionProfiles:",p.copy("pickIMElutionProfiles:", true));
      combined.insert("aggregation:",          p.copy("aggregation:", true));
      return combined;
    }
    return Param();
  }

#ifdef WITH_OPENTIMS
  BrukerTimsFile::Config getBrukerConfig_()
  {
    BrukerTimsFile::Config c;
    c.calibration_tolerance = getDoubleOption_("bruker:calibration_tolerance");
    c.calibrate = (getStringOption_("bruker:calibrate") == "true");
    String mode = getStringOption_("bruker:export_mode");
    if (mode == "frame") c.export_mode = BrukerTimsFile::Config::FRAME;
    else c.export_mode = BrukerTimsFile::Config::AUTO;
    c.ms1_centroid_mz_ppm = static_cast<float>(getDoubleOption_("bruker:ms1_centroid_mz_ppm"));
    c.ms1_centroid_im_pct = static_cast<float>(getDoubleOption_("bruker:ms1_centroid_im_pct"));
    c.ms1_n_neighbors         = getIntOption_("bruker:ms1_n_neighbors");
    c.ms1_min_support         = getIntOption_("bruker:ms1_min_support");
    c.ms1_max_rt_distance_sec = getDoubleOption_("bruker:ms1_max_rt_distance_sec");
    c.ms1_centroid_max_peaks  = getIntOption_("bruker:ms1_centroid_max_peaks");
    c.dia_ms2_n_neighbors = getIntOption_("bruker:dia_ms2_n_neighbors");
    c.dia_ms2_min_support = getIntOption_("bruker:dia_ms2_min_support");
    c.dia_ms2_centroid = (getStringOption_("bruker:dia_ms2_centroid") == "true");

    using CA = BrukerTimsFile::Config::CentroidAlgo;
    auto parse_algo = [](const String& s) {
      if (s == "greedy2d")  return CA::GREEDY2D;
      if (s == "hillbased") return CA::HILL_BASED;
      return CA::OFF;
    };
    c.ms1_centroid_algo            = parse_algo(getStringOption_("bruker:ms1_centroid_algo"));
    c.ms2_centroid_algo            = parse_algo(getStringOption_("bruker:ms2_centroid_algo"));
    c.ms2_centroid_mz_ppm          = static_cast<float>(getDoubleOption_("bruker:ms2_centroid_mz_ppm"));
    c.centroid_valley_factor       = getDoubleOption_("bruker:centroid_valley_factor");
    c.ms1_centroid_min_hill_length = static_cast<Size>(getIntOption_("bruker:ms1_centroid_min_hill_length"));
    c.ms2_centroid_min_hill_length = static_cast<Size>(getIntOption_("bruker:ms2_centroid_min_hill_length"));
    c.centroid_max_scan_gap        = static_cast<Size>(getIntOption_("bruker:centroid_max_scan_gap"));
    c.expose_hill_bounds           = (getStringOption_("bruker:expose_hill_bounds") == "true");
    c.isotopic_prefilter           = (getStringOption_("bruker:isotopic_prefilter") == "true");
    c.isotopic_prefilter_tol_ppm    = getDoubleOption_("bruker:isotopic_prefilter_tol_ppm");
    return c;
  }
#endif

  // -------------------- Low-memory consumer --------------------
  // Picks only MS1 spectra; MS2 spectra pass through as raw.
  class Consumer : public MSDataWritingConsumer
  {
  public:
    Consumer(String filename, const String& method, const PeakPickerIM& pp) :
        MSDataWritingConsumer(std::move(filename)), pp_(pp), method_(method) {}

    void processSpectrum_(MapType::SpectrumType& spectrum) override
    {
      // Only pick MS1 spectra; pass MS2 through as raw
      if (spectrum.getMSLevel() != 1) return;

      if (method_ == "mobilogram")
      {
        pp_.pickIMTraces(spectrum);
      }
      else if (method_ == "cluster")
      {
        pp_.pickIMCluster(spectrum);
      }
      else if (method_ == "traces")
      {
        pp_.pickIMElutionProfiles(spectrum);
      }
    }

    void processChromatogram_(MapType::ChromatogramType&) override {}

  private:
    PeakPickerIM pp_;
    String method_;
  };

  // -------------------- Format detection consumer (reads first MS1 spectrum only) --------------------
  class FormatDetector : public Interfaces::IMSDataConsumer
  {
  public:
    IMFormat detected_format = IMFormat::NONE;

    // Exception to abort after first spectrum (efficient early exit)
    struct FirstSpectrumRead : std::exception {};

    void consumeSpectrum(SpectrumType& s) override
    {
      if (s.getMSLevel() != 1) return; // Only check MS1 spectra (consistent with in-memory path)
      detected_format = IMTypes::determineIMFormat(s);
      throw FirstSpectrumRead(); // Abort after first MS1 spectrum
    }
    void consumeChromatogram(ChromatogramType&) override {}
    void setExperimentalSettings(const ExperimentalSettings&) override {}
    void setExpectedSize(size_t, size_t) override {}
  };

  // -------------------- Passthrough consumer (copies without processing) --------------------
  class PassthroughConsumer : public MSDataWritingConsumer
  {
  public:
    PassthroughConsumer(const String& filename) : MSDataWritingConsumer(filename) {}
    void processSpectrum_(MapType::SpectrumType&) override {} // No processing
    void processChromatogram_(MapType::ChromatogramType&) override {}
  };

  // -------------------- Helper for low-memory path --------------------
  ExitCodes doLowMemAlgorithm(const String& method, PeakPickerIM& pp,
                              const String& input_file, const String& output_file,
                              bool aggregate_scans)
  {
    if (aggregate_scans)
    {
      OPENMS_LOG_WARN << "Warning: 'aggregate_across_scans' is not supported in low-memory mode "
                      << "(requires random access to neighboring spectra). "
                      << "Proceeding without aggregation." << std::endl;
    }

    MzMLFile mzml;
    mzml.setLogType(log_type_);

    // Step 1: Detect IMFormat by reading only the first spectrum (minimal I/O)
    IMFormat im_format = IMFormat::NONE;
    {
      FormatDetector detector;
      try
      {
        mzml.transform(input_file, &detector);
        // If we reach here, file has no MS1 spectra - format stays NONE
      }
      catch (const FormatDetector::FirstSpectrumRead&)
      {
        im_format = detector.detected_format;
      }
    }

    // Step 2: Validate format
    if (im_format == IMFormat::IM_SPECTRUM)
    {
      OPENMS_LOG_ERROR << "Error: Input data has single drift time per spectrum (IM_SPECTRUM format). "
                       << "PeakPickerIM requires per-peak IM arrays (IM_PEAK format)." << std::endl;
      return ILLEGAL_PARAMETERS;
    }
    if (im_format == IMFormat::NONE)
    {
      OPENMS_LOG_WARN << "Warning: Input file does not contain ion mobility data. "
                      << "No peak picking will be performed." << std::endl;
      // Pass through unchanged - create fresh MzMLFile to avoid StopWatch state issues
      MzMLFile mzml_passthrough;
      mzml_passthrough.setLogType(log_type_);
      PassthroughConsumer passthrough(output_file);
      mzml_passthrough.transform(input_file, &passthrough);
      return EXECUTION_OK;
    }

    // Step 3: Proceed with streaming processing
    // Create a fresh MzMLFile object since the previous transform (for format detection)
    // may have left internal state (StopWatch) in an invalid state due to the exception
    MzMLFile mzml_writer;
    mzml_writer.setLogType(log_type_);

    Consumer pp_consumer(output_file, method, pp);
    pp_consumer.addDataProcessing(getProcessingInfo_(DataProcessing::PEAK_PICKING));
    mzml_writer.transform(input_file, &pp_consumer);

    return EXECUTION_OK;
  }

  ExitCodes main_(int, const char**) override
  {
    const String input_file  = getStringOption_("in");
    const String output_file = getStringOption_("out");
    const String process_opt = getStringOption_("processOption");
    const String method      = getStringOption_("method");
    const bool aggregate_scans = getFlag_("aggregate_across_scans");

    // Collect algorithm parameters from 'algorithm:' We strip and pass the remaining keys directly to PeakPickerIM.
    Param algo = getParam_().copy("algorithm:",true);

    PeakPickerIM picker;
    picker.setParameters(algo);

    // Detect input file type
    FileTypes::Type in_type = FileHandler::getType(input_file);

#ifdef WITH_OPENTIMS
    if (in_type == FileTypes::BRUKER_TDF)
    {
      if (process_opt == "lowmemory")
      {
        OPENMS_LOG_WARN << "Warning: 'lowmemory' processing is not yet supported for Bruker .d files. "
                        << "Data will be loaded fully into memory." << std::endl;
      }

      auto bruker_config = getBrukerConfig_();
      BrukerTimsFile tims_file;
      tims_file.setLogType(log_type_);

      PeakMap exp;
      tims_file.load(input_file, exp, bruker_config);

      // If built-in IM centroiding was enabled, BrukerTimsFile already produced
      // IM_CENTROIDED spectra — skip PeakPickerIM and write directly.
      bool builtin_centroiding = (bruker_config.ms1_centroid_mz_ppm > 0.0f
                                  && bruker_config.ms1_centroid_im_pct > 0.0f);
      if (builtin_centroiding)
      {
        OPENMS_LOG_INFO << "Built-in Bruker IM centroiding was applied during .d loading "
                        << "(ms1_centroid_mz_ppm=" << bruker_config.ms1_centroid_mz_ppm
                        << ", ms1_centroid_im_pct=" << bruker_config.ms1_centroid_im_pct
                        << "). Skipping PeakPickerIM algorithm." << std::endl;
        addDataProcessing_(exp, getProcessingInfo_(DataProcessing::PEAK_PICKING));
        MzMLFile().store(output_file, exp);
        return EXECUTION_OK;
      }

      // Check MS1 spectra for IM format
      IMFormat im_format = IMTypes::determineIMFormat(exp, 1);
      if (im_format == IMFormat::NONE)
      {
        OPENMS_LOG_WARN << "Warning: Input file does not contain ion mobility data. "
                        << "No peak picking will be performed." << std::endl;
        MzMLFile().store(output_file, exp);
        return EXECUTION_OK;
      }
      if (im_format == IMFormat::IM_SPECTRUM)
      {
        OPENMS_LOG_ERROR << "Error: Input data has single drift time per spectrum (IM_SPECTRUM format). "
                         << "PeakPickerIM requires per-peak IM arrays (IM_PEAK format). "
                         << "Try using bruker:export_mode=frame." << std::endl;
        return ILLEGAL_PARAMETERS;
      }

      std::exception_ptr first_error = nullptr;
#pragma omp parallel for
      for (SignedSize i = 0; i < static_cast<SignedSize>(exp.size()); ++i)
      {
        try
        {
          MSSpectrum& spectrum = exp[static_cast<Size>(i)];
          // Skip already-centroided spectra (e.g., DIA MS2 with bruker:dia_ms2_centroid=true)
          if (spectrum.getIMPeakType() == IMPeakType::IM_CENTROIDED) continue;
          if (method == "mobilogram")       picker.pickIMTraces(spectrum);
          else if (method == "cluster")     picker.pickIMCluster(spectrum);
          else if (method == "traces")      picker.pickIMElutionProfiles(spectrum);
        }
        catch (...)
        {
#pragma omp critical
          { if (!first_error) first_error = std::current_exception(); }
        }
      }
      if (first_error) std::rethrow_exception(first_error);

      addDataProcessing_(exp, getProcessingInfo_(DataProcessing::PEAK_PICKING));
      MzMLFile().store(output_file, exp);
      return EXECUTION_OK;
    }
#endif

    if (process_opt == "lowmemory")
    {
      return doLowMemAlgorithm(method, picker, input_file, output_file,
                               aggregate_scans);
    }
    else
    {
      PeakMap exp;
      MzMLFile mzml;
      mzml.load(input_file, exp);

      // Check MS1 spectra for IM format (PeakPickerIM works on per-peak IM data in MS1 frames)
      IMFormat im_format = IMTypes::determineIMFormat(exp, 1);
      if (im_format == IMFormat::NONE)
      {
        OPENMS_LOG_WARN << "Warning: Input file does not contain ion mobility data. "
                        << "No peak picking will be performed." << std::endl;
        mzml.store(output_file, exp);
        return EXECUTION_OK;
      }
      if (im_format == IMFormat::IM_SPECTRUM)
      {
        OPENMS_LOG_ERROR << "Error: Input data has single drift time per spectrum (IM_SPECTRUM format). "
                         << "PeakPickerIM requires per-peak IM arrays (IM_PEAK format)." << std::endl;
        return ILLEGAL_PARAMETERS;
      }

      if (aggregate_scans && method == "mobilogram")
      {
        // Parallel aggregation + peak picking into a separate output experiment.
        // exp stays immutable (raw data) so aggregateSpectrum_ always reads
        // original neighbors. Each thread writes to a unique index in picked_exp.
        OPENMS_LOG_INFO << "Aggregating signal across neighboring scans..." << std::endl;

        PeakMap picked_exp;
        picked_exp.resize(exp.size());

#pragma omp parallel
        {
          PeakPickerIM thread_picker;
          thread_picker.setParameters(algo);

#pragma omp for schedule(dynamic, 1)
          for (SignedSize i = 0; i < static_cast<SignedSize>(exp.size()); ++i)
          {
            Size idx = static_cast<Size>(i);

            if (exp[idx].getMSLevel() == 1 && exp[idx].containsIMData())
            {
              // Aggregate neighbors from immutable exp, pick, store in picked_exp
              MSSpectrum aggregated;
              aggregateSpectrum_(exp, idx, thread_picker, aggregated);
              thread_picker.pickIMTraces(aggregated);
              picked_exp[idx] = std::move(aggregated);
            }
            else
            {
              // Non-MS1: pass through as raw
              picked_exp[idx] = exp[idx];
            }
          }
        }

        exp = std::move(picked_exp);
      }
      else
      {
        std::exception_ptr first_error = nullptr;
#pragma omp parallel for
        for (SignedSize i = 0; i < static_cast<SignedSize>(exp.size()); ++i)
        {
          try
          {
            MSSpectrum& spectrum = exp[static_cast<Size>(i)];
            if (method == "mobilogram")       picker.pickIMTraces(spectrum);
            else if (method == "cluster")     picker.pickIMCluster(spectrum);
            else if (method == "traces")      picker.pickIMElutionProfiles(spectrum);
          }
          catch (...)
          {
#pragma omp critical
            { if (!first_error) first_error = std::current_exception(); }
          }
        }
        if (first_error) std::rethrow_exception(first_error);
      }

      // Annotate processing info (same as low-memory path)
      addDataProcessing_(exp, getProcessingInfo_(DataProcessing::PEAK_PICKING));

      mzml.store(output_file, exp);
      return EXECUTION_OK;
    }
  }
};

int main(int argc, const char** argv)
{
  TOPPPeakPickerIM tool;
  return tool.main(argc, argv);
}

/// @endcond

