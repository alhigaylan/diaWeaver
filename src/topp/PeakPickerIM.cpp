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

#ifdef WITH_OPENTIMS
#include <OpenMS/FORMAT/BrukerTimsFile.h>
#endif

#include <cmath>

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
      "Both this and ms1_centroid_im_pct must be > 0 to enable. Suggested value: 5.0. "
      "When enabled, this replaces the PeakPickerIM algorithm for MS1 frames.", false, true);
    setMinFloat_("bruker:ms1_centroid_mz_ppm", 0.0);
    registerDoubleOption_("bruker:ms1_centroid_im_pct", "<float>", 0.0,
      "MS1 frame IM-centroiding ion mobility tolerance in percent. Both this and ms1_centroid_mz_ppm "
      "must be > 0 to enable. Suggested value: 3.0.", false, true);
    setMinFloat_("bruker:ms1_centroid_im_pct", 0.0);
    registerIntOption_("bruker:dia_ms2_n_neighbors", "<int>", 0,
      "DIA MS2 frame aggregation: number of adjacent frames on each side to sum per SWATH window. "
      "0 = disabled (raw export), 1 = 3-frame sum, 2 = 5-frame sum. "
      "Boosts signal by summing intensity across neighboring RT frames, then removes isolated noise.", false, true);
    setMinInt_("bruker:dia_ms2_n_neighbors", 0);
    registerIntOption_("bruker:dia_ms2_min_support", "<int>", 1,
      "DIA MS2 denoising: minimum occupied neighbor cells in a 3x3 (m/z x IM) grid to keep a point "
      "(center cell excluded from count). Applied after frame aggregation. Only effective when dia_ms2_n_neighbors > 0.", false, true);
    setMinInt_("bruker:dia_ms2_min_support", 1);
    registerStringOption_("bruker:dia_ms2_centroid", "<toggle>", "false",
      "Apply 2D Gaussian smoothing + local maxima peak picking to the denoised DIA MS2 grid. "
      "Produces IM_CENTROIDED spectra with sub-bin (m/z, IM) precision. Only effective when dia_ms2_n_neighbors > 0.", false, true);
    setValidStrings_("bruker:dia_ms2_centroid", {"true", "false"});
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
    c.dia_ms2_n_neighbors = getIntOption_("bruker:dia_ms2_n_neighbors");
    c.dia_ms2_min_support = getIntOption_("bruker:dia_ms2_min_support");
    c.dia_ms2_centroid = (getStringOption_("bruker:dia_ms2_centroid") == "true");
    return c;
  }
#endif

  // -------------------- Low-memory consumer --------------------
  class Consumer : public MSDataWritingConsumer
  {
  public:
    Consumer(String filename, const String& method, const PeakPickerIM& pp) :
        MSDataWritingConsumer(std::move(filename)), pp_(pp), method_(method) {}

    void processSpectrum_(MapType::SpectrumType& spectrum) override
    {
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

    // Step 1: Detect IMFormat by reading only the first MS1 spectrum (minimal I/O)
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
      // Pass through unchanged
      PassthroughConsumer passthrough(output_file);
      mzml.transform(input_file, &passthrough);
      return EXECUTION_OK;
    }

    // Step 3: Proceed with streaming processing
    Consumer pp_consumer(output_file, method, pp);
    pp_consumer.addDataProcessing(getProcessingInfo_(DataProcessing::PEAK_PICKING));
    mzml.transform(input_file, &pp_consumer);

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
      return doLowMemAlgorithm(method, picker, input_file, output_file, aggregate_scans);
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

        std::exception_ptr first_error = nullptr;
#pragma omp parallel
        {
          PeakPickerIM thread_picker;
          thread_picker.setParameters(algo);

#pragma omp for schedule(dynamic, 1)
          for (SignedSize i = 0; i < static_cast<SignedSize>(exp.size()); ++i)
          {
            try
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
            catch (...)
            {
#pragma omp critical
              { if (!first_error) first_error = std::current_exception(); }
            }
          }
        }
        if (first_error) std::rethrow_exception(first_error);

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

            if (method == "mobilogram")
            {
              picker.pickIMTraces(spectrum);
            }
            else if (method == "cluster")
            {
              picker.pickIMCluster(spectrum);
            }
            else if (method == "traces")
            {
              picker.pickIMElutionProfiles(spectrum);
            }
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
