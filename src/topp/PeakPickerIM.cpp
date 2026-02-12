// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Author: Mohammed Alhigaylan $
// $Maintainer: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerIM.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>

#include <cmath>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_PeakPickerIM PeakPickerIM

@brief A tool for peak detection in the ion mobility dimension for mzML files.

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
The input mzML file should contain ion mobility data in concatenated format
(where each spectrum contains an ion mobility float data array).

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
      TOPPBase("PeakPickerIM", "Applies PeakPickerIM to an mzML file", false)
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
    registerInputFile_("in", "<file>", "", "Input mzML file");
    setValidFormats_("in", { "mzML" });

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

  // -------------------- Format detection consumer (reads first spectrum only) --------------------
  class FormatDetector : public Interfaces::IMSDataConsumer
  {
  public:
    IMFormat detected_format = IMFormat::NONE;

    // Exception to abort after first spectrum (efficient early exit)
    struct FirstSpectrumRead : std::exception {};

    void consumeSpectrum(SpectrumType& s) override
    {
      detected_format = IMTypes::determineIMFormat(s);
      throw FirstSpectrumRead(); // Abort after reading first spectrum
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
        // If we reach here, file has no spectra - format stays NONE
      }
      catch (const FormatDetector::FirstSpectrumRead&)
      {
        im_format = detector.detected_format;
      }
    }

    // Step 2: Validate format
    if (im_format == IMFormat::CENTROIDED)
    {
      OPENMS_LOG_ERROR << "Error: Input file contains ion mobility data that is already centroided. "
                       << "PeakPickerIM expects raw (concatenated) IM data. "
                       << "Re-picking already centroided data is not supported." << std::endl;
      return ILLEGAL_PARAMETERS;
    }
    if (im_format == IMFormat::MULTIPLE_SPECTRA)
    {
      OPENMS_LOG_ERROR << "Error: Input file contains ion mobility data in MULTIPLE_SPECTRA format "
                       << "(one spectrum per IM frame). PeakPickerIM expects raw (concatenated) IM data "
                       << "where each spectrum contains an ion mobility float data array. "
                       << "This format is not supported." << std::endl;
      return ILLEGAL_PARAMETERS;
    }
    if (im_format == IMFormat::MIXED)
    {
      OPENMS_LOG_ERROR << "Error: Input file contains mixed ion mobility formats "
                       << "(both CONCATENATED and MULTIPLE_SPECTRA). PeakPickerIM expects raw (concatenated) IM data "
                       << "where each spectrum contains an ion mobility float data array. "
                       << "Mixed formats are not supported." << std::endl;
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

      // Check if input contains centroided IM data (error) or no IM data (warning)
      IMFormat im_format = IMTypes::determineIMFormat(exp);
      if (im_format == IMFormat::CENTROIDED)
      {
        OPENMS_LOG_ERROR << "Error: Input file contains ion mobility data that is already centroided. "
                         << "PeakPickerIM expects raw (concatenated) IM data. "
                         << "Re-picking already centroided data is not supported." << std::endl;
        return ILLEGAL_PARAMETERS;
      }
      if (im_format == IMFormat::NONE)
      {
        OPENMS_LOG_WARN << "Warning: Input file does not contain ion mobility data. "
                        << "No peak picking will be performed." << std::endl;
        mzml.store(output_file, exp);
        return EXECUTION_OK;
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
#pragma omp parallel for
        for (SignedSize i = 0; i < static_cast<SignedSize>(exp.size()); ++i)
        {
          MSSpectrum& spectrum = exp[static_cast<Size>(i)];

          // Only pick MS1 spectra; pass MS2 through as raw
          if (spectrum.getMSLevel() != 1) continue;

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

