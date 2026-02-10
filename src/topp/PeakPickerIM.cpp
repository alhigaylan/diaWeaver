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
#include <deque>

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

  // -------------------- Aggregating low-memory consumer --------------------
  // Buffers MS1 spectra in a sliding window to enable cross-scan aggregation
  // before peak picking. MS2 spectra are passed through immediately.
  class AggregatingConsumer : public MSDataWritingConsumer
  {
  public:
    AggregatingConsumer(String filename, const String& method, PeakPickerIM& pp,
                        double rt_fwhm, double cutoff) :
        MSDataWritingConsumer(std::move(filename)), pp_(pp), method_(method),
        rt_fwhm_(rt_fwhm), cutoff_(cutoff)
    {
      factor_ = -4.0 * std::log(2.0) / (rt_fwhm_ * rt_fwhm_);
      // Compute max RT difference where weight >= cutoff
      // weight = exp(factor * rt_diff^2) >= cutoff
      // factor * rt_diff^2 >= ln(cutoff)
      // rt_diff^2 <= ln(cutoff) / factor  (factor is negative)
      // rt_diff <= sqrt(ln(cutoff) / factor)
      max_rt_diff_ = std::sqrt(std::log(cutoff_) / factor_);
      OPENMS_LOG_INFO << "AggregatingConsumer: max RT difference for aggregation: "
                      << max_rt_diff_ << " seconds\n";
    }

    // Override consumeSpectrum to control buffering and writing directly
    void consumeSpectrum(MapType::SpectrumType& spectrum) override
    {
      if (spectrum.getMSLevel() != 1)
      {
        // Non-MS1 spectra: process and write immediately via parent
        applyPeakPicking(spectrum);
        MSDataWritingConsumer::consumeSpectrum(spectrum);
        return;
      }

      // Check for IM data on first MS1 spectrum
      if (!checked_im_data_)
      {
        checked_im_data_ = true;
        if (!spectrum.containsIMData())
        {
          OPENMS_LOG_WARN << "AggregatingConsumer: First MS1 spectrum has no IM data. "
                          << "Aggregation will be skipped.\n";
          has_im_data_ = false;
        }
        else
        {
          has_im_data_ = true;
        }
      }

      // If no IM data, just process and write immediately
      if (!has_im_data_)
      {
        applyPeakPicking(spectrum);
        MSDataWritingConsumer::consumeSpectrum(spectrum);
        return;
      }

      // Add MS1 spectrum to buffer (do NOT write yet)
      ms1_buffer_.push_back(spectrum);

      // Try to flush ready spectra
      flushReady();
    }

    void processSpectrum_(MapType::SpectrumType&) override
    {
      // Not used - we handle everything in consumeSpectrum
    }

    void processChromatogram_(MapType::ChromatogramType&) override {}

    // Override the consumer to flush at the end
    ~AggregatingConsumer() override
    {
      // Flush all remaining buffered spectra
      flushAll();
      OPENMS_LOG_INFO << "AggregatingConsumer: Finished. Wrote " << spectra_written_
                      << " aggregated MS1 spectra.\n";
    }

  private:
    void applyPeakPicking(MapType::SpectrumType& spectrum)
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

    void flushReady()
    {
      // We can output spectrum at index i when we've seen all spectra that could
      // contribute to its aggregation (i.e., spectra within max_rt_diff_)
      while (ms1_buffer_.size() >= 2)
      {
        double first_rt = ms1_buffer_.front().getRT();
        double last_rt = ms1_buffer_.back().getRT();

        // If the last spectrum is far enough past the first,
        // we have all neighbors for the first spectrum
        if (last_rt - first_rt > max_rt_diff_)
        {
          outputAggregatedSpectrum(0);
          ms1_buffer_.pop_front();
        }
        else
        {
          break; // Need more spectra before we can output
        }
      }
    }

    void flushAll()
    {
      // Output all remaining buffered spectra
      while (!ms1_buffer_.empty())
      {
        outputAggregatedSpectrum(0);
        ms1_buffer_.pop_front();
      }
    }

    void outputAggregatedSpectrum(Size center_idx)
    {
      if (center_idx >= ms1_buffer_.size()) return;

      MSSpectrum& center_spectrum = ms1_buffer_[center_idx];
      double center_rt = center_spectrum.getRT();

      // Collect neighbors and their weights
      std::vector<MSSpectrum> spectra_to_aggregate;
      std::vector<double> weights;

      for (Size i = 0; i < ms1_buffer_.size(); ++i)
      {
        double rt_diff = ms1_buffer_[i].getRT() - center_rt;
        double weight = std::exp(factor_ * rt_diff * rt_diff);

        if (weight >= cutoff_ || i == center_idx)
        {
          spectra_to_aggregate.push_back(ms1_buffer_[i]);
          weights.push_back(weight);
        }
      }

      // Normalize weights
      double sum_weights = 0.0;
      for (double w : weights) sum_weights += w;
      for (double& w : weights) w /= sum_weights;

      OPENMS_LOG_DEBUG << "AggregatingConsumer: Outputting spectrum at RT=" << center_rt
                       << " with " << spectra_to_aggregate.size() << " neighbors aggregated.\n";

      // Aggregate
      MSSpectrum aggregated;
      pp_.aggregateScans(spectra_to_aggregate, weights, aggregated);

      // Apply peak picking to the aggregated spectrum
      applyPeakPicking(aggregated);

      // Write the aggregated spectrum directly (bypass processSpectrum_)
      doWriteSpectrum_(aggregated);
      ++spectra_written_;
    }

    // Direct write without going through processSpectrum_
    void doWriteSpectrum_(MapType::SpectrumType& spectrum)
    {
      // Call the parent's consumeSpectrum but since processSpectrum_ is empty,
      // it will just write the spectrum
      MSDataWritingConsumer::consumeSpectrum(spectrum);
    }

    PeakPickerIM& pp_;
    String method_;
    double rt_fwhm_;
    double cutoff_;
    double factor_;
    double max_rt_diff_;
    std::deque<MSSpectrum> ms1_buffer_;
    bool checked_im_data_ = false;
    bool has_im_data_ = true;
    Size spectra_written_ = 0;
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
                              bool aggregate_scans, double rt_fwhm, double cutoff)
  {
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
      // Pass through unchanged
      PassthroughConsumer passthrough(output_file);
      mzml.transform(input_file, &passthrough);
      return EXECUTION_OK;
    }

    // Step 3: Proceed with streaming processing
    if (aggregate_scans)
    {
      OPENMS_LOG_INFO << "Low-memory mode with aggregation: using sliding window buffer." << std::endl;
      AggregatingConsumer agg_consumer(output_file, method, pp, rt_fwhm, cutoff);
      agg_consumer.addDataProcessing(getProcessingInfo_(DataProcessing::PEAK_PICKING));
      mzml.transform(input_file, &agg_consumer);
    }
    else
    {
      Consumer pp_consumer(output_file, method, pp);
      pp_consumer.addDataProcessing(getProcessingInfo_(DataProcessing::PEAK_PICKING));
      mzml.transform(input_file, &pp_consumer);
    }
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
      // Get aggregation parameters for low-memory mode
      double rt_fwhm = algo.getValue("aggregation:rt_FWHM");
      double cutoff = algo.getValue("aggregation:cutoff");
      return doLowMemAlgorithm(method, picker, input_file, output_file,
                               aggregate_scans, rt_fwhm, cutoff);
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

      // Aggregate signal across scans if requested (improves S/N before peak picking)
      if (aggregate_scans)
      {
        OPENMS_LOG_INFO << "Aggregating signal across neighboring scans..." << std::endl;
        picker.pickExperimentWithAggregation(exp);
      }

#pragma omp parallel for
      for (SignedSize i = 0; i < static_cast<SignedSize>(exp.size()); ++i)
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

