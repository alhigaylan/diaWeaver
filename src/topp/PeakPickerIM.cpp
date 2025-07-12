// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Author: Mohammed Alhigaylan $
// $Maintainer: Timo Sachsenberg $
// -------------------------------------------------------------------------------------------------------------------------------------------

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerIM.h>

using namespace OpenMS;
using namespace std;

class TOPPPeakPickerIM : public TOPPBase
{
public:
  TOPPPeakPickerIM() : TOPPBase("PeakPickerIM", "Applies PeakPickerIM to an mzML file", false) {}

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input mzML file");
    setValidFormats_("in", { "mzML" });

    registerOutputFile_("out", "<file>", "", "Output mzML file");
    setValidFormats_("out", { "mzML" });

    registerStringOption_("processOption", "<name>", "inmemory", "Whether to load all data and process them in-memory or whether to process the data on the fly (lowmemory) without loading the whole file into memory first", false, true);
    setValidStrings_("processOption", { "inmemory", "lowmemory" } );

    registerStringOption_("method", "<name>", "", "Method to pick peaks in IM dimension", false, true);
    setValidStrings_("method", { "mobilogram", "cluster", "traces" } );

    registerFlag_("merge_neighbors", "If set, merge peaks from neighboring spectra (n-1 and n+1) with half intensity");
    registerFlag_("include_unclaimed", "If set, include unclaimed raw peaks into the centroided output");

  }

    /**
    @brief Helper class for the Low Memory peak-picking
  */
  class Consumer : public MSDataWritingConsumer
  {
    public:

    Consumer(String filename, const String& method, const PeakPickerIM& pp, bool add_unclaimed) :
      MSDataWritingConsumer(std::move(filename)), pp_(pp), method_(method), add_unclaimed_(add_unclaimed) {}

    void processSpectrum_(MapType::SpectrumType& spectrum) override
    {
      if (method_ == "mobilogram")
      {
        pp_.pickIMTraces(spectrum, add_unclaimed_);
      }
      else if (method_ == "cluster")
      {
        PeakPickerIM::pickIMCluster(spectrum, 50.0, 0.08);
      }
      else if (method_ == "traces")
      {    
        PeakPickerIM::pickIMElutionProfiles(spectrum, 50.0);
      }
    }

    void processChromatogram_(MapType::ChromatogramType & ) override
    {
    }

  private:
    PeakPickerIM pp_;
    String method_ = "mobilogram";
    bool add_unclaimed_;
  };


  ExitCodes doLowMemAlgorithm(const String method, const PeakPickerIM& pp, const String& input_file, const String& output_file, bool add_unclaimed)
  {
    ///////////////////////////////////
    // Create the consumer object, add data processing
    ///////////////////////////////////
    Consumer pp_consumer(output_file, method, pp, add_unclaimed);
    pp_consumer.addDataProcessing(getProcessingInfo_(DataProcessing::PEAK_PICKING));

    ///////////////////////////////////
    // Create new MSDataReader and set our consumer
    ///////////////////////////////////
    MzMLFile mz_data_file;
    mz_data_file.setLogType(log_type_);
    mz_data_file.transform(input_file, &pp_consumer);

    return EXECUTION_OK;
  }

  ExitCodes main_(int, const char**) override
  {
    // Get input and output file paths
    String input_file = getStringOption_("in");
    String output_file = getStringOption_("out");
    String process_option = getStringOption_("processOption");
    String method = getStringOption_("method");

    // Retrieve user parameter
    bool merge_neighbors = getFlag_("merge_neighbors");
    bool add_unclaimed = getFlag_("include_unclaimed");

    PeakPickerIM picker;
    if (process_option == "lowmemory")
    {
      return doLowMemAlgorithm(method, picker, input_file, output_file, add_unclaimed); // TODO: needs parallelization
    }
    else
    {
      // Load input mzML file
      PeakMap exp;
      MzMLFile mzml;
      mzml.load(input_file, exp);

      // --------- To boost signal-to-noise, we will borrow peaks
      // from n-1 and n+1 spectra and half their intensities.
      // Cache modified spectra
      PeakMap modified_exp = exp;
      #pragma omp parallel for
      for (Int64 i = 0; i < static_cast<Int64>(exp.size()); ++i)
      {
        MSSpectrum& current = modified_exp[i];
        const MSSpectrum& original = exp[i];

        // Start with a deep copy of the original spectrum (incl. peaks + meta)
        current = original;

        auto append_neighbor_peaks = [&](const MSSpectrum& neighbor, double intensity_scale)
        {
          const auto& neighbor_peaks = neighbor;
          const auto& neighbor_fda = neighbor.getFloatDataArrays();
          auto& current_fda = current.getFloatDataArrays();

          Size num_arrays = current_fda.size();

          for (Size j = 0; j < neighbor_peaks.size(); ++j)
          {
            Peak1D p = neighbor_peaks[j];
            p.setIntensity(p.getIntensity() * intensity_scale);
            current.push_back(p);

            // Copy over float metadata as-is (no scaling)
            for (Size k = 0; k < num_arrays; ++k)
            {
              if (k < neighbor_fda.size() && j < neighbor_fda[k].size())
              {
                current_fda[k].push_back(neighbor_fda[k][j]);
              }
              else
              {
                throw std::runtime_error("ERROR! FloatDataArray size mismatch when appending neighbor peaks.");
              }
            }
          }
        };
        if (merge_neighbors)
        {
          if (i >= 1) append_neighbor_peaks(exp[i - 1], 0.5);
          if (i + 1 < static_cast<Int64>(exp.size())) append_neighbor_peaks(exp[i + 1], 0.5);
          // merge across 5 frames
          //if (i >= 2) append_neighbor_peaks(exp[i - 2], 0.25);
          //if (i + 2 < static_cast<Int64>(exp.size())) append_neighbor_peaks(exp[i + 2], 0.25);
        }
        current.sortByPosition(); // optional, to preserve m/z order
      }

      // ---- sanity check. Print peaks from raw spectrum and modified spectrum ----
      /*
      const MSSpectrum& original_spec = exp[50];
      const MSSpectrum& modified_spec = modified_exp[50];

      std::cout << "\n===== ORIGINAL SPECTRUM @ index 50 =====" << std::endl;
      const auto& orig_fda = original_spec.getFloatDataArrays();
      for (Size i = 0; i < original_spec.size(); ++i)
      {
        std::cout << "mz: " << original_spec[i].getMZ()
                  << ", intensity: " << original_spec[i].getIntensity();

        std::cout << ", float[0]: " << orig_fda[0][i] << std::endl;
      }

      std::cout << "\n===== MODIFIED SPECTRUM @ index 50 =====" << std::endl;
      const auto& mod_fda = modified_spec.getFloatDataArrays();
      for (Size i = 0; i < modified_spec.size(); ++i)
      {
        std::cout << "mz: " << modified_spec[i].getMZ()
                  << ", intensity: " << modified_spec[i].getIntensity();

        std::cout << ", float[0]: " << mod_fda[0][i] << std::endl;
      }
       */


      // Process each spectrum with PeakPickerIM
      #pragma omp parallel for
      for (Int64 i = 0; i != modified_exp.size(); i++)
      {
        MSSpectrum& spectrum = modified_exp[i];
        OPENMS_LOG_DEBUG << "Processing MS" << spectrum.getMSLevel() << " spectrum with " 
          << spectrum.size() << " peaks in the IM frame." << std::endl;
        if (method == "mobilogram")
        {
          picker.pickIMTraces(spectrum, add_unclaimed);
        }
        else if (method == "cluster")
        {
          PeakPickerIM::pickIMCluster(spectrum, 100.0, 0.08);
        }
        else if (method == "traces")
        {    
          PeakPickerIM::pickIMElutionProfiles(spectrum, 100.0);
        }
        OPENMS_LOG_DEBUG << "Processed spectrum has " << spectrum.size() << " centroided IM peaks." << std::endl;
      }
      // Save output mzML file
      OPENMS_LOG_DEBUG << "Saving output mzML file: " << output_file << std::endl;
      mzml.store(output_file, modified_exp);

      return EXECUTION_OK;
    }
  }
};

int main(int argc, const char** argv)
{
  TOPPPeakPickerIM tool;
  return tool.main(argc, argv);
}

