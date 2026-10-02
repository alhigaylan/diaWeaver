// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Author: Timo Sachsenberg, Mohammed Alhigaylan $
// $Maintainer: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerIM.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/PROCESSING/SMOOTHING/GaussFilter.h>
#include <OpenMS/PROCESSING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/PROCESSING/RESAMPLING/LinearResamplerAlign.h>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>
#include <OpenMS/MATH/MISC/SplineBisection.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>
#include <OpenMS/KERNEL/SpectrumHelper.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include <OpenMS/FEATUREFINDER/ElutionPeakDetection.h>
#include <iostream>
#include <deque>
#include <cmath>
#include <algorithm>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <map>
#include <utility>


using namespace std;

namespace OpenMS
{

    // Function to compute the lower and upper m/z bounds based on ppm tolerance
    std::pair<double, double> PeakPickerIM::ppmBounds(double mz, double ppm)
    {
      ppm = ppm / 1e6;
      double delta_mz = (ppm * mz) / 2.0;

      double low = mz - delta_mz;
      double high = mz + delta_mz;

      return std::make_pair(low, high);
    }

    // PRECONDITION: input_spectrum is sorted by m/z
    // This function sums peaks if they are nearly identical
    // OpenMS represents TimsTOF data in MSSpectrum() objects as one-array.
    // Example: There could be multiple 500.0 m/z peaks with different ion mobility values.
    // Example2: extracted mobilogram could have multiple 0.88 1/k values from different m/z peaks.
    // Peak picking (such as HiRes) will not work properly if there are multiple y measurements at a given x position.
    // Note: does not clear the output_spectrum but add peaks to it (required for fast padding)
    void PeakPickerIM::sumFrame_(const MSSpectrum& input_spectrum,
                                 MSSpectrum& output_spectrum,
                                 double tolerance,
                                 bool use_ppm)
    {
      OPENMS_PRECONDITION(input_spectrum.isSorted(), "Spectrum must be sorted by m/z before summing peaks.");
      
      if (input_spectrum.empty()) return;

      double current_mz = input_spectrum[0].getMZ();
      double current_intensity = input_spectrum[0].getIntensity();

      for (Size i = 1; i < input_spectrum.size(); ++i)
      {
        double next_mz = input_spectrum[i].getMZ();
        double next_intensity = input_spectrum[i].getIntensity();

        double delta_mz = std::abs(next_mz - current_mz);
        bool within_tolerance = use_ppm
                                  ? ((delta_mz / current_mz) * 1e6 <= tolerance)
                                  : (delta_mz <= tolerance);

        if (within_tolerance)
        {
          current_intensity += next_intensity;
        }
        else // new peak is outside of tolerance window
        {
          output_spectrum.emplace_back(current_mz, current_intensity);
          current_mz = next_mz;
          current_intensity = next_intensity;
        }
      }
      output_spectrum.emplace_back(current_mz, current_intensity);
    }

    // We use peak FWHM (from PeakPickerHiRes) to extract ion mobility traces.
    // Given a picked m/z peak, we extract a Mobilogram from the raw peaks within its FWHM.
    // To enable recomputing of m/z center after ion mobility peak picking, we tack raw m/z peak values
    // in FloatDataArrays().

    std::vector<Mobilogram> PeakPickerIM::extractIonMobilityTraces(
      const MSSpectrum& picked_spectrum,
      const MSSpectrum& raw_spectrum)
    {
      const auto& float_data_arrays = picked_spectrum.getFloatDataArrays();

      // Find FWHM array in picked_spectrum
      const MSSpectrum::FloatDataArray* fwhm_array = nullptr;

      for (const auto& array : float_data_arrays)
      {
        if (array.getName() == "FWHM_ppm")
        {
          fwhm_array = &array;
          break;
        }
      }

      if (!fwhm_array)
      {
        OPENMS_LOG_WARN << "FWHM data array not found!\n";
        return {};
      }

      if (fwhm_array->size() != picked_spectrum.size())
      {
        OPENMS_LOG_WARN << "Size mismatch between FWHM array and picked peaks!\n";
        return {};
      }
      // Get the Ion Mobility array index from raw_spectrum
      if (!raw_spectrum.containsIMData())
      {
        OPENMS_LOG_WARN << "No ion mobility data found in raw_spectrum.\n";
        return {};
      }
      const auto [im_data_index, im_unit] = raw_spectrum.getIMData();
      const auto& ion_mobility_array = raw_spectrum.getFloatDataArrays()[im_data_index];

      // One Mobilogram per picked m/z peak
      std::vector<Mobilogram> mobility_traces;

      for (size_t i = 0; i < picked_spectrum.size(); ++i)
      {
        double picked_mz = picked_spectrum[i].getMZ();
        double fwhm_ppm = (*fwhm_array)[i];

        auto bounds = ppmBounds(picked_mz, fwhm_ppm);
        double lower_bound = bounds.first;
        double upper_bound = bounds.second;

        SignedSize center_idx = raw_spectrum.findNearest(picked_mz);

        if (center_idx == -1)
        {
          OPENMS_LOG_WARN << "No raw peaks found near picked m/z: " << picked_mz << '\n';
          mobility_traces.emplace_back();
          continue;
        }

        // Indices of the raw peaks within the FWHM window
        std::vector<Size> trace_indices;

        // Expand left
        SignedSize left_idx = center_idx;
        while (left_idx >= 0 && raw_spectrum[left_idx].getMZ() >= lower_bound)
        {
          trace_indices.push_back(left_idx);
          --left_idx;
        }

        // Expand right
        SignedSize right_idx = center_idx + 1;
        while (right_idx < static_cast<SignedSize>(raw_spectrum.size()) &&
               raw_spectrum[right_idx].getMZ() <= upper_bound)
        {
          trace_indices.push_back(right_idx);
          ++right_idx;
        }

        // Mobilogram::sortByPosition does not reorder data arrays, so sort the indices by ion mobility instead
        std::stable_sort(trace_indices.begin(), trace_indices.end(),
                         [&ion_mobility_array](Size a, Size b) { return ion_mobility_array[a] < ion_mobility_array[b]; });

        Mobilogram trace;
        Mobilogram::FloatDataArray raw_mz_array;
        raw_mz_array.setName("raw_mz");
        for (const Size idx : trace_indices)
        {
          trace.emplace_back(ion_mobility_array[idx], raw_spectrum[idx].getIntensity());
          raw_mz_array.push_back(raw_spectrum[idx].getMZ());
        }
        trace.getFloatDataArrays().push_back(std::move(raw_mz_array));

        mobility_traces.push_back(std::move(trace));
      }

      return mobility_traces;
    }

    // Function to compute m/z centers from mobilogram_traces and picked_traces
    MSSpectrum PeakPickerIM::computeCentroids_(const vector<Mobilogram>& mobilogram_traces,
                              const vector<Mobilogram>& picked_traces)
    {
      MSSpectrum centroided_frame;

      // Create float data arrays to house ion mobility data and peaks FWHM
      MSSpectrum::FloatDataArray ion_mobility_array;
      ion_mobility_array.setName(Constants::UserParam::ION_MOBILITY_CENTROID);

      MSSpectrum::FloatDataArray ion_mobility_fwhm;
      ion_mobility_fwhm.setName("IM FWHM");

      MSSpectrum::FloatDataArray mz_fwhm_array;
      mz_fwhm_array.setName("MZ FWHM");

      // Loop over picked traces and their corresponding raw mobilogram traces
      for (size_t i = 0; i < picked_traces.size(); ++i)
      {
        // std::cout << "Looping through picked_trace that has .. " << picked_traces[i].size() << '\n';
        const Mobilogram& picked_trace = picked_traces[i];
        const Mobilogram& raw_trace = mobilogram_traces[i];

        const auto& picked_float_arrays = picked_trace.getFloatDataArrays();

        if (picked_float_arrays.empty())
        {
          OPENMS_LOG_WARN << "No IM FWHM array found for picked_trace " << i << "!\n";
          continue;
        }

        // Assuming the first FloatDataArray holds the ion mobility peak FWHM values
        const auto& fwhm_array = picked_float_arrays[0];

        if (fwhm_array.size() != picked_trace.size())
        {
          OPENMS_LOG_WARN << "FWHM array size mismatch with picked_trace size!\n";
          continue;
        }

        // Get the FloatDataArrays from raw_trace (assumed to hold the raw m/z values)
        const auto& raw_float_arrays = raw_trace.getFloatDataArrays();

        if (raw_float_arrays.empty())
        {
          OPENMS_LOG_WARN << "No raw m/z peaks found for raw_trace " << i << "!\n";
          continue;
        }

        // Assume the first array holds the raw m/z values
        const auto& raw_mz_values = raw_float_arrays[0];

        if (raw_mz_values.size() != raw_trace.size())
        {
          OPENMS_LOG_WARN << "raw_mz_values size mismatch with raw_trace size!\n";
          continue;
        }


        // Create reusable objects outside the loop to reduce memory allocations
        MSSpectrum raw_peaks_within_bounds;
        MSSpectrum raw_mz_peaks;
        vector<double> mz_values;
        vector<double> intensity_values;
        vector<size_t> indices;
        vector<double> sorted_mz;
        vector<double> sorted_intensity;
        
        // Iterate through picked peaks in this trace
        for (Size j = 0; j < picked_trace.size(); ++j)
        {
          double centroid_im = picked_trace[j].getMobility();
          double fwhm = fwhm_array[j];

          double im_lower = centroid_im - (fwhm / 2.0);
          double im_upper = centroid_im + (fwhm / 2.0);

          // Use findNearest() to get the index of the closest peak in the raw mobilogram trace
          SignedSize center_idx = raw_trace.findNearest(centroid_im);

          if (center_idx == -1)
          {
            OPENMS_LOG_WARN << "Could not find nearest peak to centroid_im in raw_trace!\n";
            continue;
          }

          // Clear the spectrum for reuse
          raw_peaks_within_bounds.clear(false);

          // --- Expand Left ---
          SignedSize left_idx = center_idx;
          while (left_idx >= 0 && raw_trace[left_idx].getMobility() >= im_lower)
          {
            Peak1D new_peak;
            new_peak.setMZ(raw_mz_values[left_idx]);                      // m/z from FloatDataArray
            new_peak.setIntensity(raw_trace[left_idx].getIntensity());    // intensity from raw_trace
            raw_peaks_within_bounds.push_back(new_peak);

            --left_idx;
          }

          // --- Expand Right ---
          SignedSize right_idx = center_idx + 1;
          while (right_idx < static_cast<SignedSize>(raw_trace.size()) &&
                 raw_trace[right_idx].getMobility() <= im_upper)
          {
            Peak1D new_peak;
            new_peak.setMZ(raw_mz_values[right_idx]);
            new_peak.setIntensity(raw_trace[right_idx].getIntensity());
            raw_peaks_within_bounds.push_back(new_peak);

            ++right_idx;
          }


          // If we only retrieved one raw peak, pass it over to centroided_frame as is
          // Resampling and smoothing the raw data distorts the intensity values.
          // We recompute the m/z peak maxima and intensity using spline
          if (raw_peaks_within_bounds.size() == 1)
          {
            const Peak1D& single_peak = raw_peaks_within_bounds.front();

            // Add it directly to centroided_frame
            centroided_frame.push_back(single_peak);

            // Push corresponding ion mobility and FWHM arrays
            ion_mobility_array.push_back(centroid_im);
            ion_mobility_fwhm.push_back(fwhm);
            mz_fwhm_array.push_back(0.0);

            // Skip the rest of the loop and move on to the next picked_trace peak
            continue;
          }
        
          // Sort by m/z for sumFrame_() which requires sorted input.
          // Data is typically unsorted because peaks were collected by walking through IM-sorted
          // data (expand-left/right), and m/z has no correlation with IM order within FWHM window.
          raw_peaks_within_bounds.sortByPosition();

          // Clear the spectrum for reuse
          raw_mz_peaks.clear(false);
          sumFrame_(raw_peaks_within_bounds, raw_mz_peaks, sum_tolerance_mz_, true);
          if (raw_mz_peaks.empty())
          {
            OPENMS_LOG_DEBUG << "No data in raw_mz_peaks for picked IM peak " << j << "!\n";
            continue;
          }

          // Summing peaks could result in spectrum.size() == 1 which causes error.
          // in this case, simply sum the intensity values
          if (raw_mz_peaks.size() == 1)
          {
            centroided_frame.push_back(raw_mz_peaks[0]);
            ion_mobility_array.push_back(centroid_im);
            ion_mobility_fwhm.push_back(fwhm);
            mz_fwhm_array.push_back(0.0);

            continue;
          }

          // Clear and reuse vectors for spline data
          mz_values.clear();
          intensity_values.clear();
          
          // Reserve space for efficiency
          mz_values.reserve(raw_mz_peaks.size());
          intensity_values.reserve(raw_mz_peaks.size());

          // Initialize sorting flag
          bool is_sorted = true;

          // Populate the vectors and check if sorted at the same time
          for (size_t i = 0; i < raw_mz_peaks.size(); ++i)
          {
            double current_mz = raw_mz_peaks[i].getMZ();
            double current_intensity = raw_mz_peaks[i].getIntensity();
            
            // Check if still sorted (compare with previous value if not the first element)
            if (i > 0 && current_mz < mz_values.back())
            {
              is_sorted = false;
            }
            
            mz_values.push_back(current_mz);
            intensity_values.push_back(current_intensity);
          }

          // Sort vectors if needed (safety check - data should already be sorted from sumFrame_)
          // CubicSpline2d requires sorted x-coordinates
          if (!is_sorted)
          {
            // Reuse indices vector
            indices.resize(mz_values.size());
            for (size_t i = 0; i < indices.size(); ++i)
            {
              indices[i] = i;
            }
            
            // Sort indices based on m/z values
            std::sort(indices.begin(), indices.end(),
                    [&mz_values](size_t i1, size_t i2) {
                      return mz_values[i1] < mz_values[i2];
                    });
            
            // Reuse sorted vectors
            sorted_mz.resize(mz_values.size());
            sorted_intensity.resize(intensity_values.size());
            
            // Reorder both vectors using the sorted indices
            for (size_t i = 0; i < indices.size(); ++i)
            {
              sorted_mz[i] = mz_values[indices[i]];
              sorted_intensity[i] = intensity_values[indices[i]];
            }
            
            // Replace the original vectors with the sorted ones
            mz_values = std::move(sorted_mz);
            intensity_values = std::move(sorted_intensity);
          }

          // Initialize spline with the two vectors
          CubicSpline2d spline(mz_values, intensity_values);

          // Define boundaries
          const double left_bound = mz_values.front();
          const double right_bound = mz_values.back();

          // Find maximum via spline bisection
          double apex_mz = (left_bound + right_bound) / 2.0;
          double apex_intensity = 0.0;

          const double max_search_threshold = 1e-6;

          Math::spline_bisection(spline, left_bound, right_bound, apex_mz, apex_intensity, max_search_threshold);


          // FWHM calculation (same binary search as before)
          double half_height = apex_intensity / 2.0;
          const double fwhm_search_threshold = 0.01 * half_height;

          // ---- Left side search ----
          double mz_left = left_bound;
          double mz_center = apex_mz;
          double int_mid = 0.0;
          double mz_mid = mz_left;

          if (spline.eval(mz_left) > half_height)
          {
            mz_mid = mz_left;
          }
          else
          {
            do
            {
              mz_mid = (mz_left + mz_center) / 2.0;
              int_mid = spline.eval(mz_mid);

              if (int_mid < half_height)
              {
                mz_left = mz_mid;
              }
              else
              {
                mz_center = mz_mid;
              }
            } while (std::fabs(int_mid - half_height) > fwhm_search_threshold);
          }
          double fwhm_left_mz = mz_mid;

          // ---- Right side search ----
          double mz_right = right_bound;
          mz_center = apex_mz;

          if (spline.eval(mz_right) > half_height)
          {
            mz_mid = mz_right;
          }
          else
          {
            do
            {
              mz_mid = (mz_right + mz_center) / 2.0;
              int_mid = spline.eval(mz_mid);

              if (int_mid < half_height)
              {
                mz_right = mz_mid;
              }
              else
              {
                mz_center = mz_mid;
              }

            } while (std::fabs(int_mid - half_height) > fwhm_search_threshold);
          }
          double fwhm_right_mz = mz_mid;

          // ---- FWHM result ----
          double mz_fwhm = fwhm_right_mz - fwhm_left_mz;


          centroided_frame.emplace_back(apex_mz, apex_intensity);
          ion_mobility_array.push_back(centroid_im);
          ion_mobility_fwhm.push_back(fwhm);
          mz_fwhm_array.push_back(mz_fwhm);
        }

      }

      auto& centroided_frame_fda = centroided_frame.getFloatDataArrays();
      centroided_frame_fda.push_back(std::move(ion_mobility_array));
      centroided_frame_fda.push_back(std::move(ion_mobility_fwhm));
      centroided_frame_fda.push_back(std::move(mz_fwhm_array));
      centroided_frame.sortByPosition();

      return centroided_frame;
    }

    void removeAllFloatDataArraysExcept(OpenMS::MSSpectrum& spectrum, const String& keep_name)
    {
      auto& float_arrays = spectrum.getFloatDataArrays();
  
      // Use remove_if to move all elements to remove to the end
      auto new_end = std::remove_if(float_arrays.begin(), float_arrays.end(),
                               [&keep_name](const MSSpectrum::FloatDataArray& array) {
                                 return array.getName() != keep_name;  // Remove if NOT the one we want to keep
                               });
  
      // Erase the removed elements
      float_arrays.erase(new_end, float_arrays.end());
    }

    PeakPickerIM::PeakPickerIM()
        : DefaultParamHandler("PeakPickerIM")
    {
      // --- PickIMTraces parameters ---
      defaults_.setValue("pickIMTraces:sum_tolerance_mz",          1.0,  "Tolerance for summing adjacent m/z peaks (ppm)");
      defaults_.setValue("pickIMTraces:mobilogram_sampling_grid", 0.01, "Grid spacing for linear resampling of mobilograms (in 1/K0 units).");
      defaults_.setValue("pickIMTraces:gauss_ppm_tolerance",      5.0,  "Gaussian smoothing m/z tolerance in ppm");
      defaults_.setValue("pickIMTraces:sgolay_frame_length",     5,     "Savitzky-Golay smoothing frame length");
      defaults_.setValue("pickIMTraces:sgolay_polynomial_order", 3,     "Savitzky-Golay smoothing polynomial order");
#if 0
      // --- PickIMCluster parameters ---
      defaults_.setValue("pickIMCluster:ppm_tolerance_cluster", 50.0, "m/z tolerance in ppm for clustering");
      defaults_.setValue("pickIMCluster:im_tolerance_cluster", 0.1, "Ion mobility tolerance for clustering (in 1/K0 units). For CCS data, use larger values (e.g., 10-20).");
      // --- PickIMElutionProfiles parameters ---
      defaults_.setValue("pickIMElutionProfiles:ppm_tolerance_elution", 50.0, "Mass trace m/z tolerance in ppm");
#endif
      // --- Aggregation parameters ---
      defaults_.setValue("aggregation:rt_FWHM", 5.0, "Full width at half maximum for Gaussian weighting of scans (in seconds)");
      defaults_.setValue("aggregation:cutoff", 0.01, "Weight cutoff below which spectra are not included in aggregation");

      defaultsToParam_();  // copies defaults_ into param_
      updateMembers_();    // caches into member variables
    }
    void PeakPickerIM::updateMembers_()
    {
      sum_tolerance_mz_         = (double)param_.getValue("pickIMTraces:sum_tolerance_mz");
      mobilogram_sampling_grid_ = (double)param_.getValue("pickIMTraces:mobilogram_sampling_grid");
      gauss_ppm_tolerance_      = (double)param_.getValue("pickIMTraces:gauss_ppm_tolerance");
      sgolay_frame_length_   = (int)param_.getValue("pickIMTraces:sgolay_frame_length");
      sgolay_polynomial_order_= (int)param_.getValue("pickIMTraces:sgolay_polynomial_order");

      Param gauss_params;
      gauss_params.setValue("ppm_tolerance", gauss_ppm_tolerance_);
      gauss_params.setValue("use_ppm_tolerance", "true");
      gauss_filter_.setParameters(gauss_params);

      Param picker_mz_p;
      picker_mz_p.setValue("signal_to_noise", 0.0);
      picker_mz_p.setValue("report_FWHM", "true");
      picker_mz_p.setValue("report_FWHM_unit", "relative");
      picker_mz_p.setValue("allow_missing_flank", "true");
      picker_mz_.setParameters(picker_mz_p);

      Param resampler_param;
      resampler_param.setValue("spacing", mobilogram_sampling_grid_);
      resampler_param.setValue("ppm", "false");
      lin_resampler_.setParameters(resampler_param);

      Param sgolay_params;
      sgolay_params.setValue("frame_length", sgolay_frame_length_);
      sgolay_params.setValue("polynomial_order", sgolay_polynomial_order_);
      sgolay_filter_.setParameters(sgolay_params);

      Param picker_im_p;
      picker_im_p.setValue("signal_to_noise", 0.0);
      picker_im_p.setValue("spacing_difference_gap", 0.0);
      picker_im_p.setValue("spacing_difference", 0.0);
      picker_im_p.setValue("missing", 0);
      picker_im_p.setValue("report_FWHM", "true");
      picker_im_p.setValue("report_FWHM_unit", "absolute");
      picker_im_.setParameters(picker_im_p);

#if 0
      ppm_tolerance_cluster_ = (double)param_.getValue("pickIMCluster:ppm_tolerance_cluster");
      im_tolerance_cluster_ = (double)param_.getValue("pickIMCluster:im_tolerance_cluster");

      ppm_tolerance_elution_ = (double)param_.getValue("pickIMElutionProfiles:ppm_tolerance_elution");
#endif

      aggregation_rt_fwhm_ = (double)param_.getValue("aggregation:rt_FWHM");
      aggregation_cutoff_ = (double)param_.getValue("aggregation:cutoff");
    }

    namespace
    {
      /**
       * @brief Helper function to validate that a spectrum contains IM data in the correct format for peak picking
       *
       * @param[in] spectrum The spectrum to validate
       * @return true if the spectrum should be processed (has concatenated IM data)
       * @return false if the spectrum should be skipped (no IM data)
       * @throws Exception::InvalidValue if the data is already centroided, UNKNOWN, or unhandled format
       */
      bool validateIMFormatForPicking(const MSSpectrum& spectrum)
      {
        IMFormat format = IMTypes::determineIMFormat(spectrum);
        switch (format)
        {
          case IMFormat::NONE:
            return false;
          case IMFormat::IM_PEAK:
          {
            // Check IM peak type -- reject already-centroided data
            IMPeakType peak_type = spectrum.getIMPeakType();
            if (peak_type == IMPeakType::IM_CENTROIDED)
            {
              throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                "Ion mobility data is already centroided. PeakPickerIM expects raw (profile) IM data. "
                "Re-picking already centroided data is not supported.",
                imPeakTypeToString(peak_type));
            }
            OPENMS_LOG_DEBUG << "Processing IM_PEAK data (profile).\n";
            return true;
          }
          case IMFormat::UNKNOWN:
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "IMFormat is UNKNOWN. Call IMTypes::determineIMFormat() first.", "");
          default:
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Unsupported IMFormat for picking.", imFormatToString(format));
        }
      }
    }

    void PeakPickerIM::pickIMTraces(MSSpectrum& spectrum)
    {
      // Only process MS1 spectra; non-MS1 spectra are passed through unchanged
      if (spectrum.getMSLevel() != 1)
      {
        return;
      }

      // Validate IM format - returns false if we should skip processing
      if (!validateIMFormatForPicking(spectrum))
      {
        return;
      }
      
      
      // Spectrum is in IM_PEAK IM format. Now sort by m/z to prepare for m/z peak picking
      spectrum.sortByPosition();

      // ************************************************* PART I *****************************************************
      // ------------------------------------------ mass-to-charge peak picking -------------------------------------
      // ------------------------------------------ Step a1: Sum m/z peaks ------------------------------------------
      // project all timsTOF peaks into the m/z axis using sumFrame_
      // The ppm tolerance is a dynamic way of testing m/z floats being almost identical. The raw intensity is summed.
      MSSpectrum summed_spectrum;
      sumFrame_(spectrum, summed_spectrum, sum_tolerance_mz_, true);

      // ------------------------------------------ step 2a: smooth ------------------------------------------
      // Apply gaussian smoothing to the peaks projected into the m/z axis. This facilitates peak picking
      // in the m/z dimension and subseqent mobilogram extraction for each picked m/z peak.
      gauss_filter_.filter(summed_spectrum);

      // ------------------------------------------ step 3a: m/z Peak Picking ------------------------------------------
      // Pick peaks in the m/z axis and toggle reporting peak width at half max (FWHM)
      // we will use the FWHM of each picked m/z peak to extract mobilograms.
      MSSpectrum picked_spectrum;
      picker_mz_.pick(summed_spectrum, picked_spectrum);
      if (picked_spectrum.empty())
      {
        OPENMS_LOG_WARN << "No m/z peaks picked. Returning empty spectrum.\n";
        // Preserve metadata (RT, MS level, precursor info) but clear peak data.
        // Mirrors the successful picking path (lines below spectrum = centroided_frame)
        // and PeakPickerHiRes::pick(), both of which retain spectrum metadata even
        // when no peaks are found.
        MSSpectrum empty_frame;
        copySpectrumMeta(spectrum, empty_frame);
        // Add an empty IM centroid array so that this spectrum is consistent
        // with successfully picked spectra, which always carry ION_MOBILITY_CENTROID.
        // MassTraceDetection requires this array to be uniformly present or absent
        // across all spectra in the experiment.
        MSSpectrum::FloatDataArray empty_im_array;
        empty_im_array.setName(Constants::UserParam::ION_MOBILITY_CENTROID);
        empty_frame.getFloatDataArrays().push_back(std::move(empty_im_array));
        empty_frame.setType(SpectrumSettings::SpectrumType::CENTROID);
        empty_frame.setIMPeakType(IMPeakType::IM_CENTROIDED);
        spectrum = std::move(empty_frame);
        return;
      }

      // ------------------------------------------ step 4a: Extraction mobilograms ------------------------------------------
      // Using m/z peaks FWHM, we iteratively extract ion mobility traces (mobilograms) from the raw spectrum.
      // To rescue weak signal in the extracted mobilograms, we use linear resampling.
      // For linear resampling, it is recommended to use a sampling rate equal or higher than the raw sampling rate.
      // We dynamically determine the raw sampling rate from well-populated extracted mobilograms
      // (currently we have this hard-coded as +20 raw peaks in a mobilogram to be considered well-populated).

      auto mobilogram_traces = PeakPickerIM::extractIonMobilityTraces(picked_spectrum, spectrum);

      // ************************************************* PART II *****************************************************
      // ------------------------------------------ Ion mobility peak picking ------------------------------------------

      vector<Mobilogram> picked_traces;
      // Remove empty traces that can occur when no raw peaks are found within the FWHM window
      // of a picked m/z peak during extractIonMobilityTraces()
      mobilogram_traces.erase(
        std::remove_if(mobilogram_traces.begin(), mobilogram_traces.end(),
                   [](const auto& trace) { return trace.empty(); }),
        mobilogram_traces.end());

      // apply PeakPickerHiRes to pick ion mobility peaks.
      // PeakPickerHiRes can be applied to chromatograms. We reasoned the same set of parameters ideal for
      // chromatograms is also applicable for mobilograms.
      // Each raw mobilogram contains a float data array with raw m/z values.
      // We will use the ion mobility peak FWHM to define min/max ion mobility boundary
      // and recompute the m/z center based on the ion mobility peak.
      for (const Mobilogram& trace : mobilogram_traces)
      {
        // ------------------------------------------ part 2b: smooth and resample --------------------------------
        // Prepare mobilograms for SGolay smoothing.
        // To avoid edge effects, pad both edges with zero-intensity points at distance
        // (sgolay_frame_length_ - 1) / 2 * mobilogram_sampling_grid_ from the data boundary.
        // After linear resampling this produces (sgolay_frame_length_ - 1) / 2 zero-baseline
        // nodes on each side, giving SGolay a valid window at the mobilogram edges.
        double im_start = trace.front().getMobility();
        double im_end = trace.back().getMobility();
        int padding_points = static_cast<int>(std::ceil((sgolay_frame_length_ - 1) / 2.0));

        Mobilogram padded_trace;
        padded_trace.reserve(trace.size() + 2);
        padded_trace.emplace_back(im_start - padding_points * mobilogram_sampling_grid_, 0.0);
        for (const auto& peak : trace) padded_trace.push_back(peak);
        padded_trace.emplace_back(im_end + padding_points * mobilogram_sampling_grid_, 0.0);

        // Linear resample onto uniform grid; LinearResamplerAlign handles duplicate
        // IM positions in the raw trace correctly (distributes intensity to bracketing nodes).
        lin_resampler_.raster(padded_trace);
        // SGolay smooth prior to peak picking
        sgolay_filter_.filter(padded_trace);

        Mobilogram picked_trace;
        std::vector<PeakPickerHiRes::PeakBoundary> boundaries;
        picker_im_.pick(padded_trace, picked_trace, boundaries, true);
        picked_traces.push_back(std::move(picked_trace));
      }

      // Recompute m/z centers and output centroided frame
      MSSpectrum centroided_frame = computeCentroids_(mobilogram_traces, picked_traces);

      // Remove extra FDAs (IM FWHM, MZ FWHM)
      removeAllFloatDataArraysExcept(centroided_frame, Constants::UserParam::ION_MOBILITY_CENTROID);

      // Copy only SpectrumSettings from the input into the centroided result
      static_cast<SpectrumSettings&>(centroided_frame) = static_cast<const SpectrumSettings&>(spectrum);
      centroided_frame.setMSLevel(spectrum.getMSLevel());
      centroided_frame.setName(spectrum.getName());
      centroided_frame.setRT(spectrum.getRT());
      removeAllFloatDataArraysExcept(centroided_frame, Constants::UserParam::ION_MOBILITY_CENTROID);
      centroided_frame.setIMFormat(IMFormat::IM_PEAK);
      centroided_frame.setIMPeakType(IMPeakType::IM_CENTROIDED);
      spectrum = std::move(centroided_frame);

    }

#if 0
    void PeakPickerIM::pickIMCluster(OpenMS::MSSpectrum& spectrum) const
    {
      if (spectrum.empty()) return;

      // Validate IM format - returns false if we should skip processing
      if (!validateIMFormatForPicking(spectrum))
      {
        return;
      }

      // Get IM data array
      if (!spectrum.containsIMData())
      {
        OPENMS_LOG_WARN << "No ion mobility data found in spectrum.\n";
        return;
      }
      const auto [im_data_index, im_unit] = spectrum.getIMData();
      auto& im_data = spectrum.getFloatDataArrays()[im_data_index];

      // Warn if CCS data with small tolerance (only once per PeakPickerIM instance)
      if (im_unit == DriftTimeUnit::CCS && im_tolerance_cluster_ < 1.0 && !ccs_warning_shown_)
      {
        OPENMS_LOG_WARN << "Warning: Ion mobility data is in CCS units (square angstroms), but im_tolerance_cluster"
                        << " = " << im_tolerance_cluster_ << " appears to be set for 1/K0 data. "
                        << "For CCS data, consider using larger values (e.g., 10-20 for clustering, 1.0 for summing)." << '\n';
        ccs_warning_shown_ = true;
      }

      struct Point {
        double mz;
        double im;
        float intensity;
        OpenMS::Size original_index;

        Point(double mz_val, double im_val, float int_val, OpenMS::Size idx) :
            mz(mz_val), im(im_val), intensity(int_val), original_index(idx) {}
      };

      // Convert peaks to Points (same as before)
      std::vector<Point> points;
      points.reserve(spectrum.size());
      for (OpenMS::Size i = 0; i < spectrum.size(); ++i) {
        const auto& peak = spectrum[i];
        points.emplace_back(peak.getMZ(), im_data[i], peak.getIntensity(), i);
      }

      // --- Setup for Clustering ---

      // 1. Create m/z-sorted indices (needed for cluster expansion)
      std::vector<OpenMS::Size> mz_sorted_indices(points.size());
      std::iota(mz_sorted_indices.begin(), mz_sorted_indices.end(), 0);
      std::sort(mz_sorted_indices.begin(), mz_sorted_indices.end(),
                [&](OpenMS::Size a, OpenMS::Size b) {
                  return points[a].mz < points[b].mz;
                });

      // Create reverse lookup: original_index -> mz_sorted_position (needed for cluster expansion)
      std::vector<OpenMS::Size> original_to_sorted_pos(points.size());
      for (OpenMS::Size i = 0; i < points.size(); ++i) {
        original_to_sorted_pos[mz_sorted_indices[i]] = i;
      }

      // 2. ***OPTIMIZATION: Create intensity-sorted indices for seed picking***
      std::vector<OpenMS::Size> intensity_sorted_indices(points.size());
      std::iota(intensity_sorted_indices.begin(), intensity_sorted_indices.end(), 0);
      std::sort(intensity_sorted_indices.begin(), intensity_sorted_indices.end(),
                [&](OpenMS::Size a, OpenMS::Size b) {
                  // Sort DESCENDING by intensity
                  return points[a].intensity > points[b].intensity;
                });

      // 3. Keep track of used points
      std::vector<bool> used(points.size(), false);
      OpenMS::Size num_used = 0;

      // Store results temporarily
      std::vector<Point> averaged_points;
      averaged_points.reserve(points.size());

      double ppm_factor = ppm_tolerance_cluster_ * 1e-6;

      // --- Main Clustering Loop ---
      // Iterate through peaks in descending order of intensity to find seeds
      for (OpenMS::Size current_intensity_rank = 0; current_intensity_rank < points.size(); ++current_intensity_rank)
      {
        // 4. Get the original index of the next potential seed (highest intensity first)
        OpenMS::Size seed_original_idx = intensity_sorted_indices[current_intensity_rank];

        // 5. Check if this peak has already been used (part of a previous cluster)
        if (used[seed_original_idx]) {
          continue; // Skip to the next highest intensity peak
        }

        // --- Found an unused seed, start clustering ---

        // 6. Initialize the cluster with the seed
        std::vector<OpenMS::Size> current_cluster_indices;
        current_cluster_indices.push_back(seed_original_idx);
        used[seed_original_idx] = true;
        num_used++;

        double cluster_mz_min = points[seed_original_idx].mz;
        double cluster_mz_max = points[seed_original_idx].mz;
        double cluster_im_min = points[seed_original_idx].im;
        double cluster_im_max = points[seed_original_idx].im;

        // 7. Expand the cluster using m/z sorted neighbors (same logic as before)
        OpenMS::Size seed_sorted_pos = original_to_sorted_pos[seed_original_idx];
        OpenMS::SignedSize left_idx = static_cast<OpenMS::SignedSize>(seed_sorted_pos) - 1;
        OpenMS::SignedSize right_idx = static_cast<OpenMS::SignedSize>(seed_sorted_pos) + 1;

        bool changed = true;
        while (changed) {
          changed = false;

          // Check left neighbor
          while (left_idx >= 0) {
            OpenMS::Size candidate_original_idx = mz_sorted_indices[left_idx];
            if (!used[candidate_original_idx]) {
              const auto& candidate_point = points[candidate_original_idx];
              double potential_mz_min = std::min(cluster_mz_min, candidate_point.mz);
              double potential_mz_max = std::max(cluster_mz_max, candidate_point.mz); // Fixed: Use max for the max value
              double potential_im_min = std::min(cluster_im_min, candidate_point.im);
              double potential_im_max = std::max(cluster_im_max, candidate_point.im);

              // Fixed m/z tolerance calculation
              bool mz_ok = (potential_mz_max - potential_mz_min) <= (potential_mz_min * ppm_factor);
              bool im_ok = (potential_im_max - potential_im_min) <= im_tolerance_cluster_;

              if (mz_ok && im_ok) {
                current_cluster_indices.push_back(candidate_original_idx);
                used[candidate_original_idx] = true;
                num_used++;
                cluster_mz_min = potential_mz_min;
                cluster_mz_max = potential_mz_max; // Fixed: Update max value
                cluster_im_min = potential_im_min;
                cluster_im_max = potential_im_max;
                left_idx--;
                changed = true;
                break; // Added point, break inner while to re-evaluate from outer changed loop
              } else {
                left_idx = -1; // Stop checking left this round
                break;
              }
            } else {
              left_idx--; // Skip used point
            }
          }

          // Check right neighbor
          while (right_idx < static_cast<OpenMS::SignedSize>(points.size())) {
            OpenMS::Size candidate_original_idx = mz_sorted_indices[right_idx];
            if (!used[candidate_original_idx]) {
              const auto& candidate_point = points[candidate_original_idx];
              double potential_mz_min = std::min(cluster_mz_min, candidate_point.mz); // Fixed: Use min
              double potential_mz_max = std::max(cluster_mz_max, candidate_point.mz);
              double potential_im_min = std::min(cluster_im_min, candidate_point.im);
              double potential_im_max = std::max(cluster_im_max, candidate_point.im);

              // Fixed m/z tolerance calculation
              bool mz_ok = (potential_mz_max - potential_mz_min) <= (potential_mz_min * ppm_factor);
              bool im_ok = (potential_im_max - potential_im_min) <= im_tolerance_cluster_;

              if (mz_ok && im_ok) {
                current_cluster_indices.push_back(candidate_original_idx);
                used[candidate_original_idx] = true;
                num_used++;
                cluster_mz_min = potential_mz_min; // Fixed: Update min value
                cluster_mz_max = potential_mz_max;
                cluster_im_min = potential_im_min;
                cluster_im_max = potential_im_max;
                right_idx++;
                changed = true;
                break; // Added point, break inner while to re-evaluate from outer changed loop
              } else {
                right_idx = points.size(); // Stop checking right this round
                break;
              }
            } else {
              right_idx++; // Skip used point
            }
          }
        } // End cluster expansion (while changed)

        // 8. Finalize the current cluster (same logic as before)
        if (!current_cluster_indices.empty()) {
          double sum_intensity = 0.0;
          double sum_mz_intensity = 0.0;
          double sum_im_intensity = 0.0;

          for (OpenMS::Size original_idx : current_cluster_indices) {
            const auto& p = points[original_idx];
            sum_intensity += p.intensity;
            sum_mz_intensity += p.mz * p.intensity;
            sum_im_intensity += p.im * p.intensity;
          }

          if (sum_intensity > std::numeric_limits<double>::epsilon()) {
            averaged_points.emplace_back(
              sum_mz_intensity / sum_intensity,
              sum_im_intensity / sum_intensity,
              static_cast<float>(sum_intensity),
              0 // Original index is meaningless for averaged points
            );
          }
        }

        // Optimization: Check if all points are processed
        if (num_used == points.size()) {
          break; // Exit the main loop early
        }

      } // End main loop (for intensity_sorted_indices)

      // 9. Update spectrum (same logic as before)
      // Clear DataArrays that won't be valid after averaging (StringDataArrays and IntegerDataArrays
      // don't correspond to the new averaged peaks, so they must be cleared to maintain consistency)
      spectrum.getStringDataArrays().clear();
      spectrum.getIntegerDataArrays().clear();

      spectrum.resize(averaged_points.size());
      spectrum.shrink_to_fit();
      im_data.resize(averaged_points.size());
      im_data.shrink_to_fit();

      for (size_t i = 0; i != averaged_points.size(); ++i)
      {
        const auto& p = averaged_points[i];
        spectrum[i].setMZ(p.mz);
        spectrum[i].setIntensity(p.intensity);
        im_data[i] = p.im;
      }

      spectrum.sortByPosition();
      spectrum.updateRanges();
      // ensure the output IM array is updated
      spectrum.getFloatDataArrays()[im_data_index].setName(Constants::UserParam::ION_MOBILITY_CENTROID);
      spectrum.setIMFormat(IMFormat::IM_PEAK);
      spectrum.setIMPeakType(IMPeakType::IM_CENTROIDED);
      removeAllFloatDataArraysExcept(spectrum, Constants::UserParam::ION_MOBILITY_CENTROID);
    } // End of pickIMCluster function

    void PeakPickerIM::pickIMElutionProfiles(MSSpectrum& input) const
    {
      if (input.empty()) return;

      // Validate IM format - returns false if we should skip processing
      if (!validateIMFormatForPicking(input))
      {
        return;
      }

      // Get IM data array
      if (!input.containsIMData())
      {
        OPENMS_LOG_WARN << "No ion mobility data found in spectrum.\n";
        return;
      }
      const auto [im_data_index, im_unit] = input.getIMData();
      auto& im_data = input.getFloatDataArrays()[im_data_index];
      // convert to MSExperiment and set drift time as RT
      MSExperiment frame_as_spectra = IMDataConverter::reshapeIMFrameToMany(input);
      for (auto& s : frame_as_spectra)
      {
        s.setRT(s.getDriftTime());
        s.setDriftTime(-1);
        s.setMSLevel(1);
      }

#ifdef DEBUG_IM_PICKER
  // write out IM frame as RT/MZ for debugging purposes to test algorithm that yet don't support the IM dimension
  MzMLFile().store("debug" + String(input.getRT()) + ".mzML", frame_as_spectra);
#endif

      if (frame_as_spectra.size() <= 3 ) return;

      // detect mass traces in IM frame
      MassTraceDetection mte;
      Param param = mte.getParameters();
      // disable most filter criteria
      param.setValue("min_trace_length", -1.0);
      param.setValue("max_trace_length", -1.0);
      param.setValue("noise_threshold_int", 0.1); // only ignore 0 peaks
      param.setValue("chrom_peak_snr", 0.0);
      param.setValue("reestimate_mt_sd", "false");
      param.setValue("mass_error_ppm", ppm_tolerance_elution_);
      param.setValue("trace_termination_criterion", "outlier");
      param.setValue("trace_termination_outlier", 1);

      mte.setLogType(ProgressLogger::NONE);
      mte.setParameters(param);
      vector<MassTrace> output_mt;
      mte.run(frame_as_spectra, output_mt);

      ElutionPeakDetection epd;
      param = epd.getParameters();
      param.setValue("chrom_fwhm", 0.01);
      param.setValue("chrom_peak_snr", 0.0);
      param.setValue("width_filtering", "off");
      param.setValue("min_fwhm", -1.0);
      param.setValue("max_fwhm", 1e6);
      param.setValue("masstrace_snr_filtering", "false");
      epd.setParameters(param);

      std::vector<MassTrace> split_mtraces;
      epd.detectPeaks(output_mt, split_mtraces);
      output_mt.clear();

      // copy mass traces centroids back to peaks
      // Clear DataArrays that won't be valid after peak detection (StringDataArrays and IntegerDataArrays
      // don't correspond to the new detected peaks, so they must be cleared to maintain consistency)
      input.getStringDataArrays().clear();
      input.getIntegerDataArrays().clear();

      input.resize(split_mtraces.size());
      input.shrink_to_fit();

      im_data.resize(split_mtraces.size());
      im_data.shrink_to_fit();

      for (Size i = 0; i < split_mtraces.size(); ++i)
      {
        const MassTrace& mt = split_mtraces[i];
        input[i].setMZ(mt.getCentroidMZ());
        // Use computeIntensitySum() instead of getIntensity() because getIntensity() depends on
        // FWHM indices which may not be set for short traces (returns 0 if fwhm_start_idx_ == fwhm_end_idx_ == 0)
        input[i].setIntensity(mt.computeIntensitySum());
        im_data[i] = mt.getCentroidRT(); // IM
      }
      input.sortByPosition();
      input.updateRanges();
      // ensure the output im name is updated
      input.getFloatDataArrays()[im_data_index].setName(Constants::UserParam::ION_MOBILITY_CENTROID);
      input.setIMFormat(IMFormat::IM_PEAK);
      input.setIMPeakType(IMPeakType::IM_CENTROIDED);
      removeAllFloatDataArraysExcept(input, Constants::UserParam::ION_MOBILITY_CENTROID);
    }
#endif

    void PeakPickerIM::aggregateScans(const std::vector<MSSpectrum>& spectra,
                                       const std::vector<double>& weights,
                                       MSSpectrum& aggregated_spectrum) const
    {
      std::vector<const MSSpectrum*> spectrum_ptrs;
      spectrum_ptrs.reserve(spectra.size());
      for (const MSSpectrum& spec : spectra) spectrum_ptrs.push_back(&spec);
      aggregateScans(spectrum_ptrs, weights, aggregated_spectrum);
    }

    void PeakPickerIM::aggregateScans(const std::vector<const MSSpectrum*>& spectra,
                                       const std::vector<double>& weights,
                                       MSSpectrum& aggregated_spectrum) const
    {
      aggregated_spectrum.clear(true);

      if (spectra.empty())
      {
        OPENMS_LOG_WARN << "aggregateScans: No spectra provided for aggregation.\n";
        return;
      }

      if (spectra.size() != weights.size())
      {
        OPENMS_LOG_WARN << "aggregateScans: Spectra and weights size mismatch.\n";
        return;
      }

      // Collect the non-empty input spectra with their IM arrays (empty spectra contribute no peaks)
      struct Input
      {
        const MSSpectrum* spec;
        const MSSpectrum::FloatDataArray* im;
        double weight;
      };
      std::vector<Input> inputs;
      Size total_peaks = 0;

      // m/z-sorted copies of inputs that are not sorted (reserved up front so pointers stay valid)
      std::vector<MSSpectrum> sorted_copies;
      sorted_copies.reserve(spectra.size());

      MSSpectrum::FloatDataArray aggregated_im_array;

      for (Size spec_idx = 0; spec_idx < spectra.size(); ++spec_idx)
      {
        const MSSpectrum* spec_ptr = spectra[spec_idx];
        if (spec_ptr->empty()) continue;

        const Size im_data_index = spec_ptr->getIMData().first;

        // Verify IM array size matches spectrum size
        if (spec_ptr->getFloatDataArrays()[im_data_index].size() != spec_ptr->size())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Ion mobility array size (" + String(spec_ptr->getFloatDataArrays()[im_data_index].size()) +
              ") does not match spectrum size (" + String(spec_ptr->size()) + ") at RT " + String(spec_ptr->getRT()),
              String(spec_ptr->getFloatDataArrays()[im_data_index].size()) + " != " + String(spec_ptr->size()));
        }

        // The merge below requires m/z-sorted inputs: sort a copy of an unsorted input
        if (!spec_ptr->isSorted())
        {
          sorted_copies.push_back(*spec_ptr);
          sorted_copies.back().sortByPosition();
          spec_ptr = &sorted_copies.back();
        }
        const MSSpectrum& spec = *spec_ptr;

        const auto& im_array = spec.getFloatDataArrays()[im_data_index];
        aggregated_im_array.setName(im_array.getName());

        inputs.push_back({&spec, &im_array, weights[spec_idx]});
        total_peaks += spec.size();
      }

      if (total_peaks == 0)
      {
        OPENMS_LOG_WARN << "aggregateScans: No peaks collected after aggregation.\n";
        return;
      }

      aggregated_spectrum.reserve(total_peaks);
      aggregated_im_array.reserve(total_peaks);

      // k-way merge of the sorted inputs. On equal m/z the earlier input wins, which reproduces the
      // order of a stable sort of the concatenated inputs.
      std::vector<Size> pos(inputs.size(), 0);
      for (Size n = 0; n < total_peaks; ++n)
      {
        Size best = inputs.size();
        for (Size k = 0; k < inputs.size(); ++k)
        {
          if (pos[k] == inputs[k].spec->size()) continue; // this input is exhausted
          if (best == inputs.size() ||
              (*inputs[k].spec)[pos[k]].getMZ() < (*inputs[best].spec)[pos[best]].getMZ())
          {
            best = k;
          }
        }
        const Input& in = inputs[best];
        const Size i = pos[best]++;

        Peak1D weighted_peak;
        weighted_peak.setMZ((*in.spec)[i].getMZ());
        weighted_peak.setIntensity((*in.spec)[i].getIntensity() * in.weight);
        aggregated_spectrum.push_back(weighted_peak);
        aggregated_im_array.push_back((*in.im)[i]);
      }

      // Attach the ion mobility array
      aggregated_spectrum.getFloatDataArrays().push_back(std::move(aggregated_im_array));

      // Set metadata from the center spectrum (index 0)
      const MSSpectrum& center_spec = *spectra[0];
      static_cast<SpectrumSettings&>(aggregated_spectrum) = static_cast<const SpectrumSettings&>(center_spec);
      aggregated_spectrum.setMSLevel(center_spec.getMSLevel());
      aggregated_spectrum.setName(center_spec.getName());
      aggregated_spectrum.setRT(center_spec.getRT());
      aggregated_spectrum.setIMFormat(center_spec.getIMFormat());
      aggregated_spectrum.setDriftTimeUnit(center_spec.getDriftTimeUnit());
    }

    void PeakPickerIM::pickExperimentWithAggregation(MSExperiment& exp)
    {
      if (exp.empty())
      {
        OPENMS_LOG_WARN << "pickExperimentWithAggregation: Empty experiment provided.\n";
        return;
      }

      // Gaussian weighting parameters
      double fwhm = aggregation_rt_fwhm_;
      double factor = -4.0 * std::log(2.0) / (fwhm * fwhm);
      double cutoff = aggregation_cutoff_;

      OPENMS_LOG_INFO << "pickExperimentWithAggregation: Using rt_FWHM=" << fwhm
                      << "s, cutoff=" << cutoff << "\n";

      // Build list of MS1 spectrum indices
      std::vector<Size> ms1_indices;
      for (Size i = 0; i < exp.size(); ++i)
      {
        if (exp[i].getMSLevel() == 1)
        {
          ms1_indices.push_back(i);
        }
      }

      OPENMS_LOG_INFO << "pickExperimentWithAggregation: Found " << ms1_indices.size()
                      << " MS1 spectra out of " << exp.size() << " total spectra.\n";

      if (ms1_indices.empty())
      {
        OPENMS_LOG_WARN << "pickExperimentWithAggregation: No MS1 spectra found.\n";
        return;
      }

      // Log RT range of MS1 spectra
      double rt_min = exp[ms1_indices.front()].getRT();
      double rt_max = exp[ms1_indices.back()].getRT();
      OPENMS_LOG_INFO << "pickExperimentWithAggregation: MS1 RT range: [" << rt_min << ", " << rt_max << "] seconds.\n";

      // Check if first MS1 spectrum has IM data - required for PeakPickerIM
      if (!exp[ms1_indices[0]].containsIMData())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Ion mobility data is not detected. PeakPickerIM is designed to pick peaks from ion mobility containing data.",
            "No IM data in first MS1 spectrum");
      }


      // Build AggregationBlocks: for each spectrum, collect neighbors with Gaussian weights
      AggregationBlocks aggregation_blocks;

      for (Size idx = 0; idx < ms1_indices.size(); ++idx)
      {
        Size center_exp_idx = ms1_indices[idx];
        double center_rt = exp[center_exp_idx].getRT();

        std::vector<std::pair<Size, double>> neighbors_with_weights;

        // Go forward from center
        for (Size j = idx; j < ms1_indices.size(); ++j)
        {
          double rt_diff = exp[ms1_indices[j]].getRT() - center_rt;
          double weight = std::exp(factor * rt_diff * rt_diff);

          if (weight < cutoff && j != idx)
          {
            break;
          }

          neighbors_with_weights.emplace_back(ms1_indices[j], weight);
        }

        // Go backward from center
        for (SignedSize j = static_cast<SignedSize>(idx) - 1; j >= 0; --j)
        {
          double rt_diff = exp[ms1_indices[j]].getRT() - center_rt;
          double weight = std::exp(factor * rt_diff * rt_diff);

          if (weight < cutoff)
          {
            break;
          }

          neighbors_with_weights.emplace_back(ms1_indices[j], weight);
        }

        // Normalize weights to sum to 1
        double sum_weights = 0.0;
        for (const auto& p : neighbors_with_weights)
        {
          sum_weights += p.second;
        }
        for (auto& p : neighbors_with_weights)
        {
          p.second /= sum_weights;
        }

        aggregation_blocks[center_exp_idx] = std::move(neighbors_with_weights);
      }

      // Log aggregation statistics
      Size total_neighbors = 0;
      Size max_neighbors = 0;
      Size min_neighbors = std::numeric_limits<Size>::max();
      for (const auto& [center_idx, neighbors] : aggregation_blocks)
      {
        total_neighbors += neighbors.size();
        max_neighbors = std::max(max_neighbors, neighbors.size());
        min_neighbors = std::min(min_neighbors, neighbors.size());
      }
      double avg_neighbors = aggregation_blocks.empty() ? 0.0 : static_cast<double>(total_neighbors) / aggregation_blocks.size();
      OPENMS_LOG_INFO << "pickExperimentWithAggregation: Aggregation blocks built. "
                      << "Avg neighbors per spectrum: " << avg_neighbors
                      << ", min: " << min_neighbors << ", max: " << max_neighbors << "\n";

      if (max_neighbors == 1)
      {
        OPENMS_LOG_WARN << "pickExperimentWithAggregation: Each spectrum only has 1 neighbor (itself). "
                        << "No actual aggregation will occur. Consider increasing rt_FWHM parameter.\n";
      }

      // Process each spectrum with its aggregation block.
      // Important: collect ALL results before writing any back to exp.
      // Writing back immediately would cause cascading: later iterations would
      // read already-aggregated (bloated) neighbors instead of original spectra,
      // leading to progressive memory growth.
      std::vector<std::pair<Size, MSSpectrum>> aggregated_results;
      aggregated_results.reserve(aggregation_blocks.size());

      for (const auto& [center_idx, neighbors] : aggregation_blocks)
      {
        // Collect spectra and weights (all reads from unmodified exp)
        std::vector<MSSpectrum> spectra_to_aggregate;
        std::vector<double> weights;
        spectra_to_aggregate.reserve(neighbors.size());
        weights.reserve(neighbors.size());

        for (const auto& [spec_idx, weight] : neighbors)
        {
          spectra_to_aggregate.push_back(exp[spec_idx]);
          weights.push_back(weight);
        }

        // Aggregate with weights
        MSSpectrum aggregated;
        aggregateScans(spectra_to_aggregate, weights, aggregated);

        if (!aggregated.empty())
        {
          aggregated_results.emplace_back(center_idx, std::move(aggregated));
        }
      }

      // Write all results back at once (no cascading possible)
      for (auto& [idx, spec] : aggregated_results)
      {
        exp[idx] = std::move(spec);
      }
    }

} // namespace OpenMS
