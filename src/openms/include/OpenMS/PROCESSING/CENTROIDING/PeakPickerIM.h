// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Author: Timo Sachsenberg, Mohammed Alhigaylan $
// $Maintainer: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <map>
#include <vector>

namespace OpenMS
{
  class MSExperiment;
  /**
    @brief Peak picking algorithm for ion mobility data
    
    This class provides three specialized methods for peak picking in ion mobility (IM) data:
    
    1. **pickIMTraces**: Mobilogram-based peak picking that extracts ion mobility traces
       from raw IM data and performs centroiding on the extracted mobilograms. This method
       processes IM data by analyzing intensity profiles along the ion mobility dimension.
    
    2. **pickIMCluster**: Clustering-based peak picking that groups peaks close in both
       m/z and ion mobility space. Peaks within specified m/z (ppm) and IM tolerances are
       averaged together using intensity-weighted averaging, reducing an IM frame to a
       single spectrum with representative peak positions.
    
    3. **pickIMElutionProfiles**: Elution profile-based peak picking that extracts peaks
       based on their elution characteristics across the ion mobility dimension. This method
       uses m/z tolerance (ppm) to identify and pick peaks from IM elution profiles.

  @ingroup PeakPicking
      */
  class OPENMS_DLLAPI PeakPickerIM : public DefaultParamHandler
  {
  public:
    /// Default constructor initializing parameters with default values.
    PeakPickerIM();

    /// Destructor.
    ~PeakPickerIM() override = default;

    /**
   * @brief Centroids ion mobility data by iteratively extracting mobilograms for each m/z peak centroid
   *
   * This function processes an MS spectrum containing ion mobility data (geared towards TimsTOF data)
   * Peaks in a given MS spectrum are projected to the m/z axis and centroided.
   * Then, the mobilogram of each m/z peak centroid is retrieved using the m/z peak FWHM.
   * Peak picking algorithm is applied to the mobilogram to resolve isobaric species with different ion mobility measurement.
   *
   * @param[in,out] spectrum Spectrum containing ion mobility data in its FloatDataArrays
   */
    void pickIMTraces(MSSpectrum& spectrum);

    /// Sets the parameters for peak picking.
    using DefaultParamHandler::setParameters;
    using DefaultParamHandler::getParameters;


    /**
     * @brief Converts an ion mobility frame to a single spectrum with averaged IM values
     *
     * This function takes an MS spectrum containing ion mobility data and reduces it to
     * a single spectrum where peaks that are close in both m/z and ion mobility space
     * are averaged together. The averaging is intensity-weighted for both m/z and ion
     * mobility values.
     *
     * The algorithm processes peaks sequentially and groups them based on two criteria:
     * 1. m/z tolerance: peaks must be within the configured ppm tolerance of each other
     * 2. ion mobility tolerance: the range of IM values must not exceed the configured tolerance
     *
     * Uses parameters pickIMCluster:ppm_tolerance_cluster and pickIMCluster:im_tolerance_cluster.
     *
     * @param[in,out] spec Spectrum containing ion mobility data in its FloatDataArrays
     *
     * @throws Exception::MissingInformation if input spectrum lacks ion mobility data
     *
     * @note The input spectrum should contain ion mobility data in its FloatDataArrays.
     *       The output spectrum will contain averaged peaks with their corresponding
     *       intensity-weighted average ion mobility values.
     *
     * Example:
     * @code
     * MSSpectrum spectrum;  // spectrum with IM FloatDataArrays
     * PeakPickerIM picker;
     * picker.pickIMCluster(spectrum);
     * @endcode
     */
    void pickIMCluster(MSSpectrum& spec) const;

    /**
     * @brief Picks ion mobility elution profiles from the given spectrum using eluting profiles.
     *
     * This function processes an MS spectrum containing ion mobility data and
     * extracts IM elution profiles based on the configured ppm tolerance.
     *
     * @param[in,out] input Spectrum containing ion mobility data in its FloatDataArrays
     */
    void pickIMElutionProfiles(MSSpectrum& input) const;

    /// Blocks of spectra to aggregate: maps master spectrum index to vector of (spectrum index, weight) pairs
    typedef std::map<Size, std::vector<std::pair<Size, double>>> AggregationBlocks;

    /**
     * @brief Aggregates peaks across multiple adjacent scans with Gaussian weighting.
     *
     * This function takes a vector of spectra with corresponding weights and combines
     * their peaks into a single spectrum. Peak intensities are multiplied by their
     * respective weights before aggregation. This allows for Gaussian-weighted signal
     * boosting where closer scans contribute more than distant ones.
     *
     * @param[in] spectra Vector of spectra to aggregate (must contain ion mobility data).
     *                    The center spectrum must be at index 0.
     * @param[in] weights Vector of weights corresponding to each spectrum (must sum to 1).
     * @param[out] aggregated_spectrum Output spectrum containing weighted peaks from input spectra.
     *                                 Ion mobility data is preserved in FloatDataArrays.
     *
     * @note The center spectrum (index 0) is used for metadata (RT, MS level, name, etc.).
     * @note Weights should be normalized to sum to 1 for proper intensity scaling.
     * @note All input spectra must contain ion mobility data.
     */
    void aggregateScans(const std::vector<MSSpectrum>& spectra,
                        const std::vector<double>& weights,
                        MSSpectrum& aggregated_spectrum) const;

    /**
     * @brief Aggregates adjacent scans in an experiment using Gaussian-weighted signal boosting.
     *
     * For each MS1 spectrum in the experiment, this method aggregates adjacent scans using
     * Gaussian weights based on RT distance. This approach boosts signal-to-noise by combining
     * peaks from multiple scans with weights that decrease with distance from the center scan.
     * The aggregated spectra can then be passed to peak picking methods (e.g., pickIMTraces).
     *
     * Uses parameters:
     * - aggregation:rt_FWHM: Full width at half maximum for Gaussian weighting (in seconds)
     * - aggregation:cutoff: Weight threshold below which spectra are not included
     *
     * @param[in,out] exp The experiment to process. Each MS1 spectrum will be replaced with
     *                    an aggregated version combining signal from neighboring scans.
     *
     * @throws Exception::InvalidValue if no ion mobility data is detected in the experiment.
     */
    void pickExperimentWithAggregation(MSExperiment& exp);

  protected:
    void updateMembers_() override;

  private:
    /**
     * @brief Sum up the intensity of data points with nearly identical float values.
     *
     * By default, this function assumes the tolerance provided is in parts per million.
     * But it can be adjusted to use absolute value tolerance.
     * @param[in] input_spectrum Sorted raw spectrum with duplicate peaks due to scan merging or presence of ion mobility data.
     * @param[out] output_spectrum Output spectrum containing the summed peaks.
     * @param[in] tolerance Mass tolerance between peaks
     * @param[in] use_ppm Whether to use parts per million tolerance. If set to False, absolute tolerance will be used.
     */
    void sumFrame_(const MSSpectrum& input_spectrum, MSSpectrum& output_spectrum, double tolerance = 0.01, bool use_ppm = true);

    /// Compute lower and upper m/z bounds based on ppm
    std::pair<double, double> ppmBounds(double mz, double ppm);

    /// Extract ion mobility traces as MSSpectra from the raw TimsTOF frame
    /// Ion mobility is temporarily written in place of m/z inside Peak1D object.
    /// raw m/z values are allocated to float data arrays with the label 'raw_mz'
    std::pair<std::vector<MSSpectrum>, std::vector<bool>> extractIonMobilityTraces(
      const MSSpectrum& picked_spectrum,
      const MSSpectrum& raw_spectrum);

    /// compute m/z and ion mobility centers for picked traces. Returns centroided spectrum.
    MSSpectrum computeCentroids_(const std::vector<MSSpectrum>& mobilogram_traces,
                              const std::vector<MSSpectrum>& picked_traces);

    /// This function takes in a boolean vector of unclaimed raw peaks from the method PickIMTraces
    /// and append them to the centroided spectrum. It uses PickIMCluster to group 'leftover' raw peaks.
    /// In  other words, a different peak picking approach is used for low intensity raw peaks not sufficient
    /// for ion mobilogram extraction.
    void Add_unclaimedPeaks(MSSpectrum& centroided_frame, const MSSpectrum& raw_frame, const std::vector<bool>& claimed) const;


    double sum_tolerance_mz_{1.0};
    double gauss_ppm_tolerance_{5.0};
    double mobilogram_sampling_grid_{0.01};
    int sgolay_frame_length_{5};
    int sgolay_polynomial_order_{3};
    bool include_unclaimed_{false};

    double ppm_tolerance_cluster_{50.0};
    double im_tolerance_cluster_{0.1};

    double ppm_tolerance_elution_{50.0};

    double aggregation_rt_fwhm_{1.0};   ///< Gaussian FWHM for RT-based weighting (in seconds)
    double aggregation_cutoff_{0.01};   ///< Weight cutoff for including spectra in aggregation

    /// Flag to track if CCS tolerance warning has been shown (mutable for const methods)
    mutable bool ccs_warning_shown_{false};
  };
} // namespace OpenMS