// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Author: Timo Sachsenberg, Mohammed Alhigaylan $
// $Maintainer: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/Mobilogram.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/PROCESSING/RESAMPLING/LinearResamplerAlign.h>
#include <OpenMS/PROCESSING/SMOOTHING/GaussFilter.h>
#include <OpenMS/PROCESSING/SMOOTHING/SavitzkyGolayFilter.h>

#include <map>
#include <vector>

namespace OpenMS
{
  class MSExperiment;
  /**
    @brief Peak picking algorithm for ion mobility data
    
    **pickIMTraces**: Mobilogram-based peak picking that extracts ion mobility traces
    from raw IM data and performs centroiding on the extracted mobilograms. This method
    processes IM data by analyzing intensity profiles along the ion mobility dimension.

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


#if 0
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
#endif

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
     * @brief Same as above, but reads the spectra in place (no copies) and merges them instead of
     * concatenating and sorting.
     *
     * Inputs that are not sorted by m/z are sorted (as a copy) first. The merge gives the same peak
     * order as a stable sort of the concatenated inputs (peaks with equal m/z keep input order).
     */
    void aggregateScans(const std::vector<const MSSpectrum*>& spectra,
                        const std::vector<double>& weights,
                        MSSpectrum& aggregated_spectrum) const;

    /**
     * @brief Copies the peaks of all given spectra into one spectrum, unscaled and merged by m/z; no weights are applied.
     *
     * The center spectrum must be at index 0 and provides the metadata. All input spectra must contain ion mobility data.
     */
    void aggregateScans(const std::vector<const MSSpectrum*>& spectra,
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
     * - aggregation:frames_per_side: if >= 0, copy the peaks of this many frames on each side into the centre frame, unscaled; no weights, rt_FWHM and cutoff are then not used
     *
     * @param[in,out] exp The experiment to process. Each MS1 spectrum will be replaced with
     *                    an aggregated version combining signal from neighboring scans.
     *
     * @throws Exception::InvalidValue if no ion mobility data is detected in the experiment.
     */
    void pickExperimentWithAggregation(MSExperiment& exp);

    /**
     * @brief Frames whose peaks are copied into a centre frame when aggregating by a number of frames
     * (aggregation:frames_per_side >= 0), instead of by Gaussian weights.
     *
     * The centre frame comes first, then up to frames_per_side following frames, then up to
     * frames_per_side preceding frames (fewer near the ends of the experiment). Frames are taken by position.
     *
     * @param[in] n_frames number of frames that are aggregated together (one MS level, one window)
     * @param[in] centre index of the centre frame, must be smaller than n_frames
     * @param[in] frames_per_side frames pooled on each side of the centre frame
     * @return indices of the frames, the centre frame first
     */
    static std::vector<Size> neighbouringFrames(Size n_frames, Size centre, Size frames_per_side);

    /**
     * @brief Estimates the spacing of the raw m/z lattice in sqrt(m/z) (one m/z column).
     *
     * TimsTOF m/z values lie on a lattice in sqrt(m/z) that is constant within an acquisition file but
     * differs between files. The step is the dominant smallest gap between neighbouring distinct
     * m/z values; gaps below 1e-6 (duplicates across scans) are ignored.
     * @param[in] spectra Raw spectra sorted by m/z, ideally several frames spread over the file
     * @return the step, or 0 if it could not be estimated
     */
    static double estimateMzGridStep(const std::vector<const MSSpectrum*>& spectra);

    /// Sets the m/z grid step (sqrt(m/z)) needed by the neighbor filter; call once per file with estimateMzGridStep()
    void setMzGridStep(double step) { mz_grid_step_ = step; }

    /// True if the m/z neighbor filter is enabled (pickIMTraces:neighbor_filter)
    bool neighborFilterEnabled() const { return neighbor_filter_; }

    /// True if m/z peaks are summed per raw m/z lattice column (pickIMTraces:mz_binning = lattice); needs setMzGridStep()
    bool latticeBinningEnabled() const { return mz_lattice_binning_; }

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
     *
     * With pickIMTraces:mz_binning = lattice the tolerance is ignored: peaks are summed if they fall into the same
     * column of the raw m/z lattice, i.e. their sqrt(m/z) differs from the first peak of the group by less than half
     * of mz_grid_step_ (the lattice offset is not needed).
     */
    void sumFrame_(const MSSpectrum& input_spectrum, MSSpectrum& output_spectrum, double tolerance = 0.01, bool use_ppm = true);

    /// Compute lower and upper m/z bounds based on ppm
    std::pair<double, double> ppmBounds(double mz, double ppm);

    /// Extract ion mobility traces as Mobilograms from the raw TimsTOF frame
    /// raw m/z values are allocated to float data arrays with the label 'raw_mz'
    std::vector<Mobilogram> extractIonMobilityTraces(
      const MSSpectrum& picked_spectrum,
      const MSSpectrum& raw_spectrum);

    /// compute m/z and ion mobility centers for picked traces. Returns centroided spectrum.
    MSSpectrum computeCentroids_(const std::vector<Mobilogram>& mobilogram_traces,
                              const std::vector<Mobilogram>& picked_traces,
                              const MSSpectrum& raw_spectrum);

    /**
     * @brief m/z-column neighbour support of the raw peaks of a picked peak box.
     *
     * Raw timsTOF m/z values lie on a lattice in sqrt(m/z) with spacing mz_grid_step_; one lattice step is
     * called a column. Each raw peak is assigned the column index round((sqrt(mz) - sqrt(mz_min)) / mz_grid_step_)
     * relative to the lowest raw peak of the box.
     *
     * The box is supported (returns true) if either
     *  - one column holds at least singleton_column_min_peaks_ raw peaks, or
     *  - two occupied columns are at most neighbor_max_column_gap_ columns apart.
     * The two rules are independent: the per-column count only matters for the first, and the neighbours
     * need not be symmetric about the peak centre. The box is never trimmed or moved; it is kept or rejected whole.
     * @return false if the box has no neighbour support and should not be reported
     */
    bool hasNeighborSupport_(const MSSpectrum& peaks) const;


    double sum_tolerance_mz_{1.0};
    double gauss_ppm_tolerance_{5.0};
    bool gauss_exact_kernel_{false};
    double mobilogram_sampling_grid_{0.01};
    int sgolay_frame_length_{5};
    int sgolay_polynomial_order_{3};
    Size min_raw_peaks_{1};
    double min_signal_to_noise_{0.0};
    double noise_window_mz_{200.0};
    double noise_window_im_{0.04};
    bool neighbor_filter_{false};
    bool mz_lattice_binning_{false};
    int singleton_column_min_peaks_{2};
    int neighbor_max_column_gap_{1};
    double mz_grid_step_{0.0};   ///< spacing of the raw m/z lattice in sqrt(m/z); see setMzGridStep()

    /// Filters and pickers used by pickIMTraces(), configured once in updateMembers_()
    /// instead of per spectrum (their construction and parameter setup dominated runtime)
    GaussFilter gauss_filter_;               ///< m/z smoothing of the summed frame
    PeakPickerHiRes picker_mz_;              ///< m/z peak picking (reports FWHM in ppm)
    LinearResamplerAlign lin_resampler_;     ///< mobilogram resampling onto a uniform grid
    SavitzkyGolayFilter sgolay_filter_;      ///< mobilogram smoothing
    PeakPickerHiRes picker_im_;              ///< ion mobility peak picking (reports absolute FWHM)

#if 0
    double ppm_tolerance_cluster_{50.0};
    double im_tolerance_cluster_{0.1};

    double ppm_tolerance_elution_{50.0};
#endif

    double aggregation_rt_fwhm_{1.0};   ///< Gaussian FWHM for RT-based weighting (in seconds)
    double aggregation_cutoff_{0.01};   ///< Weight cutoff for including spectra in aggregation
    int aggregation_frames_per_side_{-1}; ///< Frames pooled on each side of the centre frame; negative = use aggregation_rt_fwhm_

#if 0
    /// Flag to track if CCS tolerance warning has been shown (mutable for const methods)
    mutable bool ccs_warning_shown_{false};
#endif
  };
} // namespace OpenMS