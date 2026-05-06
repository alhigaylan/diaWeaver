// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Mohammed Alhigaylan $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MasstraceCorrelator.h>

namespace OpenMS
{

  /**
   @brief Clusters fragment ion mass traces with precursor mass traces to generate pseudo MS/MS spectra.

   This class correlates precursor elution profiles with fragment ion mass traces
   based on their chromatographic similarity (Pearson correlation). It implements
   the ETISEQ-like algorithm for DIA/SWATH data analysis.

   The algorithm works as follows:
   1. For each precursor within the SWATH isolation window, find candidate fragment traces
   2. Score precursor-fragment pairs using Pearson correlation of their elution profiles
   3. Assign fragments to precursors if correlation exceeds threshold
   4. Generate pseudo MS/MS spectra from the assignments

   This class supports two input modes:
   - ConsensusMap-based: for use with MassTraceExtractor output (consensusXML format)
   - FeatureMap/MassTrace-based: for use with FeatureFinderMetabo output

   Reference:
   ETISEQ -- an algorithm for automated elution time ion sequencing of concurrently
   fragmented peptides for mass spectrometry-based proteomics
   BMC Bioinformatics 2009, 10:244 doi:10.1186/1471-2105-10-244

   @ingroup OpenSWATH
   */
  class OPENMS_DLLAPI ClusterMassTracesByPrecursor :
    public DefaultParamHandler,
    public ProgressLogger
  {

  public:

    /// Default constructor
    ClusterMassTracesByPrecursor();

    /// Destructor
    ~ClusterMassTracesByPrecursor() override;

    /**
     * @brief Generate pseudo spectra from ConsensusMap inputs (MS1 and MS2 mass traces)
     *
     * This method takes mass traces stored in ConsensusMap format and generates
     * pseudo MS/MS spectra by correlating precursor and fragment elution profiles.
     *
     * @param[in] ms1_traces MS1 precursor mass traces in ConsensusMap format
     * @param[in] ms2_traces MS2 fragment mass traces in ConsensusMap format
     * @param[in] swath_lower Lower m/z bound of the DIA/SWATH isolation window
     * @param[in] swath_upper Upper m/z bound of the DIA/SWATH isolation window
     * @param[out] pseudo_spectra Output pseudo MS/MS spectra
     */
    void run(const ConsensusMap& ms1_traces,
             const ConsensusMap& ms2_traces,
             double swath_lower,
             double swath_upper,
             MSExperiment& pseudo_spectra);

    /**
     * @brief Generate pseudo spectra from FeatureMap and MassTrace inputs
     *
     * This method takes MS1 features (from FeatureFinderMetabo) and MS2 mass traces
     * and generates pseudo MS/MS spectra by correlating precursor and fragment
     * elution profiles.
     *
     * @param[in] ms1_features MS1 peptide features from FeatureFinderMetabo (monoisotopic)
     * @param[in] ms1_traces MS1 mass traces (for raw intensity data used in correlation)
     * @param[in] ms2_traces MS2 fragment mass traces
     * @param[in] swath_lower Lower m/z bound of the DIA/SWATH isolation window
     * @param[in] swath_upper Upper m/z bound of the DIA/SWATH isolation window
     * @param[out] pseudo_spectra Output pseudo MS/MS spectra
     */
    void run(const FeatureMap& ms1_features,
             const std::vector<MassTrace>& ms1_traces,
             const std::vector<MassTrace>& ms2_traces,
             double swath_lower,
             double swath_upper,
             MSExperiment& pseudo_spectra);

  protected:

    /// Update members after parameter change
    void updateMembers_() override;

    /**
     * @brief Core clustering algorithm using pre-computed elution profiles
     *
     * This is the internal implementation that works with cached elution profiles
     * and metadata. Both public run() methods prepare the data and call this method.
     *
     * @param[in] precursor_profiles Elution profiles (RT/intensity pairs) for each precursor
     * @param[in] precursor_mz m/z values for each precursor
     * @param[in] precursor_rt RT values for each precursor
     * @param[in] precursor_im Ion mobility values for each precursor (0 if not available)
     * @param[in] precursor_charge Charge states for each precursor
     * @param[in] precursor_intensity Intensities for each precursor
     * @param[in] fragment_profiles Elution profiles (RT/intensity pairs) for each fragment
     * @param[in] fragment_mz m/z values for each fragment
     * @param[in] fragment_rt RT values for each fragment
     * @param[in] fragment_im Ion mobility values for each fragment (0 if not available)
     * @param[in] fragment_intensity Intensities for each fragment
     * @param[in] swath_lower Lower m/z bound of the DIA/SWATH isolation window
     * @param[in] swath_upper Upper m/z bound of the DIA/SWATH isolation window
     * @param[in] has_im_data Flag indicating whether ion mobility data is present in the input
     * @param[out] pseudo_spectra Output pseudo MS/MS spectra
     */
    void clusterAndCreateSpectra_(
        const std::vector<MasstraceCorrelator::MasstracePointsType>& precursor_profiles,
        const std::vector<double>& precursor_mz,
        const std::vector<double>& precursor_rt,
        const std::vector<double>& precursor_im,
        const std::vector<int>& precursor_charge,
        const std::vector<double>& precursor_intensity,
        const std::vector<MasstraceCorrelator::MasstracePointsType>& fragment_profiles,
        const std::vector<double>& fragment_mz,
        const std::vector<double>& fragment_rt,
        const std::vector<double>& fragment_im,
        const std::vector<double>& fragment_intensity,
        double swath_lower,
        double swath_upper,
        bool has_im_data,
        MSExperiment& pseudo_spectra);

  private:

    /// Replace profile intensities with EPD-smoothed values from ConsensusFeature MetaValues.
    /// Warns and leaves profiles raw if smoothed_intensities MetaValue is absent.
    void applySmoothedIntensities_(std::vector<MasstraceCorrelator::MasstracePointsType>& profiles,
                                   const ConsensusMap& traces,
                                   const String& label);

    /// Minimum Pearson correlation to match elution profiles
    double min_pearson_correlation_;

    /// Maximum lag (in spectra) for correlation
    int max_lag_;

    /// Maximum RT apex difference (in seconds)
    double max_rt_apex_difference_;

    /// Minimum number of ions to report a spectrum
    Size min_nr_ions_;

    /// Ion mobility tolerance for matching precursors to fragments
    double im_tolerance_;

    /// Whether to assign unassigned fragments to all matching precursors
    bool assign_unassigned_to_all_;

    /// RT tolerance for matching up mass trace points in correlation
    double rt_tolerance_;

    /// Maximum number of precursors a fragment can be assigned to
    Size nr_precursors_per_fragment_;

    /// Weight for Pearson correlation in combined ranking score
    double pearson_weight_;

    /// Weight for delta RT in combined ranking score
    double delta_rt_weight_;

    /// Weight for delta ion mobility in combined ranking score
    double delta_im_weight_;

    /// Maximum number of fragment ions per output spectrum (0 = no limit)
    Size max_nr_ions_;

    /// Whether to rank assignments by combined score (true) or Pearson only (false)
    bool use_combined_scores_;

    /// Whether to output per-fragment scores as FloatDataArrays
    bool output_fragment_scores_;

    /// Whether to apply Savitzky-Golay smoothing to MS1 elution profiles before correlation
    bool smooth_ms1_;

    /// Whether to apply Savitzky-Golay smoothing to MS2 elution profiles before correlation
    bool smooth_ms2_;
  };

} // namespace OpenMS
