// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Mohammed Alhigaylan $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FEATUREFINDER/FeatureFindingMetabo.h>

namespace OpenMS
{

  /**
    @brief Method for the assembly of mass traces belonging to the same isotope
    pattern, i.e., that are compatible in retention times, mass-to-charge ratios,
    and isotope abundances. Tailored for peptide feature finding.

    In @ref FeatureFindingPeptide, mass traces detected by the @ref
    MassTraceDetection method and afterwards split into individual
    chromatographic peaks by the @ref ElutionPeakDetection method are assembled
    to composite features if they are compatible with respect to RTs, m/z ratios,
    and isotopic intensities. To this end, feature hypotheses are formulated
    exhaustively based on the set of mass traces detected within a local RT and
    m/z region. These feature hypotheses are scored by their similarity to real
    peptide isotope patterns. The score is derived from independent models for
    retention time shifts and m/z differences between isotopic mass traces.
    Hypotheses with correct or false isotopic abundances are distinguished by a
    SVM model. Mass traces that could not be assembled or low-intensity
    peptides with only a monoisotopic mass trace to observe are left in the
    resulting @ref FeatureMap as singletons with the undefined charge state of 0.

    @htmlinclude OpenMS_FeatureFindingPeptide.parameters

    @ingroup Quantitation
  */
  class OPENMS_DLLAPI FeatureFindingPeptide :
    public DefaultParamHandler,
    public ProgressLogger
  {
public:
    /// Default constructor
    FeatureFindingPeptide();

    /// Default destructor
    ~FeatureFindingPeptide() override;

    /// main method of FeatureFindingPeptide
    void run(std::vector<MassTrace>& input_mtraces, FeatureMap& output_featmap, std::vector<std::vector< OpenMS::MSChromatogram > >& output_chromatograms);

protected:
    void updateMembers_() override;

private:
    /**
     * Calculate the maximal and minimal mass defects of isotopes for a given set of elements.
     *
     * @param[in] alphabet   chemical alphabet (elements which are expected to be present)
     * @param[in] peakOffset integer distance between isotope peak and monoisotopic peak (minimum: 1)
     * @return an interval which should contain the isotopic peak. This interval is relative to the monoisotopic peak.
     */
    Range getTheoreticIsotopicMassWindow_(const std::vector<Element const *>& alphabet, int peakOffset) const;

    /**
     * @brief Check whether a neutral peptide mass falls within the expected peptide mass defect filter.
     * This function was implemented in Tsou et al. (DIA-UMPIRE) https://doi.org/10.1002/pmic.201500526
     * and adapted from Toumi et al. https://doi.org/10.1021/pr100291q
     *
     * Peptides have a characteristic mass defect pattern that can be described by two linear
     * boundaries as a function of nominal mass. This filter, adapted from DIA-Umpire, rejects
     * mass traces whose fractional mass falls outside the filter defined by:
     *
     *   upper = frac(0.00052738 * mass + 0.066015 + d)
     *   lower = frac(0.00042565 * mass + 0.00038210 - d)
     *
     * where frac(x) = x - floor(x) and d is a user-configurable tolerance offset.
     * increasing d is needed for modified peptides.
     * @param[in] neutral_mass Neutral monoisotopic peptide mass (in Da)
     * @param[in] d            Tolerance offset (default 0.1). Increase this for modified peptides.
     * @return true if the mass defect is within the expected peptide corridor
     */
    bool isMassDefectValid_(double neutral_mass, double d) const;

    /** @brief Computes the cosine similarity between two vectors
     *
     * The cosine similarity (or cosine distance) is the cosine of the angle
     * between two vectors or the normalized dot product of two vectors.
     *
     * See also https://en.wikipedia.org/wiki/Cosine_similarity
     *
     * @param[in] vec1 First vector
     * @param[in] vec2 Second vector
    */
    double computeCosineSim_(const std::vector<double>& vec1, const std::vector<double>& vec2) const;

    /** @brief Perform mass to charge scoring of two multiple mass traces
     *
     * Scores two mass traces based on the m/z and the hypothesis that one
     * trace is an isotopic trace of the other one. The isotopic position
     * (which trace it is) and the charge for the hypothesis are given as
     * additional parameters.
     * The scoring is described in Kenar et al., and is based on a random
     * sample of 115 000 compounds drawn from a comprehensive set of 24 million
     * putative sum formulas, of which the isotopic distribution was accurately
     * calculated. Thus, a theoretical mu and sigma are calculated as:
     *
     * mu = 1.000857 * j + 0.001091 u
     * sigma = 0.0016633 j * 0.0004751
     *
     * where j is the isotopic peak considered. A similarity score based on
     * agreement with the model is then computed.
     *
     * Reference: Kenar et al., doi: 10.1074/mcp.M113.031278
     *
     * An alternative scoring was added which test if isotope m/z distances lie in an expected m/z window.
     * This window is computed from a given set of elements.
     *
     * @param[in] mt1 First mass trace
     * @param[in] mt2 Second mass trace
     * @param[in] isotopic_position Isotopic position
     * @param[in] charge Charge
    */
    double scoreMZ_(const MassTrace& mt1, const MassTrace& mt2, Size isotopic_position, Size charge) const;

    /**
     * @brief score isotope m/z distance based on the expected m/z distances using C13-C12 or Kenar method
     * @param[in] iso_pos Isotopic position
     * @param[in] charge Charge
     * @param[in] diff_mz Mass-to-charge difference
     * @param[in] mt_variances Mass trace variances
     * @return Score value
     */
    double scoreMZByExpectedMean_(Size iso_pos, Size charge, const double diff_mz, double mt_variances) const;

    /**
     * @brief score isotope m/z distance based on an expected isotope window which was calculated from a set of expected elements
     * @param[in] charge Charge
     * @param[in] diff_mz Mass-to-charge difference
     * @param[in] mt_variances m/z variance between the two mass traces which are compared
     * @param[in] isotope_window Isotope window
     * @return Score value
     */
    double scoreMZByExpectedRange_(Size charge, const double diff_mz, double mt_variances, Range isotope_window) const;

    /** @brief Unified RT similarity scoring combining profile alignment, FWHM overlap,
     *  cosine similarity, Pearson correlation and normalised cross-correlation.
     *
     *  The function:
     *  1. Aligns the full elution profiles using a tolerance-based merge (0.1 s),
     *     preferring smoothed intensities when available.
     *  2. Checks that the FWHM windows of the two traces overlap by at least 70 %
     *     (adapted from FeatureFindingMetabo). Returns 0 immediately if they do not.
     *  3. Filters on Pearson correlation (rt_min_pearson_correlation_).
     *  4. Computes normalised cross-correlation (up to rt_max_lag_ scans of shift) on the aligned profiles.
     *
     * @param[in] mt1 First mass trace (monoisotopic)
     * @param[in] mt2 Second mass trace (candidate isotope)
     */
    std::pair<double, double> scoreRT_(const MassTrace& mt1, const MassTrace& mt2) const;

    /** @brief Perform intensity scoring using the averagine model (for peptides only)
     *
     * Compare the isotopic intensity distribution with the theoretical one
     * expected for peptides, using the averagine model. Compute the cosine
     * similarity between the two values.
     *
     * @param[in] intensities Intensity values
     * @param[in] molecular_weight Molecular weight
    */
    double computeAveragineSimScore_(const std::vector<double>& intensities, const double& molecular_weight) const;

    /** @brief Identify groupings of mass traces based on a set of reasonable candidates
     *
     * Takes a set of reasonable candidates for mass trace grouping and checks
     * all combinations of charge and isotopic positions on the candidates. It
     * is assumed that candidates[0] is the monoisotopic trace.
     *
     * The resulting possible groupings are appended to output_hypotheses.
     *
     * @param[in] candidates Candidate mass traces
     * @param[in] total_intensity Total intensity
     * @param[out] output_hypotheses Output feature hypotheses
    */
    void findLocalFeatures_(const std::vector<const MassTrace*>& candidates, std::vector<FeatureHypothesis>& output_hypotheses) const;

    /// parameter stuff
    double local_rt_range_;
    double local_im_range_;
    double local_mz_range_;
    Size charge_lower_bound_;
    Size charge_upper_bound_;
    double chrom_fwhm_;

    bool report_summed_ints_;
    bool enable_RT_filtering_;
    bool use_smoothed_intensities_;
    bool report_smoothed_intensities_;

    bool report_convex_hulls_;
    bool report_chromatograms_;

    bool remove_single_traces_;
    bool overlapping_features_;
    double hypothesis_score_quantile_;

    bool enable_mass_defect_filtering_;
    double mass_defect_offset_;
    Size minimum_isotopes_nr_;

    double rt_peak_overlap_threshold_;
    double rt_min_pearson_correlation_;
    int rt_max_lag_;
  };

}
