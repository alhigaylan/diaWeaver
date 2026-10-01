// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Mohammed Alhigaylan $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/FeatureFindingPeptide.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <OpenMS/CONCEPT/ParallelFor.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/OPENSWATHALGO/ALGO/Scoring.h>
#include <OpenMS/SYSTEM/File.h>

#include <boost/dynamic_bitset.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

// #define FFM_DEBUG

#include <numeric>
#include <optional>
#include <unordered_set>

namespace OpenMS
{
  FeatureFindingPeptide::FeatureFindingPeptide() :
    DefaultParamHandler("FeatureFindingPeptide"), ProgressLogger()
  {
    defaults_.setValue("local_rt_range", 5.0, "RT range where to look for coeluting mass traces", {"advanced"});
    defaults_.setValue("local_im_range", 0.02, "IM range where to look for coeluting mass traces", {"advanced"});
    defaults_.setValue("local_mz_range", 3.0, "MZ range where to look for isotopic mass traces", {"advanced"});
    defaults_.setValue("charge_lower_bound", 1, "Lowest charge state to consider");
    defaults_.setValue("charge_upper_bound", 4, "Highest charge state to consider");
    defaults_.setValue("chrom_fwhm", 5.0, "Expected chromatographic peak width (in seconds).");
    defaults_.setValue("report_summed_ints", "false", "Set to true for a feature intensity summed up over all traces rather than using monoisotopic trace intensity alone.", {"advanced"});
    defaults_.setValidStrings("report_summed_ints", {"false","true"});
    defaults_.setValue("enable_RT_filtering", "true", "Require sufficient overlap in RT while assembling mass traces. Disable for direct injection data..");
    defaults_.setValidStrings("enable_RT_filtering", {"false","true"});

    defaults_.setValue("use_smoothed_intensities", "true", "Use Savitzky-Golay smoothed intensities (produced by ElutionPeakDetection) instead of raw intensities.", {"advanced"});
    defaults_.setValidStrings("use_smoothed_intensities", {"false","true"});
    defaults_.setValue("report_smoothed_intensities", "true", "Report smoothed intensities (only if use_smoothed_intensities is true).", {"advanced"});
    defaults_.setValidStrings("report_smoothed_intensities", {"false","true"});

    defaults_.setValue("report_convex_hulls", "false", "Augment each reported feature with the convex hull of the underlying mass traces (increases featureXML file size considerably).");
    defaults_.setValidStrings("report_convex_hulls", {"false","true"});

    defaults_.setValue("report_chromatograms", "false", "Adds Chromatogram for each reported feature (Output in mzml).");
    defaults_.setValidStrings("report_chromatograms", {"false","true"});

    defaults_.setValue("remove_single_traces", "false", "Remove unassembled traces (single traces).");
    defaults_.setValidStrings("remove_single_traces", {"false","true"});

    defaults_.setValue("rt_peak_overlap_threshold", 0.3, "Minimum FWHM overlap proportion required between two mass traces to be considered co-eluting. Computed as the intersection of the two FWHM windows divided by the longer FWHM. Range [0, 1].");
    defaults_.setMinFloat("rt_peak_overlap_threshold", 0.0);
    defaults_.setMaxFloat("rt_peak_overlap_threshold", 1.0);

    defaults_.setValue("rt_min_pearson_correlation", 0.7, "Minimum Pearson correlation required between two mass trace elution profiles before cross-correlation is computed. Pairs below this threshold are rejected without the more expensive XCorr calculation.");
    defaults_.setMinFloat("rt_min_pearson_correlation", 0.0);
    defaults_.setMaxFloat("rt_min_pearson_correlation", 1.0);

    defaults_.setValue("rt_max_lag", 5, "Maximum lag (in number of scans) allowed when computing the normalized cross-correlation between two isotopic elution profiles. A value of 5 permits isotope traces shifted by up to 5 scans relative to the monoisotopic trace. Usually should match local_rt_range");
    defaults_.setMinInt("rt_max_lag", 0);

    defaults_.setValue("mass_defect_filtering", "true", "Filter feature hypotheses by peptide mass defect boundaries (adapted from DIA-Umpire). Rejects features whose mass defect falls outside the linear boundary defined for peptides.");
    defaults_.setValidStrings("mass_defect_filtering", {"false","true"});
    defaults_.setValue("mass_defect_offset", 0.1, "Mass defect tolerance offset (in Da, applied to the fractional mass) for the peptide mass defect filter. Increasing defect tolerance is recommended for modified peptides", {"advanced"});
    defaults_.setMinFloat("mass_defect_offset", 0.0);

    defaults_.setValue("minimum_isotopes_nr", 2, "Minimum number of isotopic mass traces required for a feature hypothesis to be reported. Must be at least 2 (monoisotopic + one isotope trace).");
    defaults_.setMinInt("minimum_isotopes_nr", 2);

    defaultsToParam_();

    this->setLogType(CMD);
  }

  FeatureFindingPeptide::~FeatureFindingPeptide() = default;

  void FeatureFindingPeptide::updateMembers_()
  {
    local_rt_range_ = (double)param_.getValue("local_rt_range");
    local_im_range_ = (double)param_.getValue("local_im_range");
    local_mz_range_ = (double)param_.getValue("local_mz_range");
    chrom_fwhm_ = (double)param_.getValue("chrom_fwhm");

    charge_lower_bound_ = (Size)param_.getValue("charge_lower_bound");
    charge_upper_bound_ = (Size)param_.getValue("charge_upper_bound");

    report_summed_ints_ = param_.getValue("report_summed_ints").toBool();
    enable_RT_filtering_ = param_.getValue("enable_RT_filtering").toBool();

    use_smoothed_intensities_ = param_.getValue("use_smoothed_intensities").toBool();
    bool use_smoothed = param_.getValue("use_smoothed_intensities").toBool();
    bool report_smoothed = param_.getValue("report_smoothed_intensities").toBool();
    if (report_smoothed && !use_smoothed) {
      OPENMS_LOG_WARN << "Warning: 'report_smoothed_intensities' is set to true, but 'use_smoothed_intensities' is false. Ignoring 'report_smoothed_intensities'.\n";
      report_smoothed = false;
    }
    use_smoothed_intensities_ = use_smoothed;
    report_smoothed_intensities_ = report_smoothed;

    report_convex_hulls_ = param_.getValue("report_convex_hulls").toBool();
    report_chromatograms_ = param_.getValue("report_chromatograms").toBool();

    remove_single_traces_ = param_.getValue("remove_single_traces").toBool();

    rt_peak_overlap_threshold_ = (double)param_.getValue("rt_peak_overlap_threshold");
    rt_min_pearson_correlation_ = (double)param_.getValue("rt_min_pearson_correlation");
    rt_max_lag_ = (int)param_.getValue("rt_max_lag");
    enable_mass_defect_filtering_ = param_.getValue("mass_defect_filtering").toBool();
    mass_defect_offset_ = (double)param_.getValue("mass_defect_offset");
    minimum_isotopes_nr_ = static_cast<Size>((int)param_.getValue("minimum_isotopes_nr"));
  }


  double FeatureFindingPeptide::computeAveragineSimScore_(const std::vector<double>& hypo_ints, const double& mol_weight) const
  {
    CoarseIsotopePatternGenerator solver(hypo_ints.size());
    auto isodist = solver.estimateFromPeptideMonoWeight(mol_weight);
    // isodist.renormalize();

    IsotopeDistribution::ContainerType averagine_dist = isodist.getContainer();
    double max_int(0.0), theo_max_int(0.0);
    for (Size i = 0; i < hypo_ints.size(); ++i)
    {
      if (hypo_ints[i] > max_int)
      {
        max_int = hypo_ints[i];
      }

      if (averagine_dist[i].getIntensity() > theo_max_int)
      {
        theo_max_int = averagine_dist[i].getIntensity();
      }
    }

    // compute normalized intensities
    std::vector<double> averagine_ratios, hypo_isos;
    for (Size i = 0; i < hypo_ints.size(); ++i)
    {
      averagine_ratios.push_back(averagine_dist[i].getIntensity() / theo_max_int);
      hypo_isos.push_back(hypo_ints[i] / max_int);
    }

    double iso_score = computeCosineSim_(averagine_ratios, hypo_isos);
    return iso_score;
  }

  double FeatureFindingPeptide::scoreMZ_(const MassTrace& tr1, const MassTrace& tr2, Size iso_pos, Size charge) const
  {
    double diff_mz(std::fabs(tr2.getCentroidMZ() - tr1.getCentroidMZ()));

    double mt_sigma1(tr1.getCentroidSD());
    double mt_sigma2(tr2.getCentroidSD());
    double mt_variances(std::exp(2 * std::log(mt_sigma1)) + std::exp(2 * std::log(mt_sigma2)));

    return scoreMZByExpectedMean_(iso_pos, charge, diff_mz, mt_variances);
  }

  double FeatureFindingPeptide::scoreMZByExpectedMean_(Size iso_pos, Size charge, const double diff_mz, double mt_variances) const
  {
    // Use C13-C12 mass difference as expected isotope spacing for peptides
    const double mu = (Constants::C13C12_MASSDIFF_U * iso_pos) / charge;
    const double sd = (0.0016633 * iso_pos - 0.0004751) / charge;

    double sigma_mult(3.0);
    double mz_score(0.0);

    //standard deviation including the estimated isotope deviation
    double score_sigma(std::sqrt(std::exp(2 * std::log(sd)) + mt_variances));

    // std::cout << std::setprecision(15) << "old " << score_sigma_old << " new " << score_sigma << '\n';

    if ((diff_mz < mu + sigma_mult * score_sigma) && (diff_mz > mu - sigma_mult * score_sigma))
    {
      double tmp_exponent((diff_mz - mu) / score_sigma);
      mz_score = std::exp(-0.5 * tmp_exponent * tmp_exponent);
    }
    return mz_score;
  }

  double FeatureFindingPeptide::scoreMZByExpectedRange_(Size charge, const double diff_mz, double mt_variances, Range isotope_window) const
  {
    //This isotope picking using m/z differences of elements' isotopes is based on the approach used in SIRIUS
    double sigma_mult(3.0);
    double mz_score(0.0);

    //standard deviation of m/z distance between the 2 mass traces
    double mt_sigma(std::sqrt(mt_variances));

    double max_allowed_deviation = mt_sigma * sigma_mult;

    double lbound = isotope_window.left_boundary / charge;
    double rbound = isotope_window.right_boundary / charge;

    if ((diff_mz < rbound) && (diff_mz > lbound))
    {
      //isotope masstrace lies in the expected range
      mz_score = 1.0;
    }
    else if ((diff_mz < rbound + max_allowed_deviation) && (diff_mz > lbound - max_allowed_deviation))
    {
      //score only the m/z difference which cannot explained by the elements m/z ranges
      double tmp_exponent;
      if (diff_mz < lbound)
      {
        tmp_exponent = (lbound - diff_mz) / mt_sigma;
      }
      else
      {
        tmp_exponent = (diff_mz - rbound) / mt_sigma;
      }
      mz_score = std::exp(-0.5 * tmp_exponent * tmp_exponent);
    }
    //else mz_score stays 0

    return mz_score;
  }

  std::pair<double, double> FeatureFindingPeptide::scoreRT_(const MassTrace& tr1, const MassTrace& tr2) const
  {
    if (!enable_RT_filtering_) return {1.0, 1.0};

    // Align elution profiles using smoothed intensities
    // Tolerance-based merge (mindiff = 0.1 s) pairs scans within 0.1 s of each
    // other and zero-pads non-overlapping ends. This avoids the exact-RT-equality
    // assumption that FeatureFindingMetabo relied on.
    const std::vector<double>& sm1 = tr1.getSmoothedIntensities();
    const std::vector<double>& sm2 = tr2.getSmoothedIntensities();
    const bool use_sm1 = !sm1.empty();
    const bool use_sm2 = !sm2.empty();

    const double mindiff = 0.1;
    std::vector<double> vec1, vec2;

    Size i = 0, j = 0;
    const Size n1 = tr1.getSize(), n2 = tr2.getSize();
    while (i < n1 && j < n2)
    {
      const double rt1 = tr1[i].getRT();
      const double rt2 = tr2[j].getRT();
      if (std::fabs(rt1 - rt2) < mindiff)
      {
        vec1.push_back(use_sm1 ? sm1[i] : tr1[i].getIntensity());
        vec2.push_back(use_sm2 ? sm2[j] : tr2[j].getIntensity());
        ++i; ++j;
      }
      else if (rt1 < rt2)
      {
        vec1.push_back(use_sm1 ? sm1[i] : tr1[i].getIntensity());
        vec2.push_back(0.0);
        ++i;
      }
      else
      {
        vec1.push_back(0.0);
        vec2.push_back(use_sm2 ? sm2[j] : tr2[j].getIntensity());
        ++j;
      }
    }
    while (i < n1) { vec1.push_back(use_sm1 ? sm1[i] : tr1[i].getIntensity()); vec2.push_back(0.0); ++i; }
    while (j < n2) { vec1.push_back(0.0); vec2.push_back(use_sm2 ? sm2[j] : tr2[j].getIntensity()); ++j; }

    // FWHM overlap check (adapted from FeatureFindingMetabo)
    // Reject pairs whose FWHM windows share less than 70 % of the longer FWHM.
    const std::pair<Size, Size> fwhm1 = tr1.getFWHMborders();
    const std::pair<Size, Size> fwhm2 = tr2.getFWHMborders();

    const double fwhm1_start = tr1[fwhm1.first].getRT();
    const double fwhm1_end   = tr1[fwhm1.second].getRT();
    const double fwhm2_start = tr2[fwhm2.first].getRT();
    const double fwhm2_end   = tr2[fwhm2.second].getRT();

    const double overlap    = std::max(0.0, std::min(fwhm1_end, fwhm2_end) - std::max(fwhm1_start, fwhm2_start));
    const double max_length = std::max(fwhm1_end - fwhm1_start, fwhm2_end - fwhm2_start);

    if (max_length == 0.0) return {0.0, 0.0};

    const double proportion = overlap / max_length;

    if (proportion < rt_peak_overlap_threshold_)
    {
      return {0.0, 0.0};
    }

    // Pearson correlation
    const double pearson = Math::pearsonCorrelationCoefficient(
        vec1.begin(), vec1.end(), vec2.begin(), vec2.end());
    if (pearson < rt_min_pearson_correlation_) return {0.0, proportion};

    // Normalised cross-correlation: evaluates shape similarity at multiple lags,
    // subsumes cosine similarity (which is NCC at lag=0).
    OpenSwath::Scoring::XCorrArrayType xcorr =
        OpenSwath::Scoring::normalizedCrossCorrelation(vec1, vec2, rt_max_lag_, 1);
    const double xcorr_score = OpenSwath::Scoring::xcorrArrayGetMaxPeak(xcorr)->second;

    return {xcorr_score, proportion};
  }

  Range FeatureFindingPeptide::getTheoreticIsotopicMassWindow_(const std::vector<Element const *>& alphabet, int peakOffset) const
  {
    if (peakOffset < 1)
    {
      throw std::invalid_argument("Expect a peak offset of at least 1");
    }
    double minmz = std::numeric_limits<double>::infinity();
    double maxmz = -std::numeric_limits<double>::infinity();

    for (const Element* e : alphabet) {
      IsotopeDistribution iso = e->getIsotopeDistribution();
      for (unsigned int k = 1; k < iso.size(); ++k) {
        const double mz_mono = iso[0].getMZ();
        const double mz_iso = iso[k].getMZ();

        const int integer_mz_mono =  (int)round(mz_mono);
        const int integer_mz_iso =  (int)round(mz_iso);
        const int i = integer_mz_iso - integer_mz_mono;

        if (i > peakOffset) break;
        const double mz_diff_iso_mono = mz_iso - mz_mono;
        double diff = mz_diff_iso_mono - i;
        diff *= (peakOffset / i);
        minmz = std::min(minmz, diff);
        maxmz = std::max(maxmz, diff);
      }
    }

    Range range = Range();
    range.left_boundary = peakOffset + minmz;
    range.right_boundary = peakOffset + maxmz;
    return range;
  }

  double FeatureFindingPeptide::computeCosineSim_(const std::vector<double>& x, const std::vector<double>& y) const
  {
    if (x.size() != y.size())
    {
      return 0.0;
    }

    double mixed_sum(0.0);
    double x_squared_sum(0.0);
    double y_squared_sum(0.0);

    for (Size i = 0; i < x.size(); ++i)
    {
      mixed_sum += x[i] * y[i];
      x_squared_sum += x[i] * x[i];
      y_squared_sum += y[i] * y[i];
    }

    double denom(std::sqrt(x_squared_sum) * std::sqrt(y_squared_sum));
    return (denom > 0.0) ? mixed_sum / denom : 0.0;
  }


  bool FeatureFindingPeptide::isMassDefectValid_(double neutral_mass, double d) const
  {
    // Fractional part (mass defect) helper: frac(x) = x - floor(x)
    auto massDefect = [](double x) { return x - std::floor(x); };

    // Linear boundaries fitted to the peptide mass defect filter (from DIA-Umpire; Tsou et al.): https://doi.org/10.1002/pmic.201500526
    // which was adapted from Toumi et al. https://doi.org/10.1021/pr100291q
    const double u = massDefect(0.00052738 * neutral_mass + 0.066015  + d);
    const double l = massDefect(0.00042565 * neutral_mass + 0.00038210 - d);
    const double defect = massDefect(neutral_mass);

    if (u > l)
    {
      // Normal (non-wrapping) interval
      return defect >= l && defect <= u;
    }
    // filter wraps across the 0/1 boundary (e.g. l=0.9, u=0.1)
    return defect >= l || defect <= u;
  }

  void FeatureFindingPeptide::findLocalFeatures_(const std::vector<const MassTrace*>& candidates, std::vector<FeatureHypothesis>& output_hypotheses) const
  {
    // single Mass trace hypothesis — no pair score available, quality is 0
    FeatureHypothesis tmp_hypo;
    tmp_hypo.addMassTrace(*candidates[0]);
    tmp_hypo.setScore(0.0);

    output_hypotheses.push_back(tmp_hypo);

    // scoreRT_ only depends on the trace pair (not on charge or isotope position) and is costly,
    // so compute it at most once per candidate
    std::vector<std::optional<std::pair<double, double>>> rt_scores(candidates.size());

    for (Size charge = charge_lower_bound_; charge <= charge_upper_bound_; ++charge)
    {
      // Reject this charge hypothesis immediately if the neutral mass falls outside
      // the expected peptide mass defect filter
      if (enable_mass_defect_filtering_)
      {
        const double neutral_mass = charge * (candidates[0]->getCentroidMZ() - Constants::PROTON_MASS_U);
        if (!isMassDefectValid_(neutral_mass, mass_defect_offset_))
        {
          continue;
        }
      }

      FeatureHypothesis fh_tmp;
      fh_tmp.addMassTrace(*candidates[0]);

      Size last_iso_idx(0);
      Size iso_pos_max(static_cast<Size>(std::floor(charge * local_mz_range_)));

      // Accumulators for mean individual scores across all accepted iso_pos pairs.
      double acc_rt(0.0), acc_mz(0.0), acc_int(0.0), acc_overlap(0.0), acc_pair(0.0);
      Size acc_count(0);

      for (Size iso_pos = 1; iso_pos <= iso_pos_max; ++iso_pos)
      {
        // Find mass trace that best agrees with current hypothesis of charge and isotopic position
        // Add a new score: trace count in a given hypothesis.
        // Quadratic trace count score: rewards larger hypotheses. keeping range [0,1].
        // Ideally, this will be biased against +1 charged peptides and small hypothesis (~2 isotopes)
        double best_so_far(0.0);
        Size best_idx(0);
        double best_rt_score(0.0), best_mz_score(0.0), best_int_score(0.0), best_overlap_score(0.0);
        for (Size mt_idx = last_iso_idx + 1; mt_idx < candidates.size(); ++mt_idx)
        {
#ifdef FFM_DEBUG
          std::cout << "scoring " << candidates[0]->getLabel() << " " << candidates[0]->getCentroidMZ() <<
            " with " << candidates[mt_idx]->getLabel() << " " << candidates[mt_idx]->getCentroidMZ() << '\n';
#endif
          double mz_score(scoreMZ_(*candidates[0], *candidates[mt_idx], iso_pos, charge));
          // a pair with mz_score 0 has total_pair_score 0, so it can never become the best match
          if (mz_score <= 0.0) continue;

          if (!rt_scores[mt_idx]) rt_scores[mt_idx] = scoreRT_(*candidates[0], *candidates[mt_idx]);
          const auto [rt_score, overlap_score] = *rt_scores[mt_idx];

          // initialize int score
          double int_score(0.0);
          // double int_score((candidates[0]->getIntensity(use_smoothed_intensities_))/total_weight + (candidates[mt_idx]->getIntensity(use_smoothed_intensities_))/total_weight);
          // the pair score is only used if rt_score and mz_score are positive, so skip the costly isotope pattern otherwise
          if (rt_score > 0.0 && mz_score > 0.0)
          {
            std::vector<double> tmp_ints(fh_tmp.getAllIntensities());
            tmp_ints.push_back(candidates[mt_idx]->getIntensity(use_smoothed_intensities_));
            int_score = computeAveragineSimScore_(tmp_ints, charge * (candidates[mt_idx]->getCentroidMZ() - Constants::PROTON_MASS_U));
          }

#ifdef FFM_DEBUG
          std::cout << fh_tmp.getLabel() << "_" << candidates[mt_idx]->getLabel() <<
            "\t" << "ch: " << charge << " isotope: " << iso_pos << " rt: " <<
            rt_score << "mz: " << mz_score << "int: " << int_score << '\n';
#endif

          double total_pair_score(0.0);
          if (rt_score > 0.0 && mz_score > 0.0 && int_score > 0.0)
          {
            const double trace_count_score = (fh_tmp.getSize() + 1 > 2) ? 1.0 : 0.0;
            total_pair_score = (rt_score + mz_score + int_score + trace_count_score) / 4.0;
          }
          if (total_pair_score > best_so_far)
          {
            best_so_far = total_pair_score;
            best_idx = mt_idx;
            best_rt_score = rt_score;
            best_mz_score = mz_score;
            best_int_score = int_score;
            best_overlap_score = overlap_score;
          }
        } // end mt_idx

        // Store mass trace that best agrees with current hypothesis of charge
        // and isotopic position
        if (best_so_far > 0.0)
        {
          ++acc_count;
          acc_rt      += best_rt_score;
          acc_mz      += best_mz_score;
          acc_int     += best_int_score;
          acc_overlap += best_overlap_score;
          acc_pair    += best_so_far;

          fh_tmp.addMassTrace(*candidates[best_idx]);
          fh_tmp.setScore(acc_pair / acc_count);
          fh_tmp.setScoreRT(acc_rt / acc_count);
          fh_tmp.setScoreMZ(acc_mz / acc_count);
          fh_tmp.setScoreInt(acc_int / acc_count);
          fh_tmp.setScoreOverlap(acc_overlap / acc_count);
          fh_tmp.setCharge(charge);
          last_iso_idx = best_idx;

          output_hypotheses.push_back(fh_tmp);
        }
        else
        {
          break;
        }
      } // end for iso_pos

#ifdef FFM_DEBUG
      std::cout << "best found for ch " << charge << ":" << fh_tmp.getLabel() << " score: " << fh_tmp.getScore() << '\n';
#endif
    } // end for charge
  } // end of findLocalFeatures_(...)

  void FeatureFindingPeptide::run(std::vector<MassTrace>& input_mtraces, FeatureMap& output_featmap, std::vector<std::vector< OpenMS::MSChromatogram > >& output_chromatograms, Size num_threads)
  {

    output_featmap.clear();
    output_chromatograms.clear();

    if (input_mtraces.empty())
    {
      return;
    }

    // mass traces must be sorted by their centroid MZ
    std::sort(input_mtraces.begin(), input_mtraces.end(), CmpMassTraceByMZ());

    this->startProgress(0, input_mtraces.size(), "assembling mass traces to features");

    // *********************************************************** //
    // Step 2 Iterate through all mass traces to find likely matches
    // and generate isotopic / charge hypotheses
    // *********************************************************** //

    // One bucket per input trace, so each worker only ever writes its own bucket (no
    // shared-vector contention) and the buckets can be concatenated afterward in a fixed,
    // deterministic order -- unlike a shared vector filled via mutex-guarded push_back,
    // whose order otherwise follows thread completion timing.
    std::vector<std::vector<FeatureHypothesis>> hypos_per_trace(input_mtraces.size());

    auto process_trace = [&](Size i)
    {
      std::vector<const MassTrace*> local_traces;
      double ref_trace_mz(input_mtraces[i].getCentroidMZ());
      double ref_trace_rt(input_mtraces[i].getCentroidRT());
      double ref_trace_im(input_mtraces[i].getCentroidIM());

      local_traces.push_back(&input_mtraces[i]);

      for (Size ext_idx = i + 1; ext_idx < input_mtraces.size(); ++ext_idx)
      {
        // traces are sorted by m/z, so we can break when we leave the allowed window
        double diff_mz = std::fabs(input_mtraces[ext_idx].getCentroidMZ() - ref_trace_mz);
        if (diff_mz > local_mz_range_)
        {
          break;
        }
        double diff_rt = std::fabs(input_mtraces[ext_idx].getCentroidRT() - ref_trace_rt);
        double diff_im = std::fabs(input_mtraces[ext_idx].getCentroidIM() - ref_trace_im);
        if (diff_rt <= local_rt_range_ && diff_im < local_im_range_)
        {
          // std::cout << " accepted!\n";
          local_traces.push_back(&input_mtraces[ext_idx]);
        }
      }
      findLocalFeatures_(local_traces, hypos_per_trace[i]);
    };

    if (num_threads > 0)
    {
      // Explicit thread budget: this is a caller (diaWeaver) running inside its own
      // already-active nested parallelism, where relying on a nested OpenMP region's
      // implicit barrier proved unreliable (see diaWeaver.cpp's nested-parallelism data
      // race investigation). Use std::thread workers instead, whose join() is a
      // plain-C++-guaranteed synchronization point independent of the OpenMP runtime.
      OpenMS::parallelFor(input_mtraces.size(), num_threads, process_trace);
    }
    else
    {
      Size progress(0);
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (SignedSize i = 0; i < (SignedSize)input_mtraces.size(); ++i)
      {
        IF_MASTERTHREAD this->setProgress(progress);
#ifdef _OPENMP
#pragma omp atomic
#endif
        ++progress;

        process_trace(static_cast<Size>(i));
      }
    }
    this->endProgress();

    // Concatenate in input-trace order (fixed regardless of thread completion order)
    std::vector<FeatureHypothesis> feat_hypos;
    for (auto& local_hypos : hypos_per_trace)
    {
      for (auto& fh : local_hypos) feat_hypos.push_back(std::move(fh));
    }

    // sort feature candidates by their score (descending)
    std::sort(feat_hypos.begin(), feat_hypos.end(), CmpHypothesesByScore());

    // Remove hypotheses that don't meet the minimum isotope trace count.
    // Single-trace (charge-0) hypotheses are exempted when remove_single_traces_ is false:
    // minimum_isotopes_nr_ is a quality gate for assembled multi-isotope features, not for
    // unassembled single traces whose retention is controlled by remove_single_traces_.
    feat_hypos.erase(
      std::remove_if(feat_hypos.begin(), feat_hypos.end(),
        [this](const FeatureHypothesis& fh) {
          if (!remove_single_traces_ && fh.getCharge() == 0) return false;
          return fh.getSize() < minimum_isotopes_nr_;
        }),
      feat_hypos.end());

#ifdef FFM_DEBUG
    std::cout << "size of hypotheses: " << feat_hypos.size() << '\n';
    // output all hypotheses:
    for (Size hypo_idx = 0; hypo_idx < feat_hypos.size(); ++ hypo_idx)
    {
      std::cout << feat_hypos[hypo_idx].getLabel() << " ch: " << feat_hypos[hypo_idx].getCharge() <<
        " score: " << feat_hypos[hypo_idx].getScore() << '\n';
    }
#endif

    // *********************************************************** //
    // Step 3 Iterate through all hypotheses, starting with the highest
    // scoring one. Accept them if they do not contain traces that have
    // already been used by a higher scoring hypothesis.
    // *********************************************************** //

    // A trace claimed by an already-accepted hypothesis may not be reused by a later,
    // lower-scoring one. Traces are identified by address (they all live in input_mtraces, which
    // is not modified from here on); this matches identifying them by label, as trace labels from
    // MassTraceDetection/ElutionPeakDetection are unique.
    std::unordered_set<const MassTrace*> claimed_traces;

    for (Size hypo_idx = 0; hypo_idx < feat_hypos.size(); ++hypo_idx)
    {
      const std::vector<const MassTrace*>& traces = feat_hypos[hypo_idx].getMassTraces();

      bool collision = false;
      for (const MassTrace* trace : traces)
      {
        if (claimed_traces.count(trace))
        {
          collision = true;
          break;
        }
      }
      if (collision) continue;

      for (const MassTrace* trace : traces) claimed_traces.insert(trace);

      // filter out single traces if option is set
      if (remove_single_traces_ && feat_hypos[hypo_idx].getCharge() == 0)
      {
        continue;
      }

      //
      // Now accept hypothesis
      //

      Feature f;
      f.setRT(feat_hypos[hypo_idx].getCentroidRT());
      f.setMZ(feat_hypos[hypo_idx].getCentroidMZ());

      if (report_summed_ints_)
      {
        // f.setIntensity(feat_hypos[hypo_idx].getSummedFeatureIntensity(report_smoothed_intensities_));
        f.setIntensity(feat_hypos[hypo_idx].getSummedFeatureIntensity(use_smoothed_intensities_));
      }
      else
      {
        //f.setIntensity(feat_hypos[hypo_idx].getMonoisotopicFeatureIntensity(report_smoothed_intensities_));
        f.setIntensity(feat_hypos[hypo_idx].getMonoisotopicFeatureIntensity(use_smoothed_intensities_));
      }

      f.setWidth(feat_hypos[hypo_idx].getFWHM());
      f.setCharge(feat_hypos[hypo_idx].getCharge());
      f.setMetaValue(3, feat_hypos[hypo_idx].getLabel());
      //f.setMetaValue("max_height", feat_hypos[hypo_idx].getMaxIntensity(report_smoothed_intensities_));
      f.setMetaValue("max_height", feat_hypos[hypo_idx].getMaxIntensity(use_smoothed_intensities_));

      // store isotope intensities
      //std::vector<double> all_ints(feat_hypos[hypo_idx].getAllIntensities(report_smoothed_intensities_));
      std::vector<double> all_ints(feat_hypos[hypo_idx].getAllIntensities(use_smoothed_intensities_));
      f.setMetaValue(Constants::UserParam::NUM_OF_MASSTRACES, all_ints.size());
      if (report_convex_hulls_) f.setConvexHulls(feat_hypos[hypo_idx].getConvexHulls());
      f.setOverallQuality(feat_hypos[hypo_idx].getScore());
      f.setMetaValue("score_rt", feat_hypos[hypo_idx].getScoreRT());
      f.setMetaValue("score_mz", feat_hypos[hypo_idx].getScoreMZ());
      f.setMetaValue("score_int", feat_hypos[hypo_idx].getScoreInt());
      f.setMetaValue("score_overlap", feat_hypos[hypo_idx].getScoreOverlap());
      f.setMetaValue("masstrace_intensity", all_ints);
      f.setMetaValue("masstrace_centroid_rt", feat_hypos[hypo_idx].getAllCentroidRT());
      f.setMetaValue("masstrace_centroid_mz", feat_hypos[hypo_idx].getAllCentroidMZ());
      f.setMetaValue("masstrace_centroid_im", feat_hypos[hypo_idx].getAllCentroidIM());
      f.setMetaValue("isotope_distances", feat_hypos[hypo_idx].getIsotopeDistances());
      f.applyMemberFunction(&UniqueIdInterface::setUniqueId);
      output_featmap.push_back(std::move(f));
      const Feature& added = output_featmap.back();

      if (report_chromatograms_ && added.getIntensity() != 0)
      {
        output_chromatograms.push_back(feat_hypos[hypo_idx].getChromatograms(added.getUniqueId()));
      }


    }
    output_featmap.setUniqueId(UniqueIdGenerator::getUniqueId());

    // Sort by m/z, same as FeatureMap::sortByMZ(), but with an explicit, deterministic
    // tie-break (ion mobility, then intensity, then label) for features that share the
    // same m/z -- typically isobaric, ion-mobility-distinct species. Plain sortByMZ()'s
    // std::sort has no tie-break, so its result for m/z-tied features depends on
    // whatever order they arrived in beforehand, which can vary run to run when
    // feat_hypos above was built by concurrent worker threads.
    // The tie-break keys are read from meta values, so compute them once per feature instead of
    // in every comparison; sort feature indices by the keys, then reorder the features.
    struct SortKey
    {
      double mz;
      double im;
      double intensity;
      String label;
    };
    std::vector<SortKey> keys;
    keys.reserve(output_featmap.size());
    for (const Feature& f : output_featmap)
    {
      std::vector<double> ims = f.getMetaValue("masstrace_centroid_im");
      keys.push_back({f.getMZ(), ims.empty() ? 0.0 : ims[0], f.getIntensity(), String(f.getMetaValue("label"))});
    }

    std::vector<Size> order(output_featmap.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
      [&keys](Size ia, Size ib)
      {
        const SortKey& a = keys[ia];
        const SortKey& b = keys[ib];
        if (a.mz != b.mz) return a.mz < b.mz;
        if (a.im != b.im) return a.im < b.im;
        if (a.intensity != b.intensity) return a.intensity < b.intensity;
        return a.label < b.label;
      });

    std::vector<Feature> sorted_features;
    sorted_features.reserve(output_featmap.size());
    for (Size idx : order)
    {
      sorted_features.push_back(std::move(output_featmap[idx]));
    }
    for (Size i = 0; i < sorted_features.size(); ++i)
    {
      output_featmap[i] = std::move(sorted_features[i]);
    }
  } // end of FeatureFindingPeptide::run

}
