// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Mohammed Alhigaylan $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/FeatureFindingPeptide.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/OPENSWATHALGO/ALGO/Scoring.h>
#include <OpenMS/SYSTEM/File.h>

#include <boost/dynamic_bitset.hpp>

// #define FFM_DEBUG

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

    defaults_.setValue("use_smoothed_intensities", "true", "Use LOWESS intensities instead of raw intensities.", {"advanced"});
    defaults_.setValidStrings("use_smoothed_intensities", {"false","true"});
    defaults_.setValue("report_smoothed_intensities", "true", "Report smoothed intensities (only if use_smoothed_intensities is true).", {"advanced"});
    defaults_.setValidStrings("report_smoothed_intensities", {"false","true"});

    defaults_.setValue("report_convex_hulls", "false", "Augment each reported feature with the convex hull of the underlying mass traces (increases featureXML file size considerably).");
    defaults_.setValidStrings("report_convex_hulls", {"false","true"});

    defaults_.setValue("report_chromatograms", "false", "Adds Chromatogram for each reported feature (Output in mzml).");
    defaults_.setValidStrings("report_chromatograms", {"false","true"});

    defaults_.setValue("remove_single_traces", "false", "Remove unassembled traces (single traces).");
    defaults_.setValidStrings("remove_single_traces", {"false","true"});

    defaults_.setValue("overlapping_features", "false", "Allow mass traces to be reused to explain lower scoring features. Recommended for peptides.");
    defaults_.setValidStrings("overlapping_features", {"false","true"});

    defaults_.setValue("rt_min_pearson_correlation", 0.7, "Minimum Pearson correlation required between two mass trace elution profiles before cross-correlation is computed. Pairs below this threshold are rejected without the more expensive XCorr calculation.");
    defaults_.setMinFloat("rt_min_pearson_correlation", 0.0);
    defaults_.setMaxFloat("rt_min_pearson_correlation", 1.0);

    defaults_.setValue("rt_max_lag", 5, "Maximum lag (in number of scans) allowed when computing the normalized cross-correlation between two isotopic elution profiles. A value of 5 permits isotope traces shifted by up to 5 scans relative to the monoisotopic trace. Usually should match local_rt_range");
    defaults_.setMinInt("rt_max_lag", 0);

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

    rt_min_pearson_correlation_ = (double)param_.getValue("rt_min_pearson_correlation");
    rt_max_lag_ = (int)param_.getValue("rt_max_lag");
    rt_max_lag_ = (int)param_.getValue("rt_max_lag");
    overlapping_features_ = param_.getValue("overlapping_features").toBool();
  }


  double FeatureFindingPeptide::computeAveragineSimScore_(const std::vector<double>& hypo_ints, const double& mol_weight) const
  {
    CoarseIsotopePatternGenerator solver(hypo_ints.size());
    auto isodist = solver.estimateFromPeptideWeight(mol_weight);
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

  double FeatureFindingPeptide::scoreRT_(const MassTrace& tr1, const MassTrace& tr2) const
  {
    // return success if this filter is disabled
    if (!enable_RT_filtering_) return 1.0;

    // Use Savitzky-Golay smoothed intensities from ElutionPeakDetection if
    // available; fall back to raw intensities if smoothing was not run
    // (e.g. when EPD is disabled).
    const std::vector<double>& sm1 = tr1.getSmoothedIntensities();
    const std::vector<double>& sm2 = tr2.getSmoothedIntensities();
    const bool use_sm1 = !sm1.empty();
    const bool use_sm2 = !sm2.empty();

    // Build RT-aligned intensity vectors using a tolerance-based merge
    // (handles missing data points by filling in with zeros).
    // Non-overlapping ends are zero-padded.
    // mindiff is hard-coded as 0.1 s. In FFMetabo, exact double equality
    // was used since we work on the same map.
    // But this is more robust for same-map traces and matches OpenSWATH/diaWeaver behavior.
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

    // compute pearson correlation. If a correlation is below a cutoff,
    // quit calculation and filter this pair as unlikely.
    const double pearson = Math::pearsonCorrelationCoefficient(
        vec1.begin(), vec1.end(), vec2.begin(), vec2.end());
    if (pearson < rt_min_pearson_correlation_) return 0.0;

    // normalized cross correlation to identify traces that are slightly shifted due to noise
    // maxdelay --> how many offsets (in units of scans) to consider. Usually, this should match local_rt_range
    // lag --> I believe this is the step size when considering max lag. Hard-coding this as 1 means it will
    // evaluate every possible scan offset within +/- maxdelay.
    // Since this function is borrowed from OpenSWATH, lag parameter may be helpful in XIC with 100s of data points.
    // but here, mass traces are roughly around 15 ~ 30 points.
    OpenSwath::Scoring::XCorrArrayType xcorr =
        OpenSwath::Scoring::normalizedCrossCorrelation(vec1, vec2, rt_max_lag_, 1);
    OpenSwath::Scoring::XCorrArrayType::const_iterator best =
        OpenSwath::Scoring::xcorrArrayGetMaxPeak(xcorr);

    return best->second;
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


  void FeatureFindingPeptide::findLocalFeatures_(const std::vector<const MassTrace*>& candidates, const double total_intensity, std::vector<FeatureHypothesis>& output_hypotheses) const
  {
    // single Mass trace hypothesis
    FeatureHypothesis tmp_hypo;
    tmp_hypo.addMassTrace(*candidates[0]);
    tmp_hypo.setScore((candidates[0]->getIntensity(use_smoothed_intensities_)) / total_intensity);

#ifdef _OPENMP
#pragma omp critical (OPENMS_FFMetabo_output_hypos)
#endif
    {
      // pushing back to shared vector needs to be synchronized
      output_hypotheses.push_back(tmp_hypo);
    }

    for (Size charge = charge_lower_bound_; charge <= charge_upper_bound_; ++charge)
    {
      FeatureHypothesis fh_tmp;
      fh_tmp.addMassTrace(*candidates[0]);
      fh_tmp.setScore((candidates[0]->getIntensity(use_smoothed_intensities_)) / total_intensity);

      // double mono_iso_rt(candidates[0]->getCentroidRT());
      // double mono_iso_mz(candidates[0]->getCentroidMZ());
      // double mono_iso_int(candidates[0]->computePeakArea());

      Size last_iso_idx(0);
      Size iso_pos_max(static_cast<Size>(std::floor(charge * local_mz_range_)));
      for (Size iso_pos = 1; iso_pos <= iso_pos_max; ++iso_pos)
      {

        // Constants
        //const double neutron_mass = 1.0033548378;  // Da
        //const double iso_ppm_tolerance = 20.0; // ppm tolerance

        // Find mass trace that best agrees with current hypothesis of charge
        // and isotopic position
        double best_so_far(0.0);
        Size best_idx(0);
        for (Size mt_idx = last_iso_idx + 1; mt_idx < candidates.size(); ++mt_idx)
        {

          /*
          // before computing any of the scores -- implement a strict mz pattern cutoff
          double mono_mz = candidates[0]->getCentroidMZ();
          double cand_mz = candidates[mt_idx]->getCentroidMZ();

          // Compute expected theoretical m/z shift for this isotope position and charge
          double expected_shift = (neutron_mass * iso_pos) / charge;

          // Actual observed m/z shift
          double observed_shift = std::fabs(cand_mz - mono_mz);
          double ppm_error = std::fabs(observed_shift - expected_shift) / expected_shift * 1e6;
          // skip loop if isotope position is outside of expected range
          if (ppm_error > iso_ppm_tolerance) continue;
           */

          // double tmp_iso_rt(candidates[mt_idx]->getCentroidRT());
          // double tmp_iso_mz(candidates[mt_idx]->getCentroidMZ());
          // double tmp_iso_int(candidates[mt_idx]->computePeakArea());

#ifdef FFM_DEBUG
          std::cout << "scoring " << candidates[0]->getLabel() << " " << candidates[0]->getCentroidMZ() <<
            " with " << candidates[mt_idx]->getLabel() << " " << candidates[mt_idx]->getCentroidMZ() << '\n';
#endif

          // Score current mass trace candidates against hypothesis
          // currently, if pearson correlation score is below rt_min_pearson_correlation,
          // scoreRT_ will return 0.
          double rt_score(scoreRT_(*candidates[0], *candidates[mt_idx]));
          double mz_score(scoreMZ_(*candidates[0], *candidates[mt_idx], iso_pos, charge));

          // disable intensity scoring for now...
          double int_score(1.0);
          // double int_score((candidates[0]->getIntensity(use_smoothed_intensities_))/total_weight + (candidates[mt_idx]->getIntensity(use_smoothed_intensities_))/total_weight);

          {
            std::vector<double> tmp_ints(fh_tmp.getAllIntensities());
            tmp_ints.push_back(candidates[mt_idx]->getIntensity(use_smoothed_intensities_));
            int_score = computeAveragineSimScore_(tmp_ints, candidates[mt_idx]->getCentroidMZ() * charge);
          }

#ifdef FFM_DEBUG
          std::cout << fh_tmp.getLabel() << "_" << candidates[mt_idx]->getLabel() <<
            "\t" << "ch: " << charge << " isopos: " << iso_pos << " rt: " <<
            rt_score << "mz: " << mz_score << "int: " << int_score << '\n';
#endif

          double total_pair_score(0.0);
          if (rt_score > 0.0 && mz_score > 0.0 && int_score > 0.0)
          {
            total_pair_score = std::exp(std::log(rt_score) + log(mz_score) + log(int_score));
          }
          if (total_pair_score > best_so_far)
          {
            best_so_far = total_pair_score;
            best_idx = mt_idx;
          }
        } // end mt_idx

        // Store mass trace that best agrees with current hypothesis of charge
        // and isotopic position
        if (best_so_far > 0.0)
        {
          fh_tmp.addMassTrace(*candidates[best_idx]);
          double weighted_score(((candidates[best_idx]->getIntensity(use_smoothed_intensities_)) * best_so_far) / total_intensity);

          fh_tmp.setScore(fh_tmp.getScore() + weighted_score);
          fh_tmp.setCharge(charge);
          last_iso_idx = best_idx;

#ifdef _OPENMP
#pragma omp critical (OPENMS_FFMetabo_output_hypos)
#endif
          {
            // pushing back to shared vector needs to be synchronized
            output_hypotheses.push_back(fh_tmp);
          }
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

  void FeatureFindingPeptide::run(std::vector<MassTrace>& input_mtraces, FeatureMap& output_featmap, std::vector<std::vector< OpenMS::MSChromatogram > >& output_chromatograms)
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

    double total_intensity(0.0);
    for (Size i = 0; i < input_mtraces.size(); ++i)
    {
      total_intensity += input_mtraces[i].getIntensity(use_smoothed_intensities_);
    }

    // *********************************************************** //
    // Step 2 Iterate through all mass traces to find likely matches
    // and generate isotopic / charge hypotheses
    // *********************************************************** //

    std::vector<FeatureHypothesis> feat_hypos;
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
      findLocalFeatures_(local_traces, total_intensity, feat_hypos);
    }
    this->endProgress();

    // sort feature candidates by their score
    std::sort(feat_hypos.begin(), feat_hypos.end(), CmpHypothesesByScore());

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

    // ---------- allow trace collision if it uses different isotopes / charge ----
    std::multimap<String, int> trace_excl_map;

    for (Size hypo_idx = 0; hypo_idx < feat_hypos.size(); ++hypo_idx)
    {
      const std::vector<String>& labels = feat_hypos[hypo_idx].getLabels();
      int current_charge = feat_hypos[hypo_idx].getCharge();

      bool trace_coll = false;

      if (overlapping_features_)
      {
        // Only check the *first* label, but allow multiple charges
        const String& first_label = labels[0];

        auto range = trace_excl_map.equal_range(first_label);
        for (auto it = range.first; it != range.second; ++it)
        {
          if (it->second == current_charge)
          {
            // Found same label + same charge → collision
            trace_coll = true;
            break;
          }
        }
      }
      else
      {
        // Old behavior: check all labels regardless of charge
        for (Size lab_idx = 0; lab_idx < labels.size(); ++lab_idx)
        {
          auto range = trace_excl_map.equal_range(labels[lab_idx]);
          if (range.first != range.second)  // label already used
          {
            trace_coll = true;
            break;
          }
        }
      }

      // skip if collision
      if (trace_coll)
      {
        continue;
      }

      // Accept hypothesis → mark its traces
      if (overlapping_features_)
      {
        // mark only the *first* label + charge
        trace_excl_map.emplace(labels[0], current_charge);
      }
      else
      {
        // mark all labels (but still allow multiple charges)
        for (const auto& lab : labels)
        {
          trace_excl_map.emplace(lab, current_charge);
        }
      }

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
      f.setMetaValue("masstrace_intensity", all_ints);
      f.setMetaValue("masstrace_centroid_rt", feat_hypos[hypo_idx].getAllCentroidRT());
      f.setMetaValue("masstrace_centroid_mz", feat_hypos[hypo_idx].getAllCentroidMZ());
      f.setMetaValue("masstrace_centroid_im", feat_hypos[hypo_idx].getAllCentroidIM());
      f.setMetaValue("isotope_distances", feat_hypos[hypo_idx].getIsotopeDistances());
      f.applyMemberFunction(&UniqueIdInterface::setUniqueId);
      output_featmap.push_back(f);

      if (report_chromatograms_ && f.getIntensity() != 0)
      {
        output_chromatograms.push_back(feat_hypos[hypo_idx].getChromatograms(f.getUniqueId()));
      }


      if (overlapping_features_)
      {
        // We only consider the monoisotope and charge combination as used.
        trace_excl_map.emplace(labels[0], current_charge);
      }
      else
      {
        // old behavior: mark all traces in the feature as used and prevent them from forming new hypotheses
        for (Size lab_idx = 0; lab_idx < labels.size(); ++lab_idx)
        {
          trace_excl_map.emplace(labels[lab_idx], current_charge);
        }
      }

    }
    output_featmap.setUniqueId(UniqueIdGenerator::getUniqueId());
    output_featmap.sortByMZ();
  } // end of FeatureFindingPeptide::run

}
