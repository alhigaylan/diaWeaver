// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Mohammed Alhigaylan $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/ClusterMassTracesByPrecursor.h>
#include <OpenMS/PROCESSING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>

#include <algorithm>
#include <map>
#include <set>
#include <cmath>

namespace OpenMS
{

  ClusterMassTracesByPrecursor::ClusterMassTracesByPrecursor() :
    DefaultParamHandler("ClusterMassTracesByPrecursor"),
    ProgressLogger(),
    min_pearson_correlation_(0.7),
    max_lag_(1),
    max_rt_apex_difference_(5.0),
    min_nr_ions_(3),
    im_tolerance_(0.02),
    assign_unassigned_to_all_(false),
    rt_tolerance_(2.0),
    nr_precursors_per_fragment_(25),
    pearson_weight_(1.0),
    delta_rt_weight_(1.0),
    delta_im_weight_(1.0)
  {
    defaults_.setValue("min_pearson_correlation", 0.7,
      "Minimal Pearson correlation score to match elution profiles to each other.");
    defaults_.setMinFloat("min_pearson_correlation", 0.0);
    defaults_.setMaxFloat("min_pearson_correlation", 1.0);

    defaults_.setValue("max_lag", 1,
      "Maximal lag (e.g., by how many spectra the peak may be shifted at most). "
      "This parameter depends on chromatographic setup but 1-3 is usually sensible.");
    defaults_.setMinInt("max_lag", 0);

    defaults_.setValue("max_rt_apex_difference", 5.0,
      "Maximal difference of the apex in retention time (in seconds). "
      "Profiles further away will not be considered.");
    defaults_.setMinFloat("max_rt_apex_difference", 0.0);

    defaults_.setValue("min_nr_ions", 3,
      "Minimal number of fragment ions to report a pseudo spectrum.");
    defaults_.setMinInt("min_nr_ions", 1);

    defaults_.setValue("im_tolerance", 0.02,
      "Ion mobility tolerance for matching precursors to fragments. Set to 0 to disable IM filtering.");
    defaults_.setMinFloat("im_tolerance", 0.0);

    defaults_.setValue("assign_unassigned_to_all", "false",
      "Assign unassigned MS2 fragments to all precursors within RT range.");
    defaults_.setValidStrings("assign_unassigned_to_all", {"true", "false"});

    defaults_.setValue("rt_tolerance", 2.0,
      "RT tolerance (in seconds) for matching up mass trace points during correlation. "
      "Points within this tolerance are considered equal in RT.");
    defaults_.setMinFloat("rt_tolerance", 0.0);

    defaults_.setValue("nr_precursors_per_fragment", 25,
      "Maximum number of precursors a fragment can be assigned to. "
      "If a fragment correlates with more precursors, only the top N with highest combined scores are kept.");
    defaults_.setMinInt("nr_precursors_per_fragment", 1);

    defaults_.setValue("pearson_weight", 1.0,
      "Weight for the Pearson correlation component in the combined score used to rank precursor-fragment assignments.");
    defaults_.setMinFloat("pearson_weight", 0.0);

    defaults_.setValue("delta_rt_weight", 1.0,
      "Weight for the delta RT component in the combined score used to rank precursor-fragment assignments.");
    defaults_.setMinFloat("delta_rt_weight", 0.0);

    defaults_.setValue("delta_im_weight", 1.0,
      "Weight for the delta ion mobility component in the combined score used to rank precursor-fragment assignments.");
    defaults_.setMinFloat("delta_im_weight", 0.0);

    defaults_.setValue("output_fragment_scores", "false",
      "If true, output per-fragment scores (pearson_score, xcorr_lag, xcorr_lag_intensity, fragment_rt, fragment_ion_mobility, combined_score) as FloatDataArrays in the output spectra.");
    defaults_.setValidStrings("output_fragment_scores", {"false", "true"});

    defaultsToParam_();
  }

  ClusterMassTracesByPrecursor::~ClusterMassTracesByPrecursor() = default;

  void ClusterMassTracesByPrecursor::updateMembers_()
  {
    min_pearson_correlation_ = param_.getValue("min_pearson_correlation");
    max_lag_ = param_.getValue("max_lag");
    max_rt_apex_difference_ = param_.getValue("max_rt_apex_difference");
    min_nr_ions_ = static_cast<Size>(static_cast<int>(param_.getValue("min_nr_ions")));
    im_tolerance_ = param_.getValue("im_tolerance");
    assign_unassigned_to_all_ = param_.getValue("assign_unassigned_to_all").toBool();
    rt_tolerance_ = param_.getValue("rt_tolerance");
    nr_precursors_per_fragment_ = static_cast<Size>(static_cast<int>(param_.getValue("nr_precursors_per_fragment")));
    pearson_weight_ = param_.getValue("pearson_weight");
    delta_rt_weight_ = param_.getValue("delta_rt_weight");
    delta_im_weight_ = param_.getValue("delta_im_weight");
    output_fragment_scores_ = param_.getValue("output_fragment_scores").toBool();
  }

  void ClusterMassTracesByPrecursor::run(
      const ConsensusMap& ms1_traces,
      const ConsensusMap& ms2_traces,
      double swath_lower,
      double swath_upper,
      MSExperiment& pseudo_spectra)
  {
    if (ms1_traces.empty() || ms2_traces.empty())
    {
      return;
    }

    MasstraceCorrelator mtcorr;

    // Create cache for MS1 traces
    std::vector<MasstraceCorrelator::MasstracePointsType> precursor_profiles;
    std::vector<std::pair<double, double>> max_intensities_ms1;
    std::vector<double> precursor_rt;
    mtcorr.createConsensusMapCache(ms1_traces, precursor_profiles, max_intensities_ms1, precursor_rt);

    // Create cache for MS2 traces
    std::vector<MasstraceCorrelator::MasstracePointsType> fragment_profiles;
    std::vector<std::pair<double, double>> max_intensities_ms2;
    std::vector<double> fragment_rt;
    mtcorr.createConsensusMapCache(ms2_traces, fragment_profiles, max_intensities_ms2, fragment_rt);

    // Apply Savitzky-Golay smoothing to MS1 traces
    SavitzkyGolayFilter sg_filter;
    Param sg_params = sg_filter.getParameters();
    sg_params.setValue("frame_length", 5);
    sg_params.setValue("polynomial_order", 3);
    sg_filter.setParameters(sg_params);

    for (auto& trace : precursor_profiles)
    {
      MSSpectrum spec;
      for (const auto& p : trace)
      {
        Peak1D peak;
        peak.setMZ(p.first);  // Store RT as m/z for filtering
        peak.setIntensity(p.second);
        spec.push_back(peak);
      }

      sg_filter.filter(spec);

      for (Size i = 0; i < trace.size(); ++i)
      {
        trace[i].second = spec[i].getIntensity();
      }
    }

    // Check if IM data exists in the input maps
    bool has_im_data = false;
    if (!ms1_traces.empty() && ms1_traces[0].metaValueExists("Ion Mobility Centroid"))
    {
      has_im_data = true;
    }

    // Cache m/z values
    std::vector<double> precursor_mz;
    std::vector<double> precursor_im;
    std::vector<int> precursor_charge;
    std::vector<double> precursor_intensity;

    for (Size i = 0; i < ms1_traces.size(); ++i)
    {
      precursor_mz.push_back(ms1_traces[i].getMZ());
      precursor_charge.push_back(ms1_traces[i].getCharge());
      precursor_intensity.push_back(ms1_traces[i].getIntensity());

      if (has_im_data)
      {
        precursor_im.push_back(ms1_traces[i].getMetaValue("Ion Mobility Centroid"));
      }
      else
      {
        precursor_im.push_back(0.0);
      }
    }

    // Cache MS2 m/z, IM, and intensity values
    std::vector<double> fragment_mz;
    std::vector<double> fragment_im;
    std::vector<double> fragment_intensity;

    for (Size i = 0; i < ms2_traces.size(); ++i)
    {
      fragment_mz.push_back(ms2_traces[i].getMZ());
      fragment_intensity.push_back(ms2_traces[i].getIntensity());

      if (has_im_data)
      {
        fragment_im.push_back(ms2_traces[i].getMetaValue("Ion Mobility Centroid"));
      }
      else
      {
        fragment_im.push_back(0.0);
      }
    }

    // Run the core clustering algorithm
    clusterAndCreateSpectra_(
        precursor_profiles, precursor_mz, precursor_rt, precursor_im, precursor_charge, precursor_intensity,
        fragment_profiles, fragment_mz, fragment_rt, fragment_im, fragment_intensity,
        swath_lower, swath_upper, has_im_data, pseudo_spectra);
  }

  void ClusterMassTracesByPrecursor::run(
      const FeatureMap& ms1_features,
      const std::vector<MassTrace>& ms1_traces,
      const std::vector<MassTrace>& ms2_traces,
      double swath_lower,
      double swath_upper,
      MSExperiment& pseudo_spectra)
  {
    if (ms1_features.empty() || ms2_traces.empty())
    {
      return;
    }

    // Build lookup map from trace label to MassTrace for O(log N) access
    std::map<String, const MassTrace*> trace_lookup;
    for (const auto& tr : ms1_traces)
    {
      trace_lookup[tr.getLabel()] = &tr;
    }

    // Check if IM data exists in the input maps
    bool has_im_data = false;
    if (!ms1_features.empty())
    {
      const Feature& first_feature = ms1_features[0];
      if (first_feature.metaValueExists(Constants::UserParam::ION_MOBILITY_CENTROID) ||
          first_feature.metaValueExists("masstrace_centroid_im"))
      {
        has_im_data = true;
      }
    }
    if (!has_im_data && !ms2_traces.empty())
    {
      // Also check MS2 traces for IM data
      has_im_data = (ms2_traces[0].getCentroidIM() != 0.0);
    }

    // Extract elution profiles from MS1 features using their underlying mass traces
    std::vector<MasstraceCorrelator::MasstracePointsType> precursor_profiles;
    std::vector<double> precursor_rt;
    std::vector<double> precursor_mz;
    std::vector<double> precursor_im;
    std::vector<int> precursor_charge;
    std::vector<double> precursor_intensity;

    for (Size i = 0; i < ms1_features.size(); ++i)
    {
      const Feature& f = ms1_features[i];

      // Use monoisotopic m/z from the feature
      precursor_mz.push_back(f.getMZ());
      precursor_rt.push_back(f.getRT());
      precursor_charge.push_back(f.getCharge());
      precursor_intensity.push_back(f.getIntensity());

      // Get ion mobility if available
      if (has_im_data)
      {
        if (f.metaValueExists(Constants::UserParam::ION_MOBILITY_CENTROID))
        {
          precursor_im.push_back(f.getMetaValue(Constants::UserParam::ION_MOBILITY_CENTROID));
        }
        else if (f.metaValueExists("masstrace_centroid_im"))
        {
          std::vector<double> im_values = f.getMetaValue("masstrace_centroid_im");
          precursor_im.push_back(im_values.empty() ? 0.0 : im_values[0]);
        }
        else
        {
          precursor_im.push_back(0.0);
        }
      }
      else
      {
        precursor_im.push_back(0.0);
      }

      // Extract elution profile from the monoisotopic mass trace (real intensity values)
      MasstraceCorrelator::MasstracePointsType points;

      // Get the feature label and extract monoisotopic trace label (first part before '_')
      if (f.metaValueExists("label"))
      {
        String feat_label = f.getMetaValue("label");
        StringList label_tokens;
        feat_label.split("_", label_tokens);

        if (!label_tokens.empty())
        {
          String mono_label = label_tokens[0];
          auto it = trace_lookup.find(mono_label);
          if (it != trace_lookup.end())
          {
            const MassTrace& mono_trace = *(it->second);
            // Extract real RT/intensity pairs from the mass trace
            for (const Peak2D& peak : mono_trace)
            {
              points.push_back(std::make_pair(peak.getRT(), peak.getIntensity()));
            }
          }
        }
      }

      // Fallback: create single point at feature apex if no trace found
      if (points.empty())
      {
        points.push_back(std::make_pair(f.getRT(), f.getIntensity()));
        OPENMS_LOG_DEBUG << "Warning: No mass trace found for feature at m/z " << f.getMZ()
                         << ", using apex only" << std::endl;
      }

      precursor_profiles.push_back(points);
    }

    // Build cache for MS2 traces (elution profiles, m/z, RT, IM)
    std::vector<MasstraceCorrelator::MasstracePointsType> fragment_profiles;
    std::vector<double> fragment_rt;
    std::vector<double> fragment_mz;
    std::vector<double> fragment_im;
    std::vector<double> fragment_intensity;

    for (const auto& trace : ms2_traces)
    {
      // Extract elution profile (RT/intensity pairs)
      MasstraceCorrelator::MasstracePointsType points;
      for (const Peak2D& peak : trace)
      {
        points.push_back(std::make_pair(peak.getRT(), peak.getIntensity()));
      }
      fragment_profiles.push_back(points);

      // Cache centroid values
      fragment_rt.push_back(trace.getCentroidRT());
      fragment_mz.push_back(trace.getCentroidMZ());
      fragment_im.push_back(has_im_data ? trace.getCentroidIM() : 0.0);
      fragment_intensity.push_back(trace.getIntensity(false));
    }

    // Run the core clustering algorithm
    clusterAndCreateSpectra_(
        precursor_profiles, precursor_mz, precursor_rt, precursor_im, precursor_charge, precursor_intensity,
        fragment_profiles, fragment_mz, fragment_rt, fragment_im, fragment_intensity,
        swath_lower, swath_upper, has_im_data, pseudo_spectra);
  }

  void ClusterMassTracesByPrecursor::clusterAndCreateSpectra_(
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
      MSExperiment& pseudo_spectra)
  {
    MasstraceCorrelator mtcorr;

    // Stores all scores produced by scoreHullpoints for a precursor-fragment pair
    struct HullScore
    {
      int index;           // precursor or fragment index
      double pearson;      // Pearson correlation at zero lag
      int lag;             // lag at max cross-correlation
      double lag_intensity; // normalized cross-correlation at optimal lag
      double delta_rt;     // absolute RT difference between precursor and fragment
      double delta_im;     // absolute IM difference between precursor and fragment
    };

    // Track fragment assignments: for each fragment, store list of scored precursor matches
    std::vector<std::vector<HullScore>> fragment_assignments(fragment_profiles.size());

    // Assignment map: precursor index -> list of scored fragment matches
    std::map<int, std::vector<HullScore>> assignment_map;

    // -----------------------------------
    // Step 1 - assign fragment mass traces to precursors
    // -----------------------------------
    startProgress(0, precursor_profiles.size(), "Assigning fragments to precursors");
    for (Size i = 0; i < precursor_profiles.size(); ++i)
    {
      setProgress(i);

      // Only consider precursors within the SWATH window
      if (precursor_mz[i] < swath_lower || precursor_mz[i] > swath_upper)
      {
        continue;
      }

      assignment_map[i].clear();
      double current_rt = precursor_rt[i];

      for (Size j = 0; j < fragment_profiles.size(); ++j)
      {
        // Check RT distance
        if (std::fabs(current_rt - fragment_rt[j]) > max_rt_apex_difference_)
        {
          continue;
        }

        // Check ion mobility tolerance (if IM data is present)
        if (has_im_data && im_tolerance_ > 0)
        {
          if (std::fabs(precursor_im[i] - fragment_im[j]) > im_tolerance_)
          {
            continue;
          }
        }

        // Exclude fragments within precursor isolation window
        if (fragment_mz[j] >= swath_lower && fragment_mz[j] <= swath_upper)
        {
          continue;
        }

        // Score the precursor elution profile against the fragment mass trace
        int lag;
        double lag_intensity, pearson_score;
        mtcorr.scoreHullpoints(precursor_profiles[i], fragment_profiles[j],
                               lag, lag_intensity, pearson_score,
                               min_pearson_correlation_, max_lag_, rt_tolerance_);

        if (pearson_score > min_pearson_correlation_ && lag >= -max_lag_ && lag <= max_lag_)
        {
          // Collect all valid precursor-fragment assignments
          double delta_rt = std::fabs(precursor_rt[i] - fragment_rt[j]);
          double delta_im = has_im_data ? std::fabs(precursor_im[i] - fragment_im[j]) : 0.0;
          fragment_assignments[j].push_back({static_cast<int>(i), pearson_score, lag, lag_intensity, delta_rt, delta_im});
          assignment_map[i].push_back({static_cast<int>(j), pearson_score, lag, lag_intensity, delta_rt, delta_im});
        }
      }
    }
    endProgress();

    // -----------------------------------
    // Filter fragments assigned to too many precursors - keep only top N precursor per fragment using a combined score
    // -----------------------------------
    // Combined score = w_pearson * normalized_pearson + w_rt * normalized_rt + w_im * normalized_im
    // normalized_pearson: (pearson - min_pearson_correlation_) / (1 - min_pearson_correlation_)  -> [0, 1]
    // normalized_rt:      1 - (delta_rt / max_rt_apex_difference_)                               -> [0, 1]
    // normalized_im:      1 - (delta_im / im_tolerance_)                                         -> [0, 1]
    double pearson_range = 1.0 - min_pearson_correlation_;
    Size cnt_fragments_filtered = 0;
    for (Size j = 0; j < fragment_assignments.size(); ++j)
    {
      if (fragment_assignments[j].size() > nr_precursors_per_fragment_)
      {
        cnt_fragments_filtered++;

        // Sort by combined score (descending)
        std::sort(fragment_assignments[j].begin(), fragment_assignments[j].end(),
                  [&](const HullScore& a, const HullScore& b)
                  {
                    auto combined = [&](const HullScore& hs)
                    {
                      double norm_pearson = (hs.pearson - min_pearson_correlation_) / pearson_range;
                      double norm_rt = 1.0 - (hs.delta_rt / max_rt_apex_difference_);
                      double norm_im = (has_im_data && im_tolerance_ > 0) ? 1.0 - (hs.delta_im / im_tolerance_) : 1.0;
                      return (pearson_weight_ * norm_pearson) + (delta_rt_weight_ * norm_rt) + (delta_im_weight_ * norm_im);
                    };
                    return combined(a) > combined(b);
                  });

        // Keep only top N precursors
        fragment_assignments[j].resize(nr_precursors_per_fragment_);
      }
    }

    if (cnt_fragments_filtered > 0)
    {
      OPENMS_LOG_INFO << "Filtered " << cnt_fragments_filtered << " fragments that were assigned to more than "
                      << nr_precursors_per_fragment_ << " precursors (kept top " << nr_precursors_per_fragment_
                      << " by combined score)" << std::endl;
    }

    // -----------------------------------
    // Rebuild assignment_map from filtered fragment_assignments
    // -----------------------------------
    assignment_map.clear();
    for (Size j = 0; j < fragment_assignments.size(); ++j)
    {
      for (const auto& hs : fragment_assignments[j])
      {
        assignment_map[hs.index].push_back({static_cast<int>(j), hs.pearson, hs.lag, hs.lag_intensity, hs.delta_rt, hs.delta_im});
      }
    }

    // Stats
    Size cnt_ms2_used = 0;
    Size cnt_ms1_used = 0;
    for (Size i = 0; i < precursor_profiles.size(); ++i)
    {
      if (!assignment_map[i].empty())
      {
        cnt_ms1_used++;
      }
    }
    for (Size i = 0; i < fragment_assignments.size(); ++i)
    {
      if (!fragment_assignments[i].empty())
      {
        cnt_ms2_used++;
      }
    }

    OPENMS_LOG_INFO << "Assigned " << cnt_ms2_used << " (out of " << fragment_profiles.size()
                    << ") MS2 fragments to " << cnt_ms1_used << " (out of " << precursor_profiles.size()
                    << ") precursors" << std::endl;

    // -----------------------------------
    // Step 2 - assign the unused fragment ions (if requested)
    // -----------------------------------
    int cnt_unassigned = 0;
    if (assign_unassigned_to_all_)
    {
      startProgress(0, fragment_profiles.size(), "Assigning unused fragments");
      for (Size j = 0; j < fragment_profiles.size(); ++j)
      {
        setProgress(j);
        if (!fragment_assignments[j].empty())
        {
          continue;
        }
        cnt_unassigned++;

        // Find suitable precursors to assign these fragments to
        for (Size i = 0; i < precursor_profiles.size(); ++i)
        {
          if (precursor_mz[i] < swath_lower || precursor_mz[i] > swath_upper)
          {
            continue;
          }
          if (assignment_map[i].empty())
          {
            continue;
          }
          if (std::fabs(precursor_rt[i] - fragment_rt[j]) > max_rt_apex_difference_)
          {
            continue;
          }

          // Assign to all matching precursors with score 0 (will be filtered first if too many)
          double delta_rt = std::fabs(precursor_rt[i] - fragment_rt[j]);
          double delta_im = has_im_data ? std::fabs(precursor_im[i] - fragment_im[j]) : 0.0;
          assignment_map[i].push_back({static_cast<int>(j), 0.0, 0, 0.0, delta_rt, delta_im});
        }
      }
      endProgress();

      OPENMS_LOG_INFO << "There were " << cnt_unassigned << " (out of " << fragment_profiles.size()
                      << ") unused fragment ions assigned to all precursors within RT range." << std::endl;
    }

    // -----------------------------------
    // Step 3 - create spectra from assignments
    // -----------------------------------
    int cnt_spectra = 0;
    startProgress(0, precursor_profiles.size(), "Creating pseudo spectra");
    for (Size i = 0; i < precursor_profiles.size(); ++i)
    {
      setProgress(i);

      if (precursor_mz[i] < swath_lower || precursor_mz[i] > swath_upper)
      {
        continue;
      }

      if (assignment_map[i].size() < min_nr_ions_)
      {
        continue;
      }

      MSSpectrum spectrum;
      spectrum.setRT(precursor_rt[i]);
      spectrum.setMSLevel(2);
      spectrum.setType(SpectrumSettings::SpectrumType::CENTROID);

      // Set ion mobility if available
      if (has_im_data)
      {
        spectrum.setDriftTime(precursor_im[i]);
      }

      // Set precursor information
      Precursor p;
      p.setMZ(precursor_mz[i]);
      p.setCharge(precursor_charge[i]);
      p.setIntensity(precursor_intensity[i]);
      if (has_im_data)
      {
        p.setDriftTime(precursor_im[i]);
        p.setDriftTimeUnit(DriftTimeUnit::VSSC);
      }
      spectrum.setPrecursors({p});

      // Get assignments for this precursor
      auto& assignments = assignment_map[i];

      // If more than 500 assignments to one precursor, sort by combined score and keep top 500
      if (assignments.size() > 500)
      {
        std::sort(assignments.begin(), assignments.end(),
                  [&](const HullScore& a, const HullScore& b)
                  {
                    auto combined = [&](const HullScore& hs)
                    {
                      double norm_pearson = (hs.pearson - min_pearson_correlation_) / pearson_range;
                      double norm_rt = 1.0 - (hs.delta_rt / max_rt_apex_difference_);
                      double norm_im = (has_im_data && im_tolerance_ > 0) ? 1.0 - (hs.delta_im / im_tolerance_) : 1.0;
                      return (pearson_weight_ * norm_pearson) + (delta_rt_weight_ * norm_rt) + (delta_im_weight_ * norm_im);
                    };
                    return combined(a) > combined(b);
                  });
        assignments.resize(500);
      }

      if (output_fragment_scores_)
      {
        spectrum.getFloatDataArrays().resize(6);
        spectrum.getFloatDataArrays()[0].setName("pearson_score");
        spectrum.getFloatDataArrays()[1].setName("xcorr_lag");
        spectrum.getFloatDataArrays()[2].setName("xcorr_lag_intensity");
        spectrum.getFloatDataArrays()[3].setName("fragment_rt");
        spectrum.getFloatDataArrays()[4].setName("fragment_ion_mobility");
        spectrum.getFloatDataArrays()[5].setName("combined_score");
      }

      // Add fragment ions with their scores, RT, and ion mobility
      for (const auto& hs : assignments)
      {
        Peak1D peak;
        peak.setMZ(fragment_mz[hs.index]);
        peak.setIntensity(fragment_intensity[hs.index]);
        spectrum.push_back(peak);

        if (output_fragment_scores_)
        {
          double norm_pearson = (pearson_range > 0) ? (hs.pearson - min_pearson_correlation_) / pearson_range : 1.0;
          double norm_rt = 1.0 - (hs.delta_rt / max_rt_apex_difference_);
          double norm_im = (has_im_data && im_tolerance_ > 0) ? 1.0 - (hs.delta_im / im_tolerance_) : 1.0;
          double combined = (pearson_weight_ * norm_pearson) + (delta_rt_weight_ * norm_rt) + (delta_im_weight_ * norm_im);
          spectrum.getFloatDataArrays()[0].push_back(hs.pearson);
          spectrum.getFloatDataArrays()[1].push_back(static_cast<float>(hs.lag));
          spectrum.getFloatDataArrays()[2].push_back(hs.lag_intensity);
          spectrum.getFloatDataArrays()[3].push_back(fragment_rt[hs.index]);
          spectrum.getFloatDataArrays()[4].push_back(fragment_im[hs.index]);
          spectrum.getFloatDataArrays()[5].push_back(combined);
        }
      }

      if (spectrum.size() >= min_nr_ions_)
      {
        pseudo_spectra.addSpectrum(spectrum);
        cnt_spectra++;
      }
    }
    endProgress();

    OPENMS_LOG_INFO << "Created " << cnt_spectra << " pseudo spectra with more than "
                    << min_nr_ions_ << " fragment ions." << std::endl;
  }

} // namespace OpenMS