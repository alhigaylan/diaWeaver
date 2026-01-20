// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Mohammed Alhigaylan $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/KERNEL/ConvexHull2D.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/APPLICATIONS/diaWeaver.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerIM.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include <OpenMS/FEATUREFINDER/ElutionPeakDetection.h>
#include <OpenMS/FEATUREFINDER/FeatureFindingMetabo.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MasstraceCorrelator.h>
#include <OpenMS/PROCESSING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <cmath>
#include <set>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace OpenMS;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_diaWeaver diaWeaver

@brief Splits a DIA mzML file into per-window MS1 and MS2 mzML files with feature detection and pseudo spectra generation.

This tool extracts DIA windows from ion mobility DIA data and writes separate mzML files
for each precursor isolation window. It processes both MS1 and MS2 spectra, filtering
peaks by m/z and ion mobility ranges defined by the DIA acquisition windows.

The tool applies peak picking to all extracted spectra using either PeakPickerIM (for
ion mobility data) or PeakPickerHiRes (for non-IM data). After peak picking:
- FeatureFinderMetabo runs on MS1 data to detect monoisotopic peptide features (featureXML output)
- FeatureFinderMetabo runs on unfragmented precursor data to detect precursor features (featureXML output)
- MassTraceExtractor runs on MS2 data to extract fragment mass traces (consensusXML output)
- ClusterMassTraces correlates MS1 features (monoisotopic) with MS2 fragment traces to create pseudo spectra (mzML output)

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_diaWeaver.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_diaWeaver.html
*/

class TOPPDiaWeaver :
  public TOPPBase
{
public:
  TOPPDiaWeaver() :
    TOPPBase(
      "diaWeaver",
      "Split a DIA mzML into per-window MS1 and MS2 mzML files",
      false)
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_(
      "in",
      "<file>",
      "",
      "Input DIA mzML file",
      true);
    setValidFormats_("in", {"mzML"});

    registerOutputPrefix_(
      "out",
      "<prefix>",
      "",
      "Output directory prefix (default: <input>_diaWindows)");

    registerFlag_(
      "save_unfragmented_precursors",
      "If set, save peaks within precursor isolation window in MS2 and apply precursor detection algorithm");

    registerSubsection_("PeakPickerIM", "Parameters for ion mobility peak picking (used when input has IM data)");

    registerSubsection_("PeakPickerHiRes", "Parameters for high-resolution peak picking (used when input has no IM data)");

    registerSubsection_("FeatureFinderMetabo", "Parameters for FeatureFinderMetabo algorithm (precursor detection on MS1 and unfragmented precursors)");

    registerSubsection_("MassTraceExtractor", "Parameters for MassTraceExtractor algorithm (mass trace detection on MS2 data)");

    registerSubsection_("ClusterMassTraces", "Parameters for clustering mass traces to create pseudo spectra");

    registerIntOption_(
      "threads",
      "<n>",
      1,
      "Total number of threads to use for processing",
      false);
    setMinInt_("threads", 1);

    registerIntOption_(
      "threads_outer_loop",
      "<n>",
      -1,
      "Number of threads for the outer loop (over DIA windows). Remaining threads are used for "
      "inner loop (peak picking within each window). Set to -1 to use all threads in outer loop only (no nested parallelism). "
      "Example: with 24 total threads and 4 outer threads, each window gets 6 threads for peak picking.",
      false);
  }

  Param getSubsectionDefaults_(const String& name) const override
  {
    if (name == "PeakPickerIM")
    {
      PeakPickerIM ppim;
      return ppim.getDefaults();
    }
    if (name == "PeakPickerHiRes")
    {
      PeakPickerHiRes pphr;
      return pphr.getDefaults();
    }
    if (name == "FeatureFinderMetabo")
    {
      Param combined;

      // Common parameters for all FFM sub-algorithms
      Param p_com;
      p_com.setValue("noise_threshold_int", 10.0, "Intensity threshold below which peaks are regarded as noise.");
      p_com.setValue("chrom_peak_snr", 3.0, "Minimum signal-to-noise a mass trace should have.");
      p_com.setValue("chrom_fwhm", 5.0, "Expected chromatographic peak width (in seconds).");
      combined.insert("common:", p_com);
      combined.setSectionDescription("common", "Common parameters for all other subsections");

      // MassTraceDetection parameters
      Param p_mtd = MassTraceDetection().getDefaults();
      p_mtd.remove("noise_threshold_int");
      p_mtd.remove("chrom_peak_snr");
      combined.insert("mtd:", p_mtd);
      combined.setSectionDescription("mtd", "Mass Trace Detection parameters");

      // ElutionPeakDetection parameters
      Param p_epd;
      p_epd.setValue("enabled", "true", "Enable splitting of isobaric mass traces by chromatographic peak detection. Disable for direct injection.");
      p_epd.setValidStrings("enabled", {"true", "false"});
      p_epd.insert("", ElutionPeakDetection().getDefaults());
      p_epd.remove("chrom_peak_snr");
      p_epd.remove("chrom_fwhm");
      combined.insert("epd:", p_epd);
      combined.setSectionDescription("epd", "Elution Profile Detection (to separate isobaric Mass Traces by elution time).");

      // FeatureFindingMetabo parameters
      Param p_ffm = FeatureFindingMetabo().getDefaults();
      p_ffm.remove("chrom_fwhm");
      p_ffm.remove("report_chromatograms");
      combined.insert("ffm:", p_ffm);
      combined.setSectionDescription("ffm", "FeatureFinder parameters (assembling mass traces to charged features)");

      return combined;
    }
    if (name == "MassTraceExtractor")
    {
      Param combined;

      // Common parameters
      Param p_com;
      p_com.setValue("noise_threshold_int", 10.0, "Intensity threshold below which peaks are regarded as noise.");
      p_com.setValue("chrom_peak_snr", 3.0, "Minimum signal-to-noise a mass trace should have.");
      p_com.setValue("chrom_fwhm", 5.0, "Expected chromatographic peak width (in seconds).");
      combined.insert("common:", p_com);
      combined.setSectionDescription("common", "Common parameters for all other subsections");

      // MassTraceDetection parameters
      Param p_mtd = MassTraceDetection().getDefaults();
      p_mtd.remove("noise_threshold_int");
      p_mtd.remove("chrom_peak_snr");
      combined.insert("mtd:", p_mtd);
      combined.setSectionDescription("mtd", "Mass Trace Detection parameters");

      // ElutionPeakDetection parameters
      Param p_epd;
      p_epd.setValue("enabled", "true", "Enable splitting of isobaric mass traces by chromatographic peak detection.");
      p_epd.setValidStrings("enabled", {"true", "false"});
      p_epd.insert("", ElutionPeakDetection().getDefaults());
      p_epd.remove("chrom_peak_snr");
      p_epd.remove("chrom_fwhm");
      combined.insert("epd:", p_epd);
      combined.setSectionDescription("epd", "Elution Profile Detection (to separate isobaric Mass Traces by elution time).");

      return combined;
    }
    if (name == "ClusterMassTraces")
    {
      Param p;
      p.setValue("min_pearson_correlation", 0.7, "Minimal pearson correlation score to match elution profiles to each other.");
      p.setMinFloat("min_pearson_correlation", 0.0);
      p.setMaxFloat("min_pearson_correlation", 1.0);

      p.setValue("max_lag", 1, "Maximal lag (e.g. by how many spectra the peak may be shifted at most).");
      p.setMinInt("max_lag", 0);

      p.setValue("min_nr_ions", 3, "Minimal number of ions to report a spectrum.");
      p.setMinInt("min_nr_ions", 1);

      p.setValue("max_rt_apex_difference", 5.0, "Maximal difference of the apex in retention time (in seconds).");
      p.setMinFloat("max_rt_apex_difference", 0.0);

      p.setValue("im_tolerance", 0.02, "Ion mobility tolerance for matching precursors to fragments.");
      p.setMinFloat("im_tolerance", 0.0);

      return p;
    }
    return Param();
  }

  /**
   * @brief Run FeatureFinderMetabo pipeline on a centroided MSExperiment
   * @param[in,out] ms_peakmap Input centroided peak map (will be sorted)
   * @param[in] common_param Common parameters for FFM algorithms
   * @param[in] mtd_param MassTraceDetection parameters
   * @param[in] epd_param ElutionPeakDetection parameters
   * @param[in] ffm_param FeatureFindingMetabo parameters
   * @param[out] feat_map Output feature map
   * @return True on success, false on error
   */
  bool runFeatureFinderMetabo_(MSExperiment& ms_peakmap,
                               const Param& common_param,
                               Param mtd_param,
                               Param epd_param,
                               Param ffm_param,
                               FeatureMap& feat_map)
  {
    if (ms_peakmap.empty())
    {
      return true; // Nothing to process
    }

    // Ensure spectra are sorted by m/z
    ms_peakmap.sortSpectra(true);

    std::vector<MassTrace> m_traces;

    // Configure and run mass trace detection
    MassTraceDetection mtdet;
    mtd_param.insert("", common_param);
    mtd_param.remove("chrom_fwhm");
    mtdet.setParameters(mtd_param);
    mtdet.run(ms_peakmap, m_traces);

    if (m_traces.empty())
    {
      OPENMS_LOG_INFO << "No mass traces detected." << std::endl;
      return true;
    }

    // Configure and run elution peak detection
    std::vector<MassTrace> m_traces_final;
    if (epd_param.getValue("enabled").toBool())
    {
      std::vector<MassTrace> split_mtraces;
      epd_param.remove("enabled");
      epd_param.insert("", common_param);
      epd_param.remove("noise_threshold_int");
      ElutionPeakDetection epdet;
      epdet.setParameters(epd_param);
      epdet.detectPeaks(m_traces, split_mtraces);
      if (epdet.getParameters().getValue("width_filtering") == "auto")
      {
        m_traces_final.clear();
        epdet.filterByPeakWidth(split_mtraces, m_traces_final);
      }
      else
      {
        m_traces_final = split_mtraces;
      }
    }
    else
    {
      m_traces_final = m_traces;
      for (Size i = 0; i < m_traces_final.size(); ++i)
      {
        m_traces_final[i].estimateFWHM(false);
      }
      if (ffm_param.getValue("use_smoothed_intensities").toBool())
      {
        OPENMS_LOG_WARN << "Without EPD, smoothing is not supported. Setting 'use_smoothed_intensities' to false!" << std::endl;
        ffm_param.setValue("use_smoothed_intensities", "false");
      }
    }

    // Configure and run feature finding
    ffm_param.insert("", common_param);
    ffm_param.remove("noise_threshold_int");
    ffm_param.remove("chrom_peak_snr");
    ffm_param.setValue("report_chromatograms", "false");

    std::vector<std::vector<MSChromatogram>> feat_chromatograms;
    FeatureFindingMetabo ffmet;
    ffmet.setParameters(ffm_param);
    ffmet.run(m_traces_final, feat_map, feat_chromatograms);

    // Filter features with zero intensity
    auto intensity_zero = [](Feature& f) { return f.getIntensity() == 0; };
    feat_map.erase(std::remove_if(feat_map.begin(), feat_map.end(), intensity_zero), feat_map.end());

    OPENMS_LOG_INFO << "FFMetabo: " << m_traces_final.size() << " traces -> "
                    << feat_map.size() << " features" << std::endl;

    return true;
  }

  /**
   * @brief Run MassTraceExtractor pipeline on a centroided MSExperiment
   * @param[in,out] ms_peakmap Input centroided peak map (will be sorted)
   * @param[in] common_param Common parameters for MTE algorithms
   * @param[in] mtd_param MassTraceDetection parameters
   * @param[in] epd_param ElutionPeakDetection parameters
   * @param[out] consensus_map Output consensus map with mass traces
   * @return True on success, false on error
   */
  bool runMassTraceExtractor_(MSExperiment& ms_peakmap,
                              const Param& common_param,
                              Param mtd_param,
                              Param epd_param,
                              ConsensusMap& consensus_map)
  {
    if (ms_peakmap.empty())
    {
      return true; // Nothing to process
    }

    // Ensure spectra are sorted by m/z
    ms_peakmap.sortSpectra(true);

    std::vector<MassTrace> m_traces;

    // Configure and run mass trace detection
    MassTraceDetection mtdet;
    mtd_param.insert("", common_param);
    mtd_param.remove("chrom_fwhm");
    mtdet.setParameters(mtd_param);
    mtdet.run(ms_peakmap, m_traces);

    if (m_traces.empty())
    {
      OPENMS_LOG_INFO << "No mass traces detected." << std::endl;
      return true;
    }

    // Configure and run elution peak detection if enabled
    std::vector<MassTrace> m_traces_final;
    bool use_epd = epd_param.getValue("enabled").toBool();

    if (use_epd)
    {
      std::vector<MassTrace> split_mtraces;
      epd_param.remove("enabled");
      epd_param.insert("", common_param);
      epd_param.remove("noise_threshold_int");
      ElutionPeakDetection epdet;
      epdet.setParameters(epd_param);
      epdet.detectPeaks(m_traces, split_mtraces);
      if (epdet.getParameters().getValue("width_filtering") == "auto")
      {
        m_traces_final.clear();
        epdet.filterByPeakWidth(split_mtraces, m_traces_final);
      }
      else
      {
        m_traces_final = split_mtraces;
      }
    }
    else
    {
      m_traces_final = m_traces;
    }

    // Convert mass traces to ConsensusMap
    for (Size i = 0; i < m_traces_final.size(); ++i)
    {
      if (m_traces_final[i].getSize() == 0)
      {
        continue;
      }
      ConsensusFeature fcons;
      int k = 0;
      for (const Peak2D& mss : m_traces_final[i])
      {
        FeatureHandle fhandle;
        fhandle.setRT(mss.getRT());
        fhandle.setMZ(mss.getMZ());
        fhandle.setIntensity(mss.getIntensity());
        fhandle.setUniqueId(++k);
        fcons.insert(fhandle);
      }

      fcons.setMetaValue(3, m_traces_final[i].getLabel());
      fcons.setCharge(0);
      fcons.setWidth(m_traces_final[i].estimateFWHM(use_epd));
      fcons.setQuality(1 - (1.0 / m_traces_final[i].getSize()));

      fcons.setRT(m_traces_final[i].getCentroidRT());
      fcons.setMZ(m_traces_final[i].getCentroidMZ());
      fcons.setIntensity(m_traces_final[i].getIntensity(false));

      // Attach mz peak FWHM if available
      if (m_traces_final[i].fwhm_mz_avg > 0)
      {
        fcons.setMetaValue(Constants::UserParam::FWHM_MZ_AVG, m_traces_final[i].fwhm_mz_avg);
      }
      // Annotate with ion mobility centroid if available
      if (m_traces_final[i].getCentroidIM() > 0)
      {
        fcons.setMetaValue(Constants::UserParam::ION_MOBILITY_CENTROID, m_traces_final[i].getCentroidIM());
      }
      // Add ion mobility peak FWHM if available
      if (m_traces_final[i].fwhm_im_avg > 0)
      {
        fcons.setMetaValue(Constants::UserParam::FWHM_IM_AVG, m_traces_final[i].fwhm_im_avg);
      }

      consensus_map.push_back(fcons);
    }

    consensus_map.applyMemberFunction(&UniqueIdInterface::setUniqueId);
    consensus_map.setUniqueId();

    OPENMS_LOG_INFO << "MassTraceExtractor: " << m_traces.size() << " traces -> "
                    << m_traces_final.size() << " final traces -> "
                    << consensus_map.size() << " consensus features" << std::endl;

    return true;
  }

  /**
   * @brief Cluster MS1 features (precursors) with MS2 mass traces to create pseudo spectra
   * @param[in] ms1_features MS1 peptide features from FeatureFinderMetabo (monoisotopic)
   * @param[in] ms2_traces MS2 mass traces (fragments)
   * @param[in] swath_lower Lower m/z bound of the DIA window
   * @param[in] swath_upper Upper m/z bound of the DIA window
   * @param[in] cluster_param Clustering parameters
   * @param[out] pseudo_spectra Output pseudo spectra
   * @return True on success
   */
  bool clusterMassTraces_(const FeatureMap& ms1_features,
                          const ConsensusMap& ms2_traces,
                          double swath_lower,
                          double swath_upper,
                          const Param& cluster_param,
                          MSExperiment& pseudo_spectra)
  {
    if (ms1_features.empty() || ms2_traces.empty())
    {
      return true; // Nothing to cluster
    }

    // Get parameters
    double min_pscore = cluster_param.getValue("min_pearson_correlation");
    int max_lag = cluster_param.getValue("max_lag");
    double rt_max_distance = cluster_param.getValue("max_rt_apex_difference");
    Size min_nr_ions = (Size)((int)cluster_param.getValue("min_nr_ions"));
    double im_tolerance = cluster_param.getValue("im_tolerance");
    double mindiff = 2.0; // RT tolerance for correlation

    MasstraceCorrelator mtcorr;

    // Extract elution profiles from MS1 features (using convex hulls)
    std::vector<MasstraceCorrelator::MasstracePointsType> feature_points_ms1;
    std::vector<double> rt_cache_ms1;
    std::vector<double> mz_cache_ms1;
    std::vector<double> im_cache_ms1;
    std::vector<int> charge_cache_ms1;

    for (Size i = 0; i < ms1_features.size(); ++i)
    {
      const Feature& f = ms1_features[i];

      // Use monoisotopic m/z from the feature
      mz_cache_ms1.push_back(f.getMZ());
      rt_cache_ms1.push_back(f.getRT());
      charge_cache_ms1.push_back(f.getCharge());

      // Get ion mobility if available (check both possible meta value names)
      if (f.metaValueExists(Constants::UserParam::ION_MOBILITY_CENTROID))
      {
        im_cache_ms1.push_back(f.getMetaValue(Constants::UserParam::ION_MOBILITY_CENTROID));
      }
      else if (f.metaValueExists("masstrace_centroid_im"))
      {
        // Use first isotope's IM centroid
        std::vector<double> im_values = f.getMetaValue("masstrace_centroid_im");
        im_cache_ms1.push_back(im_values.empty() ? 0.0 : im_values[0]);
      }
      else
      {
        im_cache_ms1.push_back(0.0);
      }

      // Extract elution profile from mass trace meta values (stored by FeatureFindingMetabo)
      MasstraceCorrelator::MasstracePointsType points;
      if (f.metaValueExists("masstrace_centroid_rt") && f.metaValueExists("masstrace_intensity"))
      {
        // Get the RT and intensity arrays for the monoisotopic trace (first isotope)
        std::vector<double> rt_values = f.getMetaValue("masstrace_centroid_rt");
        std::vector<double> int_values = f.getMetaValue("masstrace_intensity");

        // Use the monoisotopic trace (index 0) if available
        // Note: These are centroid values per isotope, not the full elution profile
        // For full profile, we need to look at convex hull + reconstruct from mass trace
        if (!rt_values.empty() && !int_values.empty())
        {
          // FFM stores one centroid RT/intensity per isotope trace, not the full profile
          // Use feature apex as single point for now
          points.push_back(std::make_pair(f.getRT(), f.getIntensity()));
        }
      }

      // Try to extract from convex hull if we have one (hull points are at actual RT positions)
      if (points.empty() && !f.getConvexHulls().empty())
      {
        const ConvexHull2D& hull = f.getConvexHulls()[0]; // Monoisotopic trace hull
        const auto& hull_points = hull.getHullPoints();

        // Convex hull stores (RT, m/z) - we need intensity
        // Get unique RT values from hull and use feature intensity as approximation
        std::set<double> unique_rts;
        for (const auto& pt : hull_points)
        {
          unique_rts.insert(pt[0]); // RT is first coordinate
        }

        // Create points at each RT with interpolated intensity (simplified: use feature intensity)
        double peak_rt = f.getRT();
        double peak_int = f.getIntensity();
        double width = f.getWidth() > 0 ? f.getWidth() : 10.0; // FWHM or default

        for (double rt : unique_rts)
        {
          // Gaussian approximation of elution profile
          double rt_diff = rt - peak_rt;
          double sigma = width / 2.355; // FWHM to sigma
          double intensity = peak_int * std::exp(-0.5 * (rt_diff / sigma) * (rt_diff / sigma));
          points.push_back(std::make_pair(rt, intensity));
        }

        // Sort by RT
        std::sort(points.begin(), points.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });
      }

      // Fallback: create single point at feature apex
      if (points.empty())
      {
        points.push_back(std::make_pair(f.getRT(), f.getIntensity()));
      }
      feature_points_ms1.push_back(points);
    }

    // Cache data structures for MS2 traces
    std::vector<MasstraceCorrelator::MasstracePointsType> feature_points_ms2;
    std::vector<std::pair<double, double>> max_intensities_ms2;
    std::vector<double> rt_cache_ms2;
    mtcorr.createConsensusMapCache(ms2_traces, feature_points_ms2, max_intensities_ms2, rt_cache_ms2);

    // Cache MS2 m/z and IM values
    std::vector<double> mz_cache_ms2;
    std::vector<double> im_cache_ms2;
    for (Size i = 0; i < ms2_traces.size(); ++i)
    {
      mz_cache_ms2.push_back(ms2_traces[i].getMZ());
      if (ms2_traces[i].metaValueExists(Constants::UserParam::ION_MOBILITY_CENTROID))
      {
        im_cache_ms2.push_back(ms2_traces[i].getMetaValue(Constants::UserParam::ION_MOBILITY_CENTROID));
      }
      else
      {
        im_cache_ms2.push_back(0.0);
      }
    }

    // Assignment map: MS1 feature index -> list of MS2 trace indices
    std::map<int, std::vector<int>> ms1_assignment_map;

    // Assign fragment mass traces to precursor features
    for (Size i = 0; i < ms1_features.size(); ++i)
    {
      // Only consider MS1 features within the SWATH window as potential precursors
      if (mz_cache_ms1[i] < swath_lower || mz_cache_ms1[i] > swath_upper) continue;

      ms1_assignment_map[i].clear();
      double current_rt = rt_cache_ms1[i];

      for (Size j = 0; j < ms2_traces.size(); ++j)
      {
        // Check RT distance
        if (fabs(current_rt - rt_cache_ms2[j]) > rt_max_distance) continue;

        // Check ion mobility tolerance (if both have IM data)
        if (im_tolerance > 0 && im_cache_ms1[i] > 0 && im_cache_ms2[j] > 0)
        {
          if (fabs(im_cache_ms1[i] - im_cache_ms2[j]) > im_tolerance) continue;
        }

        // Exclude fragments within precursor isolation window
        if (mz_cache_ms2[j] >= swath_lower && mz_cache_ms2[j] <= swath_upper) continue;

        // Score the MS1 feature elution profile against the MS2 mass trace
        int lag;
        double lag_intensity, pearson_score;
        mtcorr.scoreHullpoints(feature_points_ms1[i], feature_points_ms2[j],
                               lag, lag_intensity, pearson_score, min_pscore, max_lag, mindiff);

        if (pearson_score > min_pscore && lag >= -max_lag && lag <= max_lag)
        {
          ms1_assignment_map[i].push_back(j);
        }
      }

      // Only keep assignments with enough ions
      if (ms1_assignment_map[i].size() < min_nr_ions)
      {
        ms1_assignment_map[i].clear();
      }
    }

    // Create pseudo spectra from assignments
    for (Size i = 0; i < ms1_features.size(); ++i)
    {
      if (mz_cache_ms1[i] < swath_lower || mz_cache_ms1[i] > swath_upper) continue;
      if (ms1_assignment_map[i].size() < min_nr_ions) continue;

      MSSpectrum spectrum;
      spectrum.setRT(ms1_features[i].getRT());
      spectrum.setMSLevel(2);
      spectrum.setType(SpectrumSettings::SpectrumType::CENTROID);

      // Set ion mobility if available
      if (im_cache_ms1[i] > 0)
      {
        spectrum.setDriftTime(im_cache_ms1[i]);
      }

      // Set precursor with monoisotopic m/z
      Precursor p;
      p.setMZ(ms1_features[i].getMZ());  // Monoisotopic m/z from FFM
      p.setCharge(charge_cache_ms1[i]);
      p.setIntensity(ms1_features[i].getIntensity());
      if (im_cache_ms1[i] > 0)
      {
        p.setDriftTime(im_cache_ms1[i]);
        p.setDriftTimeUnit(DriftTimeUnit::VSSC);
      }
      spectrum.setPrecursors({p});

      // Add fragment ions
      for (int ms2_idx : ms1_assignment_map[i])
      {
        Peak1D peak;
        peak.setMZ(ms2_traces[ms2_idx].getMZ());
        peak.setIntensity(ms2_traces[ms2_idx].getIntensity());
        spectrum.push_back(peak);
      }

      if (spectrum.size() >= min_nr_ions)
      {
        pseudo_spectra.addSpectrum(spectrum);
      }
    }

    OPENMS_LOG_INFO << "ClusterMassTraces: Created " << pseudo_spectra.size()
                    << " pseudo spectra from " << ms1_features.size() << " MS1 features and "
                    << ms2_traces.size() << " MS2 traces" << std::endl;

    return true;
  }

  ExitCodes main_(int, const char**) override
  {
    const String in = getStringOption_("in");
    const bool save_precursors = getFlag_("save_unfragmented_precursors");
    const Param ppim_params = getParam_().copy("PeakPickerIM:", true);
    const Param pphr_params = getParam_().copy("PeakPickerHiRes:", true);

    // FeatureFinderMetabo parameters (for MS1 and precursor data)
    const Param ffm_common_param = getParam_().copy("FeatureFinderMetabo:common:", true);
    Param ffm_mtd_param = getParam_().copy("FeatureFinderMetabo:mtd:", true);
    Param ffm_epd_param = getParam_().copy("FeatureFinderMetabo:epd:", true);
    Param ffm_ffm_param = getParam_().copy("FeatureFinderMetabo:ffm:", true);

    // MassTraceExtractor parameters (for MS1 and MS2 data)
    const Param mte_common_param = getParam_().copy("MassTraceExtractor:common:", true);
    Param mte_mtd_param = getParam_().copy("MassTraceExtractor:mtd:", true);
    Param mte_epd_param = getParam_().copy("MassTraceExtractor:epd:", true);

    // ClusterMassTraces parameters (for pseudo spectra generation)
    const Param cluster_param = getParam_().copy("ClusterMassTraces:", true);

#ifdef _OPENMP
    const int num_threads = getIntOption_("threads");
    const int threads_outer_loop = getIntOption_("threads_outer_loop");
#endif

    String out = getStringOption_("out");
    if (out.empty())
    {
      out = File::path(in) + "/" +
            File::basename(in) + "_diaWindows";
    }
    File::makeDir(out);

    // ------------------------------
    // Step 1: Use OnDiscMSExperiment for memory-efficient metadata access
    // This determines DIA windows and IM info without loading all peak data
    // ------------------------------
    OPENMS_LOG_INFO << "Opening file for metadata access..." << std::endl;
    OnDiscMSExperiment on_disc;
    if (!on_disc.openFile(in))
    {
      OPENMS_LOG_ERROR << "Failed to open file as indexed mzML." << std::endl;
      return INPUT_FILE_NOT_FOUND;
    }

    // Determine DIA windows from MS2 metadata (efficient - no peak data loaded)
    DiaWeaver::WindowMap windows;
    DiaWeaver::determineWindows(on_disc, windows);

    // Determine IM info (loads only representative spectra)
    DiaWeaver::IMInfo im_info = DiaWeaver::determineIMInfo(on_disc, windows);

    // Convert map to vector for OpenMP indexed access
    std::vector<std::pair<DiaWeaver::DIAWindow, std::vector<Size>>> window_vec(
      windows.begin(), windows.end());
    const Size total_windows = window_vec.size();

    // Note: OnDiscMSExperiment provides memory-efficient on-demand spectrum loading
    // Each thread will get its own copy via firstprivate (creates separate file handles)

    OPENMS_LOG_INFO << "Processing " << total_windows << " DIA windows";
#ifdef _OPENMP
    OPENMS_LOG_INFO << " using " << num_threads << " thread(s)";
#endif
    OPENMS_LOG_INFO << " with on-disc random access..." << std::endl;

    if (im_info.available)
    {
      OPENMS_LOG_INFO << "Ion mobility data detected. Using PeakPickerIM (mobilogram method)." << std::endl;
    }
    else
    {
      OPENMS_LOG_INFO << "No ion mobility data detected. Using PeakPickerHiRes." << std::endl;
    }

    OPENMS_LOG_INFO << "FeatureFinderMetabo will run on MS1 and precursor data." << std::endl;
    OPENMS_LOG_INFO << "MassTraceExtractor will run on MS2 data." << std::endl;
    OPENMS_LOG_INFO << "ClusterMassTraces will create pseudo spectra from MS1 features and MS2 traces." << std::endl;

    // ------------------------------
    // Step 3: Process windows in parallel with nested parallelism
    // Outer loop: over DIA windows
    // Inner loop: over spectra within each window (for peak picking)
    // ------------------------------
    Size processed = 0;

#ifdef _OPENMP
    // Store total number of threads available
    const int total_threads = num_threads;

    // Calculate outer and inner thread counts
    int outer_threads = total_threads;
    int inner_threads = 1;

    if (threads_outer_loop > 0)
    {
      // User specified nested parallelism
      outer_threads = std::min(threads_outer_loop, total_threads);
      inner_threads = std::max(1, total_threads / outer_threads);
      omp_set_nested(1);
      omp_set_dynamic(0);
      OPENMS_LOG_INFO << "Using nested parallelism: " << outer_threads << " outer threads x "
                      << inner_threads << " inner threads for peak picking." << std::endl;
    }
    else
    {
      OPENMS_LOG_INFO << "Using " << outer_threads << " threads for window processing (no nested parallelism)." << std::endl;
    }

    omp_set_num_threads(outer_threads);
#pragma omp parallel for schedule(dynamic, 1) firstprivate(on_disc)
#endif
    for (SignedSize idx = 0; idx < static_cast<SignedSize>(total_windows); ++idx)
    {
      const DiaWeaver::DIAWindow& w = window_vec[idx].first;
      const std::vector<Size>& indices = window_vec[idx].second;

      // Thread-local buffers
      MzMLFile mzml;
      MSExperiment ms2_exp;
      MSExperiment ms1_exp;
      MSExperiment precursor_exp;
      ConsensusMap ms2_traces;  // MS2 mass traces for clustering

      // Build human-readable filename: mz<lower>-<upper>_im<lower>-<upper>.mzML
      String fname_base = "mz" + String(static_cast<int>(std::round(w.lower_mz))) +
                          "-" + String(static_cast<int>(std::round(w.upper_mz)));
      if (w.hasIonMobility())
      {
        String im_lower = String::number(w.lower_im, 2);
        String im_upper = String::number(w.upper_im, 2);
        im_lower.substitute(".", "p");
        im_upper.substitute(".", "p");
        fname_base += "_im" + im_lower + "-" + im_upper;
      }
      fname_base += ".mzML";

      // Extract MS2 (on-demand from disk - each thread has its own file handle)
      DiaWeaver::extractSingleMS2Window(on_disc, w, indices, im_info, ms2_exp,
                                         save_precursors ? &precursor_exp : nullptr);

      // Apply peak picking to MS2 spectra (inner parallel loop)
      // Use parallel region to create thread-private pickers (more efficient than per-iteration)
#ifdef _OPENMP
#pragma omp parallel num_threads(inner_threads)
      {
        PeakPickerIM picker_im;
        PeakPickerHiRes picker_hr;
        if (im_info.available)
        {
          picker_im.setParameters(ppim_params);
        }
        else
        {
          picker_hr.setParameters(pphr_params);
        }

#pragma omp for schedule(dynamic, 1)
        for (SignedSize s = 0; s < static_cast<SignedSize>(ms2_exp.size()); ++s)
        {
          if (im_info.available)
          {
            picker_im.pickIMTraces(ms2_exp[s]);
          }
          else
          {
            MSSpectrum picked;
            picker_hr.pick(ms2_exp[s], picked);
            ms2_exp[s] = std::move(picked);
          }
        }
      }
#else
      for (SignedSize s = 0; s < static_cast<SignedSize>(ms2_exp.size()); ++s)
      {
        PeakPickerIM picker_im;
        PeakPickerHiRes picker_hr;
        if (im_info.available)
        {
          picker_im.setParameters(ppim_params);
          picker_im.pickIMTraces(ms2_exp[s]);
        }
        else
        {
          picker_hr.setParameters(pphr_params);
          MSSpectrum picked;
          picker_hr.pick(ms2_exp[s], picked);
          ms2_exp[s] = std::move(picked);
        }
      }
#endif

      if (!ms2_exp.empty())
      {
        mzml.store(out + "/ms2_" + fname_base, ms2_exp);

        // Run MassTraceExtractor on MS2 data
        // Create local copies of parameters since they get modified
        Param mte_mtd_copy = mte_mtd_param;
        Param mte_epd_copy = mte_epd_param;
        if (runMassTraceExtractor_(ms2_exp, mte_common_param, mte_mtd_copy, mte_epd_copy, ms2_traces))
        {
          if (!ms2_traces.empty())
          {
            String traces_fname = fname_base;
            traces_fname.substitute(".mzML", ".consensusXML");
            FileHandler().storeConsensusFeatures(out + "/ms2_traces_" + traces_fname, ms2_traces, {FileTypes::CONSENSUSXML});
          }
        }
      }

      // Apply peak picking to precursors (inner parallel loop)
      if (save_precursors && !precursor_exp.empty())
      {
#ifdef _OPENMP
#pragma omp parallel num_threads(inner_threads)
        {
          PeakPickerIM picker_im;
          PeakPickerHiRes picker_hr;
          if (im_info.available)
          {
            picker_im.setParameters(ppim_params);
          }
          else
          {
            picker_hr.setParameters(pphr_params);
          }

#pragma omp for schedule(dynamic, 1)
          for (SignedSize s = 0; s < static_cast<SignedSize>(precursor_exp.size()); ++s)
          {
            if (im_info.available)
            {
              picker_im.pickIMTraces(precursor_exp[s]);
            }
            else
            {
              MSSpectrum picked;
              picker_hr.pick(precursor_exp[s], picked);
              precursor_exp[s] = std::move(picked);
            }
          }
        }
#else
        for (SignedSize s = 0; s < static_cast<SignedSize>(precursor_exp.size()); ++s)
        {
          PeakPickerIM picker_im;
          PeakPickerHiRes picker_hr;
          if (im_info.available)
          {
            picker_im.setParameters(ppim_params);
            picker_im.pickIMTraces(precursor_exp[s]);
          }
          else
          {
            picker_hr.setParameters(pphr_params);
            MSSpectrum picked;
            picker_hr.pick(precursor_exp[s], picked);
            precursor_exp[s] = std::move(picked);
          }
        }
#endif
        mzml.store(out + "/precursor_" + fname_base, precursor_exp);

        // Run FeatureFinderMetabo on precursor data
        FeatureMap precursor_features;
        // Create local copies of parameters since they get modified
        Param mtd_copy = ffm_mtd_param;
        Param epd_copy = ffm_epd_param;
        Param ffm_copy = ffm_ffm_param;
        if (runFeatureFinderMetabo_(precursor_exp, ffm_common_param, mtd_copy, epd_copy, ffm_copy, precursor_features))
        {
          if (!precursor_features.empty())
          {
            precursor_features.ensureUniqueId();
            String feature_fname = fname_base;
            feature_fname.substitute(".mzML", ".featureXML");
            FileHandler().storeFeatures(out + "/precursor_features_" + feature_fname, precursor_features, {FileTypes::FEATUREXML});
          }
        }
      }

      // Extract MS1 (on-demand from disk - each thread has its own file handle)
      DiaWeaver::extractSingleMS1Window(on_disc, w, im_info, ms1_exp);

      // Apply peak picking to MS1 spectra (inner parallel loop)
#ifdef _OPENMP
#pragma omp parallel num_threads(inner_threads)
      {
        PeakPickerIM picker_im;
        PeakPickerHiRes picker_hr;
        if (im_info.available)
        {
          picker_im.setParameters(ppim_params);
        }
        else
        {
          picker_hr.setParameters(pphr_params);
        }

#pragma omp for schedule(dynamic, 1)
        for (SignedSize s = 0; s < static_cast<SignedSize>(ms1_exp.size()); ++s)
        {
          if (im_info.available)
          {
            picker_im.pickIMTraces(ms1_exp[s]);
          }
          else
          {
            MSSpectrum picked;
            picker_hr.pick(ms1_exp[s], picked);
            ms1_exp[s] = std::move(picked);
          }
        }
      }
#else
      for (SignedSize s = 0; s < static_cast<SignedSize>(ms1_exp.size()); ++s)
      {
        PeakPickerIM picker_im;
        PeakPickerHiRes picker_hr;
        if (im_info.available)
        {
          picker_im.setParameters(ppim_params);
          picker_im.pickIMTraces(ms1_exp[s]);
        }
        else
        {
          picker_hr.setParameters(pphr_params);
          MSSpectrum picked;
          picker_hr.pick(ms1_exp[s], picked);
          ms1_exp[s] = std::move(picked);
        }
      }
#endif

      if (!ms1_exp.empty())
      {
        mzml.store(out + "/ms1_" + fname_base, ms1_exp);

        // Run FeatureFinderMetabo on MS1 data to detect peptide features (monoisotopic)
        FeatureMap ms1_features;
        // Create local copies of parameters since they get modified
        Param mtd_copy = ffm_mtd_param;
        Param epd_copy = ffm_epd_param;
        Param ffm_copy = ffm_ffm_param;
        if (runFeatureFinderMetabo_(ms1_exp, ffm_common_param, mtd_copy, epd_copy, ffm_copy, ms1_features))
        {
          if (!ms1_features.empty())
          {
            ms1_features.ensureUniqueId();
            String feature_fname = fname_base;
            feature_fname.substitute(".mzML", ".featureXML");
            FileHandler().storeFeatures(out + "/ms1_features_" + feature_fname, ms1_features, {FileTypes::FEATUREXML});

            // Cluster MS1 features with MS2 traces to create pseudo spectra
            if (!ms2_traces.empty())
            {
              MSExperiment pseudo_spectra;
              if (clusterMassTraces_(ms1_features, ms2_traces, w.lower_mz, w.upper_mz, cluster_param, pseudo_spectra))
              {
                if (!pseudo_spectra.empty())
                {
                  String pseudo_fname = fname_base;
                  pseudo_fname.substitute(".mzML", "_pseudo.mzML");
                  MzMLFile().store(out + "/pseudo_spectra_" + pseudo_fname, pseudo_spectra);
                }
              }
            }
          }
        }
      }

      // Progress logging (thread-safe)
#ifdef _OPENMP
#pragma omp critical (progress_log)
#endif
      {
        ++processed;
        OPENMS_LOG_INFO << "Processed window " << processed << "/" << total_windows
                        << " (m/z: " << w.lower_mz << "-" << w.upper_mz << ")" << std::endl;
      }
    }

#ifdef _OPENMP
    // Restore total thread count if nested parallelism was used
    if (threads_outer_loop > 0)
    {
      omp_set_num_threads(total_threads);
    }
#endif

    OPENMS_LOG_INFO << "Finished processing all windows." << std::endl;

    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPDiaWeaver tool;
  return tool.main(argc, argv);
}
