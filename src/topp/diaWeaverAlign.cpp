// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Mohammed Alhigaylan $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/APPLICATIONS/diaWeaver.h>
#include <OpenMS/APPLICATIONS/diaWeaverAlign.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <fstream>
#include <unordered_map>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace OpenMS;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_diaWeaverAlign diaWeaverAlign

@brief Aligns diaWeaver pseudo-MS2 spectra across multiple runs against a
peak-picked experimental DIA file.

For each DIA window detected in the query file the tool:
1. Collects pseudo-MS2 spectra whose precursor m/z falls within that window
   from every input file.
2. Builds a fragment index (DiaWeaverAlign) and a precursor index.
3. Scores every query MS2 spectrum in the window against the fragment index
   and accumulates per-spectrum-id score traces.
4. Groups spectrum-ids with similar precursor m/z, close apex RT, and
   sufficient fragment overlap into FeatureGroups across source files.

Output is a tab-separated file with one row per feature-group member.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_diaWeaverAlign.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_diaWeaverAlign.html
*/

class TOPPDiaWeaverAlign :
  public TOPPBase
{
public:
  TOPPDiaWeaverAlign() :
    TOPPBase(
      "diaWeaverAlign",
      "Align diaWeaver pseudo-MS2 spectra across runs using fragment-index scoring",
      false)
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFileList_(
      "in",
      "<files>",
      {},
      "Pseudo-MS2 mzML files produced by diaWeaver (one per replicate run). "
      "Each file may contain spectra from all DIA windows.",
      true);
    setValidFormats_("in", {"mzML"});

    registerInputFile_(
      "in_query",
      "<file>",
      "",
      "Peak-picked experimental DIA mzML file. Its MS2 spectra are scored "
      "against the fragment index built from '-in'.",
      true);
    setValidFormats_("in_query", {"mzML"});

    registerOutputFile_(
      "out",
      "<file>",
      "",
      "Output TSV with one row per feature-group member.",
      true);

    registerDoubleOption_(
      "rt_tolerance",
      "<seconds>",
      10.0,
      "Maximum apex RT difference (seconds) between two spectrum-ids to be "
      "considered for grouping.",
      false);
    setMinFloat_("rt_tolerance", 0.0);

    registerIntOption_(
      "min_fragment_overlap",
      "<n>",
      5,
      "Minimum number of shared fragment bins required to group two spectrum-ids.",
      false);
    setMinInt_("min_fragment_overlap", 1);

    registerIntOption_(
      "min_matched_peaks",
      "<n>",
      2,
      "Minimum matched-peaks count for a query spectrum to contribute to a score trace.",
      false);
    setMinInt_("min_matched_peaks", 1);

    registerSubsection_(
      "DiaWeaverAlign",
      "Parameters for the DiaWeaverAlign fragment-index algorithm");

    registerIntOption_(
      "threads",
      "<n>",
      1,
      "Number of threads used by matchExperiment.",
      false);
    setMinInt_("threads", 1);
  }

  Param getSubsectionDefaults_(const String& section) const override
  {
    if (section == "DiaWeaverAlign")
    {
      DiaWeaverAlign aligner;
      return aligner.getDefaults();
    }
    return Param();
  }

  ExitCodes main_(int /*argc*/, const char** /*argv*/) override
  {
    // ------------------------------------------------------------------
    // Read options
    // ------------------------------------------------------------------
    const StringList in_files    = getStringList_("in");
    const String     in_query    = getStringOption_("in_query");
    const String     out_file    = getStringOption_("out");
    const double     rt_tol      = getDoubleOption_("rt_tolerance");
    const uint32_t   min_frags   = static_cast<uint32_t>(getIntOption_("min_fragment_overlap"));
    const uint32_t   min_peaks   = static_cast<uint32_t>(getIntOption_("min_matched_peaks"));

#ifdef _OPENMP
    omp_set_num_threads(getIntOption_("threads"));
#endif

    // ------------------------------------------------------------------
    // Load query experiment and detect DIA windows
    // ------------------------------------------------------------------
    OPENMS_LOG_INFO << "Loading query: " << in_query << std::endl;
    MSExperiment query_exp;
    MzMLFile().load(in_query, query_exp);

    // The query peak-picked file has MS level 1 (set by PeakPickerIM), so we
    // use the local helper rather than DiaWeaver::determineWindows which
    // would skip all spectra with getMSLevel() != 2.
    DiaWeaver::WindowMap query_windows;
    determineWindowsPeakPicked_(query_exp, query_windows);
    OPENMS_LOG_INFO << "Detected " << query_windows.size()
                    << " DIA window(s) in query." << std::endl;

    // ------------------------------------------------------------------
    // Load all pseudo-MS2 input files
    // ------------------------------------------------------------------
    OPENMS_LOG_INFO << "Loading " << in_files.size() << " pseudo-MS2 file(s)..." << std::endl;
    std::vector<MSExperiment> pseudo_exps(in_files.size());
    for (Size i = 0; i < in_files.size(); ++i)
    {
      OPENMS_LOG_INFO << "  [" << (i + 1) << "/" << in_files.size() << "] "
                      << in_files[i] << std::endl;
      MzMLFile().load(in_files[i], pseudo_exps[i]);
    }

    // ------------------------------------------------------------------
    // Set up DiaWeaverAlign with user parameters
    // ------------------------------------------------------------------
    DiaWeaverAlign aligner;
    Param aligner_param = getParam_().copy("DiaWeaverAlign:", true);
    aligner.setParameters(aligner_param);

    // ------------------------------------------------------------------
    // Open output TSV
    // ------------------------------------------------------------------
    std::ofstream tsv(out_file.c_str());
    if (!tsv.is_open())
    {
      OPENMS_LOG_ERROR << "Cannot open output file: " << out_file << std::endl;
      return CANNOT_WRITE_OUTPUT_FILE;
    }

    tsv << "group_id\tn_members\tmean_apex_rt_s\tn_shared_frags"
        << "\twindow_lower_mz\twindow_upper_mz"
        << "\tspectrum_id\tsource_file\tnative_id"
        << "\tprecursor_mz\tprecursor_charge"
        << "\tpseudo_rt_s\tquery_apex_rt_s\n";

    // ------------------------------------------------------------------
    // Process windows
    // ------------------------------------------------------------------
    uint32_t global_group_id = 0;

    // If no DIA windows were detected (e.g. non-DIA query), fall back to
    // global processing using an empty window sentinel.
    if (query_windows.empty())
    {
      OPENMS_LOG_WARN << "No DIA windows detected; processing all spectra globally." << std::endl;

      aligner.buildIndex(pseudo_exps, in_files);
      aligner.buildPrecursorIndex();

      auto traces = aligner.matchExperiment(query_exp, min_peaks);
      auto groups = aligner.groupFeatures(traces, rt_tol, min_frags);

      writeGroups_(tsv, groups, traces, aligner, global_group_id,
                   0.0, 0.0);  // window bounds unknown
      return EXECUTION_OK;
    }

    for (const auto& [window, q_indices] : query_windows)
    {
      OPENMS_LOG_INFO << "Window m/z [" << window.lower_mz << ", " << window.upper_mz << "]"
                      << "  query spectra: " << q_indices.size() << std::endl;

      // ---- Collect query spectra for this window ----
      // q_indices already contains only spectra belonging to this window
      // (as determined by determineWindowsPeakPicked_). They are MS level 1
      // peak-picked MS2 spectra; no additional MS-level filter is needed.
      MSExperiment window_query;
      for (Size idx : q_indices)
        window_query.addSpectrum(query_exp[idx]);

      if (window_query.empty())
      {
        OPENMS_LOG_WARN << "  No MS2 spectra for window ["
                        << window.lower_mz << ", " << window.upper_mz << "]. Skipping." << std::endl;
        continue;
      }

      // ---- Partition pseudo spectra by precursor m/z within this window ----
      std::vector<MSExperiment> window_pseudo(in_files.size());
      for (Size i = 0; i < pseudo_exps.size(); ++i)
      {
        for (const MSSpectrum& spec : pseudo_exps[i])
        {
          if (spec.getPrecursors().empty()) continue;
          const double prec_mz = spec.getPrecursors()[0].getMZ();
          if (prec_mz >= window.lower_mz && prec_mz <= window.upper_mz)
            window_pseudo[i].addSpectrum(spec);
        }
      }

      // Count how many pseudo spectra we have for this window
      Size total_pseudo = 0;
      for (const auto& exp : window_pseudo)
        total_pseudo += exp.size();

      OPENMS_LOG_INFO << "  Pseudo spectra for this window: " << total_pseudo << std::endl;

      if (total_pseudo == 0)
      {
        OPENMS_LOG_WARN << "  No pseudo spectra for window ["
                        << window.lower_mz << ", " << window.upper_mz << "]. Skipping." << std::endl;
        continue;
      }

      // ---- Build index, match, group ----
      aligner.buildIndex(window_pseudo, in_files);
      aligner.buildPrecursorIndex();

      const auto traces = aligner.matchExperiment(window_query, min_peaks);

      OPENMS_LOG_INFO << "  Score traces above threshold: " << traces.size() << std::endl;

      const auto groups = aligner.groupFeatures(traces, rt_tol, min_frags);

      OPENMS_LOG_INFO << "  Feature groups (>=2 members): " << groups.size() << std::endl;

      writeGroups_(tsv, groups, traces, aligner, global_group_id,
                   window.lower_mz, window.upper_mz);
    }

    OPENMS_LOG_INFO << "Written " << global_group_id << " feature group(s) to " << out_file << std::endl;

    return EXECUTION_OK;
  }

private:

  /// Detect DIA windows from a peak-picked DIA experiment whose spectra carry
  /// MS level 1 (set by PeakPickerIM).  Works identically to
  /// DiaWeaver::determineWindows() but accepts any MS level and selects
  /// spectra by precursor presence instead.
  static void determineWindowsPeakPicked_(
    const MSExperiment& exp,
    DiaWeaver::WindowMap& window_map)
  {
    window_map.clear();
    std::vector<DiaWeaver::DIAWindow> known_windows;

    for (Size i = 0; i < exp.size(); ++i)
    {
      const MSSpectrum& spec = exp[i];
      if (spec.getPrecursors().empty()) continue;   // skip genuine MS1 scans

      const Precursor& p = spec.getPrecursors()[0];
      if (p.getIsolationWindowLowerOffset() == 0 &&
          p.getIsolationWindowUpperOffset() == 0)
        continue;   // no isolation window info — not a DIA spectrum

      DiaWeaver::DIAWindow candidate;
      candidate.center_mz = p.getMZ();
      candidate.lower_mz  = p.getMZ() - p.getIsolationWindowLowerOffset();
      candidate.upper_mz  = p.getMZ() + p.getIsolationWindowUpperOffset();

      if (spec.metaValueExists("ion mobility lower limit"))
        candidate.lower_im = (double)spec.getMetaValue("ion mobility lower limit");
      if (spec.metaValueExists("ion mobility upper limit"))
        candidate.upper_im = (double)spec.getMetaValue("ion mobility upper limit");

      bool found = false;
      for (const DiaWeaver::DIAWindow& known : known_windows)
      {
        if (candidate.isEqual(known))
        {
          window_map[known].push_back(i);
          found = true;
          break;
        }
      }
      if (!found)
      {
        known_windows.push_back(candidate);
        window_map[candidate].push_back(i);
      }
    }
  }

  /// Write feature groups for one window to the TSV.
  static void writeGroups_(
    std::ofstream& tsv,
    const std::vector<DiaWeaverAlign::FeatureGroup>& groups,
    const std::vector<DiaWeaverAlign::ScoreTrace>& traces,
    const DiaWeaverAlign& aligner,
    uint32_t& global_group_id,
    double window_lower_mz,
    double window_upper_mz)
  {
    // Build spectrum_id → apex_rt lookup from score traces.
    std::unordered_map<uint32_t, double> apex_rt_map;
    apex_rt_map.reserve(traces.size());
    for (const DiaWeaverAlign::ScoreTrace& t : traces)
      apex_rt_map[t.spectrum_id] = t.apex_rt;

    for (const DiaWeaverAlign::FeatureGroup& group : groups)
    {
      for (uint32_t sid : group.spectrum_ids)
      {
        const DiaWeaverAlign::SpectrumEntry& se = aligner.getSpectrumEntry(sid);
        const auto it = apex_rt_map.find(sid);
        const double query_apex_rt = (it != apex_rt_map.end()) ? it->second : -1.0;

        tsv << global_group_id
            << '\t' << group.spectrum_ids.size()
            << '\t' << group.mean_apex_rt
            << '\t' << group.shared_fragments.size()
            << '\t' << window_lower_mz
            << '\t' << window_upper_mz
            << '\t' << sid
            << '\t' << aligner.getSourceFile(se.source_file_idx)
            << '\t' << se.native_id
            << '\t' << se.precursor_mz
            << '\t' << se.precursor_charge
            << '\t' << se.retention_time
            << '\t' << query_apex_rt
            << '\n';
      }
      ++global_group_id;
    }
  }
};


int main(int argc, const char** argv)
{
  TOPPDiaWeaverAlign tool;
  return tool.main(argc, argv);
}

/// @cond TOPP_DOCS
