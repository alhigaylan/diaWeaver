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

#include <cstdio>
#include <filesystem>
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
peak-picked experimental DIA file using a streaming index architecture.

The tool operates in three phases:

**Phase 1 – Streaming accumulation.**
The peak-picked query mzML is loaded once to detect DIA isolation windows.
One per-window raw index file pair is opened in a working directory.  Each
pseudo-MS2 input file is then loaded one at a time; every spectrum is
dispatched to the raw index file of the matching window and the file is
unloaded.  Peak RAM during this phase is dominated by one loaded pseudo file.

**Phase 2 – Finalization.**
For each window the raw index files are read back, sorted, and assembled
into a complete fragment+precursor CSR index which is saved as a single
.dwaindex binary.  The raw files are deleted immediately afterwards.

**Phase 3 – Alignment.**
Each .dwaindex is loaded in turn, query spectra for that window are
extracted, fragment scoring and feature grouping are performed, and results
are written to the output TSV.

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
      "Align diaWeaver pseudo-MS2 spectra across runs using a streaming fragment-index",
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
      "Peak-picked experimental DIA mzML file (MS level 1, produced by "
      "PeakPickerIM). Its MS2 spectra are scored against the fragment index "
      "built from '-in'.",
      true);
    setValidFormats_("in_query", {"mzML"});

    registerOutputFile_(
      "out",
      "<file>",
      "",
      "Output TSV with one row per feature-group member.",
      true);

    registerStringOption_(
      "workdir",
      "<directory>",
      "",
      "Working directory for intermediate per-window index files (.dwaf, .dwam, "
      ".dwaindex). The directory is created if it does not exist. Intermediate "
      "files are deleted after use; only .dwaindex files remain.",
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
    const StringList in_files  = getStringList_("in");
    const String     in_query  = getStringOption_("in_query");
    const String     out_file  = getStringOption_("out");
    const String     workdir   = getStringOption_("workdir");
    const double     rt_tol    = getDoubleOption_("rt_tolerance");
    const uint32_t   min_frags = static_cast<uint32_t>(getIntOption_("min_fragment_overlap"));
    const uint32_t   min_peaks = static_cast<uint32_t>(getIntOption_("min_matched_peaks"));

#ifdef _OPENMP
    omp_set_num_threads(getIntOption_("threads"));
#endif

    // Create workdir if needed
    {
      namespace fs = std::filesystem;
      std::error_code ec;
      fs::create_directories(workdir.c_str(), ec);
      if (ec)
      {
        OPENMS_LOG_ERROR << "Cannot create workdir '" << workdir
                         << "': " << ec.message() << std::endl;
        return CANNOT_WRITE_OUTPUT_FILE;
      }
    }

    // ------------------------------------------------------------------
    // Load query experiment and detect DIA windows
    // ------------------------------------------------------------------
    OPENMS_LOG_INFO << "Loading query: " << in_query << std::endl;
    MSExperiment query_exp;
    MzMLFile().load(in_query, query_exp);

    // The query peak-picked file has MS level 1 (set by PeakPickerIM), so we
    // use the local helper rather than DiaWeaver::determineWindows which
    // would skip all spectra with getMSLevel() != 2.
    DiaWeaver::WindowMap query_windows_map;
    determineWindowsPeakPicked_(query_exp, query_windows_map);
    OPENMS_LOG_INFO << "Detected " << query_windows_map.size()
                    << " DIA window(s) in query." << std::endl;

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

    uint32_t global_group_id = 0;

    Param aligner_param = getParam_().copy("DiaWeaverAlign:", true);

    // ------------------------------------------------------------------
    // No-window fallback: process all spectra globally using in-memory index.
    // This covers non-DIA or single-window query files.
    // ------------------------------------------------------------------
    if (query_windows_map.empty())
    {
      OPENMS_LOG_WARN << "No DIA windows detected; processing all spectra globally." << std::endl;

      OPENMS_LOG_INFO << "Loading " << in_files.size() << " pseudo-MS2 file(s)..." << std::endl;
      std::vector<MSExperiment> pseudo_exps(in_files.size());
      for (Size i = 0; i < in_files.size(); ++i)
        MzMLFile().load(in_files[i], pseudo_exps[i]);

      DiaWeaverAlign aligner;
      aligner.setParameters(aligner_param);
      aligner.buildIndex(pseudo_exps, in_files);
      aligner.buildPrecursorIndex();

      const auto traces = aligner.matchExperiment(query_exp, min_peaks);
      const auto groups = aligner.groupFeatures(traces, rt_tol, min_frags);
      writeGroups_(tsv, groups, traces, aligner, global_group_id, 0.0, 0.0);
      return EXECUTION_OK;
    }

    // Build parallel window / query-index vectors from the map
    // (preserves the insertion order guaranteed by the underlying std::map-like type)
    std::vector<DiaWeaver::DIAWindow> windows;
    std::vector<std::vector<Size>>    window_q_indices;
    windows.reserve(query_windows_map.size());
    window_q_indices.reserve(query_windows_map.size());
    for (const auto& [win, q_idx] : query_windows_map)
    {
      windows.push_back(win);
      window_q_indices.push_back(q_idx);
    }
    const Size n_windows = windows.size();

    // ================================================================
    // PHASE 1: Stream pseudo spectra into per-window raw index files.
    //
    // One DiaWeaverAlign per window keeps a FILE* pair open.
    // Pseudo files are loaded one at a time and freed after streaming.
    // ================================================================
    OPENMS_LOG_INFO << "Phase 1: streaming " << in_files.size()
                    << " pseudo file(s) into " << n_windows << " window index(es)." << std::endl;

    std::vector<DiaWeaverAlign> aligners(n_windows);
    for (Size w = 0; w < n_windows; ++w)
    {
      aligners[w].setParameters(aligner_param);
      aligners[w].openStream(
        workdir + "/window_" + String(w) + ".dwaf",
        workdir + "/window_" + String(w) + ".dwam");
    }

    for (Size file_idx = 0; file_idx < in_files.size(); ++file_idx)
    {
      OPENMS_LOG_INFO << "  [" << (file_idx + 1) << "/" << in_files.size() << "] "
                      << in_files[file_idx] << std::endl;

      MSExperiment pseudo_exp;
      MzMLFile().load(in_files[file_idx], pseudo_exp);

      for (const MSSpectrum& spec : pseudo_exp)
      {
        if (spec.getMSLevel() != 2 || spec.empty() || spec.getPrecursors().empty()) continue;
        const double prec_mz = spec.getPrecursors()[0].getMZ();

        for (Size w = 0; w < n_windows; ++w)
        {
          if (prec_mz >= windows[w].lower_mz && prec_mz <= windows[w].upper_mz)
          {
            aligners[w].appendSpectrumToStream(spec, file_idx);
            break;
          }
        }
      }
      // pseudo_exp destructs here, freeing peak data
    }

    for (Size w = 0; w < n_windows; ++w)
      aligners[w].closeStream();

    // ================================================================
    // PHASE 2: Finalize each window's index and save to a .dwaindex file.
    //
    // Processing is serial: one window's sort+CSR lives in RAM at a time.
    // Raw streaming files are deleted after the index is saved.
    // ================================================================
    OPENMS_LOG_INFO << "Phase 2: finalizing " << n_windows << " window index(es)." << std::endl;

    for (Size w = 0; w < n_windows; ++w)
    {
      const String frags_path = workdir + "/window_" + String(w) + ".dwaf";
      const String meta_path  = workdir + "/window_" + String(w) + ".dwam";
      const String index_path = workdir + "/window_" + String(w) + ".dwaindex";

      OPENMS_LOG_INFO << "  Window " << w
                      << " [" << windows[w].lower_mz << ", " << windows[w].upper_mz << "]"
                      << std::endl;

      aligners[w].finalizeFromStream(frags_path, meta_path, in_files);
      aligners[w].buildPrecursorIndex();
      aligners[w].saveIndex(index_path);

      // Delete raw streaming files immediately to free disk space
      std::remove(frags_path.c_str());
      std::remove(meta_path.c_str());

      // Free in-memory index (replace with empty default object)
      aligners[w] = DiaWeaverAlign();
    }

    // ================================================================
    // PHASE 3: Load each index, score query spectra, group, write TSV.
    //
    // Only one .dwaindex is in RAM at a time.
    // ================================================================
    OPENMS_LOG_INFO << "Phase 3: aligning " << n_windows << " window(s)." << std::endl;

    for (Size w = 0; w < n_windows; ++w)
    {
      const String index_path = workdir + "/window_" + String(w) + ".dwaindex";
      const DiaWeaver::DIAWindow& window   = windows[w];
      const std::vector<Size>&    q_idxs   = window_q_indices[w];

      DiaWeaverAlign aligner;
      aligner.loadIndex(index_path);

      // Collect query spectra for this window
      MSExperiment window_query;
      for (Size idx : q_idxs)
        window_query.addSpectrum(query_exp[idx]);

      OPENMS_LOG_INFO << "  Window " << w
                      << " [" << window.lower_mz << ", " << window.upper_mz << "]"
                      << "  pseudo=" << aligner.getSpectrumCount()
                      << "  query=" << window_query.size() << std::endl;

      if (window_query.empty() || aligner.getSpectrumCount() == 0)
      {
        OPENMS_LOG_WARN << "    Skipping (empty window)." << std::endl;
        continue;
      }

      const auto traces = aligner.matchExperiment(window_query, min_peaks);
      OPENMS_LOG_INFO << "    Score traces: " << traces.size() << std::endl;

      const auto groups = aligner.groupFeatures(traces, rt_tol, min_frags);
      OPENMS_LOG_INFO << "    Feature groups (>=2 members): " << groups.size() << std::endl;

      writeGroups_(tsv, groups, traces, aligner, global_group_id,
                   window.lower_mz, window.upper_mz);
    }

    OPENMS_LOG_INFO << "Written " << global_group_id
                    << " feature group(s) to " << out_file << std::endl;

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
