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
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
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

    registerDoubleOption_(
      "im_tolerance",
      "<Vs/cm^2>",
      0.01,
      "Maximum drift time difference (Vs/cm^2) between two spectrum-ids to be "
      "considered the same precursor identity. Only applied when IM data is present.",
      false);
    setMinFloat_("im_tolerance", 0.0);

    registerDoubleOption_(
      "precursor_ppm_tolerance",
      "<ppm>",
      20.0,
      "Maximum precursor m/z difference (ppm) between two spectrum-ids to be "
      "considered for grouping.",
      false);
    setMinFloat_("precursor_ppm_tolerance", 0.0);

    registerIntOption_(
      "min_fragment_overlap",
      "<n>",
      5,
      "Minimum number of shared fragment bins required to group two spectrum-ids.",
      false);
    setMinInt_("min_fragment_overlap", 1);

    registerIntOption_(
      "isotope_error_tol",
      "<n>",
      3,
      "Maximum number of isotope errors to tolerate when matching precursor m/z across runs. "
      "For each pair, the m/z difference is checked against k * 1.003355 / charge for "
      "k in {0, ±1, …, ±isotope_error_tol}. Set to 0 to disable isotope correction.",
      false);
    setMinInt_("isotope_error_tol", 0);

    registerDoubleOption_(
      "min_overlap_similarity",
      "<fraction>",
      0.7,
      "Minimum Overlap Coefficient (|A∩B| / min(|A|,|B|) of apex fragment bins) required "
      "for two members to stay in the same sub-group during complete-linkage re-clustering. "
      "Using the smaller set as denominator removes the bias from same-run fingerprint inflation. "
      "Re-clustering is only triggered for groups that exceed the number of source files and "
      "show heterogeneity in m/z, ion mobility, and RT simultaneously.",
      false);
    setMinFloat_("min_overlap_similarity", 0.0);
    setMaxFloat_("min_overlap_similarity", 1.0);

    registerDoubleOption_(
      "min_within_file_jaccard",
      "<fraction>",
      0.7,
      "Minimum Jaccard similarity (|A∩B| / |A∪B|) required for two spectrum_ids from the "
      "same source file to be merged into a single file-cluster before cross-file Overlap "
      "Coefficient clustering. Spectrum_ids below this threshold are treated as different "
      "compounds and kept as separate file-clusters.",
      false);
    setMinFloat_("min_within_file_jaccard", 0.0);
    setMaxFloat_("min_within_file_jaccard", 1.0);

    registerIntOption_(
      "singleton_min_frags",
      "<n>",
      50,
      "Minimum apex matched-peak count for a spectrum_id with no cross-file partner to be "
      "retained as a singleton group. Set to 0 to discard all singletons.",
      false);
    setMinInt_("singleton_min_frags", 0);

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
    const double     im_tol    = getDoubleOption_("im_tolerance");
    const double     prec_ppm  = getDoubleOption_("precursor_ppm_tolerance");
    const uint32_t   min_frags    = static_cast<uint32_t>(getIntOption_("min_fragment_overlap"));
    const uint32_t   iso_err_tol  = static_cast<uint32_t>(getIntOption_("isotope_error_tol"));
    const uint32_t   min_peaks    = static_cast<uint32_t>(getIntOption_("min_matched_peaks"));
    const double     min_overlap         = getDoubleOption_("min_overlap_similarity");
    const double     min_file_jaccard    = getDoubleOption_("min_within_file_jaccard");
    const uint32_t   singleton_min_frags = static_cast<uint32_t>(getIntOption_("singleton_min_frags"));

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

    tsv << "group_id\tn_members\tmean_apex_rt_s\tn_shared_frags\tn_quant_frags\tmin_internal_overlap"
        << "\twindow_lower_mz\twindow_upper_mz\twindow_lower_im\twindow_upper_im"
        << "\tdecoy_window_lower_mz\tdecoy_window_upper_mz\tdecoy_window_lower_im\tdecoy_window_upper_im"
        << "\tspectrum_id\tsource_file\tnative_id"
        << "\tprecursor_mz\tprecursor_charge\tprecursor_drift_time"
        << "\tpseudo_rt_s\tquery_apex_rt_s\tinput_fragments\tapex_matched_peaks\tbest_decoy\tfold_changes\n";

    uint32_t global_group_id = 0;

    Param aligner_param = getParam_().copy("DiaWeaverAlign:", true);

    // ------------------------------------------------------------------
    // Require DIA windows. Without them we cannot verify that query precursors
    // belong to the correct isolation window, making any score meaningless.
    // ------------------------------------------------------------------
    if (query_windows_map.empty())
    {
      OPENMS_LOG_ERROR << "No DIA isolation windows detected in '" << in_query << "'.\n"
                       << "Every MS2 spectrum in the query file must carry a precursor "
                       << "with non-zero isolation window offsets. "
                       << "Ensure the file was produced by PeakPickerIM and retains "
                       << "the original isolation window metadata." << std::endl;
      return ILLEGAL_PARAMETERS;
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
        const double drift_t = spec.getDriftTime();
        const bool   spec_has_im = (drift_t != IMTypes::DRIFTTIME_NOT_SET);

        for (Size w = 0; w < n_windows; ++w)
        {
          // Precursor m/z must fall within the window's isolation range.
          if (prec_mz < windows[w].lower_mz || prec_mz > windows[w].upper_mz) continue;

          // If both the window and the spectrum carry ion mobility data,
          // the precursor drift time must also fall within the window's IM range.
          if (windows[w].hasIonMobility() && spec_has_im)
          {
            if (drift_t < windows[w].lower_im || drift_t > windows[w].upper_im) continue;
          }

          aligners[w].appendSpectrumToStream(spec, file_idx);
          break;
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
      aligners[w].setWindowBounds(
        windows[w].lower_mz, windows[w].upper_mz,
        windows[w].lower_im, windows[w].upper_im);
      aligners[w].saveIndex(index_path);

      // Delete raw streaming files immediately to free disk space
      std::remove(frags_path.c_str());
      std::remove(meta_path.c_str());

      // Free in-memory index (replace with empty default object)
      aligners[w] = DiaWeaverAlign();
    }

    // ================================================================
    // Pre-compute decoy window pairings.
    //
    // For each target window, find the non-adjacent window with the most
    // similar mean peak count (busyness).  This pairing is used in Phase 3
    // to estimate a per-spectrum-id FDR proxy ("best_decoy").
    // ================================================================
    std::vector<double> mean_peak_counts(n_windows);
    for (Size w = 0; w < n_windows; ++w)
      mean_peak_counts[w] = computeMeanPeakCount_(query_exp, window_q_indices[w]);

    std::vector<Size> decoy_window_idx(n_windows);
    OPENMS_LOG_INFO << "Decoy window pairings:" << std::endl;
    for (Size w = 0; w < n_windows; ++w)
    {
      decoy_window_idx[w] = findDecoyWindow_(windows, mean_peak_counts, w);
      if (decoy_window_idx[w] == std::numeric_limits<Size>::max())
      {
        OPENMS_LOG_WARN << "  Window " << w
                        << " [" << windows[w].lower_mz << ", " << windows[w].upper_mz << "]"
                        << ": no valid (non-adjacent) decoy window found." << std::endl;
      }
      else
      {
        const Size d = decoy_window_idx[w];
        double im_overlap = 0.0;
        const double target_im_width = windows[w].upper_im - windows[w].lower_im;
        if (windows[w].hasIonMobility() && windows[d].hasIonMobility() && target_im_width > 0.0)
        {
          const double raw_overlap = std::max(0.0,
            std::min(windows[w].upper_im, windows[d].upper_im) -
            std::max(windows[w].lower_im, windows[d].lower_im));
          im_overlap = raw_overlap / target_im_width;
        }
        OPENMS_LOG_INFO << "  Window " << w
                        << " [" << windows[w].lower_mz << ", " << windows[w].upper_mz << "]"
                        << " -> decoy window " << d
                        << " [" << windows[d].lower_mz << ", " << windows[d].upper_mz << "]"
                        << "  (im_overlap=" << im_overlap
                        << "  mean_peaks: " << mean_peak_counts[w]
                        << " vs " << mean_peak_counts[d] << ")" << std::endl;
      }
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

      // Sanity-check: the loaded index should describe the same window we expect.
      // A mismatch means the workdir contains stale files from a different run.
      if (std::abs(aligner.getWindowLowerMz() - window.lower_mz) > 1e-6 ||
          std::abs(aligner.getWindowUpperMz() - window.upper_mz) > 1e-6)
      {
        OPENMS_LOG_ERROR << "  Window " << w << ": index file '" << index_path
                         << "' covers precursor range ["
                         << aligner.getWindowLowerMz() << ", "
                         << aligner.getWindowUpperMz() << "] but expected ["
                         << window.lower_mz << ", " << window.upper_mz << "]. "
                         << "Delete workdir and re-run." << std::endl;
        return ILLEGAL_PARAMETERS;
      }

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

      // ------------------------------------------------------------------
      // Decoy scoring: run the same target index against query spectra from
      // a non-adjacent window with similar busyness.  Any match is a false
      // positive by construction because the query precursors do not belong
      // to this isolation window.
      // ------------------------------------------------------------------
      std::unordered_map<uint32_t, uint32_t> decoy_apex_map;
      const Size d_idx = decoy_window_idx[w];
      if (d_idx != std::numeric_limits<Size>::max())
      {
        MSExperiment decoy_query;
        for (Size idx : window_q_indices[d_idx])
          decoy_query.addSpectrum(query_exp[idx]);

        if (!decoy_query.empty())
        {
          // min_matched_peaks=1 to capture every non-zero decoy hit
          const auto decoy_traces = aligner.matchExperiment(decoy_query, 1);
          for (const DiaWeaverAlign::ScoreTrace& dt : decoy_traces)
            decoy_apex_map[dt.spectrum_id] = dt.apex_score;
          OPENMS_LOG_INFO << "    Decoy traces: " << decoy_traces.size() << std::endl;
        }
      }

      const auto groups = aligner.groupFeatures(traces, rt_tol, min_frags, im_tol, prec_ppm, iso_err_tol, min_overlap, min_file_jaccard, singleton_min_frags);
      OPENMS_LOG_INFO << "    Feature groups (>=2 members): " << groups.size() << std::endl;

      const DiaWeaver::DIAWindow* decoy_win =
          (d_idx != std::numeric_limits<Size>::max()) ? &windows[d_idx] : nullptr;
      writeGroups_(tsv, groups, traces, decoy_apex_map, aligner, global_group_id,
                   window, decoy_win);
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
  /// @param decoy_win  Pointer to the chosen decoy DIAWindow; nullptr if none available.
  static void writeGroups_(
    std::ofstream& tsv,
    const std::vector<DiaWeaverAlign::FeatureGroup>& groups,
    const std::vector<DiaWeaverAlign::ScoreTrace>& traces,
    const std::unordered_map<uint32_t, uint32_t>& decoy_apex_map,
    const DiaWeaverAlign& aligner,
    uint32_t& global_group_id,
    const DiaWeaver::DIAWindow& window,
    const DiaWeaver::DIAWindow* decoy_win)
  {
    // Helper: format an IM bound as string, or "NA" when not set.
    auto im_str = [](double im) -> std::string {
      return (im > DiaWeaver::IM_NOT_SET + 1e-9) ? std::to_string(im) : "NA";
    };

    // Pre-format window IM columns (same for every row in this window).
    const std::string win_lo_im = im_str(window.lower_im);
    const std::string win_hi_im = im_str(window.upper_im);

    const std::string decoy_lo_mz = decoy_win ? std::to_string(decoy_win->lower_mz) : "NA";
    const std::string decoy_hi_mz = decoy_win ? std::to_string(decoy_win->upper_mz) : "NA";
    const std::string decoy_lo_im = decoy_win ? im_str(decoy_win->lower_im) : "NA";
    const std::string decoy_hi_im = decoy_win ? im_str(decoy_win->upper_im) : "NA";

    // Build spectrum_id → ScoreTrace pointer lookup.
    std::unordered_map<uint32_t, const DiaWeaverAlign::ScoreTrace*> trace_map;
    trace_map.reserve(traces.size());
    for (const DiaWeaverAlign::ScoreTrace& t : traces)
      trace_map[t.spectrum_id] = &t;

    for (const DiaWeaverAlign::FeatureGroup& group : groups)
    {
      for (uint32_t sid : group.spectrum_ids)
      {
        const DiaWeaverAlign::SpectrumEntry& se = aligner.getSpectrumEntry(sid);
        const auto it = trace_map.find(sid);
        const double   query_apex_rt      = (it != trace_map.end()) ? it->second->apex_rt    : -1.0;
        const uint32_t apex_matched_peaks = (it != trace_map.end()) ? it->second->apex_score : 0;

        const auto decoy_it = decoy_apex_map.find(sid);
        const uint32_t best_decoy = (decoy_it != decoy_apex_map.end()) ? decoy_it->second : 0;

        const bool has_im = (se.drift_time >= 0.0);

        // Fold changes: one value per quantification bin.
        // Numerator  = group-level max experimental intensity (group.quantification_max_exp).
        // Denominator = this member's index log2 intensity for that bin.
        std::string fold_changes;
        if (!group.quantification_bins.empty() && it != trace_map.end())
        {
          // Build flat_bin_idx → index_log2_intensity map for this member.
          std::unordered_map<uint32_t, float> idx_log2;
          idx_log2.reserve(it->second->apex_fingerprint.size());
          for (const DiaWeaverAlign::ApexFragment& af : it->second->apex_fingerprint)
            idx_log2[af.flat_bin_idx] = af.index_log2_intensity;

          std::ostringstream oss;
          bool first = true;
          for (Size k = 0; k < group.quantification_bins.size(); ++k)
          {
            const auto jt = idx_log2.find(group.quantification_bins[k]);
            if (jt == idx_log2.end()) continue; // bin not in this member's fingerprint
            if (!first) oss << ',';
            oss << (std::log2f(std::max(1.0f, group.quantification_max_exp[k])) - jt->second);
            first = false;
          }
          fold_changes = oss.str();
          if (fold_changes.empty()) fold_changes = "NA";
        }
        else
        {
          fold_changes = "NA";
        }

        tsv << global_group_id
            << '\t' << group.spectrum_ids.size()
            << '\t' << group.mean_apex_rt
            << '\t' << group.shared_fragments.size()
            << '\t' << group.quantification_bins.size()
            << '\t' << group.min_internal_overlap
            << '\t' << window.lower_mz
            << '\t' << window.upper_mz
            << '\t' << win_lo_im
            << '\t' << win_hi_im
            << '\t' << decoy_lo_mz
            << '\t' << decoy_hi_mz
            << '\t' << decoy_lo_im
            << '\t' << decoy_hi_im
            << '\t' << sid
            << '\t' << aligner.getSourceFile(se.source_file_idx)
            << '\t' << se.native_id
            << '\t' << se.precursor_mz
            << '\t' << se.precursor_charge
            << '\t' << (has_im ? std::to_string(se.drift_time) : "NA")
            << '\t' << se.retention_time
            << '\t' << query_apex_rt
            << '\t' << se.n_indexed_fragments
            << '\t' << apex_matched_peaks
            << '\t' << best_decoy
            << '\t' << fold_changes
            << '\n';
      }
      ++global_group_id;
    }
  }

  /// Compute mean number of peaks per spectrum for a set of query spectrum indices.
  static double computeMeanPeakCount_(
    const MSExperiment& exp,
    const std::vector<Size>& q_idxs)
  {
    if (q_idxs.empty()) return 0.0;
    double total = 0.0;
    for (Size idx : q_idxs)
      total += static_cast<double>(exp[idx].size());
    return total / static_cast<double>(q_idxs.size());
  }

  /// Find the best decoy window for window target_w.
  ///
  /// Hard constraint: decoy must NOT overlap or be adjacent to the target in m/z.
  ///
  /// Ranking (two-level):
  ///   1. Primary  — maximize IM overlap with the target window.
  ///                 overlap = max(0, min(upper_im) - max(lower_im))
  ///                 Windows without IM bounds contribute 0 overlap.
  ///   2. Secondary — minimize |mean_peaks[decoy] - mean_peaks[target]|
  ///                 (closest busyness, not highest).
  ///
  /// Returns std::numeric_limits<Size>::max() if no valid decoy window exists.
  static Size findDecoyWindow_(
    const std::vector<DiaWeaver::DIAWindow>& windows,
    const std::vector<double>& mean_peak_counts,
    Size target_w)
  {
    const DiaWeaver::DIAWindow& target = windows[target_w];
    Size   best_idx         = std::numeric_limits<Size>::max();
    double best_im_overlap  = -1.0;
    double best_busy_diff   = std::numeric_limits<double>::max();

    for (Size d = 0; d < windows.size(); ++d)
    {
      if (d == target_w) continue;

      // Hard constraint: decoy must not overlap or be adjacent in m/z
      const double d_lo = windows[d].lower_mz;
      const double d_hi = windows[d].upper_mz;
      const bool non_adjacent = (d_lo > target.upper_mz + 1e-6) ||
                                (d_hi < target.lower_mz - 1e-6);
      if (!non_adjacent) continue;

      // Primary: IM overlap fraction in [0, 1] relative to the target window width.
      // 1.0 = decoy fully covers the target IM range; 0.0 = no overlap.
      // Falls back to 0.0 when either window has no IM bounds.
      double im_overlap = 0.0;
      const double target_im_width = target.upper_im - target.lower_im;
      if (target.hasIonMobility() && windows[d].hasIonMobility() && target_im_width > 0.0)
      {
        const double raw_overlap = std::max(0.0,
          std::min(target.upper_im, windows[d].upper_im) -
          std::max(target.lower_im, windows[d].lower_im));
        im_overlap = raw_overlap / target_im_width;
      }

      // Secondary: busyness closeness
      const double busy_diff = std::abs(mean_peak_counts[d] - mean_peak_counts[target_w]);

      // Accept if: more IM overlap, or equal IM overlap with closer busyness
      if (im_overlap > best_im_overlap ||
          (im_overlap == best_im_overlap && busy_diff < best_busy_diff))
      {
        best_im_overlap = im_overlap;
        best_busy_diff  = busy_diff;
        best_idx        = d;
      }
    }
    return best_idx;
  }
};


int main(int argc, const char** argv)
{
  TOPPDiaWeaverAlign tool;
  return tool.main(argc, argv);
}

/// @cond TOPP_DOCS
