// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mohammed Alhigaylan $
// $Authors: Mohammed Alhigaylan $
// --------------------------------------------------------------------------

/**
@page TOPP_diaWeaverCounter diaWeaverCounter

@brief Accounts for MS2 mass traces by mapping theoretical peptide b/y fragment ions to experimental ion traces.

For each peptide identified at 1% FDR from an external search engine (DIA-NN, Spectronaut, etc.),
this tool generates theoretical b/y fragment ions using TheoreticalSpectrumGenerator and maps them
to MS2 mass traces detected in the raw DIA data via MassTraceDetection + ElutionPeakDetection.
A 2D KD-tree (RT x m/z) accelerates the mapping, with an optional post-filter on ion mobility.

The best representative trace per fragment ion is selected by apex SNR (from ElutionPeakDetection).
When multiple peptides claim the same trace, a collision row is written to out_collisions.

Output files:
- out_traces    : per-trace accountability (orphan / unique / ambiguous)
- out_peptides  : per-peptide fragment ion coverage (optional)
- out_collisions: long-format collision table, one row per (trace, peptide, ion) (optional)

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_diaWeaverCounter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_diaWeaverCounter.html
*/

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include <OpenMS/FEATUREFINDER/ElutionPeakDetection.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/DATASTRUCTURES/KDTree.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <fstream>
#include <map>
#include <set>
#include <cmath>
#include <algorithm>

using namespace OpenMS;
using namespace std;

// ---------------------------------------------------------------------------
// KD-tree node: stores centroid RT, m/z, IM and the index into ms2_traces
// Dimensions: [0] = RT, [1] = m/z  (2D tree; IM is post-filtered separately)
// ---------------------------------------------------------------------------
struct MS2TraceNode
{
  double rt;
  double mz;
  double im;         // 0.0 when IM data are absent
  Size   trace_idx;  // index into the ms2_traces vector

  typedef double value_type;

  double operator[](size_t i) const
  {
    return i == 0 ? rt : mz;
  }
};

typedef KDTree::KDTree<2, MS2TraceNode> TraceKDTree;

// ---------------------------------------------------------------------------
// One peptide entry read from the input identification TSV
// ---------------------------------------------------------------------------
struct PeptideEntry
{
  String sequence;   // OpenMS-compatible (modified) sequence string
  int    charge  = 0;
  double rt      = 0.0;  // apex RT in seconds
  double mz      = 0.0;  // precursor m/z (reference only)
  double im      = 0.0;  // precursor IM centroid; 0.0 = not available
  String id;             // row identifier (used in output, e.g. "row 42")
};

// ---------------------------------------------------------------------------
// Best trace selected for a specific (peptide, ion) pair
// ---------------------------------------------------------------------------
struct IonBestMatch
{
  Size   best_trace_idx = 0;
  double best_snr       = -1.0;
  Size   n_traces       = 0;   // total traces that matched this ion
};

// ===========================================================================

class TOPPDiaWeaverCounter : public TOPPBase
{
public:
  TOPPDiaWeaverCounter() :
    TOPPBase("diaWeaverCounter",
             "Accounts for MS2 mass traces by mapping theoretical peptide fragment ions.",
             false)
  {}

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input DIA raw data (mzML)", true);
    setValidFormats_("in", {"mzML"});

    registerInputFile_("in_ids", "<file>", "",
                       "Peptide identifications at 1% FDR (TSV). "
                       "See documentation for expected column names.", true);
    setValidFormats_("in_ids", {"tsv"});

    registerOutputFile_("out_traces", "<file>", "",
                        "Per-trace accountability output (TSV)", true);
    setValidFormats_("out_traces", {"tsv"});

    registerOutputFile_("out_peptides", "<file>", "",
                        "Per-peptide fragment-ion coverage output (TSV). "
                        "Leave empty to skip.", false);
    setValidFormats_("out_peptides", {"tsv"});

    registerOutputFile_("out_collisions", "<file>", "",
                        "Long-format collision table: one row per (trace, peptide, ion) "
                        "for every trace claimed by more than one peptide. "
                        "Leave empty to skip.", false);
    setValidFormats_("out_collisions", {"tsv"});

    registerDoubleOption_("mz_tolerance", "<ppm>", 20.0,
                          "Fragment ion m/z tolerance (±, ppm)", false);
    setMinFloat_("mz_tolerance", 0.0);

    registerDoubleOption_("rt_tolerance", "<s>", 15.0,
                          "RT tolerance (±, seconds) around the peptide apex", false);
    setMinFloat_("rt_tolerance", 0.0);

    registerDoubleOption_("im_tolerance", "<1/K0>", 0.02,
                          "Ion mobility tolerance (±, 1/K0; applied only when IM data are present)", false);
    setMinFloat_("im_tolerance", 0.0);

    registerSubsection_("MassTraceDetection",
                        "Parameters forwarded to the internal MassTraceDetection algorithm");

    registerSubsection_("ElutionPeakDetection",
                        "Parameters forwarded to the internal ElutionPeakDetection algorithm");
  }

  Param getSubsectionDefaults_(const String& section) const override
  {
    if (section == "MassTraceDetection")
    {
      MassTraceDetection mtd;
      Param p = mtd.getDefaults();
      p.setValue("mass_error_ppm",   10.0, "Allowed mass deviation (ppm)");
      p.setValue("min_trace_length",  2.0, "Minimum trace length (seconds)");
      return p;
    }
    if (section == "ElutionPeakDetection")
    {
      ElutionPeakDetection epd;
      return epd.getDefaults();
    }
    return Param();
  }

  // -------------------------------------------------------------------------
  // TSV parser for DIA-NN library-free report.tsv (and Spectronaut equivalents)
  //
  // run_id    : stem of the input mzML (no path, no extension). Only rows
  //             whose Run or File.Name stem matches this value are loaded.
  //             Pass an empty string to disable filtering (loads all rows).
  //
  // Recognised column aliases (case-insensitive):
  //   run      : "run", "file.name"   ← used for per-run filtering
  //   sequence : "sequence", "modified.sequence", "stripped.sequence",
  //              "pep.strippedsequence"
  //   charge   : "charge", "precursor.charge", "eg.precursorcharge"
  //   rt       : "rt", "eg.apexrt", "eg.meanapexrt"   [DIA-NN: minutes → ×60]
  //   mz       : "mz", "precursor.mz", "fg.precmz"
  //   im       : "im", "ionmobility", "eg.ionmobility"
  //
  // TODO: Add Spectronaut-specific sequence notation stripping
  // TODO: Add support for idXML / pepXML inputs via IdXMLFile / PepXMLFile
  // -------------------------------------------------------------------------

  // Returns the run identifier embedded in a DIA-NN File.Name value.
  // DIA-NN stores the full raw-file path (e.g. /data/foo/run42.d); we want
  // just the stem ("run42") so it can be compared with the mzML basename.
  static String runIdFromFileName_(const String& file_name_field)
  {
    String stem = File::basename(file_name_field);  // "run42.d"  or  "run42.mzML"
    // Strip any known raw-file extension
    for (const char* ext : {".d", ".mzML", ".mzml", ".raw", ".wiff"})
    {
      if (stem.hasSuffix(ext)) { stem = stem.prefix(stem.size() - strlen(ext)); break; }
    }
    return stem;
  }

  vector<PeptideEntry> parsePeptideTSV_(const String& filename,
                                        const String& run_id) const
  {
    vector<PeptideEntry> entries;

    ifstream file(filename.c_str());
    if (!file.is_open())
    {
      OPENMS_LOG_ERROR << "[diaWeaverCounter] Cannot open peptide TSV: " << filename << "\n";
      return entries;
    }

    // --- Parse header row ---
    string raw_header;
    if (!getline(file, raw_header))
    {
      OPENMS_LOG_ERROR << "[diaWeaverCounter] Empty TSV file: " << filename << "\n";
      return entries;
    }

    vector<String> col_names;
    String(raw_header).split("\t", col_names);
    for (auto& c : col_names) c.trim();

    // Map lowercased column name → internal field tag
    const map<String, String> aliases = {
      {"run",                   "run"}, {"file.name",            "run"},
      {"sequence",              "seq"}, {"modified.sequence",    "seq"},
      {"stripped.sequence",     "seq"}, {"pep.strippedsequence", "seq"},
      {"charge",                "charge"}, {"precursor.charge",  "charge"},
      {"eg.precursorcharge",    "charge"},
      {"rt",                    "rt"},  {"eg.apexrt",            "rt"},
      {"eg.meanapexrt",         "rt"},
      {"mz",                    "mz"},  {"precursor.mz",         "mz"},
      {"fg.precmz",             "mz"},
      {"im",                    "im"},  {"ionmobility",          "im"},
      {"eg.ionmobility",        "im"}
    };

    // field tag → column index
    map<String, int> field_col;
    for (Size i = 0; i < col_names.size(); ++i)
    {
      String lower = col_names[i];
      lower.toLower();
      auto it = aliases.find(lower);
      if (it != aliases.end())
        field_col[it->second] = static_cast<int>(i);
    }

    if (field_col.find("seq")    == field_col.end() ||
        field_col.find("charge") == field_col.end() ||
        field_col.find("rt")     == field_col.end())
    {
      OPENMS_LOG_WARN << "[diaWeaverCounter] TSV is missing required columns "
                         "(sequence / charge / rt). Returning empty peptide list.\n"
                         "  Required (case-insensitive aliases):\n"
                         "    sequence : sequence | modified.sequence | stripped.sequence | "
                                        "pep.strippedsequence\n"
                         "    charge   : charge | precursor.charge | eg.precursorcharge\n"
                         "    rt       : rt | eg.apexrt | eg.meanapexrt (minutes for DIA-NN)\n";
      return entries;
    }

    const bool has_run_col   = field_col.count("run") > 0;
    const bool filter_by_run = !run_id.empty();
    // DIA-NN "Run" column holds just the stem; "File.Name" holds the full path.
    // runIdFromFileName_() normalises both to a bare stem for comparison.
    const bool run_col_is_filename =
      has_run_col &&
      col_names[field_col.at("run")].toLower() == String("file.name");

    if (filter_by_run && !has_run_col)
    {
      OPENMS_LOG_WARN << "[diaWeaverCounter] No run/file.name column found in TSV; "
                         "run filter '" << run_id << "' cannot be applied. "
                         "All rows will be loaded.\n";
    }
    else if (filter_by_run)
    {
      OPENMS_LOG_INFO << "[diaWeaverCounter] Filtering TSV to run: " << run_id << "\n";
    }

    // --- Parse data rows ---
    string raw_line;
    Size row = 1;
    while (getline(file, raw_line))
    {
      ++row;
      if (raw_line.empty()) continue;

      vector<String> fields;
      String(raw_line).split("\t", fields);

      try
      {
        // --- Run filter ---
        if (filter_by_run && has_run_col)
        {
          String row_run = fields.at(field_col.at("run")).trim();
          if (run_col_is_filename) row_run = runIdFromFileName_(row_run);
          if (row_run != run_id) continue;
        }

        PeptideEntry e;
        e.sequence = fields.at(field_col.at("seq")).trim();

        // DIA-NN reports ambiguous X residue. We cannot readily compute theoretical fragment ion masses.
        {
          bool has_bare_x = false;
          for (Size ci = 0; ci < e.sequence.size(); ++ci)
          {
            if (e.sequence[ci] == 'X' &&
                (ci + 1 >= e.sequence.size() || e.sequence[ci + 1] != '['))
            {
              has_bare_x = true;
              break;
            }
          }
          if (has_bare_x)
          {
            OPENMS_LOG_WARN << "[diaWeaverCounter] Sequence '" << e.sequence
                            << "' (row " << row << ") contains ambiguous residue X "
                               "with unknown mass. Cannot compute theoretical fragment "
                               "ions. Skipping.\n";
            continue;
          }
        }

        e.charge   = fields.at(field_col.at("charge")).trim().toInt();
        // DIA-NN reports RT in minutes; convert to seconds
        e.rt       = fields.at(field_col.at("rt")).trim().toDouble() * 60.0;
        e.mz       = field_col.count("mz") ?
                       fields.at(field_col.at("mz")).trim().toDouble() : 0.0;
        e.im       = field_col.count("im") ?
                       fields.at(field_col.at("im")).trim().toDouble() : 0.0;
        e.id       = "row" + String(row);
        entries.push_back(move(e));
      }
      catch (const Exception::BaseException& ex)
      {
        OPENMS_LOG_WARN << "[diaWeaverCounter] Skipping malformed row " << row
                        << ": " << ex.getMessage() << "\n";
      }
      catch (const out_of_range&)
      {
        OPENMS_LOG_WARN << "[diaWeaverCounter] Skipping short row " << row << "\n";
      }
    }

    OPENMS_LOG_INFO << "[diaWeaverCounter] Loaded " << entries.size()
                    << " peptide entries from " << filename << "\n";
    return entries;
  }

  // -------------------------------------------------------------------------

  ExitCodes main_(int, const char**) override
  {
    const String in_file         = getStringOption_("in");
    const String in_ids_file     = getStringOption_("in_ids");
    const String out_traces_file = getStringOption_("out_traces");
    const String out_pep_file    = getStringOption_("out_peptides");
    const String out_coll_file   = getStringOption_("out_collisions");

    const double mz_tol_ppm = getDoubleOption_("mz_tolerance");
    const double rt_tol     = getDoubleOption_("rt_tolerance");
    const double im_tol     = getDoubleOption_("im_tolerance");

    // ------------------------------------------------------------------
    // Step 1: Load mzML and filter to MS2 spectra
    // ------------------------------------------------------------------
    OPENMS_LOG_INFO << "Loading: " << in_file << "\n";
    MSExperiment full_exp;
    MzMLFile().load(in_file, full_exp);

    MSExperiment ms2_exp;
    for (const auto& spec : full_exp)
    {
      if (spec.getMSLevel() == 2)
        ms2_exp.addSpectrum(spec);
    }
    OPENMS_LOG_INFO << "MS2 spectra: " << ms2_exp.size() << "\n";

    if (ms2_exp.empty())
    {
      OPENMS_LOG_ERROR << "No MS2 spectra found in " << in_file << "\n";
      return INCOMPATIBLE_INPUT_DATA;
    }

    // ------------------------------------------------------------------
    // Step 2: Detect MS2 mass traces (MTD → EPD), matching diaWeaver's
    //         runMassTraceExtractor_ pipeline. EPD adds smoothed
    //         intensities required for apex SNR computation.
    // ------------------------------------------------------------------
    ms2_exp.sortSpectra(true);

    MassTraceDetection mtd;
    mtd.setParameters(getParam_().copy("MassTraceDetection:", true));

    vector<MassTrace> raw_traces;
    OPENMS_LOG_INFO << "Running MassTraceDetection...\n";
    mtd.run(ms2_exp, raw_traces);
    OPENMS_LOG_INFO << "MTD: " << raw_traces.size() << " raw traces.\n";

    ElutionPeakDetection epd;
    epd.setParameters(getParam_().copy("ElutionPeakDetection:", true));

    vector<MassTrace> ms2_traces;
    OPENMS_LOG_INFO << "Running ElutionPeakDetection...\n";
    epd.detectPeaks(raw_traces, ms2_traces);

    if (epd.getParameters().getValue("width_filtering") == "auto")
    {
      vector<MassTrace> filtered;
      epd.filterByPeakWidth(ms2_traces, filtered);
      ms2_traces = move(filtered);
    }

    ms2_traces.erase(
      remove_if(ms2_traces.begin(), ms2_traces.end(),
                [](const MassTrace& t) { return t.getSize() == 0; }),
      ms2_traces.end());

    OPENMS_LOG_INFO << "After EPD: " << ms2_traces.size() << " traces.\n";

    if (ms2_traces.empty())
    {
      OPENMS_LOG_WARN << "No MS2 mass traces remain after EPD. "
                         "Consider relaxing MassTraceDetection / ElutionPeakDetection parameters.\n";
    }

    // Precompute apex SNR for every trace (requires smoothed intensities from EPD)
    vector<double> trace_snr(ms2_traces.size(), 0.0);
    for (Size i = 0; i < ms2_traces.size(); ++i)
    {
      if (!ms2_traces[i].getSmoothedIntensities().empty())
        trace_snr[i] = epd.computeApexSNR(ms2_traces[i]);
    }

    // Check once whether any trace carries IM data
    bool dataset_has_im = false;
    for (const auto& mt : ms2_traces)
    {
      if (mt.containsIMData()) { dataset_has_im = true; break; }
    }

    // ------------------------------------------------------------------
    // Step 3: Build 2D KD-tree (RT x m/z) over detected MS2 traces.
    //         IM is stored in the node but used only as a post-filter.
    // ------------------------------------------------------------------
    TraceKDTree kd_tree;

    for (Size i = 0; i < ms2_traces.size(); ++i)
    {
      MS2TraceNode node;
      node.rt        = ms2_traces[i].getCentroidRT();
      node.mz        = ms2_traces[i].getCentroidMZ();
      node.im        = ms2_traces[i].containsIMData() ?
                         ms2_traces[i].getCentroidIM() : 0.0;
      node.trace_idx = i;
      kd_tree.insert(node);
    }
    kd_tree.optimise();
    OPENMS_LOG_INFO << "KD-tree built (" << kd_tree.size() << " nodes).\n";

    // ------------------------------------------------------------------
    // Step 4: Parse peptide identification TSV, filtered to this run.
    //
    // The run ID is the mzML basename without extension (e.g. "run42").
    // DIA-NN's Run column holds exactly this stem; File.Name holds the
    // full raw-file path which runIdFromFileName_() normalises to the stem.
    // ------------------------------------------------------------------
    String run_id = File::basename(in_file);   // "run42.mzML"
    for (const char* ext : {".mzML", ".mzml"})
    {
      if (run_id.hasSuffix(ext)) { run_id = run_id.prefix(run_id.size() - strlen(ext)); break; }
    }
    OPENMS_LOG_INFO << "Run ID derived from mzML: " << run_id << "\n";

    vector<PeptideEntry> peptides = parsePeptideTSV_(in_ids_file, run_id);
    OPENMS_LOG_INFO << "Peptides to map: " << peptides.size() << "\n";

    // ------------------------------------------------------------------
    // Step 5: Generate theoretical b/y ions and map to traces via KD-tree
    //
    // trace_claims[i] = list of (pep_idx, ion_name) pairs that claimed trace i.
    //                   Multiple entries from the same peptide (different ions)
    //                   and from different peptides are all recorded.
    //
    // best_ion_match[(pep_idx, ion_name)] = trace with highest apex SNR
    //                   among all traces that matched this ion.
    // ------------------------------------------------------------------
    vector<vector<pair<Size, String>>> trace_claims(ms2_traces.size());
    map<pair<Size,String>, IonBestMatch> best_ion_match;

    vector<Size>      pep_n_ions(peptides.size(), 0);
    vector<Size>      pep_matched_ions(peptides.size(), 0);      // ions with >=1 trace
    vector<Size>      pep_multi_trace_ions(peptides.size(), 0);  // ions with >1 trace
    vector<set<Size>> pep_matched_trace_set(peptides.size());    // distinct traces matched

    TheoreticalSpectrumGenerator tsg;
    Param tsg_params = tsg.getDefaults();
    tsg_params.setValue("add_metainfo", "true",  "");  // ion names needed for collision table
    tsg_params.setValue("add_a_ions",   "false", "");
    tsg_params.setValue("add_c_ions",   "false", "");
    tsg_params.setValue("add_x_ions",   "false", "");
    tsg_params.setValue("add_z_ions",   "false", "");
    tsg_params.setValue("add_losses",   "false", "");
    tsg.setParameters(tsg_params);

    for (Size pep_idx = 0; pep_idx < peptides.size(); ++pep_idx)
    {
      const PeptideEntry& pep = peptides[pep_idx];

      if (pep.charge < 1)
      {
        OPENMS_LOG_ERROR << "[diaWeaverCounter] Invalid charge " << pep.charge
                         << " for peptide '" << pep.sequence << "' (" << pep.id
                         << "). Charge must be >= 1. Aborting.\n";
        return INPUT_FILE_CORRUPT;
      }

      AASequence aa_seq;
      try
      {
        aa_seq = AASequence::fromString(pep.sequence);
      }
      catch (const Exception::BaseException& ex)
      {
        OPENMS_LOG_WARN << "[diaWeaverCounter] Cannot parse sequence '"
                        << pep.sequence << "': " << ex.getMessage() << ". Skipping.\n";
        continue;
      }

      if (aa_seq.size() < 2) continue;  // no b/y ions possible

      PeakSpectrum theo_spec;
      tsg.getSpectrum(theo_spec, aa_seq, 1, pep.charge);
      pep_n_ions[pep_idx] = theo_spec.size();

      // Ion names in StringDataArrays[0], fragment charges in IntegerDataArrays[0]
      const PeakSpectrum::StringDataArray&  ion_names    =
        theo_spec.getStringDataArrays().at(0);
      const PeakSpectrum::IntegerDataArray& frag_charges =
        theo_spec.getIntegerDataArrays().at(0);

      // Each fragment ion is searched at 4 isotope positions (M+0 .. M+3);
      // the precursor is searched at 5 (M+0 .. M+4).
      pep_n_ions[pep_idx] = theo_spec.size() * 4 + 5;

      // Helper: KD-tree query + IM post-filter → trace index list
      auto queryTraces = [&](double mz_center) -> vector<Size>
      {
        const double mz_tol_da = mz_center * mz_tol_ppm * 1e-6;
        TraceKDTree::_Region_ region;
        region._M_low_bounds[0]  = pep.rt - rt_tol;
        region._M_high_bounds[0] = pep.rt + rt_tol;
        region._M_low_bounds[1]  = mz_center - mz_tol_da;
        region._M_high_bounds[1] = mz_center + mz_tol_da;

        vector<MS2TraceNode> candidates;
        kd_tree.find_within_range(region, back_inserter(candidates));

        vector<Size> result;
        for (const auto& cand : candidates)
        {
          if (dataset_has_im && pep.im > 0.0 &&
              ms2_traces[cand.trace_idx].containsIMData())
          {
            if (std::abs(cand.im - pep.im) > im_tol) continue;
          }
          result.push_back(cand.trace_idx);
        }
        return result;
      };

      // Helper: record trace claims and update best-SNR book-keeping for one ion name
      auto recordMatches = [&](const String& ion_name, const vector<Size>& tidx_list)
      {
        auto ion_key = make_pair(pep_idx, ion_name);
        IonBestMatch& best = best_ion_match[ion_key];
        best.n_traces += tidx_list.size();

        for (Size tidx : tidx_list)
        {
          trace_claims[tidx].emplace_back(pep_idx, ion_name);
          pep_matched_trace_set[pep_idx].insert(tidx);

          if (trace_snr[tidx] > best.best_snr)
          {
            best.best_snr       = trace_snr[tidx];
            best.best_trace_idx = tidx;
          }
        }
      };

      // --- Fragment b/y ions: each isotope searched and reported independently.
      //     Heavier isotopes are only searched if the monoisotopic trace is found. ---
      for (Size ion_i = 0; ion_i < theo_spec.size(); ++ion_i)
      {
        const double mz_mono  = theo_spec[ion_i].getMZ();
        const String ion_name = ion_names[ion_i];
        const int    frag_z   = frag_charges[ion_i];

        const vector<Size> mono_tidxs = queryTraces(mz_mono);
        if (mono_tidxs.empty()) continue;

        ++pep_matched_ions[pep_idx];
        if (mono_tidxs.size() > 1) ++pep_multi_trace_ions[pep_idx];
        recordMatches(ion_name + "[M+0]", mono_tidxs);

        for (int iso = 1; iso <= 3; ++iso)
        {
          const double mz_iso      = mz_mono + iso * Constants::C13C12_MASSDIFF_U / frag_z;
          const vector<Size> tidxs = queryTraces(mz_iso);

          if (tidxs.empty()) continue;

          ++pep_matched_ions[pep_idx];
          if (tidxs.size() > 1) ++pep_multi_trace_ions[pep_idx];
          recordMatches(ion_name + "[M+" + String(iso) + "]", tidxs);
        }
      }

      // --- Precursor ion: heavier isotopes only searched if monoisotopic is found. ---
      {
        const double mz_prec_mono = aa_seq.getMZ(pep.charge);

        const vector<Size> mono_tidxs = queryTraces(mz_prec_mono);
        if (!mono_tidxs.empty())
        {
          ++pep_matched_ions[pep_idx];
          if (mono_tidxs.size() > 1) ++pep_multi_trace_ions[pep_idx];
          recordMatches("p[M+0]", mono_tidxs);

          for (int iso = 1; iso <= 4; ++iso)
          {
            const double mz_iso      = mz_prec_mono + iso * Constants::C13C12_MASSDIFF_U / pep.charge;
            const vector<Size> tidxs = queryTraces(mz_iso);

            if (tidxs.empty()) continue;

            ++pep_matched_ions[pep_idx];
            if (tidxs.size() > 1) ++pep_multi_trace_ions[pep_idx];
            recordMatches("p[M+" + String(iso) + "]", tidxs);
          }
        }
      }
    }

    // ------------------------------------------------------------------
    // Step 6: Derive per-trace and per-peptide aggregate metrics
    // ------------------------------------------------------------------

    // Number of unique peptides claiming each trace
    vector<Size> trace_n_peptides(ms2_traces.size(), 0);
    for (Size i = 0; i < ms2_traces.size(); ++i)
    {
      set<Size> unique_peps;
      for (const auto& [pi, ion] : trace_claims[i])
        unique_peps.insert(pi);
      trace_n_peptides[i] = unique_peps.size();
    }

    // Per-peptide exclusive vs shared trace counts
    vector<Size> pep_exclusive(peptides.size(), 0);
    vector<Size> pep_shared(peptides.size(), 0);
    for (Size pi = 0; pi < peptides.size(); ++pi)
    {
      for (Size tidx : pep_matched_trace_set[pi])
      {
        if (trace_n_peptides[tidx] == 1) ++pep_exclusive[pi];
        else                              ++pep_shared[pi];
      }
    }

    Size n_orphan    = 0;
    Size n_unique    = 0;
    Size n_ambiguous = 0;
    for (Size i = 0; i < ms2_traces.size(); ++i)
    {
      switch (trace_n_peptides[i])
      {
        case 0:  ++n_orphan;    break;
        case 1:  ++n_unique;    break;
        default: ++n_ambiguous; break;
      }
    }

    OPENMS_LOG_INFO
      << "\n=== diaWeaverCounter Summary ===\n"
      << "Total MS2 traces : " << ms2_traces.size() << "\n"
      << "  Orphan    (0 peptides) : " << n_orphan    << "\n"
      << "  Unique    (1 peptide)  : " << n_unique    << "\n"
      << "  Ambiguous (>1 peptide) : " << n_ambiguous << "\n"
      << "Total peptides    : " << peptides.size() << "\n";

    // ------------------------------------------------------------------
    // Step 7: Write per-trace accountability TSV
    // ------------------------------------------------------------------
    {
      ofstream ofs(out_traces_file.c_str());
      if (!ofs.is_open())
      {
        OPENMS_LOG_ERROR << "Cannot write trace output: " << out_traces_file << "\n";
        return CANNOT_WRITE_OUTPUT_FILE;
      }

      ofs << "trace_idx\tcentroid_mz\tcentroid_rt\tcentroid_im\tapex_snr\t"
             "n_claiming_peptides\tn_claiming_ions\tstatus\n";

      for (Size i = 0; i < ms2_traces.size(); ++i)
      {
        const auto& mt    = ms2_traces[i];
        const Size  n_pep = trace_n_peptides[i];

        String status;
        if      (n_pep == 0) status = "orphan";
        else if (n_pep == 1) status = "unique";
        else                 status = "ambiguous";

        ofs << i                           << "\t"
            << mt.getCentroidMZ()          << "\t"
            << mt.getCentroidRT()          << "\t"
            << (mt.containsIMData() ?
                String(mt.getCentroidIM()) : String("N/A")) << "\t"
            << trace_snr[i]                << "\t"
            << n_pep                       << "\t"
            << trace_claims[i].size()      << "\t"
            << status                      << "\n";
      }
    }

    // ------------------------------------------------------------------
    // Step 8: Write per-peptide coverage TSV (optional)
    // ------------------------------------------------------------------
    if (!out_pep_file.empty())
    {
      ofstream ofs(out_pep_file.c_str());
      if (!ofs.is_open())
      {
        OPENMS_LOG_ERROR << "Cannot write peptide output: " << out_pep_file << "\n";
        return CANNOT_WRITE_OUTPUT_FILE;
      }

      ofs << "peptide_id\tsequence\tcharge\trt\t"
             "n_theoretical_ions\tn_matched_ions\tn_multi_trace_ions\t"
             "n_exclusive_traces\tn_shared_traces\tcoverage_pct\n";

      for (Size pi = 0; pi < peptides.size(); ++pi)
      {
        const auto& pep     = peptides[pi];
        const Size  n_theo  = pep_n_ions[pi];
        const double cov    = n_theo > 0 ?
                                100.0 * static_cast<double>(pep_matched_ions[pi]) / n_theo
                                : 0.0;

        ofs << pep.id                    << "\t"
            << pep.sequence              << "\t"
            << pep.charge                << "\t"
            << pep.rt                    << "\t"
            << n_theo                    << "\t"
            << pep_matched_ions[pi]      << "\t"
            << pep_multi_trace_ions[pi]  << "\t"
            << pep_exclusive[pi]         << "\t"
            << pep_shared[pi]            << "\t"
            << cov                       << "\n";
      }
    }

    // ------------------------------------------------------------------
    // Step 9: Write collision TSV (optional)
    //
    // Long format: one row per (trace, peptide, ion) for every trace
    // claimed by more than one unique peptide.
    //
    // is_best = true  when this trace has the highest apex SNR among all
    //                 traces that matched (peptide_id, ion_name).
    // ------------------------------------------------------------------
    if (!out_coll_file.empty())
    {
      ofstream ofs(out_coll_file.c_str());
      if (!ofs.is_open())
      {
        OPENMS_LOG_ERROR << "Cannot write collision output: " << out_coll_file << "\n";
        return CANNOT_WRITE_OUTPUT_FILE;
      }

      ofs << "trace_idx\ttrace_mz\ttrace_rt\ttrace_im\tapex_snr\t"
             "n_claiming_peptides\tpeptide_id\tsequence\tion_name\tis_best\n";

      for (Size i = 0; i < ms2_traces.size(); ++i)
      {
        if (trace_n_peptides[i] <= 1) continue;  // only collision traces

        const auto& mt     = ms2_traces[i];
        const String im_str = mt.containsIMData() ?
                                String(mt.getCentroidIM()) : String("N/A");

        for (const auto& [pi, ion_name] : trace_claims[i])
        {
          const auto ion_key = make_pair(pi, ion_name);
          const bool is_best = best_ion_match.count(ion_key) &&
                               best_ion_match.at(ion_key).best_trace_idx == i;

          ofs << i                           << "\t"
              << mt.getCentroidMZ()          << "\t"
              << mt.getCentroidRT()          << "\t"
              << im_str                      << "\t"
              << trace_snr[i]                << "\t"
              << trace_n_peptides[i]         << "\t"
              << peptides[pi].id             << "\t"
              << peptides[pi].sequence       << "\t"
              << ion_name                    << "\t"
              << (is_best ? "true" : "false") << "\n";
        }
      }
    }

    return EXECUTION_OK;
  }
};

// ---------------------------------------------------------------------------

int main(int argc, const char** argv)
{
  TOPPDiaWeaverCounter tool;
  return tool.main(argc, argv);
}
