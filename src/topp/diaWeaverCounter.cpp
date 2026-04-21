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
to MS2 mass traces detected in the raw DIA data. A 2D KD-tree (RT x m/z) accelerates the mapping,
with an optional post-filter on ion mobility when IM data are present.

Output metrics:
- Per-trace: number of peptides claiming each trace (orphan / unique / ambiguous)
- Per-peptide: number of theoretical ions matched to traces and fragment coverage (%)
- Console summary of trace accounting statistics

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
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/DATASTRUCTURES/KDTree.h>

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
  }

  Param getSubsectionDefaults_(const String& section) const override
  {
    if (section == "MassTraceDetection")
    {
      MassTraceDetection mtd;
      Param p = mtd.getDefaults();
      p.setValue("mass_error_ppm",  10.0, "Allowed mass deviation (ppm)");
      p.setValue("min_trace_length", 2.0, "Minimum trace length (seconds)");
      return p;
    }
    return Param();
  }

  // -------------------------------------------------------------------------
  // TSV placeholder parser
  //
  // Recognised column aliases (case-insensitive):
  //   sequence : "sequence", "modified.sequence", "stripped.sequence",
  //              "pep.strippedsequence"
  //   charge   : "charge", "precursor.charge", "eg.precursorcharge"
  //   rt       : "rt", "eg.apexrt", "eg.meanapexrt"           [seconds]
  //   mz       : "mz", "precursor.mz", "fg.precmz"
  //   im       : "im", "ionmobility", "eg.ionmobility"
  //
  // TODO: Add DIA-NN-specific RT unit detection (DIA-NN reports RT in minutes)
  // TODO: Add Spectronaut-specific sequence notation stripping
  // TODO: Add support for idXML / pepXML inputs via IdXMLFile / PepXMLFile
  // -------------------------------------------------------------------------
  vector<PeptideEntry> parsePeptideTSV_(const String& filename) const
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
      {"sequence",              "seq"}, {"modified.sequence", "seq"},
      {"stripped.sequence",     "seq"}, {"pep.strippedsequence", "seq"},
      {"charge",                "charge"}, {"precursor.charge", "charge"},
      {"eg.precursorcharge",    "charge"},
      {"rt",                    "rt"}, {"eg.apexrt", "rt"}, {"eg.meanapexrt", "rt"},
      {"mz",                    "mz"}, {"precursor.mz", "mz"}, {"fg.precmz", "mz"},
      {"im",                    "im"}, {"ionmobility", "im"}, {"eg.ionmobility", "im"}
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
                         "    rt       : rt | eg.apexrt | eg.meanapexrt (seconds)\n";
      return entries;
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
        PeptideEntry e;
        e.sequence = fields.at(field_col.at("seq")).trim();
        e.charge   = fields.at(field_col.at("charge")).trim().toInt();
        e.rt       = fields.at(field_col.at("rt")).trim().toDouble();
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
    const String in_file          = getStringOption_("in");
    const String in_ids_file      = getStringOption_("in_ids");
    const String out_traces_file  = getStringOption_("out_traces");
    const String out_pep_file     = getStringOption_("out_peptides");

    const double mz_tol_ppm    = getDoubleOption_("mz_tolerance");
    const double rt_tol        = getDoubleOption_("rt_tolerance");
    const double im_tol        = getDoubleOption_("im_tolerance");

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
    // Step 2: Detect MS2 mass traces
    // ------------------------------------------------------------------
    MassTraceDetection mtd;
    mtd.setParameters(getParam_().copy("MassTraceDetection:", true));

    vector<MassTrace> ms2_traces;
    OPENMS_LOG_INFO << "Running MassTraceDetection on MS2 spectra...\n";
    mtd.run(ms2_exp, ms2_traces);
    OPENMS_LOG_INFO << "Detected " << ms2_traces.size() << " MS2 mass traces.\n";

    if (ms2_traces.empty())
    {
      OPENMS_LOG_WARN << "No MS2 mass traces detected. "
                         "Consider relaxing MassTraceDetection parameters.\n";
    }

    // Check once whether any trace carries IM data
    bool dataset_has_im = false;
    for (const auto& mt : ms2_traces)
    {
      if (mt.containsIMData()) { dataset_has_im = true; break; }
    }

    // ------------------------------------------------------------------
    // Step 3: Build 2D KD-tree (RT x m/z) over detected MS2 traces
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
    // Step 4: Parse peptide identification TSV (placeholder)
    // ------------------------------------------------------------------
    vector<PeptideEntry> peptides = parsePeptideTSV_(in_ids_file);
    OPENMS_LOG_INFO << "Peptides to map: " << peptides.size() << "\n";

    // ------------------------------------------------------------------
    // Step 5: Generate theoretical b/y ions and map to traces via KD-tree
    //
    // trace_claims[i] = set of peptide indices that have at least one
    //                   theoretical ion matching trace i within tolerances.
    // ------------------------------------------------------------------
    vector<set<Size>> trace_claims(ms2_traces.size());
    vector<Size>      pep_n_ions(peptides.size(), 0);       // total theoretical ions
    vector<Size>      pep_matched_traces(peptides.size(), 0); // distinct traces matched

    TheoreticalSpectrumGenerator tsg;
    // Use b and y ions only (defaults); no neutral losses for the baseline
    Param tsg_params = tsg.getDefaults();
    tsg_params.setValue("add_a_ions",   "false", "");
    tsg_params.setValue("add_c_ions",   "false", "");
    tsg_params.setValue("add_x_ions",   "false", "");
    tsg_params.setValue("add_z_ions",   "false", "");
    tsg_params.setValue("add_losses",   "false", "");
    tsg.setParameters(tsg_params);

    for (Size pep_idx = 0; pep_idx < peptides.size(); ++pep_idx)
    {
      const PeptideEntry& pep = peptides[pep_idx];

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

      // Need at least 2 residues for any b/y ion
      if (aa_seq.size() < 2) continue;

      // Fragment charges 1 .. precursor_charge (same charge as precursor is valid)
      const int eff_max_z = max(1, pep.charge);

      PeakSpectrum theo_spec;
      tsg.getSpectrum(theo_spec, aa_seq, 1, eff_max_z);
      pep_n_ions[pep_idx] = theo_spec.size();

      set<Size> matched_for_pep;

      for (Size ion_i = 0; ion_i < theo_spec.size(); ++ion_i)
      {
        const double mz_theo   = theo_spec[ion_i].getMZ();
        const double mz_tol_da = mz_theo * mz_tol_ppm * 1e-6;

        // Rectangular KD-tree region query: RT x m/z
        TraceKDTree::_Region_ region;
        region._M_low_bounds[0]  = pep.rt - rt_tol;
        region._M_high_bounds[0] = pep.rt + rt_tol;
        region._M_low_bounds[1]  = mz_theo - mz_tol_da;
        region._M_high_bounds[1] = mz_theo + mz_tol_da;

        vector<MS2TraceNode> candidates;
        kd_tree.find_within_range(region, back_inserter(candidates));

        for (const auto& cand : candidates)
        {
          // IM post-filter: only applied when both peptide and trace carry IM
          if (dataset_has_im && pep.im > 0.0 &&
              ms2_traces[cand.trace_idx].containsIMData())
          {
            if (std::abs(cand.im - pep.im) > im_tol) continue;
          }

          matched_for_pep.insert(cand.trace_idx);
          trace_claims[cand.trace_idx].insert(pep_idx);
        }
      }

      pep_matched_traces[pep_idx] = matched_for_pep.size();
    }

    // ------------------------------------------------------------------
    // Step 6: Aggregate trace-level metrics
    // ------------------------------------------------------------------
    Size n_orphan    = 0;
    Size n_unique    = 0;
    Size n_ambiguous = 0;

    for (const auto& claims : trace_claims)
    {
      switch (claims.size())
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

      ofs << "trace_idx\tcentroid_mz\tcentroid_rt\tcentroid_im\t"
             "n_claiming_peptides\tstatus\tclaiming_peptide_ids\n";

      for (Size i = 0; i < ms2_traces.size(); ++i)
      {
        const auto& mt    = ms2_traces[i];
        const Size  n_pep = trace_claims[i].size();

        String status;
        if      (n_pep == 0) status = "orphan";
        else if (n_pep == 1) status = "unique";
        else                 status = "ambiguous";

        // Collect comma-separated claiming peptide IDs
        String pep_ids;
        for (const Size pi : trace_claims[i])
        {
          if (!pep_ids.empty()) pep_ids += ",";
          pep_ids += peptides[pi].id;
        }

        ofs << i                         << "\t"
            << mt.getCentroidMZ()        << "\t"
            << mt.getCentroidRT()        << "\t"
            << (mt.containsIMData() ?
                String(mt.getCentroidIM()) : String("N/A")) << "\t"
            << n_pep                     << "\t"
            << status                    << "\t"
            << pep_ids                   << "\n";
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
             "n_theoretical_ions\tn_matched_traces\tcoverage_pct\n";

      for (Size pi = 0; pi < peptides.size(); ++pi)
      {
        const auto& pep      = peptides[pi];
        const Size  n_theo   = pep_n_ions[pi];
        const Size  n_match  = pep_matched_traces[pi];
        const double cov_pct = n_theo > 0 ?
                                 100.0 * static_cast<double>(n_match) / n_theo : 0.0;

        ofs << pep.id       << "\t"
            << pep.sequence << "\t"
            << pep.charge   << "\t"
            << pep.rt       << "\t"
            << n_theo       << "\t"
            << n_match      << "\t"
            << cov_pct      << "\n";
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
