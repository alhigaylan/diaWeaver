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
#include <OpenMS/APPLICATIONS/diaWeaver.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include <OpenMS/FEATUREFINDER/ElutionPeakDetection.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerIM.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
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
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace OpenMS;

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

  // Open-search fields (populated only by parseMSFraggerOpenSearchTSV_)
  String         bare_sequence;           // unmodified peptide sequence
  std::map<int,double> assigned_mods;          // positional delta mods from Assigned Modifications
  double         delta_mass = 0.0;        // open-search delta mass
  std::vector<int>    localization_candidates; // 1-indexed residue positions from Best Positions
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

    registerOutputFile_("out_summary", "<file>", "",
                        "Run-level ion accounting summary (TSV): "
                        "total traces, explained/unexplained counts and percentages.", true);
    setValidFormats_("out_summary", {"tsv"});

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

    registerIntOption_("threads", "<n>", 1,
                       "Number of threads for parallel window processing.", false);
    setMinInt_("threads", 1);

    registerFlag_("aggregate_across_scans",
                  "If set, aggregate signal across neighboring scans using Gaussian weighting "
                  "before peak picking. Improves signal-to-noise for low-intensity peaks. "
                  "Only effective when IM data are present (TimsTOF).", false);

    registerSubsection_("MassTraceExtractor",
                        "Parameters for per-window MS2 mass trace extraction (mirrors diaWeaver)");

    registerSubsection_("PeakPickerIM",
                        "Parameters for ion mobility peak picking (used when input has IM data)");

    registerSubsection_("PeakPickerHiRes",
                        "Parameters for high-resolution peak picking (used when input has no IM data)");
  }

  Param getSubsectionDefaults_(const String& section) const override
  {
    if (section == "MassTraceExtractor")
    {
      Param combined;

      Param p_com;
      p_com.setValue("noise_threshold_int",       30.0,   "Intensity threshold below which peaks are regarded as noise.");
      p_com.setValue("chrom_peak_snr",             1.0,   "Minimum signal-to-noise a mass trace should have.");
      p_com.setValue("chrom_fwhm",                 3.0,   "Expected chromatographic peak width (in seconds).");
      p_com.setValue("auto_noise_threshold",      "false", "Automatically estimate noise threshold from the input map.");
      p_com.setValidStrings("auto_noise_threshold", {"true", "false"});
      p_com.setValue("noise_estimation_n_scans",  50,     "Number of scans sampled to estimate noise when auto_noise_threshold is true.");
      p_com.setValue("noise_estimation_percentile", 80.0, "Intensity percentile used to define noise level from sampled scans.");
      combined.insert("common:", p_com);
      combined.setSectionDescription("common", "Common parameters shared by MTD and EPD");

      Param p_mtd = MassTraceDetection().getDefaults();
      p_mtd.setValue("mass_error_ppm",              10.0, "Allowed mass deviation (ppm).");
      p_mtd.setValue("min_trace_length",             2.0, "Minimum expected trace length (seconds).");
      p_mtd.setValue("trace_termination_outliers",     2, "Consecutive missing-peak scans before trace termination.");
      p_mtd.setValue("reestimate_mt_sd",          "false", "Enable dynamic re-estimation of m/z variance.");
      p_mtd.setValue("quant_method",          "max_height", "Quantification method for mass traces.");
      p_mtd.setValue("impute_zeros_missing_scans", "true", "Insert zero-intensity points at scans with no peak (preserves RT grid).");
      p_mtd.remove("noise_threshold_int");
      p_mtd.remove("chrom_peak_snr");
      p_mtd.remove("auto_noise_threshold");
      p_mtd.remove("noise_estimation_n_scans");
      p_mtd.remove("noise_estimation_percentile");
      combined.insert("mtd:", p_mtd);
      combined.setSectionDescription("mtd", "MassTraceDetection parameters");

      Param p_epd;
      p_epd.setValue("enabled",         "true", "Enable splitting of isobaric mass traces by elution peak detection.");
      p_epd.setValue("width_filtering",  "off", "Filter unlikely peak widths. 'auto' uses 5th/95th percentile; 'fixed' uses min/max_fwhm.");
      p_epd.setValidStrings("enabled",         {"true", "false"});
      p_epd.setValidStrings("width_filtering", {"off", "fixed", "auto"});
      p_epd.insert("", ElutionPeakDetection().getDefaults());
      p_epd.remove("chrom_peak_snr");
      p_epd.remove("chrom_fwhm");
      combined.insert("epd:", p_epd);
      combined.setSectionDescription("epd", "ElutionPeakDetection parameters");

      return combined;
    }
    if (section == "PeakPickerIM")
    {
      PeakPickerIM pp;
      return pp.getDefaults();
    }
    if (section == "PeakPickerHiRes")
    {
      PeakPickerHiRes pp;
      return pp.getDefaults();
    }
    return Param();
  }

  // -------------------------------------------------------------------------
  // TSV parser for DIA-NN library-free report.tsv and MSFragger psm.tsv
  //
  // run_id    : stem of the input mzML (no path, no extension). Only rows
  //             whose run column stem matches this value are loaded.
  //             Pass an empty string to disable filtering (loads all rows).
  //
  // Recognised column aliases (case-insensitive):
  //   run      : "run", "file.name"          ← DIA-NN
  //              "spectrum"                   ← MSFragger (stem extracted via runIdFromSpectrum_)
  //   sequence : "modified peptide", "peptide"          ← MSFragger
  //              "modified.sequence", "stripped.sequence",
  //              "sequence", "pep.strippedsequence"      ← DIA-NN / generic
  //   charge   : "charge", "precursor.charge", "eg.precursorcharge"
  //   rt       : "retention"                 ← MSFragger [seconds, no conversion]
  //              "rt", "eg.apexrt", "eg.meanapexrt"      ← DIA-NN [minutes → ×60]
  //   mz       : "calibrated observed m/z", "observed m/z"  ← MSFragger
  //              "mz", "precursor.mz", "fg.precmz"           ← DIA-NN / generic
  //   im       : "im", "ionmobility", "eg.ionmobility"
  //
  // TODO: Add Spectronaut-specific sequence notation stripping
  // TODO: Add support for idXML / pepXML inputs via IdXMLFile / PepXMLFile
  // -------------------------------------------------------------------------

  // Returns the run identifier from a DIA-NN File.Name value (full raw-file path).
  // Strips path and known raw-file extensions to produce a bare stem.
  static String runIdFromFileName_(const String& file_name_field)
  {
    String stem = File::basename(file_name_field);  // "/data/foo/run42.d" → "run42.d"
    for (const char* ext : {".d", ".mzML", ".mzml", ".raw", ".wiff"})
    {
      if (stem.hasSuffix(ext)) { stem = stem.prefix(stem.size() - strlen(ext)); break; }
    }
    return stem;
  }

  // Returns the run identifier from an MSFragger Spectrum field value.
  // Format: "<run_name>[_diatracer].<scan>.<scan>.<charge>"
  // e.g.  "run42_diatracer.00030.00030.3" → "run42"
  static String runIdFromSpectrum_(const String& spectrum_field)
  {
    // Strip trailing .scan.scan.charge  (three dot-delimited numeric suffixes)
    String stem = spectrum_field;
    for (int i = 0; i < 3; ++i)
    {
      Size dot = stem.rfind('.');
      if (dot == String::npos) break;
      stem = stem.prefix(dot);
    }
    // Strip _diatracer suffix added by FragPipe DIA processing
    const String diatracer_suffix = "_diatracer";
    if (stem.hasSuffix(diatracer_suffix))
      stem = stem.prefix(stem.size() - diatracer_suffix.size());
    return stem;
  }

  // Parse "13C(57.0214), N-term(42.0106), ..." → {position → delta_mass}
  // pos 0 = N-terminal, pos > 0 = 1-indexed residue, pos -1 = C-terminal
  static std::map<int, double> parseAssignedMods_(const String& s)
  {
    std::map<int, double> result;
    if (String(s).trim().empty()) return result;

    std::vector<String> parts;
    String(s).split(",", parts);

    for (String& part : parts)
    {
      part.trim();
      Size open_p  = part.find('(');
      Size close_p = part.rfind(')');
      if (open_p == String::npos || close_p == String::npos || close_p <= open_p) continue;

      String prefix   = part.prefix(open_p).trim();
      String mass_str = part.substr(open_p + 1, close_p - open_p - 1).trim();

      double delta = 0.0;
      try { delta = mass_str.toDouble(); } catch (...) { continue; }

      String prefix_lower = prefix;
      prefix_lower.toLower();

      if (prefix_lower == "n-term" || prefix_lower == "nterm")
      {
        result[0] += delta;
      }
      else if (prefix_lower == "c-term" || prefix_lower == "cterm")
      {
        result[-1] += delta;
      }
      else
      {
        // "13C" → strip trailing amino acid letter(s) to get numeric position
        Size digit_end = prefix.size();
        while (digit_end > 0 && isalpha((unsigned char)prefix[digit_end - 1]))
          --digit_end;
        if (digit_end == 0) continue;

        int pos = 0;
        try { pos = prefix.prefix(digit_end).toInt(); } catch (...) { continue; }
        if (pos < 1) continue;
        result[pos] += delta;
      }
    }
    return result;
  }

  // Parse "V1;H2;G7" → all 1-indexed positions. Skips malformed entries.
  static std::vector<int> parseAllBestPositions_(const String& s)
  {
    std::vector<int> result;
    if (String(s).trim().empty()) return result;

    std::vector<String> parts;
    String(s).split(";", parts);

    for (String& part : parts)
    {
      part.trim();
      Size digit_start = 0;
      while (digit_start < part.size() && isalpha((unsigned char)part[digit_start]))
        ++digit_start;
      if (digit_start >= part.size()) continue;
      try { result.push_back(part.substr(digit_start).toInt()); } catch (...) {}
    }
    return result;
  }

  // Returns true if a b or y ion with the given name spans residue position pos (1-indexed)
  // in a peptide of length pep_len. Used to determine if a fragment ion carries the
  // open-search delta mass when the modification is localised to pos.
  static bool ionSpansPosition_(const String& ion_name, int pos, int pep_len)
  {
    if (ion_name.empty() || pos < 1 || pos > pep_len) return false;

    const char ion_type = ion_name[0];

    // Find the ion number (digits after the type character)
    Size i = 1;
    while (i < ion_name.size() && !isdigit((unsigned char)ion_name[i])) ++i;
    Size j = i;
    while (j < ion_name.size() && isdigit((unsigned char)ion_name[j])) ++j;
    if (j == i) return false;

    int ion_num = 0;
    try { ion_num = ion_name.substr(i, j - i).toInt(); } catch (...) { return false; }

    if (ion_type == 'b') return pos <= ion_num;
    if (ion_type == 'y') return pos >= pep_len - ion_num + 1;
    return false;
  }

  // Format delta mass as an ion label tag: "[+156.1150]" or "[-20.9137]"
  static String formatDeltaTag_(double delta)
  {
    char buf[32];
    snprintf(buf, sizeof(buf), "%.4f", delta);
    String tag = "[";
    if (delta >= 0) tag += "+";
    tag += buf;
    tag += "]";
    return tag;
  }

  // Build OpenMS-compatible sequence from bare peptide + positional delta mods.
  // pos 0 = N-terminal, pos > 0 = 1-indexed residue, pos -1 = C-terminal.
  static String buildOpenSearchSequence_(const String& bare_seq,
                                          const std::map<int, double>& pos_mods)
  {
    String result;

    auto it_nterm = pos_mods.find(0);
    if (it_nterm != pos_mods.end() && std::abs(it_nterm->second) > 0.001)
    {
      const double m = it_nterm->second;
      result += (m >= 0 ? "(+" : "(") + String(m) + ")";
    }

    for (Size i = 0; i < bare_seq.size(); ++i)
    {
      result += bare_seq[i];
      int pos = static_cast<int>(i) + 1;
      auto it = pos_mods.find(pos);
      if (it != pos_mods.end() && std::abs(it->second) > 0.001)
      {
        const double m = it->second;
        result += (m >= 0 ? "[+" : "[") + String(m) + "]";
      }
    }

    auto it_cterm = pos_mods.find(-1);
    if (it_cterm != pos_mods.end() && std::abs(it_cterm->second) > 0.001)
    {
      const double m = it_cterm->second;
      result += (m >= 0 ? "(+" : "(") + String(m) + ")";
    }

    return result;
  }

  std::vector<PeptideEntry> parsePeptideTSV_(const String& filename,
                                        const String& run_id) const
  {
    std::vector<PeptideEntry> entries;

    std::ifstream file(filename.c_str());
    if (!file.is_open())
    {
      OPENMS_LOG_ERROR << "[diaWeaverCounter] Cannot open peptide TSV: " << filename << "\n";
      return entries;
    }

    // --- Parse header row ---
    std::string raw_header;
    if (!std::getline(file, raw_header))
    {
      OPENMS_LOG_ERROR << "[diaWeaverCounter] Empty TSV file: " << filename << "\n";
      return entries;
    }

    std::vector<String> col_names;
    String(raw_header).split("\t", col_names);
    for (auto& c : col_names) c.trim();

    // Map lowercased column name → internal field tag
    const std::map<String, String> aliases = {
      // run identification
      {"run",                        "run"}, {"file.name",               "run"},
      {"spectrum",                   "run_spectrum"},  // MSFragger — needs runIdFromSpectrum_
      // sequence
      {"modified peptide",           "seq"}, {"peptide",                 "seq"},
      {"sequence",                   "seq"}, {"modified.sequence",       "seq"},
      {"stripped.sequence",          "seq"}, {"pep.strippedsequence",    "seq"},
      // charge
      {"charge",                     "charge"}, {"precursor.charge",     "charge"},
      {"eg.precursorcharge",         "charge"},
      // RT — DIA-NN reports minutes; MSFragger reports seconds (tagged separately)
      {"retention",                  "rt_sec"}, // MSFragger — already seconds
      {"rt",                         "rt_min"}, {"eg.apexrt",            "rt_min"},
      {"eg.meanapexrt",              "rt_min"},
      // m/z
      {"calibrated observed m/z",    "mz"},  {"observed m/z",            "mz"},
      {"mz",                         "mz"},  {"precursor.mz",            "mz"},
      {"fg.precmz",                  "mz"},
      // ion mobility
      {"im",                         "im"},  {"ionmobility",             "im"},
      {"eg.ionmobility",             "im"}
    };

    // field tag → column index
    std::map<String, int> field_col;
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
        (field_col.find("rt_sec") == field_col.end() &&
         field_col.find("rt_min") == field_col.end()))
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

    // Determine run-column source and RT unit
    const bool has_run_diann     = field_col.count("run") > 0;
    const bool has_run_spectrum  = field_col.count("run_spectrum") > 0;
    const bool has_run_col       = has_run_diann || has_run_spectrum;
    const bool filter_by_run     = !run_id.empty();
    const bool run_col_is_filename =
      has_run_diann &&
      col_names[field_col.at("run")].toLower() == String("file.name");

    // RT: MSFragger reports seconds ("retention"); DIA-NN reports minutes ("rt*")
    const bool rt_in_seconds = field_col.count("rt_sec") > 0;
    const String rt_field    = rt_in_seconds ? "rt_sec" : "rt_min";

    if (field_col.count("rt_sec") == 0 && field_col.count("rt_min") == 0)
    {
      OPENMS_LOG_WARN << "[diaWeaverCounter] No recognised RT column found. Skipping file.\n";
      return entries;
    }

    if (filter_by_run && !has_run_col)
    {
      OPENMS_LOG_WARN << "[diaWeaverCounter] No run/spectrum column found in TSV; "
                         "run filter '" << run_id << "' cannot be applied. "
                         "All rows will be loaded.\n";
    }
    else if (filter_by_run)
    {
      OPENMS_LOG_INFO << "[diaWeaverCounter] Filtering TSV to run: " << run_id << "\n";
    }

    // --- Parse data rows ---
    std::string raw_line;
    Size row = 1;
    while (std::getline(file, raw_line))
    {
      ++row;
      if (raw_line.empty()) continue;

      std::vector<String> fields;
      String(raw_line).split("\t", fields);

      try
      {
        // --- Run filter ---
        if (filter_by_run && has_run_col)
        {
          String row_run;
          if (has_run_spectrum)
          {
            row_run = runIdFromSpectrum_(fields.at(field_col.at("run_spectrum")).trim());
          }
          else
          {
            row_run = fields.at(field_col.at("run")).trim();
            if (run_col_is_filename) row_run = runIdFromFileName_(row_run);
          }
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
        // MSFragger reports RT in seconds; DIA-NN reports in minutes (×60)
        e.rt       = fields.at(field_col.at(rt_field)).trim().toDouble()
                     * (rt_in_seconds ? 1.0 : 60.0);
        e.mz       = field_col.count("mz") ?
                       fields.at(field_col.at("mz")).trim().toDouble() : 0.0;
        e.im       = field_col.count("im") ?
                       fields.at(field_col.at("im")).trim().toDouble() : 0.0;
        e.id       = "row" + String(row);
        entries.push_back(std::move(e));
      }
      catch (const Exception::BaseException& ex)
      {
        OPENMS_LOG_WARN << "[diaWeaverCounter] Skipping malformed row " << row
                        << ": " << ex.getMessage() << "\n";
      }
      catch (const std::out_of_range&)
      {
        OPENMS_LOG_WARN << "[diaWeaverCounter] Skipping short row " << row << "\n";
      }
    }

    OPENMS_LOG_INFO << "[diaWeaverCounter] Loaded " << entries.size()
                    << " peptide entries from " << filename << "\n";
    return entries;
  }

  // -------------------------------------------------------------------------
  // MSFragger open-search psm.tsv parser
  //
  // Builds the full modified sequence from:
  //   1. Bare Peptide column  (no mods)
  //   2. Assigned Modifications column (static + variable mods, as delta masses)
  //   3. Delta Mass applied at Best Positions (open-search PTM localisation)
  //
  // Rows are filtered by run (Spectrum column) and Q-value (1% FDR).
  // Rows with a significant delta mass (>0.02 Da) but no localisation site are skipped.
  // -------------------------------------------------------------------------
  std::vector<PeptideEntry> parseMSFraggerOpenSearchTSV_(const String& filename,
                                                     const String& run_id) const
  {
    std::vector<PeptideEntry> entries;

    std::ifstream file(filename.c_str());
    if (!file.is_open())
    {
      OPENMS_LOG_ERROR << "[diaWeaverCounter] Cannot open MSFragger psm.tsv: " << filename << "\n";
      return entries;
    }

    std::string raw_header;
    if (!std::getline(file, raw_header)) return entries;

    std::vector<String> col_names;
    String(raw_header).split("\t", col_names);
    for (auto& c : col_names) c.trim();

    const std::map<String, String> aliases = {
      {"spectrum",                  "run_spectrum"},
      {"peptide",                   "seq_bare"},
      {"charge",                    "charge"},
      {"retention",                 "rt_sec"},
      {"calibrated observed m/z",   "mz"},
      {"observed m/z",              "mz"},
      {"assigned modifications",    "assigned_mods"},
      {"delta mass",                "delta_mass"},
      {"best positions",            "best_positions"},
      {"im",                        "im"},
      {"ionmobility",               "im"},
    };

    std::map<String, int> field_col;
    for (Size i = 0; i < col_names.size(); ++i)
    {
      String lower = col_names[i];
      lower.toLower();
      auto it = aliases.find(lower);
      if (it != aliases.end())
        field_col[it->second] = static_cast<int>(i);
    }

    for (const String& req : {"seq_bare", "charge", "rt_sec",
                               "assigned_mods", "delta_mass", "best_positions"})
    {
      if (field_col.count(req) == 0)
      {
        OPENMS_LOG_WARN << "[diaWeaverCounter] MSFragger open-search psm.tsv missing "
                           "required column '" << req << "'. Cannot parse.\n";
        return entries;
      }
    }

    const bool filter_by_run = !run_id.empty() && field_col.count("run_spectrum") > 0;
    if (!run_id.empty() && !filter_by_run)
      OPENMS_LOG_WARN << "[diaWeaverCounter] No Spectrum column; run filter cannot be applied.\n";
    else if (filter_by_run)
      OPENMS_LOG_INFO << "[diaWeaverCounter] Filtering MSFragger open-search psm.tsv to run: "
                      << run_id << "\n";

    std::string raw_line;
    Size row = 1;
    Size n_skipped_run   = 0;
    Size n_skipped_noloc = 0;

    while (std::getline(file, raw_line))
    {
      ++row;
      if (raw_line.empty()) continue;

      std::vector<String> fields;
      String(raw_line).split("\t", fields);

      try
      {
        // Run filter
        if (filter_by_run)
        {
          const String row_run = runIdFromSpectrum_(
            fields.at(field_col.at("run_spectrum")).trim());
          if (row_run != run_id) { ++n_skipped_run; continue; }
        }

        const String bare_seq   = fields.at(field_col.at("seq_bare")).trim();
        const double delta_mass = fields.at(field_col.at("delta_mass")).trim().toDouble();

        // Build positional mod map from Assigned Modifications (static + variable)
        // Do NOT apply delta here — it is applied per-candidate during ion matching.
        const std::map<int, double> assigned_mods =
          parseAssignedMods_(fields.at(field_col.at("assigned_mods")).trim());

        // Parse all localisation candidates from Best Positions
        std::vector<int> cand_positions;
        if (std::abs(delta_mass) > 0.02)
        {
          cand_positions = parseAllBestPositions_(
            fields.at(field_col.at("best_positions")).trim());
          if (cand_positions.empty())
          {
            OPENMS_LOG_WARN << "[diaWeaverCounter] Row " << row
                            << ": significant delta mass " << delta_mass
                            << " Da but no Best Positions localisation. Skipping.\n";
            ++n_skipped_noloc;
            continue;
          }
        }

        PeptideEntry e;
        e.bare_sequence          = bare_seq;
        e.assigned_mods          = assigned_mods;
        e.delta_mass             = delta_mass;
        e.localization_candidates = cand_positions;
        // sequence = base modified sequence (assigned mods only, no delta) — used for display
        e.sequence = buildOpenSearchSequence_(bare_seq, assigned_mods);
        e.charge   = fields.at(field_col.at("charge")).trim().toInt();
        e.rt       = fields.at(field_col.at("rt_sec")).trim().toDouble();
        e.mz       = field_col.count("mz") ?
                       fields.at(field_col.at("mz")).trim().toDouble() : 0.0;
        e.im       = field_col.count("im") ?
                       fields.at(field_col.at("im")).trim().toDouble() : 0.0;
        e.id       = "row" + String(row);

        if (e.charge < 1)
        {
          OPENMS_LOG_WARN << "[diaWeaverCounter] Row " << row
                          << ": invalid charge " << e.charge << ". Skipping.\n";
          continue;
        }

        // Validate base sequence is parseable
        try { AASequence::fromString(e.sequence); }
        catch (const Exception::BaseException& ex)
        {
          OPENMS_LOG_WARN << "[diaWeaverCounter] Row " << row
                          << ": sequence '" << e.sequence
                          << "' failed AASequence validation: " << ex.getMessage()
                          << ". Skipping.\n";
          continue;
        }

        entries.push_back(std::move(e));
      }
      catch (const Exception::BaseException& ex)
      {
        OPENMS_LOG_WARN << "[diaWeaverCounter] Skipping malformed row " << row
                        << ": " << ex.getMessage() << "\n";
      }
      catch (const std::out_of_range&)
      {
        OPENMS_LOG_WARN << "[diaWeaverCounter] Skipping short row " << row << "\n";
      }
    }

    OPENMS_LOG_INFO << "[diaWeaverCounter] MSFragger open-search: loaded " << entries.size()
                    << " peptides (skipped: " << n_skipped_run   << " wrong run, "
                    << n_skipped_noloc << " unlocalized delta) from " << filename << "\n";
    return entries;
  }

  // Returns true if the TSV looks like a MSFragger open-search psm.tsv.
  // Detection: header contains "best positions", "assigned modifications", "delta mass".
  static bool isOpenSearchPSV_(const String& filename)
  {
    std::ifstream f(filename.c_str());
    std::string header;
    if (!std::getline(f, header)) return false;

    std::vector<String> cols;
    String(header).split("\t", cols);

    bool has_best_pos = false, has_assigned = false, has_delta = false;
    for (auto& c : cols)
    {
      String cl = c;
      cl.trim().toLower();
      if (cl == "best positions")          has_best_pos = true;
      if (cl == "assigned modifications")  has_assigned = true;
      if (cl == "delta mass")              has_delta    = true;
    }
    return has_best_pos && has_assigned && has_delta;
  }

  // -------------------------------------------------------------------------
  // Aggregate a single spectrum with its RT neighbors using Gaussian weighting.
  // Reads from the immutable source experiment; safe to call from a parallel loop
  // that writes picked results into a separate output experiment.
  // Mirrors diaWeaver's aggregateSpectrum_ exactly.
  // -------------------------------------------------------------------------
  static void aggregateSpectrum_(
    const MSExperiment& exp,
    Size center_idx,
    const PeakPickerIM& picker,
    MSSpectrum& out)
  {
    if (center_idx >= exp.size()) return;

    Param params = picker.getParameters();
    double fwhm   = (double)params.getValue("aggregation:rt_FWHM");
    double cutoff = (double)params.getValue("aggregation:cutoff");
    double factor = -4.0 * std::log(2.0) / (fwhm * fwhm);

    double center_rt = exp[center_idx].getRT();

    std::vector<MSSpectrum> spectra_to_aggregate;
    std::vector<double> weights;

    // Search forward (including center)
    for (Size j = center_idx; j < exp.size(); ++j)
    {
      double rt_diff = exp[j].getRT() - center_rt;
      double weight  = std::exp(factor * rt_diff * rt_diff);
      if (weight < cutoff && j != center_idx) break;
      spectra_to_aggregate.push_back(exp[j]);
      weights.push_back(weight);
    }

    // Search backward
    for (SignedSize j = static_cast<SignedSize>(center_idx) - 1; j >= 0; --j)
    {
      double rt_diff = exp[j].getRT() - center_rt;
      double weight  = std::exp(factor * rt_diff * rt_diff);
      if (weight < cutoff) break;
      spectra_to_aggregate.push_back(exp[j]);
      weights.push_back(weight);
    }

    double sum_w = 0.0;
    for (double w : weights) sum_w += w;
    for (double& w : weights) w /= sum_w;

    picker.aggregateScans(spectra_to_aggregate, weights, out);
  }

  // -------------------------------------------------------------------------
  // Mirrors diaWeaver's runMassTraceExtractor_: MTD + optional EPD.
  // Spectra in ms_peakmap must already be MSLevel=1 (extractSingleMS2Window
  // guarantees this). Traces are appended to traces_out.
  // -------------------------------------------------------------------------
  void runMassTraceExtractor_(MSExperiment& ms_peakmap,
                              const Param& common_param,
                              Param mtd_param,
                              Param epd_param,
                              std::vector<MassTrace>& traces_out)
  {
    if (ms_peakmap.empty()) return;
    ms_peakmap.sortSpectra(true);

    MassTraceDetection mtdet;
    mtd_param.insert("", common_param);
    mtd_param.remove("chrom_fwhm");
    mtdet.setParameters(mtd_param);

    std::vector<MassTrace> m_traces;
    mtdet.run(ms_peakmap, m_traces);
    if (m_traces.empty()) return;

    const bool use_epd = epd_param.getValue("enabled").toBool();
    if (use_epd)
    {
      epd_param.remove("enabled");
      epd_param.insert("", common_param);
      epd_param.remove("noise_threshold_int");
      ElutionPeakDetection epdet;
      epdet.setParameters(epd_param);

      std::vector<MassTrace> split;
      epdet.detectPeaks(m_traces, split);

      std::vector<MassTrace> final_traces;
      if (epdet.getParameters().getValue("width_filtering") == "auto")
        epdet.filterByPeakWidth(split, final_traces);
      else
        final_traces = std::move(split);

      for (auto& t : final_traces)
        if (t.getSize() > 0) traces_out.push_back(std::move(t));
    }
    else
    {
      for (auto& t : m_traces)
        if (t.getSize() > 0) traces_out.push_back(std::move(t));
    }
  }

  // -------------------------------------------------------------------------

  ExitCodes main_(int, const char**) override
  {
    const String in_file          = getStringOption_("in");
    const String in_ids_file      = getStringOption_("in_ids");
    const String out_traces_file  = getStringOption_("out_traces");
    const String out_summary_file = getStringOption_("out_summary");
    const String out_pep_file     = getStringOption_("out_peptides");
    const String out_coll_file    = getStringOption_("out_collisions");

    const double mz_tol_ppm    = getDoubleOption_("mz_tolerance");
    const double rt_tol        = getDoubleOption_("rt_tolerance");
    const double im_tol        = getDoubleOption_("im_tolerance");
    const bool   aggregate_scans = getFlag_("aggregate_across_scans");

    // ------------------------------------------------------------------
    // Step 1: Determine DIA windows via metadata scan (no full load).
    //         Uses OnDiscMSExperiment for memory-efficient access,
    //         exactly as diaWeaver does.
    // ------------------------------------------------------------------
    OPENMS_LOG_INFO << "Opening: " << in_file << "\n";
    OnDiscMSExperiment on_disc;
    if (!on_disc.openFile(in_file))
    {
      OPENMS_LOG_ERROR << "Failed to open as indexed mzML: " << in_file << "\n";
      return INPUT_FILE_NOT_FOUND;
    }

    DiaWeaver::WindowMap windows;
    DiaWeaver::determineWindows(on_disc, windows);
    const DiaWeaver::IMInfo im_info = DiaWeaver::determineIMInfo(on_disc, windows);

    OPENMS_LOG_INFO << "Found " << windows.size() << " DIA windows.\n";
    if (windows.empty())
    {
      OPENMS_LOG_ERROR << "No DIA windows found in " << in_file
                       << ". Is this a DIA mzML?\n";
      return INCOMPATIBLE_INPUT_DATA;
    }

    if (im_info.available)
    {
      OPENMS_LOG_INFO << "Ion mobility data detected. Using PeakPickerIM.\n";
      if (aggregate_scans)
        OPENMS_LOG_INFO << "Aggregation across scans enabled. "
                           "Signal will be boosted before peak picking.\n";
    }
    else
    {
      OPENMS_LOG_INFO << "No ion mobility data detected. Using PeakPickerHiRes.\n";
      if (aggregate_scans)
        OPENMS_LOG_WARN << "aggregate_across_scans is set but no IM data detected. "
                           "Aggregation will be skipped.\n";
    }

    // ------------------------------------------------------------------
    // Step 2: Parse peptide IDs before the window loop so we can assign
    //         each peptide to its window by precursor m/z.
    // ------------------------------------------------------------------
    String run_id = File::basename(in_file);
    for (const char* ext : {".mzML", ".mzml"})
    {
      if (run_id.hasSuffix(ext)) { run_id = run_id.prefix(run_id.size() - strlen(ext)); break; }
    }
    OPENMS_LOG_INFO << "Run ID derived from mzML: " << run_id << "\n";

    std::vector<PeptideEntry> peptides = isOpenSearchPSV_(in_ids_file)
      ? parseMSFraggerOpenSearchTSV_(in_ids_file, run_id)
      : parsePeptideTSV_(in_ids_file, run_id);
    OPENMS_LOG_INFO << "Peptides to map: " << peptides.size() << "\n";

    // Theoretical spectrum generator — configured once, reused across windows
    TheoreticalSpectrumGenerator tsg;
    {
      Param p = tsg.getDefaults();
      p.setValue("add_metainfo",    "true",  "");
      p.setValue("add_a_ions",      "false", "");
      p.setValue("add_c_ions",      "false", "");
      p.setValue("add_x_ions",      "false", "");
      p.setValue("add_z_ions",      "false", "");
      p.setValue("add_losses",      "true",  "");
      p.setValue("add_term_losses", "true",  "");
      tsg.setParameters(p);
    }

    // Global accumulators — grown one window at a time
    std::vector<MassTrace>  ms2_traces;
    std::vector<double>     trace_snr;
    std::vector<std::vector<std::pair<Size, String>>> trace_claims;
    std::vector<bool>       trace_within_pep_collision;

    // Per-peptide result vectors — indexed by global peptide index, filled across windows
    std::vector<Size>           pep_n_ions(peptides.size(), 0);
    std::vector<Size>           pep_matched_ions(peptides.size(), 0);
    std::vector<Size>           pep_multi_trace_ions(peptides.size(), 0);
    std::vector<std::set<Size>> pep_matched_trace_set(peptides.size());
    std::map<std::pair<Size,String>, IonBestMatch> best_ion_match;

    ElutionPeakDetection epd_snr;  // stateless; used only for computeApexSNR

    const Param mte_common = getParam_().copy("MassTraceExtractor:common:", true);
    const Param mte_mtd    = getParam_().copy("MassTraceExtractor:mtd:",    true);
    const Param mte_epd    = getParam_().copy("MassTraceExtractor:epd:",    true);
    const Param ppim_param = getParam_().copy("PeakPickerIM:",              true);
    const Param pphr_param = getParam_().copy("PeakPickerHiRes:",           true);

    // ------------------------------------------------------------------
    // Step 3: Per-window loop — parallel over windows.
    //
    //   Heavy work (thread-local, no guards):
    //     a) Extract + peak pick + MTD+EPD  → win_traces, win_snr
    //     b) Build per-window KD-tree       → local_kd  (local indices)
    //     c) Filter peptides by precursor m/z
    //     d) Map fragment ions → local traces; accumulate into thread-local
    //        win_pep_n_ions, win_pep_matched_ions, win_pep_multi_trace_ions,
    //        win_pep_trace_sets (local indices), win_best_ion_match (local best_trace_idx)
    //
    //   Accumulation (guarded by critical section):
    //     e) Compute global_offset = ms2_traces.size()
    //        Translate local→global indices, merge into global accumulators
    //
    //   Steps 6-10 run sequentially after the loop on the completed globals.
    // ------------------------------------------------------------------

    // Convert map to vector for OpenMP indexed access
    std::vector<std::pair<DiaWeaver::DIAWindow, std::vector<Size>>> window_vec(
      windows.begin(), windows.end());
    const SignedSize total_windows = static_cast<SignedSize>(window_vec.size());

#ifdef _OPENMP
    omp_set_num_threads(getIntOption_("threads"));
    OPENMS_LOG_INFO << "Processing " << total_windows << " windows with "
                    << getIntOption_("threads") << " thread(s).\n";
#endif

    bool had_error = false;  // set inside critical section if a window hits a fatal error

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) firstprivate(on_disc)
#endif
    for (SignedSize widx = 0; widx < total_windows; ++widx)
    {
#ifdef _OPENMP
      if (had_error) continue;  // propagate fatal error across parallel iterations
#endif
      const DiaWeaver::DIAWindow& w       = window_vec[widx].first;
      const std::vector<Size>&    indices = window_vec[widx].second;

      // --- 3a: Extract + peak pick + MTD+EPD (all thread-local) ---
      MSExperiment win_exp;
      DiaWeaver::extractSingleMS2Window(on_disc, w, indices, im_info, win_exp);
      if (win_exp.empty()) continue;

      if (aggregate_scans && im_info.available)
      {
        // Aggregate each spectrum with its RT neighbors (Gaussian weighting),
        // then pick IM traces.  win_exp stays immutable as the source so that
        // aggregateSpectrum_ always reads unmodified neighbors.
        MSExperiment win_exp_picked;
        win_exp_picked.resize(win_exp.size());

        PeakPickerIM picker_im;
        picker_im.setParameters(ppim_param);

        for (Size s = 0; s < win_exp.size(); ++s)
        {
          if (win_exp[s].getIMPeakType() != IMPeakType::IM_CENTROIDED)
          {
            MSSpectrum aggregated;
            aggregateSpectrum_(win_exp, s, picker_im, aggregated);
            picker_im.pickIMTraces(aggregated);
            win_exp_picked[s] = std::move(aggregated);
          }
          else
          {
            win_exp_picked[s] = win_exp[s];
          }
        }
        win_exp = std::move(win_exp_picked);
      }
      else
      {
        PeakPickerIM    picker_im;
        PeakPickerHiRes picker_hr;
        if (im_info.available)
          picker_im.setParameters(ppim_param);
        else
          picker_hr.setParameters(pphr_param);

        for (auto& spec : win_exp)
        {
          if (im_info.available)
          {
            if (spec.getIMPeakType() != IMPeakType::IM_CENTROIDED)
              picker_im.pickIMTraces(spec);
          }
          else
          {
            if (spec.getType(true) != SpectrumSettings::SpectrumType::CENTROID)
            {
              MSSpectrum picked;
              picker_hr.pick(spec, picked);
              spec = std::move(picked);
            }
          }
        }
      }

      std::vector<MassTrace> win_traces;
      {
        Param mtd_copy = mte_mtd;
        Param epd_copy = mte_epd;
        runMassTraceExtractor_(win_exp, mte_common, mtd_copy, epd_copy, win_traces);
      }
      if (win_traces.empty()) continue;

      // --- 3b: Per-window SNR + local KD-tree ---
      std::vector<double> win_snr(win_traces.size(), 0.0);
      bool win_has_im = false;
      TraceKDTree local_kd;

      for (Size i = 0; i < win_traces.size(); ++i)
      {
        if (!win_traces[i].getSmoothedIntensities().empty())
          win_snr[i] = epd_snr.computeApexSNR(win_traces[i]);
        if (win_traces[i].containsIMData()) win_has_im = true;

        MS2TraceNode node;
        node.rt        = win_traces[i].getCentroidRT();
        node.mz        = win_traces[i].getCentroidMZ();
        node.im        = win_traces[i].containsIMData() ?
                           win_traces[i].getCentroidIM() : 0.0;
        node.trace_idx = i;  // LOCAL
        local_kd.insert(node);
      }
      local_kd.optimise();

      // --- 3c: Filter peptides to this window by precursor m/z ---
      std::vector<Size> win_pep_indices;
      for (Size pi = 0; pi < peptides.size(); ++pi)
      {
        double pmz = peptides[pi].mz;
        if (pmz <= 0.0)
        {
          try { pmz = AASequence::fromString(peptides[pi].sequence).getMZ(peptides[pi].charge); }
          catch (...) { continue; }
        }
        if (pmz >= w.lower_mz && pmz <= w.upper_mz)
          win_pep_indices.push_back(pi);
      }

      // --- 3d: Thread-local claim accumulators ---
      std::vector<std::vector<std::pair<Size, String>>> win_trace_claims(win_traces.size());
      std::vector<bool>       win_within_pep_collision(win_traces.size(), false);
      std::map<Size, String>  pep_trace_first_ion;

      // Thread-local per-peptide accumulators (indexed by global pep_idx)
      std::vector<Size>           win_pep_n_ions(peptides.size(), 0);
      std::vector<Size>           win_pep_matched(peptides.size(), 0);
      std::vector<Size>           win_pep_multi(peptides.size(), 0);
      std::vector<std::set<Size>> win_pep_trace_sets(peptides.size());  // LOCAL indices
      std::map<std::pair<Size,String>, IonBestMatch> win_best_ion_match; // LOCAL best_trace_idx

      bool win_had_error = false;

      for (Size pep_idx : win_pep_indices)
      {
        if (win_had_error) break;
        const PeptideEntry& pep = peptides[pep_idx];

        if (pep.charge < 1)
        {
          OPENMS_LOG_ERROR << "[diaWeaverCounter] Invalid charge " << pep.charge
                           << " for peptide '" << pep.sequence << "' (" << pep.id
                           << "). Charge must be >= 1.\n";
          win_had_error = true;
          break;
        }

        pep_trace_first_ion.clear();

        // Query local KD-tree — returns LOCAL trace indices
        auto queryTraces = [&](double mz_center) -> std::vector<Size>
        {
          const double mz_tol_da = mz_center * mz_tol_ppm * 1e-6;
          TraceKDTree::_Region_ region;
          region._M_low_bounds[0]  = pep.rt - rt_tol;
          region._M_high_bounds[0] = pep.rt + rt_tol;
          region._M_low_bounds[1]  = mz_center - mz_tol_da;
          region._M_high_bounds[1] = mz_center + mz_tol_da;

          std::vector<MS2TraceNode> candidates;
          local_kd.find_within_range(region, back_inserter(candidates));

          std::vector<Size> result;
          for (const auto& cand : candidates)
          {
            if (win_has_im && pep.im > 0.0 &&
                win_traces[cand.trace_idx].containsIMData())
            {
              if (std::abs(cand.im - pep.im) > im_tol) continue;
            }
            result.push_back(cand.trace_idx);  // LOCAL
          }
          return result;
        };

        // Write only to thread-local structures; global translation deferred to critical section
        auto recordMatches = [&](const String& ion_label, const std::vector<Size>& local_tidx_list)
        {
          auto ion_key = std::make_pair(pep_idx, ion_label);
          IonBestMatch& best = win_best_ion_match[ion_key];
          best.n_traces += local_tidx_list.size();

          for (Size local_tidx : local_tidx_list)
          {
            win_trace_claims[local_tidx].emplace_back(pep_idx, ion_label);
            win_pep_trace_sets[pep_idx].insert(local_tidx);  // LOCAL

            if (win_snr[local_tidx] > best.best_snr)
            {
              best.best_snr       = win_snr[local_tidx];
              best.best_trace_idx = local_tidx;  // LOCAL — translated to global in critical section
            }

            auto fit = pep_trace_first_ion.find(local_tidx);
            if (fit == pep_trace_first_ion.end())
              pep_trace_first_ion[local_tidx] = ion_label;
            else if (fit->second != ion_label)
              win_within_pep_collision[local_tidx] = true;
          }
        };

        if (pep.localization_candidates.empty())
        {
          // Standard path (DIA-NN / MSFragger closed-search)
          AASequence aa_seq;
          try { aa_seq = AASequence::fromString(pep.sequence); }
          catch (const Exception::BaseException& ex)
          {
            OPENMS_LOG_WARN << "[diaWeaverCounter] Cannot parse sequence '"
                            << pep.sequence << "': " << ex.getMessage() << ". Skipping.\n";
            continue;
          }
          if (aa_seq.size() < 2) continue;

          PeakSpectrum theo_spec;
          tsg.getSpectrum(theo_spec, aa_seq, 1, pep.charge);

          const PeakSpectrum::StringDataArray&  ion_names    =
            theo_spec.getStringDataArrays().at(0);
          const PeakSpectrum::IntegerDataArray& frag_charges =
            theo_spec.getIntegerDataArrays().at(0);

          win_pep_n_ions[pep_idx] = theo_spec.size() * 4 + 5;

          for (Size ion_i = 0; ion_i < theo_spec.size(); ++ion_i)
          {
            const double mz_mono  = theo_spec[ion_i].getMZ();
            const String ion_name = ion_names[ion_i];
            const int    frag_z   = frag_charges[ion_i];

            const std::vector<Size> mono_tidxs = queryTraces(mz_mono);
            if (mono_tidxs.empty()) continue;

            ++win_pep_matched[pep_idx];
            if (mono_tidxs.size() > 1) ++win_pep_multi[pep_idx];
            recordMatches(ion_name + "[M+0]", mono_tidxs);

            for (int iso = 1; iso <= 3; ++iso)
            {
              const double mz_iso = mz_mono + iso * Constants::C13C12_MASSDIFF_U / frag_z;
              const std::vector<Size> tidxs = queryTraces(mz_iso);
              if (tidxs.empty()) continue;
              ++win_pep_matched[pep_idx];
              if (tidxs.size() > 1) ++win_pep_multi[pep_idx];
              recordMatches(ion_name + "[M+" + String(iso) + "]", tidxs);
            }
          }

          const double mz_prec_mono = aa_seq.getMZ(pep.charge);
          const std::vector<Size> prec_tidxs = queryTraces(mz_prec_mono);
          if (!prec_tidxs.empty())
          {
            ++win_pep_matched[pep_idx];
            if (prec_tidxs.size() > 1) ++win_pep_multi[pep_idx];
            recordMatches("p[M+0]", prec_tidxs);

            for (int iso = 1; iso <= 4; ++iso)
            {
              const double mz_iso = mz_prec_mono + iso * Constants::C13C12_MASSDIFF_U / pep.charge;
              const std::vector<Size> tidxs = queryTraces(mz_iso);
              if (tidxs.empty()) continue;
              ++win_pep_matched[pep_idx];
              if (tidxs.size() > 1) ++win_pep_multi[pep_idx];
              recordMatches("p[M+" + String(iso) + "]", tidxs);
            }
          }
        }
        else
        {
          // Open-search path
          const String delta_tag = formatDeltaTag_(pep.delta_mass);
          const int    pep_len   = static_cast<int>(pep.bare_sequence.size());
          std::set<String> seen_ion_labels;

          for (int cand_pos : pep.localization_candidates)
          {
            std::map<int,double> cand_mods = pep.assigned_mods;
            cand_mods[cand_pos] += pep.delta_mass;

            AASequence cand_aa_seq;
            try { cand_aa_seq = AASequence::fromString(
                    buildOpenSearchSequence_(pep.bare_sequence, cand_mods)); }
            catch (...) { continue; }
            if (cand_aa_seq.size() < 2) continue;

            PeakSpectrum cand_spec;
            tsg.getSpectrum(cand_spec, cand_aa_seq, 1, pep.charge);

            const PeakSpectrum::StringDataArray&  cand_names   =
              cand_spec.getStringDataArrays().at(0);
            const PeakSpectrum::IntegerDataArray& cand_charges =
              cand_spec.getIntegerDataArrays().at(0);

            for (Size ion_i = 0; ion_i < cand_spec.size(); ++ion_i)
            {
              const double mz_mono   = cand_spec[ion_i].getMZ();
              const String base_name = cand_names[ion_i];
              const int    frag_z    = cand_charges[ion_i];

              const bool   is_mod    = ionSpansPosition_(base_name, cand_pos, pep_len);
              const String ion_label = is_mod ? base_name + delta_tag : base_name;

              if (seen_ion_labels.count(ion_label)) continue;
              seen_ion_labels.insert(ion_label);
              win_pep_n_ions[pep_idx] += 4;

              const std::vector<Size> mono_tidxs = queryTraces(mz_mono);
              if (mono_tidxs.empty()) continue;

              ++win_pep_matched[pep_idx];
              if (mono_tidxs.size() > 1) ++win_pep_multi[pep_idx];
              recordMatches(ion_label + "[M+0]", mono_tidxs);

              for (int iso = 1; iso <= 3; ++iso)
              {
                const double mz_iso = mz_mono + iso * Constants::C13C12_MASSDIFF_U / frag_z;
                const std::vector<Size> tidxs = queryTraces(mz_iso);
                if (tidxs.empty()) continue;
                ++win_pep_matched[pep_idx];
                if (tidxs.size() > 1) ++win_pep_multi[pep_idx];
                recordMatches(ion_label + "[M+" + String(iso) + "]", tidxs);
              }
            }
          }

          // Precursor — invariant across candidates; compute from first
          {
            std::map<int,double> first_mods = pep.assigned_mods;
            first_mods[pep.localization_candidates[0]] += pep.delta_mass;
            win_pep_n_ions[pep_idx] += 5;
            try
            {
              const AASequence prec_seq = AASequence::fromString(
                buildOpenSearchSequence_(pep.bare_sequence, first_mods));
              const double mz_prec_mono = prec_seq.getMZ(pep.charge);
              const std::vector<Size> prec_tidxs = queryTraces(mz_prec_mono);
              if (!prec_tidxs.empty())
              {
                ++win_pep_matched[pep_idx];
                if (prec_tidxs.size() > 1) ++win_pep_multi[pep_idx];
                recordMatches("p[M+0]", prec_tidxs);
                for (int iso = 1; iso <= 4; ++iso)
                {
                  const double mz_iso = mz_prec_mono +
                    iso * Constants::C13C12_MASSDIFF_U / pep.charge;
                  const std::vector<Size> tidxs = queryTraces(mz_iso);
                  if (tidxs.empty()) continue;
                  ++win_pep_matched[pep_idx];
                  if (tidxs.size() > 1) ++win_pep_multi[pep_idx];
                  recordMatches("p[M+" + String(iso) + "]", tidxs);
                }
              }
            }
            catch (...) {}
          }
        }
      }  // end peptide loop

      // --- 3e: Guarded accumulation into global structures ---
#ifdef _OPENMP
#pragma omp critical (accumulate_results)
#endif
      {
        if (win_had_error)
        {
          had_error = true;
        }
        else
        {
          // Atomically claim the global offset for this window's traces
          const Size global_offset = ms2_traces.size();

          // Translate local→global in win_best_ion_match, then merge
          for (auto& [key, val] : win_best_ion_match)
          {
            val.best_trace_idx += global_offset;
            auto& g = best_ion_match[key];
            g.n_traces += val.n_traces;
            if (val.best_snr > g.best_snr)
            {
              g.best_snr       = val.best_snr;
              g.best_trace_idx = val.best_trace_idx;
            }
          }

          // Translate local→global pep_matched_trace_set and merge per-peptide scalars
          for (Size pi : win_pep_indices)
          {
            pep_n_ions[pi]           += win_pep_n_ions[pi];
            pep_matched_ions[pi]     += win_pep_matched[pi];
            pep_multi_trace_ions[pi] += win_pep_multi[pi];
            for (Size local_tidx : win_pep_trace_sets[pi])
              pep_matched_trace_set[pi].insert(global_offset + local_tidx);
          }

          // Append trace data
          for (auto& t : win_traces)  ms2_traces.push_back(std::move(t));
          for (double s : win_snr)    trace_snr.push_back(s);
          for (auto& c : win_trace_claims)         trace_claims.push_back(std::move(c));
          for (bool  f : win_within_pep_collision) trace_within_pep_collision.push_back(f);

          OPENMS_LOG_INFO << "Window [" << w.lower_mz << "-" << w.upper_mz << " m/z]: "
                          << win_traces.size() << " traces, "
                          << win_pep_indices.size() << " peptides.\n";
        }
      }  // end critical section
    }  // end window loop

    if (had_error) return INPUT_FILE_CORRUPT;

    OPENMS_LOG_INFO << "Total MS2 traces: " << ms2_traces.size() << "\n";

    // ------------------------------------------------------------------
    // Step 6: Derive per-trace and per-peptide aggregate metrics
    // ------------------------------------------------------------------

    // Number of unique peptides claiming each trace
    std::vector<Size> trace_n_peptides(ms2_traces.size(), 0);
    for (Size i = 0; i < ms2_traces.size(); ++i)
    {
      std::set<Size> unique_peps;
      for (const auto& [pi, ion] : trace_claims[i])
        unique_peps.insert(pi);
      trace_n_peptides[i] = unique_peps.size();
    }

    // Per-peptide exclusive vs shared trace counts
    std::vector<Size> pep_exclusive(peptides.size(), 0);
    std::vector<Size> pep_shared(peptides.size(), 0);
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

    const Size n_total     = ms2_traces.size();
    const Size n_explained = n_unique + n_ambiguous;

    auto pct = [](Size num, Size denom) -> double {
      return denom > 0 ? 100.0 * static_cast<double>(num) / static_cast<double>(denom) : 0.0;
    };

    OPENMS_LOG_INFO << "\n=== diaWeaverCounter Summary ===\n"
      << "Explained: " << n_explained << " / " << n_total
      << " traces (" << pct(n_explained, n_total) << " %)"
      << "  |  Orphan: " << n_orphan << " (" << pct(n_orphan, n_total) << " %)\n"
      << "See " << out_summary_file << " for full metrics.\n";

    // ------------------------------------------------------------------
    // Step 7: Write run-level summary TSV
    // ------------------------------------------------------------------
    {
      std::ofstream ofs(out_summary_file.c_str());
      if (!ofs.is_open())
      {
        OPENMS_LOG_ERROR << "Cannot write summary output: " << out_summary_file << "\n";
        return CANNOT_WRITE_OUTPUT_FILE;
      }

      ofs << "metric\tvalue\tpct_of_total_traces\n"
          << "total_ms2_traces\t"      << n_total            << "\t\n"
          << "explained_traces\t"      << n_explained        << "\t" << pct(n_explained, n_total) << "\n"
          << "unique_traces\t"         << n_unique           << "\t" << pct(n_unique,    n_total) << "\n"
          << "ambiguous_traces\t"      << n_ambiguous        << "\t" << pct(n_ambiguous, n_total) << "\n"
          << "unexplained_traces\t"    << n_orphan           << "\t" << pct(n_orphan,    n_total) << "\n"
          << "total_peptides_mapped\t" << peptides.size()    << "\t\n";
    }

    // ------------------------------------------------------------------
    // Step 8: Write per-trace accountability TSV
    // ------------------------------------------------------------------
    {
      std::ofstream ofs(out_traces_file.c_str());
      if (!ofs.is_open())
      {
        OPENMS_LOG_ERROR << "Cannot write trace output: " << out_traces_file << "\n";
        return CANNOT_WRITE_OUTPUT_FILE;
      }

      ofs << "trace_idx\tcentroid_mz\tcentroid_rt\tcentroid_im\tapex_snr\t"
             "n_claiming_peptides\tn_claiming_ions\tstatus\twithin_pep_ion_collision\n";

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
            << status                      << "\t"
            << (trace_within_pep_collision[i] ? "true" : "false") << "\n";
      }
    }

    // ------------------------------------------------------------------
    // Step 9: Write per-peptide coverage TSV (optional)
    // ------------------------------------------------------------------
    if (!out_pep_file.empty())
    {
      std::ofstream ofs(out_pep_file.c_str());
      if (!ofs.is_open())
      {
        OPENMS_LOG_ERROR << "Cannot write peptide output: " << out_pep_file << "\n";
        return CANNOT_WRITE_OUTPUT_FILE;
      }

      ofs << "peptide_id\tsequence\tcharge\trt\t"
             "n_theoretical_ions\tn_matched_ions\tn_multi_trace_ions\t"
             "n_exclusive_traces\tn_shared_traces\tcoverage_pct\t"
             "n_localization_candidates\tlocalization_candidates\n";

      for (Size pi = 0; pi < peptides.size(); ++pi)
      {
        const auto& pep     = peptides[pi];
        const Size  n_theo  = pep_n_ions[pi];
        const double cov    = n_theo > 0 ?
                                100.0 * static_cast<double>(pep_matched_ions[pi]) / n_theo
                                : 0.0;

        // Format localization candidates as "5;6;7" (numeric positions)
        String cand_str;
        for (Size ci = 0; ci < pep.localization_candidates.size(); ++ci)
        {
          if (ci > 0) cand_str += ";";
          cand_str += String(pep.localization_candidates[ci]);
        }

        ofs << pep.id                                    << "\t"
            << pep.sequence                              << "\t"
            << pep.charge                                << "\t"
            << pep.rt                                    << "\t"
            << n_theo                                    << "\t"
            << pep_matched_ions[pi]                      << "\t"
            << pep_multi_trace_ions[pi]                  << "\t"
            << pep_exclusive[pi]                         << "\t"
            << pep_shared[pi]                            << "\t"
            << cov                                       << "\t"
            << pep.localization_candidates.size()        << "\t"
            << cand_str                                  << "\n";
      }
    }

    // ------------------------------------------------------------------
    // Step 10: Write collision TSV (optional)
    //
    // Long format: one row per (trace, peptide, ion) for every trace
    // claimed by more than one unique peptide.
    //
    // is_best = true  when this trace has the highest apex SNR among all
    //                 traces that matched (peptide_id, ion_name).
    // ------------------------------------------------------------------
    if (!out_coll_file.empty())
    {
      std::ofstream ofs(out_coll_file.c_str());
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
          const auto ion_key = std::make_pair(pi, ion_name);
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
