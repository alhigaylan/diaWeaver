// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mohammed Alhigaylan $
// $Authors: Mohammed Alhigaylan $
// --------------------------------------------------------------------------

/**
@page TOPP_diaWeaverIonAccount diaWeaverIonAccount

@brief Reads MSFragger PSMs (psm.tsv), translates each peptide into an OpenMS-compatible
       modified sequence, resolves each PSM's exact source spectrum in the original
       diaWeaver pseudo spectra mzML, and greedily claims theoretical b/y fragment matches
       in MSFragger Hyperscore order so that fragment reuse across competing PSMs can be
       tracked and penalized.

<B>Sequence translation.</B> MSFragger's "Modified Peptide" column omits fixed
modifications entirely (e.g. fixed Carbamidomethyl-C never appears there), so this tool
does not use it. Instead every sequence is rebuilt from the bare "Peptide" column plus
the "Assigned Modifications" column, which lists every modification (fixed and variable)
as a raw numeric mass delta at a position -- never a PTM name or UniMod accession.
Translation uses OpenMS's delta-mass bracket notation:
  - internal residue : "X[+delta]"        e.g. "C[+57.02146]"
  - N-terminal        : "[+delta]" prefix  e.g. "[+42.0106]SEQUENCE..."
  - C-terminal        : "[+delta]." suffix on the last residue

OpenMS does not require a delta mass to match a cataloged UniMod entry: when it doesn't,
it creates an unnamed modification carrying the exact literal mass and uses it for ion
generation identically to a cataloged one. This was verified empirically (multiple
modifications combine additively; N-/C-terminal and internal deltas each affect exactly
the b/y ions that span them) before relying on it here.

<B>Fragment ion charge (max_fragment_charge, default 2).</B> Theoretical b/y ions are
generated up to min(max_fragment_charge, precursor charge), not the full precursor charge.
MSFragger's own fragment charge cap (fragger.params) is a small fixed value independent of
precursor charge, not precursor_charge-1 -- matching ion charges above that would compare
against fragment charge states MSFragger's own search never considered, inflating matched
ion counts and recomputed hyperscore for higher-charge PSMs specifically (verified: capping
at 2 instead of the full precursor charge dropped one real charge-4 PSM's Internal Initial
Hyperscore from 109.6 to 88.3, closing about 40% of its gap to MSFragger's own reported
score of 58.9).

<B>Spectrum resolution.</B> psm.tsv's "Spectrum" field encodes a 0-based scan number that
maps to the original mzML's 1-based "scan=N" native ID as scan+1 (verified empirically).
Resolution is cross-checked against RT and precursor m/z (not charge -- MSFragger can
reassign charge state, so it is not a reliable identity check) and any PSM that fails the
check is logged and skipped rather than silently mismatched.

<B>Greedy fragment claiming.</B> PSMs are sorted by their original MSFragger Hyperscore
(descending), tie-broken by pseudo spectrum precursor intensity (descending, a real
quality signal rather than an arbitrary tie-break), then by row number as a final
deterministic fallback. In that order, each PSM attempts to claim every (window_id,
fragment_id) trace its matched ions touch, via OpenMS::FragmentClaimRegistry (first
successful claimer wins; every trace can be claimed by at most one PSM as the registry's
own official owner). Because psm.tsv contains exactly one row per spectrum (FragPipe's
output_report_topN=1), there is no same-sequence/same-spectrum exemption case to handle.

The same peptide sequence can be identified across several pseudo spectra in the same RT x
IM dimension. This could be due to duplicates in peak picking, or isotopic errors. If the
PSMs use the same fragment ions (by their window_id, fragment_id), then they are forced to
compete and only best scoring PSM wins. We include an exception for modified peptides that
can co-elute with their naked peptide (e.g. oxidation). It can be difficult to determine
whether the modified peptide exists or if it is re-using fragment ions from the naked
peptide. By default, we allow the co-eluting modified peptide to escape the fragment ion
competition. Users can toggle this feature off if they wish.

<B>Peptidoform fragment sharing (allow_peptidoform_fragment_sharing, default true).</B> A
naked and a modified form of the same backbone peptide, or two differently-modified forms
of the same backbone, produce chemically identical b/y ions for any fragment that doesn't
span a position where they differ -- it is not a coincidental collision, it is the same
physical fragment. When enabled, a trace already claimed by another PSM is still counted
toward a querying PSM's recomputed score if the two share the same backbone AND have a
genuinely different full sequence (bare_sequence equal, but openms_sequence -- backbone
plus modification state -- different). That second condition is the whole gate: two PSMs
with the same backbone that are both fully unmodified always have identical full
sequences, so requiring a different full sequence already rules out the naked/naked case
with no separate "is either one modified" check needed, and it additionally rules out what
that naive check alone would miss -- the exact same modified form re-observed at a
different spectrum. Either kind of duplicate identification (naked or modified) must
compete for evidence like any other pair; only genuinely distinct peptidoforms of the same
backbone get the pass.

<B>Hyperscore recomputation.</B> log1p(dot_product) + 2*lnFactorial(i_min) +
lnFactorial(i_max, i_min+1), where i_min/i_max are the smaller/larger of the retained
b-ion and y-ion counts and dot_product is a single combined sum of retained matched
intensities (b and y pooled together, not logged separately). This is the exact formula
OpenMS::HyperScore::computeWithDetail uses internally for ProSE/diaWeaverPeptide's own
PSM scoring -- adopted here in place of an earlier from-the-paper reconstruction
(ln(Nb!) + ln(Ny!) + ln(sum_Ib) + ln(sum_Iy), b/y logged separately) after comparing both
against real MSFragger scores on this dataset: the ProSE-style combined-sum version
tracked MSFragger's actual reported hyperscore noticeably more closely (e.g. one PSM:
ProSE-style 65.86 vs MSFragger's 62.23, versus 76.41 from the separate-log version).
The hyperscore formula attempts to replicate MSFragger search engine scoring scheme. The
current formula does not exactly reproduce the hyperscore computed by MSFragger. On a
small test file, it was off by as high as ~2.50 -- that note was written against the
earlier separate-log formula and has not been re-measured against this one, which is
already known to track MSFragger more closely on the example above. A series with zero
retained matches contributes 0, not -inf, to the factorial terms (lnFactorial(0)=0), and
log1p(0)=0 when there is no retained intensity at all, so a PSM that retains nothing
legitimately scores exactly 0 -- not a special/undefined case, just the real output of
the formula. "Internal Initial Hyperscore" reports the same formula computed over *all*
of a PSM's matches before any claiming is applied, for transparency.

<B>Fragment usage registry.</B> Independent of per-PSM scoring, every (window_id,
fragment_id) trace touched by >=1 PSM is reported once: who officially won it (the
FragmentClaimRegistry owner), and who else attempted to use it and did not become the
owner (competitors) -- regardless of whether a competitor's match was later credited to
its own score via the peptidoform exemption. Intended for downstream use in building a
spectral library from unambiguous fragments for peptide quantification.

<B>Target/decoy type column.</B> Both output TSVs carry a "type" column ("target" or
"decoy"), read directly from psm.tsv's own "Is Decoy" column -- this tool never
re-derives decoy status itself. "Is Decoy" is a required column; psm.tsv missing it
aborts the run with an error, same as any other required column. In psm_accounting.tsv
it is the row's own PSM. In fragment_usage.tsv it is the officially winning PSM's status
(the FragmentClaimRegistry owner), not any competitor's. By default FragPipe/Philosopher's
report step (phi-report.print-decoys) strips decoys from psm.tsv even when some pass the
configured FDR threshold during filtering, so this column will show "target" for every
row unless that FragPipe setting was enabled for the search that produced psm.tsv.

<B>Multiple input files (-in).</B> MSFragger/FragPipe can search several mzML files in one
combined run, producing a single psm.tsv whose "Spectrum" field embeds the source run's
name (e.g. "0p5_msfragger-format.02345.02345.2"). Each "-in" file is matched to its own
subset of psm.tsv rows by a boundary-aware prefix comparison between the mzML's basename
(minus extension) and that embedded run stem: one is accepted as a match for the other if
it is a prefix of it and the very next character (if any) is non-alphanumeric. This makes
the match agnostic to whatever suffix a given naming convention appends when handing files
to MSFragger (e.g. a diaWeaver output "pseudo_spectra.mzML" matches a psm.tsv run stem of
"pseudo_spectra_msfragger" or "pseudo_spectra_msfragger-format" alike) while still
rejecting an unrelated file whose name happens to share only a partial prefix (e.g.
"cell1" does not match "cell10_msfragger"). A file whose run stem matches no psm.tsv rows,
or that cannot be opened as indexed mzML, is logged and skipped rather than aborting the
whole invocation. Fragment claiming is scoped independently per file -- a PSM from one
file never competes with a PSM from another, since (window_id, fragment_id) numbering is
local to each acquisition and not guaranteed unique across separate files. Every output
row across all files carries a leading "source_file" column (the exact "-in" path) so
downstream analysis can still separate or recombine files as needed.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_diaWeaverIonAccount.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_diaWeaverIonAccount.html
*/

#include <OpenMS/ANALYSIS/ID/FragmentClaimRegistry.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>

using namespace OpenMS;

// ---------------------------------------------------------------------------

class TOPPDiaWeaverIonAccount : public TOPPBase
{
public:
  TOPPDiaWeaverIonAccount() :
    TOPPBase("diaWeaverIonAccount",
             "Translates MSFragger PSMs (psm.tsv) into OpenMS sequences, resolves each to "
             "its exact diaWeaver pseudo spectrum, and greedily claims theoretical b/y "
             "fragment matches in Hyperscore order to track fragment reuse across PSMs.",
             false)
  {}

protected:

  // -------------------------------------------------------------------------
  // One PSM row read from psm.tsv, plus its translated OpenMS sequence.
  // -------------------------------------------------------------------------
  struct PsmEntry
  {
    Size row = 0;
    String bare_sequence;                // "Peptide" column
    std::map<int, double> mods;          // position -> delta mass; 0 = N-term, -1 = C-term, 1..len = residue (1-based)
    String openms_sequence;              // translated, OpenMS-compatible bracket notation
    int charge = 0;
    double rt = 0.0;                     // seconds ("Retention")
    double precursor_mz = 0.0;           // "Observed M/Z" (or "Calibrated Observed M/Z")
    double im = 0.0;                     // "Ion Mobility"; 0.0 if column absent
    double hyperscore = 0.0;             // "Hyperscore" column, MSFragger's original value, untouched by Percolator/MSBooster
    int psm_scan = -1;                   // 0-based scan number parsed from "Spectrum"
    String spectrum_field;
    String run_stem;                     // "Spectrum" field with the trailing .scan.scan.charge stripped
    bool is_decoy = false;               // "Is Decoy" column (required)

    bool isModified() const { return !mods.empty(); }
  };

  // -------------------------------------------------------------------------
  // A PsmEntry whose spectrum has been resolved and sanity-checked, with its
  // spectrum's peaks and fragment_trace_id/fragment_window_id arrays cached
  // for the claiming pass (avoids re-fetching from disk during the sorted
  // claiming loop).
  // -------------------------------------------------------------------------
  struct ResolvedPsm
  {
    const PsmEntry* psm = nullptr;
    Size spec_idx = 0;
    int native_scan = 0;                 // 1-based "scan=N"
    double resolved_rt = 0.0;
    double resolved_precursor_mz = 0.0;
    double resolved_im = 0.0;
    double precursor_intensity = 0.0;
    std::vector<double> peak_mz;         // ascending, per mzML convention
    std::vector<float> peak_intensity;
    std::vector<Int> frag_trace_id;
    std::vector<Int> frag_window_id;
  };

  // -------------------------------------------------------------------------
  // One matched-ion touch on a fragment trace, recorded regardless of
  // claiming outcome -- the raw material for fragment_usage.tsv.
  // -------------------------------------------------------------------------
  struct UsageTouch
  {
    String sequence;      // openms_sequence of the touching PSM
    int native_scan = 0;
    String ion_label;
  };

  // One retained-or-lost match, kept per PSM for the two-pass output build.
  struct MatchOutcome
  {
    String ion_label;
    FragmentClaimRegistry::TraceKey key = 0;
  };

  // -------------------------------------------------------------------------
  // Per-PSM claiming/recompute outcome, for the accounting output.
  // -------------------------------------------------------------------------
  struct PsmAccounting
  {
    const PsmEntry* psm = nullptr;
    Size spec_idx = 0;
    int native_scan = 0;
    double precursor_mz = 0.0;
    double rt = 0.0;
    double im = 0.0;
    double original_hyperscore = 0.0;
    double internal_initial_hyperscore = 0.0;
    double recomputed_hyperscore = 0.0;
    std::vector<MatchOutcome> retained;
    std::vector<MatchOutcome> lost;
  };

  // -------------------------------------------------------------------------
  // Parse "Assigned Modifications", e.g. "16C(57.0214), 3M(15.9949)" or
  // "N-term(42.0106)", into position -> delta mass.
  // Position convention: 0 = N-term, -1 = C-term, 1..len = 1-based residue
  // index into the BARE peptide sequence (verified against real psm.tsv rows:
  // "16C(...)" on an 18-residue peptide points at its 16th character).
  // Malformed tokens are logged and skipped individually rather than failing
  // the whole row, so one bad token doesn't silently drop a real modification
  // elsewhere in the same peptide.
  // -------------------------------------------------------------------------
  static std::map<int, double> parseAssignedMods_(const String& raw, Size row)
  {
    std::map<int, double> result;
    String s = raw;
    s.trim();
    if (s.empty()) return result;

    std::vector<String> parts;
    s.split(",", parts);
    for (String part : parts)
    {
      part.trim();
      if (part.empty()) continue;

      const Size op = part.find('(');
      const Size cp = part.rfind(')');
      if (op == String::npos || cp == String::npos || cp <= op)
      {
        OPENMS_LOG_WARN << "[diaWeaverIonAccount] Row " << row
                        << ": could not parse modification token '" << part << "'. Skipping token.\n";
        continue;
      }

      String prefix = part.prefix(op);
      prefix.trim();

      double delta = 0.0;
      try
      {
        String delta_str = part.substr(op + 1, cp - op - 1);
        delta_str.trim();
        delta = delta_str.toDouble();
      }
      catch (const Exception::BaseException&)
      {
        OPENMS_LOG_WARN << "[diaWeaverIonAccount] Row " << row
                        << ": could not parse modification mass in '" << part << "'. Skipping token.\n";
        continue;
      }

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
        Size digits_end = 0;
        while (digits_end < prefix.size() && isdigit(static_cast<unsigned char>(prefix[digits_end]))) ++digits_end;
        if (digits_end == 0)
        {
          OPENMS_LOG_WARN << "[diaWeaverIonAccount] Row " << row
                          << ": unrecognized modification token '" << part << "'. Skipping token.\n";
          continue;
        }
        int pos = 0;
        try { pos = prefix.prefix(digits_end).toInt(); }
        catch (const Exception::BaseException&)
        {
          OPENMS_LOG_WARN << "[diaWeaverIonAccount] Row " << row
                          << ": could not parse residue position in '" << part << "'. Skipping token.\n";
          continue;
        }
        if (pos < 1)
        {
          OPENMS_LOG_WARN << "[diaWeaverIonAccount] Row " << row
                          << ": invalid residue position in '" << part << "'. Skipping token.\n";
          continue;
        }
        result[pos] += delta;
      }
    }
    return result;
  }

  static String fmtDelta_(double d)
  {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%+.5f", d);
    return String(buf);
  }

  // -------------------------------------------------------------------------
  // Translate a bare peptide sequence + position->delta-mass map into an
  // OpenMS-compatible AASequence string using delta-mass bracket notation.
  // Intended to be fed to AASequence::fromString -- carries the trailing '.'
  // C-terminal marker OpenMS requires.
  //
  // NOTE: if a peptide has both a residue-specific mod AND a C-terminal mod
  // on its very last residue simultaneously, both bracket groups are emitted
  // back-to-back ("K[+a][+b]."). This combination was not empirically
  // verified against AASequence::fromString (unlike every other case this
  // tool relies on) -- it is rare and did not occur in the dataset this
  // tool was built against. A warning is logged so it is never silent.
  // -------------------------------------------------------------------------
  static String buildOpenMSSequence_(const String& bare_seq, const std::map<int, double>& mods, Size row)
  {
    String result;

    const auto it_n = mods.find(0);
    if (it_n != mods.end())
    {
      result += "[" + fmtDelta_(it_n->second) + "]";
    }

    const Size len = bare_seq.size();
    const auto it_c = mods.find(-1);

    for (Size i = 0; i < len; ++i)
    {
      const int pos = static_cast<int>(i) + 1;
      const bool is_last = (i + 1 == len);

      result += bare_seq[i];

      const auto it = mods.find(pos);
      const bool has_internal_mod = (it != mods.end());
      const bool has_cterm_mod = is_last && (it_c != mods.end());

      if (has_internal_mod && has_cterm_mod)
      {
        OPENMS_LOG_WARN << "[diaWeaverIonAccount] Row " << row
                        << ": last residue carries both a residue-specific and a C-terminal "
                        << "modification simultaneously -- this combination was not empirically "
                        << "verified against AASequence::fromString. Proceeding, but verify the "
                        << "translated sequence for this row.\n";
      }

      if (has_internal_mod)
      {
        result += "[" + fmtDelta_(it->second) + "]";
      }
      if (has_cterm_mod)
      {
        result += "[" + fmtDelta_(it_c->second) + "].";
      }
    }

    return result;
  }

  // -------------------------------------------------------------------------
  // Build a human-readable bracket-annotated slice of a peptide (display
  // only -- not intended to be re-parsed by AASequence::fromString, so no
  // trailing '.' marker is added).
  // -------------------------------------------------------------------------
  static String formatResidueSlice_(const String& bare_seq, const std::map<int, double>& mods,
                                     Size first_pos, Size last_pos, bool include_nterm, bool include_cterm)
  {
    String result;
    if (include_nterm)
    {
      const auto it = mods.find(0);
      if (it != mods.end()) result += "[" + fmtDelta_(it->second) + "]";
    }
    for (Size pos = first_pos; pos <= last_pos; ++pos)
    {
      result += bare_seq[pos - 1];
      const auto it = mods.find(static_cast<int>(pos));
      if (it != mods.end()) result += "[" + fmtDelta_(it->second) + "]";
    }
    if (include_cterm)
    {
      const auto it = mods.find(-1);
      if (it != mods.end()) result += "[" + fmtDelta_(it->second) + "]";
    }
    return result;
  }

  // -------------------------------------------------------------------------
  // Parse an ion label like "b3+" or "y13++" into (type, position, charge).
  // -------------------------------------------------------------------------
  static bool parseIonLabel_(const String& ion_label, char& type, int& position, int& charge)
  {
    if (ion_label.empty()) return false;
    type = ion_label[0];
    if (type != 'b' && type != 'y') return false;
    Size i = 1;
    while (i < ion_label.size() && isdigit(static_cast<unsigned char>(ion_label[i]))) ++i;
    if (i == 1) return false;
    try { position = ion_label.substr(1, i - 1).toInt(); }
    catch (const Exception::BaseException&) { return false; }
    charge = 0;
    for (Size j = i; j < ion_label.size(); ++j) if (ion_label[j] == '+') ++charge;
    return true;
  }

  // -------------------------------------------------------------------------
  // The literal chemical fragment for an ion: first N residues from the
  // N-terminus for a b-ion, last N residues from the C-terminus for a
  // y-ion, with only the modifications that actually fall within that
  // slice reapplied (N-term mod only for b, C-term mod only for y).
  // -------------------------------------------------------------------------
  static String fragmentSubsequence_(const String& bare_sequence, const std::map<int, double>& mods,
                                      const String& ion_label)
  {
    char type; int position, charge;
    if (!parseIonLabel_(ion_label, type, position, charge)) return ion_label; // fallback: raw label
    const Size len = bare_sequence.size();
    if (position < 1 || static_cast<Size>(position) >= len) return ion_label;
    if (type == 'b')
      return formatResidueSlice_(bare_sequence, mods, 1, static_cast<Size>(position), true, false);
    else
      return formatResidueSlice_(bare_sequence, mods, len - static_cast<Size>(position) + 1, len, false, true);
  }

  // -------------------------------------------------------------------------
  // Parse the 0-based scan number embedded in psm.tsv's "Spectrum" column,
  // e.g. "0p5_pseudo_spectra_msfragger.00001.00001.2" -> 1.
  // Returns -1 on any parse failure or if the two scan tokens disagree.
  // -------------------------------------------------------------------------
  // Parses "{run_stem}.{scan}.{scan}.{charge}" -> 0-based scan number, and writes the
  // run_stem (everything before the three trailing dot-segments) to out_run_stem.
  static int parseScanFromSpectrumField_(const String& spectrum_field, String& out_run_stem)
  {
    String rest = spectrum_field;
    String tail[3];
    for (int i = 2; i >= 0; --i)
    {
      const Size dot = rest.rfind('.');
      if (dot == String::npos) return -1;
      tail[i] = rest.substr(dot + 1);
      rest = rest.prefix(dot);
    }
    out_run_stem = rest;
    try
    {
      const int scan1 = tail[0].toInt();
      const int scan2 = tail[1].toInt();
      if (scan1 != scan2) return -1;
      return scan1;
    }
    catch (const Exception::BaseException&) { return -1; }
  }

  // True if `shorter` is a prefix of `longer` AND the character immediately following that
  // prefix in `longer` is not alphanumeric (end-of-string counts as a boundary too). This is
  // how an mzML's own basename ("pseudo_spectra") is matched against the run stem embedded
  // in psm.tsv's Spectrum column, which names the file MSFragger actually searched and can
  // carry an arbitrary appended suffix ("pseudo_spectra_msfragger",
  // "pseudo_spectra_msfragger-format", ...) rather than one specific hardcoded string.
  // The boundary check exists so "cell1" does not falsely prefix-match "cell10_msfragger".
  static bool isPrefixAtBoundary_(const String& shorter, const String& longer)
  {
    if (!longer.hasPrefix(shorter)) return false;
    if (longer.size() == shorter.size()) return true;
    const unsigned char next = static_cast<unsigned char>(longer[shorter.size()]);
    return !isalnum(next);
  }

  // Two run stems are considered the same run if either is a boundary-prefix of the other --
  // covers both "mzML basename is a prefix of the psm.tsv run stem" (the common case) and the
  // reverse, in case a future naming convention flips which side carries the extra suffix.
  static bool runStemsMatch_(const String& a, const String& b)
  {
    if (a == b) return true;
    return isPrefixAtBoundary_(a, b) || isPrefixAtBoundary_(b, a);
  }

  // -------------------------------------------------------------------------
  // Parse psm.tsv: required columns are Spectrum, Peptide, Charge, Retention,
  // Assigned Modifications, Hyperscore, and Observed M/Z (or Calibrated
  // Observed M/Z), and Is Decoy. Ion Mobility is read if present, optional
  // otherwise. Every row's translated sequence is validated against
  // AASequence::fromString immediately; rows that fail are logged and
  // dropped rather than silently carried forward with a broken sequence.
  // -------------------------------------------------------------------------
  std::vector<PsmEntry> parsePsmTsv_(const String& filename) const
  {
    std::vector<PsmEntry> entries;
    std::ifstream file(filename.c_str());
    if (!file.is_open())
    {
      OPENMS_LOG_ERROR << "[diaWeaverIonAccount] Cannot open psm.tsv: " << filename << "\n";
      return entries;
    }

    std::string raw_header;
    if (!std::getline(file, raw_header))
    {
      OPENMS_LOG_ERROR << "[diaWeaverIonAccount] Empty psm.tsv: " << filename << "\n";
      return entries;
    }
    std::vector<String> col_names;
    String(raw_header).split("\t", col_names);
    for (auto& c : col_names) c.trim();

    std::map<String, int> col_idx;
    for (Size i = 0; i < col_names.size(); ++i) col_idx[col_names[i]] = static_cast<int>(i);

    auto find_col = [&](const String& name) -> int
    {
      auto it = col_idx.find(name);
      return (it == col_idx.end()) ? -1 : it->second;
    };

    const int c_spectrum = find_col("Spectrum");
    const int c_peptide  = find_col("Peptide");
    const int c_charge   = find_col("Charge");
    const int c_rt       = find_col("Retention");
    const int c_assigned = find_col("Assigned Modifications");
    const int c_hyper    = find_col("Hyperscore");
    int c_mz = find_col("Observed M/Z");
    if (c_mz < 0) c_mz = find_col("Calibrated Observed M/Z");
    const int c_im = find_col("Ion Mobility");
    const int c_decoy = find_col("Is Decoy");

    std::vector<String> missing;
    if (c_spectrum < 0) missing.push_back("Spectrum");
    if (c_peptide  < 0) missing.push_back("Peptide");
    if (c_charge   < 0) missing.push_back("Charge");
    if (c_rt       < 0) missing.push_back("Retention");
    if (c_assigned < 0) missing.push_back("Assigned Modifications");
    if (c_hyper    < 0) missing.push_back("Hyperscore");
    if (c_mz       < 0) missing.push_back("Observed M/Z / Calibrated Observed M/Z");
    if (c_decoy    < 0) missing.push_back("Is Decoy");

    if (!missing.empty())
    {
      OPENMS_LOG_ERROR << "[diaWeaverIonAccount] psm.tsv is missing required column(s):";
      for (const String& m : missing) OPENMS_LOG_ERROR << " '" << m << "'";
      OPENMS_LOG_ERROR << "\n";
      return entries;
    }
    if (c_im < 0)
    {
      OPENMS_LOG_WARN << "[diaWeaverIonAccount] psm.tsv has no 'Ion Mobility' column. Proceeding without it.\n";
    }

    std::string raw_line;
    Size row = 1;
    Size n_translate_failed = 0;
    while (std::getline(file, raw_line))
    {
      ++row;
      if (raw_line.empty()) continue;
      std::vector<String> f;
      String(raw_line).split("\t", f);
      try
      {
        PsmEntry e;
        e.row = row;
        e.bare_sequence = f.at(static_cast<Size>(c_peptide));
        e.bare_sequence.trim();
        e.mods = parseAssignedMods_(f.at(static_cast<Size>(c_assigned)), row);
        e.openms_sequence = buildOpenMSSequence_(e.bare_sequence, e.mods, row);

        String charge_str = f.at(static_cast<Size>(c_charge)); charge_str.trim();
        e.charge = charge_str.toInt();

        String rt_str = f.at(static_cast<Size>(c_rt)); rt_str.trim();
        e.rt = rt_str.toDouble();

        String mz_str = f.at(static_cast<Size>(c_mz)); mz_str.trim();
        e.precursor_mz = mz_str.toDouble();

        String hyper_str = f.at(static_cast<Size>(c_hyper)); hyper_str.trim();
        e.hyperscore = hyper_str.toDouble();

        if (c_im >= 0 && static_cast<Size>(c_im) < f.size())
        {
          String im_str = f.at(static_cast<Size>(c_im)); im_str.trim();
          e.im = im_str.toDouble();
        }

        String decoy_str = f.at(static_cast<Size>(c_decoy)); decoy_str.trim().toLower();
        e.is_decoy = (decoy_str == "true");

        e.spectrum_field = f.at(static_cast<Size>(c_spectrum));
        e.spectrum_field.trim();
        e.psm_scan = parseScanFromSpectrumField_(e.spectrum_field, e.run_stem);

        try
        {
          AASequence::fromString(e.openms_sequence);
        }
        catch (const Exception::BaseException& ex)
        {
          ++n_translate_failed;
          OPENMS_LOG_WARN << "[diaWeaverIonAccount] Row " << row << ": translated sequence '"
                          << e.openms_sequence << "' (from Peptide '" << e.bare_sequence
                          << "', Assigned Modifications '" << f.at(static_cast<Size>(c_assigned))
                          << "') failed OpenMS validation: " << ex.getMessage() << ". Skipping PSM.\n";
          continue;
        }

        entries.push_back(std::move(e));
      }
      catch (const std::out_of_range&)
      {
        OPENMS_LOG_WARN << "[diaWeaverIonAccount] Skipping short row " << row << " in psm.tsv\n";
      }
      catch (const Exception::BaseException& ex)
      {
        OPENMS_LOG_WARN << "[diaWeaverIonAccount] Skipping malformed row " << row << ": " << ex.getMessage() << "\n";
      }
    }

    OPENMS_LOG_INFO << "[diaWeaverIonAccount] Parsed " << entries.size() << " PSMs from " << filename
                    << " (" << n_translate_failed << " dropped: sequence translation failed OpenMS validation)\n";
    return entries;
  }

  // -------------------------------------------------------------------------
  // Parameter registration
  // -------------------------------------------------------------------------
  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<file>...", StringList(), "One or more original diaWeaver pseudo "
                       "spectra mzML files (must retain fragment_trace_id/fragment_window_id data "
                       "arrays and 'scan=N' native IDs -- NOT the MSFragger-stripped search-input "
                       "copy). Each file is matched to its own subset of psm.tsv by comparing the "
                       "file's basename against the run stem embedded in psm.tsv's Spectrum column "
                       "(a boundary-aware prefix match, so an arbitrary MSFragger-added suffix such "
                       "as '_msfragger' or '_msfragger-format' is tolerated without needing to name "
                       "it explicitly). Claiming is scoped independently per file -- a PSM from one "
                       "file never competes with a PSM from another, since fragment_trace_id/"
                       "fragment_window_id numbering is local to each file and not guaranteed unique "
                       "across separate acquisitions.", true);
    setValidFormats_("in", {"mzML"});

    registerInputFile_("in_ids", "<file>", "", "MSFragger/FragPipe psm.tsv (already FDR-filtered).", true);
    setValidFormats_("in_ids", {"tsv"});

    registerOutputFile_("out_psm_accounting", "<file>", "", "Output TSV: one row per resolved PSM with its "
                        "original, internal-initial, and recomputed hyperscore, and its retained/lost "
                        "fragment ions.", true);
    setValidFormats_("out_psm_accounting", {"tsv"});

    registerOutputFile_("out_fragment_usage", "<file>", "", "Output TSV: one row per (window_id, fragment_id) "
                        "touched by >=1 matched ion, giving its official winning PSM and every competing PSM "
                        "that also attempted to use it. Leave empty to skip.", false);
    setValidFormats_("out_fragment_usage", {"tsv"});

    registerDoubleOption_("fragment_mz_tolerance", "<ppm>", 20.0,
                          "Fragment ion m/z tolerance for matching theoretical b/y ions to spectrum peaks (ppm).", false);
    setMinFloat_("fragment_mz_tolerance", 0.0);

    registerDoubleOption_("precursor_mz_tolerance", "<ppm>", 20.0,
                          "Precursor m/z tolerance for the PSM-to-spectrum identity sanity check (ppm). "
                          "This is an identity check, not a search window.", false);
    setMinFloat_("precursor_mz_tolerance", 0.0);

    registerDoubleOption_("rt_tolerance", "<s>", 0.01,
                          "Retention time tolerance for the PSM-to-spectrum identity sanity check (s). "
                          "The spectrum is resolved by exact scan number; RT is only cross-checked to catch "
                          "scan-mapping errors, not used as a search window.", false);
    setMinFloat_("rt_tolerance", 0.0);

    registerIntOption_("max_fragment_charge", "<charge>", 2,
                       "Maximum charge state for theoretical b/y fragment ions, independent of precursor "
                       "charge. MSFragger's own max_fragment_charge parameter (see fragger.params) is "
                       "typically a small fixed value (e.g. 2), not precursor_charge-1 -- generating "
                       "fragment ions up to the full precursor charge (the previous default here) matches "
                       "against ion charge states MSFragger never searched at all, inflating both matched-ion "
                       "counts and recomputed hyperscore for higher-charge PSMs. The effective cap applied is "
                       "min(max_fragment_charge, psm charge).", false);
    setMinInt_("max_fragment_charge", 1);

    registerStringOption_("allow_peptidoform_fragment_sharing", "<true/false>", "true",
                          "If true (default), a naked peptide and a co-eluting modified peptide of the same "
                          "backbone sequence (or two differently-modified forms of the same backbone) are "
                          "allowed to both count a shared fragment toward their own recomputed hyperscore, "
                          "since such fragments are chemically identical, not a coincidental collision. Two "
                          "identifications of the same backbone that are BOTH fully unmodified are never "
                          "exempted by this option -- they always compete normally.", false);
    setValidStrings_("allow_peptidoform_fragment_sharing", {"true", "false"});

    registerFlag_("neutral_losses",
                  "If set, include neutral loss ions (-H2O/-NH3) when generating theoretical b/y ions.");
  }

  // Exact copy of OpenMS::HyperScore's private logfactorial_ (ProSE's own scoring helper,
  // not otherwise reachable outside that class): log(x!) relative to log((base-1)!).
  static double lnFactorial_(int x, int base = 2)
  {
    base = std::max(base, 2);
    if (x < base - 1) return 0.0;
    return std::lgamma(static_cast<double>(x) + 1.0) - std::lgamma(static_cast<double>(base));
  }

  // -------------------------------------------------------------------------
  // main_
  // -------------------------------------------------------------------------
  // Basename minus a trailing .mzML/.mzml extension, for run-stem comparison against psm.tsv.
  static String mzmlBasenameStem_(const String& path)
  {
    String stem = File::basename(path);
    for (const char* ext : {".mzML", ".mzml"})
      if (stem.hasSuffix(ext)) { stem = stem.prefix(stem.size() - strlen(ext)); break; }
    return stem;
  }

  ExitCodes main_(int, const char**) override
  {
    const StringList in_files = getStringList_("in");
    const String in_ids = getStringOption_("in_ids");
    const String out_psm_accounting = getStringOption_("out_psm_accounting");
    const String out_fragment_usage = getStringOption_("out_fragment_usage");
    const double frag_ppm = getDoubleOption_("fragment_mz_tolerance");
    const double prec_ppm = getDoubleOption_("precursor_mz_tolerance");
    const double rt_tol = getDoubleOption_("rt_tolerance");
    const int max_fragment_charge = getIntOption_("max_fragment_charge");
    const bool allow_peptidoform_sharing = (getStringOption_("allow_peptidoform_fragment_sharing") == "true");
    const bool neutral_losses = getFlag_("neutral_losses");

    // ---- Step 1: parse psm.tsv (once; rows are routed to their matching mzML per file below) ----
    const std::vector<PsmEntry> all_psms = parsePsmTsv_(in_ids);
    if (all_psms.empty())
    {
      OPENMS_LOG_ERROR << "[diaWeaverIonAccount] No PSMs parsed from " << in_ids << ". Aborting.\n";
      return INCOMPATIBLE_INPUT_DATA;
    }

    // ---- Theoretical ion generator: identical for every file, built once ----
    TheoreticalSpectrumGenerator tsg;
    {
      Param p = tsg.getDefaults();
      p.setValue("add_metainfo", "true");
      p.setValue("add_b_ions", "true");
      p.setValue("add_y_ions", "true");
      p.setValue("add_a_ions", "false");
      p.setValue("add_c_ions", "false");
      p.setValue("add_x_ions", "false");
      p.setValue("add_z_ions", "false");
      p.setValue("add_losses", neutral_losses ? "true" : "false");
      p.setValue("add_term_losses", neutral_losses ? "true" : "false");
      tsg.setParameters(p);
    }

    // ---- Open output streams once; headers written once; rows appended per file below ----
    std::ofstream out_acc_stream(out_psm_accounting.c_str());
    if (!out_acc_stream.is_open())
    {
      OPENMS_LOG_ERROR << "[diaWeaverIonAccount] Cannot write output TSV: " << out_psm_accounting << "\n";
      return CANNOT_WRITE_OUTPUT_FILE;
    }
    out_acc_stream << "source_file\tnative_scan\tsequence\ttype\tcharge\tprecursor_mz\tretention_time\tion_mobility\t"
                      "original_hyperscore\tInternal Initial Hyperscore\trecomputed_hyperscore\t"
                      "Retained Fragments\tLost Fragments\n";

    std::ofstream out_usage_stream;
    if (!out_fragment_usage.empty())
    {
      out_usage_stream.open(out_fragment_usage.c_str());
      if (!out_usage_stream.is_open())
      {
        OPENMS_LOG_ERROR << "[diaWeaverIonAccount] Cannot write output TSV: " << out_fragment_usage << "\n";
        return CANNOT_WRITE_OUTPUT_FILE;
      }
      out_usage_stream << "source_file\twindow_id\tfragment_id\tmz\tintensity\tcategory\tWinning PSM\ttype\t"
                           "Winning Fragment Ion Annotation\tCompeting PSMs\tCompeting Fragment Ion Annotation\n";
    }

    Size total_resolved = 0, total_psms_seen = 0, total_claimed = 0;
    Size total_unique = 0, total_peptidoform = 0, total_shared = 0;

    // Distinct run stems present in psm.tsv, and which of them get covered by at least one
    // -in file below. Anything left uncovered after the loop is a silent-drop hazard -- those
    // psm.tsv rows would otherwise vanish from every output with no trace, so it is checked
    // and reported explicitly rather than left implicit.
    std::set<String> psm_run_stems;
    for (const PsmEntry& p : all_psms) psm_run_stems.insert(p.run_stem);
    std::set<String> matched_run_stems;
    std::vector<String> files_zero_match;
    std::vector<String> files_failed_open;

    for (const String& in : in_files)
    {
    const String file_run_stem = mzmlBasenameStem_(in);
    std::vector<PsmEntry> psms;
    for (const PsmEntry& p : all_psms)
      if (runStemsMatch_(file_run_stem, p.run_stem)) { psms.push_back(p); matched_run_stems.insert(p.run_stem); }

    OPENMS_LOG_INFO << "\n[diaWeaverIonAccount] === File: " << in << " (run stem '" << file_run_stem
                    << "') -- " << psms.size() << " / " << all_psms.size() << " psm.tsv rows matched ===\n";
    if (psms.empty())
    {
      OPENMS_LOG_WARN << "[diaWeaverIonAccount] WARNING: -in file '" << in << "' (run stem '" << file_run_stem
                       << "') matched NO psm.tsv rows. It will be skipped entirely -- if psm.tsv was meant to "
                          "contain identifications for it, check that its run stem is a boundary-prefix of (or "
                          "prefixed by) this file's basename.\n";
      files_zero_match.push_back(in);
      continue;
    }
    total_psms_seen += psms.size();

    // ---- Step 2: open mzML ----
    OPENMS_LOG_INFO << "[diaWeaverIonAccount] Opening: " << in << "\n";
    OnDiscMSExperiment on_disc;
    if (!on_disc.openFile(in))
    {
      OPENMS_LOG_ERROR << "[diaWeaverIonAccount] WARNING: failed to open '" << in << "' as indexed mzML -- "
                        << psms.size() << " matching psm.tsv row(s) will NOT be processed. Skipping file.\n";
      files_failed_open.push_back(in);
      continue;
    }
    const Size n_spec = on_disc.getNrSpectra();

    // ---- Step 3: resolve + sanity-check every PSM's spectrum ----
    // A diaWeaver/diaWeaverCounter-produced pseudo spectra mzML always writes native IDs as
    // "scan=" + (1-based index + 1) sequentially, so the true native scan (psm_scan + 1) maps
    // directly to spectrum index psm_scan. Verified per-PSM below rather than assumed blindly.
    std::vector<ResolvedPsm> resolved;
    resolved.reserve(psms.size());

    Size n_skip_scan_parse = 0, n_skip_scan_range = 0, n_skip_nativeid = 0;
    Size n_skip_no_precursor = 0, n_skip_rt = 0, n_skip_mz = 0, n_skip_no_frag_arrays = 0;

    for (const PsmEntry& psm : psms)
    {
      if (psm.psm_scan < 0)
      {
        ++n_skip_scan_parse;
        OPENMS_LOG_WARN << "[diaWeaverIonAccount] Row " << psm.row
                        << ": could not parse scan number from Spectrum '" << psm.spectrum_field << "'. Skipping.\n";
        continue;
      }

      const Size candidate_idx = static_cast<Size>(psm.psm_scan);
      const int expected_native_scan = psm.psm_scan + 1;
      if (candidate_idx >= n_spec)
      {
        ++n_skip_scan_range;
        OPENMS_LOG_WARN << "[diaWeaverIonAccount] Row " << psm.row << ": scan=" << expected_native_scan
                        << " is out of range for " << in << " (" << n_spec << " spectra). Skipping.\n";
        continue;
      }

      MSSpectrum spec = on_disc.getSpectrum(candidate_idx);
      const String expected_native_id = "scan=" + String(expected_native_scan);
      if (spec.getNativeID() != expected_native_id)
      {
        ++n_skip_nativeid;
        if (n_skip_nativeid == 1)
        {
          OPENMS_LOG_ERROR << "[diaWeaverIonAccount] Row " << psm.row << ": expected native ID '"
                           << expected_native_id << "' at spectrum index " << candidate_idx << ", found '"
                           << spec.getNativeID() << "'. Is --in the original pseudo_spectra.mzML "
                           << "(not the MSFragger-stripped search-input copy)?\n";
        }
        continue;
      }

      if (spec.getPrecursors().empty())
      {
        ++n_skip_no_precursor;
        OPENMS_LOG_WARN << "[diaWeaverIonAccount] Row " << psm.row << ": scan=" << expected_native_scan
                        << " has no precursor. Skipping.\n";
        continue;
      }
      const Precursor& prec = spec.getPrecursors()[0];

      if (std::abs(spec.getRT() - psm.rt) > rt_tol)
      {
        ++n_skip_rt;
        OPENMS_LOG_WARN << "[diaWeaverIonAccount] Row " << psm.row << ": RT mismatch (psm=" << psm.rt
                        << ", spectrum=" << spec.getRT() << "). Skipping.\n";
        continue;
      }

      const double mz_ppm_err = std::abs(prec.getMZ() - psm.precursor_mz) / psm.precursor_mz * 1e6;
      if (mz_ppm_err > prec_ppm)
      {
        ++n_skip_mz;
        OPENMS_LOG_WARN << "[diaWeaverIonAccount] Row " << psm.row << ": precursor m/z mismatch (psm="
                        << psm.precursor_mz << ", spectrum=" << prec.getMZ() << ", " << mz_ppm_err
                        << " ppm). Skipping.\n";
        continue;
      }

      const MSSpectrum::IntegerDataArray* trace_arr = nullptr;
      const MSSpectrum::IntegerDataArray* window_arr = nullptr;
      for (const auto& arr : spec.getIntegerDataArrays())
      {
        if (arr.getName() == "fragment_trace_id") trace_arr = &arr;
        else if (arr.getName() == "fragment_window_id") window_arr = &arr;
      }
      if (trace_arr == nullptr || window_arr == nullptr)
      {
        ++n_skip_no_frag_arrays;
        if (n_skip_no_frag_arrays == 1)
        {
          OPENMS_LOG_ERROR << "[diaWeaverIonAccount] Row " << psm.row << ": scan=" << expected_native_scan
                           << " is missing 'fragment_trace_id'/'fragment_window_id' data arrays.\n";
        }
        continue;
      }

      ResolvedPsm rp;
      rp.psm = &psm;
      rp.spec_idx = candidate_idx;
      rp.native_scan = expected_native_scan;
      rp.resolved_rt = spec.getRT();
      rp.resolved_precursor_mz = prec.getMZ();
      rp.resolved_im = spec.getDriftTime();
      rp.precursor_intensity = prec.getIntensity();
      rp.peak_mz.reserve(spec.size());
      rp.peak_intensity.reserve(spec.size());
      for (const auto& peak : spec) { rp.peak_mz.push_back(peak.getMZ()); rp.peak_intensity.push_back(peak.getIntensity()); }
      rp.frag_trace_id.assign(trace_arr->begin(), trace_arr->end());
      rp.frag_window_id.assign(window_arr->begin(), window_arr->end());
      resolved.push_back(std::move(rp));
    }

    OPENMS_LOG_INFO << "[diaWeaverIonAccount] Resolved " << resolved.size() << " / " << psms.size() << " PSMs to spectra.\n"
                    << "  Skipped -- scan unparseable       : " << n_skip_scan_parse << "\n"
                    << "  Skipped -- scan out of range      : " << n_skip_scan_range << "\n"
                    << "  Skipped -- native ID mismatch     : " << n_skip_nativeid << "\n"
                    << "  Skipped -- no precursor           : " << n_skip_no_precursor << "\n"
                    << "  Skipped -- RT mismatch            : " << n_skip_rt << "\n"
                    << "  Skipped -- precursor m/z mismatch : " << n_skip_mz << "\n"
                    << "  Skipped -- missing fragment arrays: " << n_skip_no_frag_arrays << "\n";

    if (resolved.empty())
    {
      OPENMS_LOG_WARN << "[diaWeaverIonAccount] No PSMs resolved to spectra for this file. Skipping.\n";
      continue;
    }
    total_resolved += resolved.size();

    // ---- Step 4: sort by original Hyperscore desc, precursor intensity desc, row asc ----
    std::stable_sort(resolved.begin(), resolved.end(),
      [](const ResolvedPsm& a, const ResolvedPsm& b)
      {
        if (a.psm->hyperscore != b.psm->hyperscore) return a.psm->hyperscore > b.psm->hyperscore;
        if (a.precursor_intensity != b.precursor_intensity) return a.precursor_intensity > b.precursor_intensity;
        return a.psm->row < b.psm->row;
      });

    // Sequence identity lookup: openms_sequence -> one PsmEntry with that sequence (any instance
    // works, since identical openms_sequence strings imply identical bare_sequence/mods).
    // Scoped per file: claiming and peptidoform comparisons never reach across files.
    std::unordered_map<String, const PsmEntry*> seq_registry;
    seq_registry.reserve(resolved.size());
    for (const ResolvedPsm& rp : resolved) seq_registry.emplace(rp.psm->openms_sequence, rp.psm);

    // ---- Step 6: greedy claiming pass (PASS 1) ----
    // Determines, for each PSM in priority order, which of its matches are retained
    // (genuinely unclaimed, or claimed but peptidoform-exempted) vs lost (claimed by a
    // genuinely different backbone). The final "who are ALL the winners of a lost
    // fragment" text cannot be resolved yet here -- a lower-priority peptidoform sibling
    // that also retains this same trace may not have been processed yet -- so lost
    // matches are only recorded as (ion_label, key) here and resolved into display text
    // in PASS 2 below, once `usage` is fully populated across every PSM.
    FragmentClaimRegistry registry;
    std::unordered_map<FragmentClaimRegistry::TraceKey, std::vector<UsageTouch>> usage;
    // A trace's own m/z and intensity, captured on first touch. First touch is always the
    // winner (whoever is processed first among everyone who ever touches a given key always
    // claims it -- there is no earlier-priority PSM left to touch it after), so this is
    // exactly the winning PSM's observed values, not an arbitrary or averaged one.
    std::unordered_map<FragmentClaimRegistry::TraceKey, std::pair<double, float>> fragment_peak_info;
    std::vector<PsmAccounting> accounting;
    accounting.reserve(resolved.size());

    for (const ResolvedPsm& rp : resolved)
    {
      const PsmEntry& psm = *rp.psm;
      AASequence aa_seq;
      try { aa_seq = AASequence::fromString(psm.openms_sequence); }
      catch (const Exception::BaseException&) { continue; } // already validated during parsing

      PeakSpectrum theo;
      // Cap fragment ion charge at max_fragment_charge (default 2, matching MSFragger's own
      // fixed cap in fragger.params), not the full precursor charge -- a fragment can't
      // physically exceed the precursor's charge, but MSFragger's search never considers
      // ion charges above its configured max regardless of precursor charge, so matching
      // higher would find ions MSFragger itself never searched for.
      const int theo_max_charge = std::min(max_fragment_charge, psm.charge);
      tsg.getSpectrum(theo, aa_seq, 1, theo_max_charge);
      const auto& ion_names = theo.getStringDataArrays().at(0);

      const String& pep_id = psm.openms_sequence;

      std::vector<FragmentClaimRegistry::TraceKey> claimable;
      PsmAccounting acc;
      acc.psm = &psm;
      acc.spec_idx = rp.spec_idx;
      acc.native_scan = rp.native_scan;
      acc.precursor_mz = rp.resolved_precursor_mz;
      acc.rt = rp.resolved_rt;
      acc.im = rp.resolved_im;
      acc.original_hyperscore = psm.hyperscore;

      int nb_ret = 0, ny_ret = 0, nb_all = 0, ny_all = 0;
      double dot_ret = 0.0, dot_all = 0.0; // combined b+y matched-intensity sum, ProSE convention

      for (Size ion_i = 0; ion_i < theo.size(); ++ion_i)
      {
        const double ion_mz = theo[ion_i].getMZ();
        const String& ion_label = ion_names[ion_i];
        const bool is_b = ion_label.hasPrefix("b");
        const bool is_y = ion_label.hasPrefix("y");
        if (!is_b && !is_y) continue;

        // Scan every peak within the tolerance window and keep the nearest one, not just
        // the first one lower_bound happens to land on -- multiple peaks can fall inside
        // a wide ppm window at high m/z, and taking the leftmost is not the same as
        // taking the closest. Picking the wrong one can silently substitute a much more
        // (or less) intense peak for the real match, distorting the recomputed hyperscore.
        // Verified empirically: real spectra in this dataset do have multiple candidates
        // per theoretical ion, sometimes differing in intensity by >40x.
        const double tol_da = ion_mz * frag_ppm * 1e-6;
        auto lo_it = std::lower_bound(rp.peak_mz.begin(), rp.peak_mz.end(), ion_mz - tol_da);
        Size scan_idx = static_cast<Size>(lo_it - rp.peak_mz.begin());
        Size peak_idx = rp.peak_mz.size(); // sentinel: no match found yet
        double best_abs_err = std::numeric_limits<double>::max();
        for (Size k = scan_idx; k < rp.peak_mz.size() && rp.peak_mz[k] <= ion_mz + tol_da; ++k)
        {
          const double abs_err = std::abs(rp.peak_mz[k] - ion_mz);
          if (abs_err < best_abs_err) { best_abs_err = abs_err; peak_idx = k; }
        }
        if (peak_idx >= rp.peak_mz.size()) continue; // no match within tolerance

        const Int window_id = rp.frag_window_id[peak_idx];
        const Int trace_id = rp.frag_trace_id[peak_idx];
        const auto key = FragmentClaimRegistry::makeKey(window_id, trace_id);
        const float inten = rp.peak_intensity[peak_idx];

        usage[key].push_back({pep_id, rp.native_scan, ion_label});
        fragment_peak_info.try_emplace(key, rp.peak_mz[peak_idx], inten);

        // Internal Initial Hyperscore: accumulate over every match, regardless of claim outcome.
        if (is_b) ++nb_all; else ++ny_all;
        dot_all += inten;

        const auto* rec = registry.getClaimRecord(key);
        bool retained;
        if (rec == nullptr)
        {
          retained = true;
          claimable.push_back(key);
        }
        else
        {
          bool exempted = false;
          if (allow_peptidoform_sharing)
          {
            const auto owner_it = seq_registry.find(rec->peptide_seq);
            if (owner_it != seq_registry.end())
            {
              const PsmEntry* owner = owner_it->second;
              const bool same_backbone = (owner->bare_sequence == psm.bare_sequence);
              const bool different_form = (owner->openms_sequence != psm.openms_sequence);
              // A peptidoform sibling is: same backbone, but a genuinely different form.
              // "different_form" alone already implies at least one of the pair carries a
              // modification -- two NAKED copies of the same backbone always have
              // identical full sequences, so there is no separate "not both naked" case
              // to check. This also means the same exact modified form re-observed at a
              // different spectrum (identical full sequence) is correctly NOT exempted --
              // it must compete like any other duplicate identification.
              exempted = same_backbone && different_form;
            }
          }
          retained = exempted;
        }

        if (retained)
        {
          if (is_b) ++nb_ret; else ++ny_ret;
          dot_ret += inten;
          acc.retained.push_back({ion_label, key});
        }
        else
        {
          acc.lost.push_back({ion_label, key});
        }
      }

      // ProSE's exact HyperScore formula (HyperScore::computeWithDetail, verified by reading
      // its source): log1p of a single combined dot product pooling b- and y-ion matched
      // intensities together, not two separate per-series logs. Confirmed empirically to
      // track MSFragger's real reported hyperscore much more closely than the from-the-paper
      // reconstruction we used before (scan 9801: ProSE-style 65.86 vs MSFragger's 62.23,
      // versus 76.41 from the separate-log version).
      const int i_min_ret = std::min(nb_ret, ny_ret), i_max_ret = std::max(nb_ret, ny_ret);
      acc.recomputed_hyperscore = std::log1p(dot_ret) + 2 * lnFactorial_(i_min_ret) + lnFactorial_(i_max_ret, i_min_ret + 1);
      const int i_min_all = std::min(nb_all, ny_all), i_max_all = std::max(nb_all, ny_all);
      acc.internal_initial_hyperscore = std::log1p(dot_all) + 2 * lnFactorial_(i_min_all) + lnFactorial_(i_max_all, i_min_all + 1);

      // Claim with the ORIGINAL score so claim priority stays fixed to the initial
      // Hyperscore ranking and never cascades from recomputed values mid-pass.
      if (!claimable.empty()) registry.tryClaim(claimable, pep_id, psm.hyperscore, rp.spec_idx);

      accounting.push_back(std::move(acc));
    }

    OPENMS_LOG_INFO << "[diaWeaverIonAccount] Claiming pass complete: " << registry.claimedCount()
                    << " fragment traces claimed across " << accounting.size() << " PSMs.\n";

    // -------------------------------------------------------------------------
    // Given a fragment key that some PSM lost, return every OTHER identity that
    // retained it: the registry's official owner, plus (if peptidoform sharing is
    // enabled) any other touching identity that shares the owner's backbone and
    // isn't both-naked with it. Excludes the querying PSM's own touch. Only
    // callable once `usage` is fully populated (PASS 2).
    // -------------------------------------------------------------------------
    auto resolveWinners = [&](FragmentClaimRegistry::TraceKey key, const String& querying_seq, int querying_scan)
      -> std::vector<String>
    {
      std::vector<String> result;
      const auto* rec = registry.getClaimRecord(key);
      if (rec == nullptr) return result;

      const auto owner_it = seq_registry.find(rec->peptide_seq);
      const PsmEntry* owner = (owner_it != seq_registry.end()) ? owner_it->second : nullptr;

      const auto uit = usage.find(key);
      if (uit == usage.end()) return result;

      std::map<std::pair<String, int>, std::vector<String>> grouped;
      for (const UsageTouch& t : uit->second) grouped[{t.sequence, t.native_scan}].push_back(t.ion_label);

      const int owner_scan = static_cast<int>(rec->spectrum_idx) + 1;

      for (const auto& [id, ion_labels] : grouped)
      {
        const String& seq = id.first;
        const int scan = id.second;
        if (seq == querying_seq && scan == querying_scan) continue; // exclude self

        const bool is_owner = (seq == rec->peptide_seq && scan == owner_scan);
        bool is_winner = is_owner;
        if (!is_winner && allow_peptidoform_sharing && owner != nullptr)
        {
          const auto cand_it = seq_registry.find(seq);
          if (cand_it != seq_registry.end())
          {
            const PsmEntry* cand = cand_it->second;
            const bool same_backbone = (cand->bare_sequence == owner->bare_sequence);
            const bool different_form = (seq != rec->peptide_seq); // same_backbone + different_form already implies not-both-naked
            is_winner = same_backbone && different_form;
          }
        }
        if (!is_winner) continue;

        const auto p_it = seq_registry.find(seq);
        const PsmEntry* p = (p_it != seq_registry.end()) ? p_it->second : nullptr;
        for (const String& ion_label : ion_labels)
        {
          // Scan number is included even though the winner may share this PSM's exact
          // sequence+ion text (point 2's same-backbone-different-scan collision) --
          // without it, "[VGA:b3+; VGA:b3+]" would look like a self-reference instead of
          // "lost to a different spectrum of the same peptide."
          const String subseq = p ? fragmentSubsequence_(p->bare_sequence, p->mods, ion_label) : ion_label;
          result.push_back(subseq + ":" + ion_label + "@scan" + String(scan));
        }
      }
      return result;
    };

    // ---- Step 7: write per-PSM accounting TSV (PASS 2: resolve Retained/Lost Fragments text) ----
    {
      for (const PsmAccounting& acc : accounting)
      {
        String retained_str;
        for (const MatchOutcome& mo : acc.retained)
        {
          const String subseq = fragmentSubsequence_(acc.psm->bare_sequence, acc.psm->mods, mo.ion_label);
          if (!retained_str.empty()) retained_str += ";";
          retained_str += subseq + ":" + mo.ion_label;
        }

        String lost_str;
        for (const MatchOutcome& mo : acc.lost)
        {
          const String my_subseq = fragmentSubsequence_(acc.psm->bare_sequence, acc.psm->mods, mo.ion_label);
          const std::vector<String> winners = resolveWinners(mo.key, acc.psm->openms_sequence, acc.native_scan);
          String winners_str;
          for (Size i = 0; i < winners.size(); ++i) { if (i) winners_str += ", "; winners_str += winners[i]; }
          if (!lost_str.empty()) lost_str += ", ";
          lost_str += "[" + my_subseq + ":" + mo.ion_label + "; " + winners_str + "]";
        }

        out_acc_stream << in << "\t" << acc.native_scan << "\t" << acc.psm->openms_sequence << "\t"
          << (acc.psm->is_decoy ? "decoy" : "target") << "\t" << acc.psm->charge << "\t"
          << acc.precursor_mz << "\t" << acc.rt << "\t" << acc.im << "\t"
          << acc.original_hyperscore << "\t" << acc.internal_initial_hyperscore << "\t" << acc.recomputed_hyperscore << "\t"
          << retained_str << "\t" << lost_str << "\n";
      }
      OPENMS_LOG_INFO << "[diaWeaverIonAccount] Wrote " << accounting.size() << " rows to: " << out_psm_accounting << "\n";
    }

    // ---- Step 8: write fragment usage registry TSV ----
    if (!out_fragment_usage.empty())
    {
      Size n_unique = 0, n_peptidoform = 0, n_shared = 0;
      for (const auto& [key, touches] : usage)
      {
        const Int window_id = static_cast<Int>(key >> 32);
        const Int fragment_id = static_cast<Int>(static_cast<uint32_t>(key));

        const auto* rec = registry.getClaimRecord(key);
        if (rec == nullptr) continue; // every touched key is claimed by construction; defensive only

        const auto info_it = fragment_peak_info.find(key);
        const double frag_mz = (info_it != fragment_peak_info.end()) ? info_it->second.first : 0.0;
        const float frag_intensity = (info_it != fragment_peak_info.end()) ? info_it->second.second : 0.0f;

        const int owner_scan = static_cast<int>(rec->spectrum_idx) + 1;
        const auto owner_it = seq_registry.find(rec->peptide_seq);
        const PsmEntry* owner = (owner_it != seq_registry.end()) ? owner_it->second : nullptr;

        std::map<std::pair<String, int>, std::vector<String>> grouped;
        for (const UsageTouch& t : touches) grouped[{t.sequence, t.native_scan}].push_back(t.ion_label);

        String winner_annotation;
        String competing_psms, competing_annotation;
        bool any_genuine_collision = false;
        bool first_competitor = true;

        for (const auto& [id, ion_labels] : grouped)
        {
          const String& seq = id.first;
          const int scan = id.second;
          const auto p_it = seq_registry.find(seq);
          const PsmEntry* p = (p_it != seq_registry.end()) ? p_it->second : nullptr;

          String annot;
          for (Size i = 0; i < ion_labels.size(); ++i)
          {
            const String subseq = p ? fragmentSubsequence_(p->bare_sequence, p->mods, ion_labels[i]) : ion_labels[i];
            if (i) annot += ",";
            annot += subseq + ":" + ion_labels[i];
          }

          const bool is_owner = (seq == rec->peptide_seq && scan == owner_scan);
          if (is_owner)
          {
            winner_annotation = annot;
            continue;
          }

          if (!first_competitor) { competing_psms += ";"; competing_annotation += ";"; }
          first_competitor = false;
          competing_psms += String(scan);
          competing_annotation += annot;

          bool exempted = false;
          if (allow_peptidoform_sharing && owner != nullptr && p != nullptr)
          {
            const bool same_backbone = (p->bare_sequence == owner->bare_sequence);
            const bool different_form = (p->openms_sequence != owner->openms_sequence); // implies not-both-naked
            exempted = same_backbone && different_form;
          }
          if (!exempted) any_genuine_collision = true;
        }

        String category;
        if (grouped.size() == 1) { category = "unique_fragment_ion"; ++n_unique; }
        else if (any_genuine_collision) { category = "shared_fragments"; ++n_shared; }
        else { category = "peptidoform_shared_fragments"; ++n_peptidoform; }

        out_usage_stream << in << "\t" << window_id << "\t" << fragment_id << "\t" << frag_mz << "\t" << frag_intensity << "\t"
          << category << "\t" << owner_scan << "\t" << ((owner != nullptr && owner->is_decoy) ? "decoy" : "target") << "\t"
          << winner_annotation << "\t" << competing_psms << "\t" << competing_annotation << "\n";
      }
      OPENMS_LOG_INFO << "[diaWeaverIonAccount] Wrote " << usage.size() << " rows to: " << out_fragment_usage << "\n"
                      << "  unique_fragment_ion             : " << n_unique << "\n"
                      << "  peptidoform_shared_fragments    : " << n_peptidoform << "\n"
                      << "  shared_fragments                : " << n_shared << "\n";

      total_unique += n_unique;
      total_peptidoform += n_peptidoform;
      total_shared += n_shared;
    }

    total_claimed += registry.claimedCount();
    } // end for (const String& in : in_files)

    // ---- Cross-check: any psm.tsv run stem not covered by ANY -in file? ----
    // This is the reverse direction of the per-file "matched NO psm.tsv rows" warning above:
    // a run stem present in psm.tsv but never matched by any -in file means those PSMs are
    // dropped from every output with no other trace, so it must be surfaced explicitly.
    std::vector<String> unmatched_stems;
    Size n_unmatched_psms = 0;
    for (const String& stem : psm_run_stems)
    {
      if (matched_run_stems.count(stem)) continue;
      unmatched_stems.push_back(stem);
      for (const PsmEntry& p : all_psms) if (p.run_stem == stem) ++n_unmatched_psms;
    }
    if (!unmatched_stems.empty())
    {
      OPENMS_LOG_WARN << "\n[diaWeaverIonAccount] WARNING: " << n_unmatched_psms << " psm.tsv row(s) belong to "
                       << unmatched_stems.size() << " run stem(s) not covered by ANY -in file -- these PSMs were "
                          "NOT processed and do not appear in any output:\n";
      for (const String& stem : unmatched_stems) OPENMS_LOG_WARN << "    '" << stem << "'\n";
      OPENMS_LOG_WARN << "  If these belong to mzML files you intended to include, pass them via -in.\n";
    }

    OPENMS_LOG_INFO << "\n[diaWeaverIonAccount] === Summary across " << in_files.size() << " file(s) ===\n"
                    << "  psm.tsv rows matched to a file  : " << total_psms_seen << " / " << all_psms.size() << "\n"
                    << "  psm.tsv run stems uncovered      : " << unmatched_stems.size() << " (" << n_unmatched_psms << " row(s))\n"
                    << "  -in files matching 0 psm.tsv rows: " << files_zero_match.size() << "\n"
                    << "  -in files that failed to open    : " << files_failed_open.size() << "\n"
                    << "  PSMs resolved to a spectrum     : " << total_resolved << "\n"
                    << "  fragment ions claimed            : " << total_claimed << "\n"
                    << "  unique_fragment_ion              : " << total_unique << "\n"
                    << "  peptidoform_shared_fragments     : " << total_peptidoform << "\n"
                    << "  shared_fragments                 : " << total_shared << "\n";

    return EXECUTION_OK;
  }
};

// ---------------------------------------------------------------------------

int main(int argc, const char** argv)
{
  TOPPDiaWeaverIonAccount tool;
  return tool.main(argc, argv);
}
