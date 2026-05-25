// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mohammed Alhigaylan $
// $Authors: Mohammed Alhigaylan $
// --------------------------------------------------------------------------

/**
@page TOPP_diaWeaverCounter diaWeaverCounter

@brief Generates ion-accounted pseudo spectra from DIA data using peptide identifications.

Runs the full diaWeaver algorithm (MS1 feature detection + MS2 mass trace extraction +
precursor-fragment correlation) and simultaneously maps identified peptides from an
external search engine (DIA-NN library-free report or MSFragger psm.tsv) onto the
detected MS2 mass traces via theoretical b/y ion matching.

Each pseudo spectrum produced by the correlation step is split into two flavours:

- <B>Orphan pseudo spectra</B> (out_orphan): fragment peaks that could not be assigned to
  any identified peptide at 1% FDR.  These represent unexplained signal.

- <B>Annotated pseudo spectra</B> (out_annotated): fragment peaks that were claimed by at
  least one identified peptide.  Each peak carries a StringDataArray named
  "peptide_ion_claims" with entries of the form
  "SEQUENCE/charge:ion[isotope];SEQUENCE2/charge:ion2[isotope]", enabling inspection
  of how a single pseudo spectrum may house multiple co-eluting peptides.

Optionally (out_full_spectra flag), all peaks (orphan + annotated) are written to a
third mzML identical in format to diaWeaver output.

Peptide-to-trace mapping tolerances are set via mz_tolerance (ppm), rt_tolerance (s),
and im_tolerance (1/K0).

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_diaWeaverCounter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_diaWeaverCounter.html
*/

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/APPLICATIONS/diaWeaver.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerIM.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include <OpenMS/FEATUREFINDER/ElutionPeakDetection.h>
#include <OpenMS/FEATUREFINDER/FeatureFindingPeptide.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ClusterMassTracesByPrecursor.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/FORMAT/MSNumpressCoder.h>
#include <OpenMS/METADATA/SourceFile.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/DATASTRUCTURES/KDTree.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <unordered_set>
#include <chrono>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace OpenMS;

// ---------------------------------------------------------------------------
// KD-tree for MS2 mass trace centroids (2D: RT x m/z)
// ---------------------------------------------------------------------------
struct MS2TraceNode
{
  double rt;
  double mz;
  double im;        ///< 0.0 when IM data are absent
  Size   trace_idx; ///< index into the per-window ms2_traces vector

  typedef double value_type;
  double operator[](size_t i) const { return i == 0 ? rt : mz; }
};
typedef KDTree::KDTree<2, MS2TraceNode> TraceKDTree;

// ---------------------------------------------------------------------------
// One peptide entry read from the identification TSV
// ---------------------------------------------------------------------------
struct PeptideEntry
{
  String sequence;
  int    charge     = 0;
  double rt         = 0.0;  ///< apex RT in seconds
  double mz         = 0.0;  ///< precursor m/z (reference only)
  double im         = 0.0;  ///< precursor IM centroid; 0.0 = not available
  String id;                ///< row identifier (e.g. "row42")

  // Open-search fields (populated only by parseMSFraggerOpenSearchTSV_)
  String              bare_sequence;
  std::map<int,double> assigned_mods;
  double              delta_mass             = 0.0;
  std::vector<int>    localization_candidates;
};

// ===========================================================================

class TOPPDiaWeaverCounter : public TOPPBase
{
public:
  TOPPDiaWeaverCounter() :
    TOPPBase("diaWeaverCounter",
             "Ion-accounted DIA pseudo spectra: splits diaWeaver output into orphan and annotated spectra.",
             false)
  {}

protected:

  // -------------------------------------------------------------------------
  // Gaussian-weighted RT aggregation (copied from diaWeaver)
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

    std::vector<MSSpectrum> to_aggregate;
    std::vector<double> weights;

    for (Size j = center_idx; j < exp.size(); ++j)
    {
      double rt_diff = exp[j].getRT() - center_rt;
      double w = std::exp(factor * rt_diff * rt_diff);
      if (w < cutoff && j != center_idx) break;
      to_aggregate.push_back(exp[j]);
      weights.push_back(w);
    }
    for (SignedSize j = static_cast<SignedSize>(center_idx) - 1; j >= 0; --j)
    {
      double rt_diff = exp[j].getRT() - center_rt;
      double w = std::exp(factor * rt_diff * rt_diff);
      if (w < cutoff) break;
      to_aggregate.push_back(exp[j]);
      weights.push_back(w);
    }

    double sum_w = 0.0;
    for (double w : weights) sum_w += w;
    for (double& w : weights) w /= sum_w;

    picker.aggregateScans(to_aggregate, weights, out);
  }

  // -------------------------------------------------------------------------
  // FeatureFinderPeptide pipeline (copied from diaWeaver)
  // -------------------------------------------------------------------------
  bool runFeatureFinderPeptide_(MSExperiment& ms_peakmap,
                                const Param& common_param,
                                Param mtd_param,
                                Param epd_param,
                                Param ffp_param,
                                FeatureMap& feat_map,
                                std::vector<MassTrace>& traces_out)
  {
    if (ms_peakmap.empty()) return true;
    ms_peakmap.sortSpectra(true);

    std::vector<MassTrace> m_traces;
    MassTraceDetection mtdet;
    mtd_param.insert("", common_param);
    mtd_param.remove("chrom_fwhm");
    mtdet.setParameters(mtd_param);
    mtdet.run(ms_peakmap, m_traces);
    if (m_traces.empty()) { OPENMS_LOG_INFO << "No mass traces detected." << std::endl; return true; }

    std::vector<MassTrace> m_traces_final;
    if (epd_param.getValue("enabled").toBool())
    {
      std::vector<MassTrace> split;
      epd_param.remove("enabled");
      epd_param.insert("", common_param);
      epd_param.remove("noise_threshold_int");
      ElutionPeakDetection epdet;
      epdet.setParameters(epd_param);
      epdet.detectPeaks(m_traces, split);
      if (epdet.getParameters().getValue("width_filtering") == "auto")
      {
        m_traces_final.clear();
        epdet.filterByPeakWidth(split, m_traces_final);
      }
      else m_traces_final = split;
    }
    else
    {
      m_traces_final = m_traces;
      for (Size i = 0; i < m_traces_final.size(); ++i) m_traces_final[i].estimateFWHM(false);
      if (ffp_param.getValue("use_smoothed_intensities").toBool())
      {
        OPENMS_LOG_WARN << "Without EPD, smoothing is not supported. Setting use_smoothed_intensities to false." << std::endl;
        ffp_param.setValue("use_smoothed_intensities", "false");
      }
    }

    ffp_param.insert("", common_param);
    ffp_param.remove("noise_threshold_int");
    ffp_param.remove("chrom_peak_snr");
    ffp_param.setValue("report_chromatograms", "false");

    std::vector<std::vector<MSChromatogram>> feat_chromatograms;
    FeatureFindingPeptide ffpep;
    ffpep.setParameters(ffp_param);
    ffpep.run(m_traces_final, feat_map, feat_chromatograms);

    auto intensity_zero = [](Feature& f) { return f.getIntensity() == 0; };
    feat_map.erase(std::remove_if(feat_map.begin(), feat_map.end(), intensity_zero), feat_map.end());

    traces_out = m_traces_final;
    OPENMS_LOG_INFO << "FFMetabo: " << m_traces_final.size() << " traces -> " << feat_map.size() << " features" << std::endl;
    return true;
  }

  // -------------------------------------------------------------------------
  // MassTraceExtractor pipeline (copied from diaWeaver)
  // -------------------------------------------------------------------------
  bool runMassTraceExtractor_(MSExperiment& ms_peakmap,
                              const Param& common_param,
                              Param mtd_param,
                              Param epd_param,
                              std::vector<MassTrace>& traces_out)
  {
    if (ms_peakmap.empty()) return true;
    ms_peakmap.sortSpectra(true);

    std::vector<MassTrace> m_traces;
    MassTraceDetection mtdet;
    mtd_param.insert("", common_param);
    mtd_param.remove("chrom_fwhm");
    mtdet.setParameters(mtd_param);
    mtdet.run(ms_peakmap, m_traces);
    if (m_traces.empty()) { OPENMS_LOG_INFO << "No mass traces detected." << std::endl; return true; }

    bool use_epd = epd_param.getValue("enabled").toBool();
    if (use_epd)
    {
      std::vector<MassTrace> split;
      epd_param.remove("enabled");
      epd_param.insert("", common_param);
      epd_param.remove("noise_threshold_int");
      ElutionPeakDetection epdet;
      epdet.setParameters(epd_param);
      epdet.detectPeaks(m_traces, split);
      if (epdet.getParameters().getValue("width_filtering") == "auto")
      { traces_out.clear(); epdet.filterByPeakWidth(split, traces_out); }
      else traces_out = std::move(split);
    }
    else traces_out = std::move(m_traces);

    traces_out.erase(
      std::remove_if(traces_out.begin(), traces_out.end(), [](const MassTrace& t){ return t.getSize() == 0; }),
      traces_out.end());

    OPENMS_LOG_INFO << "MassTraceExtractor: " << m_traces.size() << " traces -> " << traces_out.size() << " final" << std::endl;
    return true;
  }

  // -------------------------------------------------------------------------
  // TSV parsing helpers (adapted from ion-account diaWeaverCounter)
  // -------------------------------------------------------------------------
  static String runIdFromFileName_(const String& file_name_field)
  {
    String stem = File::basename(file_name_field);
    for (const char* ext : {".d", ".mzML", ".mzml", ".raw", ".wiff"})
      if (stem.hasSuffix(ext)) { stem = stem.prefix(stem.size() - strlen(ext)); break; }
    return stem;
  }

  static String runIdFromSpectrum_(const String& spectrum_field)
  {
    String stem = spectrum_field;
    for (int i = 0; i < 3; ++i)
    {
      Size dot = stem.rfind('.');
      if (dot == String::npos) break;
      stem = stem.prefix(dot);
    }
    const String suffix = "_diatracer";
    if (stem.hasSuffix(suffix)) stem = stem.prefix(stem.size() - suffix.size());
    return stem;
  }

  static std::map<int,double> parseAssignedMods_(const String& s)
  {
    std::map<int,double> result;
    if (String(s).trim().empty()) return result;
    std::vector<String> parts;
    String(s).split(",", parts);
    for (String& part : parts)
    {
      part.trim();
      Size op = part.find('('), cp = part.rfind(')');
      if (op == String::npos || cp == String::npos || cp <= op) continue;
      String prefix = part.prefix(op).trim();
      double delta = 0.0;
      try { delta = part.substr(op + 1, cp - op - 1).trim().toDouble(); } catch (...) { continue; }
      String pl = prefix; pl.toLower();
      if (pl == "n-term" || pl == "nterm") { result[0] += delta; }
      else if (pl == "c-term" || pl == "cterm") { result[-1] += delta; }
      else
      {
        Size de = prefix.size();
        while (de > 0 && isalpha((unsigned char)prefix[de - 1])) --de;
        if (de == 0) continue;
        int pos = 0;
        try { pos = prefix.prefix(de).toInt(); } catch (...) { continue; }
        if (pos >= 1) result[pos] += delta;
      }
    }
    return result;
  }

  static std::vector<int> parseAllBestPositions_(const String& s)
  {
    std::vector<int> result;
    if (String(s).trim().empty()) return result;
    std::vector<String> parts;
    String(s).split(";", parts);
    for (String& part : parts)
    {
      part.trim();
      Size ds = 0;
      while (ds < part.size() && isalpha((unsigned char)part[ds])) ++ds;
      if (ds >= part.size()) continue;
      try { result.push_back(part.substr(ds).toInt()); } catch (...) {}
    }
    return result;
  }

  static bool ionSpansPosition_(const String& ion_name, int pos, int pep_len)
  {
    if (ion_name.empty() || pos < 1 || pos > pep_len) return false;
    const char type = ion_name[0];
    Size i = 1;
    while (i < ion_name.size() && !isdigit((unsigned char)ion_name[i])) ++i;
    Size j = i;
    while (j < ion_name.size() && isdigit((unsigned char)ion_name[j])) ++j;
    if (j == i) return false;
    int num = 0;
    try { num = ion_name.substr(i, j - i).toInt(); } catch (...) { return false; }
    if (type == 'b') return pos <= num;
    if (type == 'y') return pos >= pep_len - num + 1;
    return false;
  }

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

  static String buildOpenSearchSequence_(const String& bare_seq, const std::map<int,double>& pos_mods)
  {
    String result;
    auto it_n = pos_mods.find(0);
    if (it_n != pos_mods.end() && std::abs(it_n->second) > 0.001)
    {
      const double m = it_n->second;
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
    auto it_c = pos_mods.find(-1);
    if (it_c != pos_mods.end() && std::abs(it_c->second) > 0.001)
    {
      const double m = it_c->second;
      result += (m >= 0 ? "(+" : "(") + String(m) + ")";
    }
    return result;
  }

  std::vector<PeptideEntry> parsePeptideTSV_(const String& filename, const String& run_id) const
  {
    std::vector<PeptideEntry> entries;
    std::ifstream file(filename.c_str());
    if (!file.is_open()) { OPENMS_LOG_ERROR << "Cannot open peptide TSV: " << filename << "\n"; return entries; }

    std::string raw_header;
    if (!std::getline(file, raw_header)) { OPENMS_LOG_ERROR << "Empty TSV: " << filename << "\n"; return entries; }

    std::vector<String> col_names;
    String(raw_header).split("\t", col_names);
    for (auto& c : col_names) c.trim();

    const std::map<String,String> aliases = {
      {"run","run"},{"file.name","run"},{"spectrum","run_spectrum"},
      {"modified peptide","seq"},{"peptide","seq"},{"sequence","seq"},
      {"modified.sequence","seq"},{"stripped.sequence","seq"},{"pep.strippedsequence","seq"},
      {"charge","charge"},{"precursor.charge","charge"},{"eg.precursorcharge","charge"},
      {"retention","rt_sec"},{"rt","rt_min"},{"eg.apexrt","rt_min"},{"eg.meanapexrt","rt_min"},
      {"calibrated observed m/z","mz"},{"observed m/z","mz"},{"mz","mz"},
      {"precursor.mz","mz"},{"fg.precmz","mz"},
      {"im","im"},{"ionmobility","im"},{"eg.ionmobility","im"}
    };

    std::map<String,int> field_col;
    for (Size i = 0; i < col_names.size(); ++i)
    {
      String lower = col_names[i]; lower.toLower();
      auto it = aliases.find(lower);
      if (it != aliases.end()) field_col[it->second] = static_cast<int>(i);
    }

    if (field_col.find("seq") == field_col.end() || field_col.find("charge") == field_col.end() ||
        (field_col.find("rt_sec") == field_col.end() && field_col.find("rt_min") == field_col.end()))
    {
      OPENMS_LOG_WARN << "[diaWeaverCounter] TSV missing required columns (sequence/charge/rt). Returning empty list.\n";
      return entries;
    }

    const bool has_run_diann    = field_col.count("run") > 0;
    const bool has_run_spectrum = field_col.count("run_spectrum") > 0;
    const bool has_run_col      = has_run_diann || has_run_spectrum;
    const bool filter_by_run    = !run_id.empty();
    const bool run_col_is_filename = has_run_diann &&
      col_names[field_col.at("run")].toLower() == String("file.name");
    const bool rt_in_seconds = field_col.count("rt_sec") > 0;
    const String rt_field    = rt_in_seconds ? "rt_sec" : "rt_min";

    if (filter_by_run && !has_run_col)
      OPENMS_LOG_WARN << "[diaWeaverCounter] No run column; run filter '" << run_id << "' cannot be applied.\n";
    else if (filter_by_run)
      OPENMS_LOG_INFO << "[diaWeaverCounter] Filtering TSV to run: " << run_id << "\n";

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
        if (filter_by_run && has_run_col)
        {
          String row_run;
          if (has_run_spectrum)
            row_run = runIdFromSpectrum_(fields.at(field_col.at("run_spectrum")).trim());
          else
          {
            row_run = fields.at(field_col.at("run")).trim();
            if (run_col_is_filename) row_run = runIdFromFileName_(row_run);
          }
          if (row_run != run_id) continue;
        }

        PeptideEntry e;
        e.sequence = fields.at(field_col.at("seq")).trim();

        bool has_x = false;
        for (Size ci = 0; ci < e.sequence.size(); ++ci)
          if (e.sequence[ci] == 'X' && (ci + 1 >= e.sequence.size() || e.sequence[ci+1] != '['))
            { has_x = true; break; }
        if (has_x) { OPENMS_LOG_WARN << "[diaWeaverCounter] Sequence '" << e.sequence << "' row " << row << " has ambiguous X. Skipping.\n"; continue; }

        e.charge = fields.at(field_col.at("charge")).trim().toInt();
        e.rt     = fields.at(field_col.at(rt_field)).trim().toDouble() * (rt_in_seconds ? 1.0 : 60.0);
        e.mz     = field_col.count("mz") ? fields.at(field_col.at("mz")).trim().toDouble() : 0.0;
        e.im     = field_col.count("im") ? fields.at(field_col.at("im")).trim().toDouble() : 0.0;
        e.id     = "row" + String(row);
        entries.push_back(std::move(e));
      }
      catch (const Exception::BaseException& ex)
        { OPENMS_LOG_WARN << "[diaWeaverCounter] Skipping malformed row " << row << ": " << ex.getMessage() << "\n"; }
      catch (const std::out_of_range&)
        { OPENMS_LOG_WARN << "[diaWeaverCounter] Skipping short row " << row << "\n"; }
    }

    OPENMS_LOG_INFO << "[diaWeaverCounter] Loaded " << entries.size()
                    << " peptide entries for run '" << run_id << "' from " << filename << "\n";
    return entries;
  }

  std::vector<PeptideEntry> parseMSFraggerOpenSearchTSV_(const String& filename, const String& run_id) const
  {
    std::vector<PeptideEntry> entries;
    std::ifstream file(filename.c_str());
    if (!file.is_open()) { OPENMS_LOG_ERROR << "Cannot open MSFragger psm.tsv: " << filename << "\n"; return entries; }

    std::string raw_header;
    if (!std::getline(file, raw_header)) return entries;
    std::vector<String> col_names;
    String(raw_header).split("\t", col_names);
    for (auto& c : col_names) c.trim();

    const std::map<String,String> aliases = {
      {"spectrum","run_spectrum"},{"peptide","seq_bare"},{"charge","charge"},
      {"retention","rt_sec"},{"calibrated observed m/z","mz"},{"observed m/z","mz"},
      {"assigned modifications","assigned_mods"},{"delta mass","delta_mass"},
      {"best positions","best_positions"},{"im","im"},{"ionmobility","im"}
    };
    std::map<String,int> field_col;
    for (Size i = 0; i < col_names.size(); ++i)
    {
      String lower = col_names[i]; lower.toLower();
      auto it = aliases.find(lower);
      if (it != aliases.end()) field_col[it->second] = static_cast<int>(i);
    }
    for (const String& req : {"seq_bare","charge","rt_sec","assigned_mods","delta_mass","best_positions"})
      if (!field_col.count(req))
      {
        OPENMS_LOG_WARN << "[diaWeaverCounter] MSFragger open-search psm.tsv missing column '" << req << "'. Cannot parse.\n";
        return entries;
      }

    const bool filter_by_run = !run_id.empty() && field_col.count("run_spectrum") > 0;
    if (!run_id.empty() && !filter_by_run)
      OPENMS_LOG_WARN << "[diaWeaverCounter] No Spectrum column; run filter not applied.\n";
    else if (filter_by_run)
      OPENMS_LOG_INFO << "[diaWeaverCounter] Filtering MSFragger open-search to run: " << run_id << "\n";

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
        if (filter_by_run)
        {
          if (runIdFromSpectrum_(fields.at(field_col.at("run_spectrum")).trim()) != run_id) continue;
        }
        const String bare_seq = fields.at(field_col.at("seq_bare")).trim();
        const double delta_mass = fields.at(field_col.at("delta_mass")).trim().toDouble();
        const std::map<int,double> assigned_mods = parseAssignedMods_(fields.at(field_col.at("assigned_mods")).trim());

        std::vector<int> cand_positions;
        if (std::abs(delta_mass) > 0.02)
        {
          cand_positions = parseAllBestPositions_(fields.at(field_col.at("best_positions")).trim());
          if (cand_positions.empty())
          { OPENMS_LOG_WARN << "[diaWeaverCounter] Row " << row << ": significant delta mass but no Best Positions. Skipping.\n"; continue; }
        }

        PeptideEntry e;
        e.bare_sequence           = bare_seq;
        e.assigned_mods           = assigned_mods;
        e.delta_mass              = delta_mass;
        e.localization_candidates = cand_positions;
        e.sequence = buildOpenSearchSequence_(bare_seq, assigned_mods);
        e.charge   = fields.at(field_col.at("charge")).trim().toInt();
        e.rt       = fields.at(field_col.at("rt_sec")).trim().toDouble();
        e.mz       = field_col.count("mz") ? fields.at(field_col.at("mz")).trim().toDouble() : 0.0;
        e.im       = field_col.count("im") ? fields.at(field_col.at("im")).trim().toDouble() : 0.0;
        e.id       = "row" + String(row);

        if (e.charge < 1) { OPENMS_LOG_WARN << "[diaWeaverCounter] Row " << row << ": invalid charge. Skipping.\n"; continue; }
        try { AASequence::fromString(e.sequence); }
        catch (const Exception::BaseException& ex)
        { OPENMS_LOG_WARN << "[diaWeaverCounter] Row " << row << ": sequence validation failed: " << ex.getMessage() << ". Skipping.\n"; continue; }

        entries.push_back(std::move(e));
      }
      catch (const Exception::BaseException& ex)
        { OPENMS_LOG_WARN << "[diaWeaverCounter] Skipping malformed row " << row << ": " << ex.getMessage() << "\n"; }
      catch (const std::out_of_range&)
        { OPENMS_LOG_WARN << "[diaWeaverCounter] Skipping short row " << row << "\n"; }
    }

    OPENMS_LOG_INFO << "[diaWeaverCounter] MSFragger open-search: loaded " << entries.size()
                    << " peptides from " << filename << "\n";
    return entries;
  }

  static bool isOpenSearchPSV_(const String& filename)
  {
    std::ifstream f(filename.c_str());
    std::string header;
    if (!std::getline(f, header)) return false;
    std::vector<String> cols;
    String(header).split("\t", cols);
    bool has_bp = false, has_am = false, has_dm = false;
    for (auto& c : cols)
    {
      String cl = c; cl.trim().toLower();
      if (cl == "best positions")          has_bp = true;
      if (cl == "assigned modifications")  has_am = true;
      if (cl == "delta mass")              has_dm = true;
    }
    return has_bp && has_am && has_dm;
  }

  // -------------------------------------------------------------------------
  // Ion accounting: build trace_claim_strings for one window's ms2_traces.
  //
  // For each trace index i, trace_claim_strings[i] is either "" (orphan)
  // or a semicolon-separated list of "SEQUENCE/charge:ion[isotope]" claims.
  // -------------------------------------------------------------------------
  void explainTraces_(
      const std::vector<MassTrace>& ms2_traces,
      const std::vector<PeptideEntry>& peptides,
      const std::vector<Size>& win_pep_indices,
      bool has_im,
      double mz_tol_ppm,
      double rt_tol,
      double im_tol,
      const TheoreticalSpectrumGenerator& tsg,
      std::vector<String>& trace_claim_strings)
  {
    trace_claim_strings.assign(ms2_traces.size(), String());

    if (win_pep_indices.empty() || ms2_traces.empty()) return;

    // Build KD-tree indexed by (RT, m/z) over this window's traces
    TraceKDTree kd;
    for (Size i = 0; i < ms2_traces.size(); ++i)
    {
      MS2TraceNode node;
      node.rt        = ms2_traces[i].getCentroidRT();
      node.mz        = ms2_traces[i].getCentroidMZ();
      node.im        = (has_im && ms2_traces[i].containsIMData()) ? ms2_traces[i].getCentroidIM() : 0.0;
      node.trace_idx = i;
      kd.insert(node);
    }
    kd.optimise();

    // Per-trace claims (list of (peptide_id_string, ion_label) before string formatting)
    std::vector<std::vector<std::pair<String,String>>> raw_claims(ms2_traces.size());

    auto queryAndRecord = [&](double mz_center, const String& ion_label, double pep_rt, double pep_im)
    {
      const double mz_da = mz_center * mz_tol_ppm * 1e-6;
      TraceKDTree::_Region_ region;
      region._M_low_bounds[0]  = pep_rt - rt_tol;
      region._M_high_bounds[0] = pep_rt + rt_tol;
      region._M_low_bounds[1]  = mz_center - mz_da;
      region._M_high_bounds[1] = mz_center + mz_da;

      std::vector<MS2TraceNode> candidates;
      kd.find_within_range(region, std::back_inserter(candidates));
      return candidates;
    };

    for (Size pi : win_pep_indices)
    {
      const PeptideEntry& pep = peptides[pi];
      const String pep_id = pep.sequence + "/" + String(pep.charge);

      auto recordHits = [&](const std::vector<MS2TraceNode>& cands, const String& ion_label)
      {
        for (const auto& cand : cands)
        {
          if (has_im && pep.im > 0.0 && ms2_traces[cand.trace_idx].containsIMData())
            if (std::abs(cand.im - pep.im) > im_tol) continue;
          raw_claims[cand.trace_idx].emplace_back(pep_id, ion_label);
        }
      };

      if (pep.localization_candidates.empty())
      {
        // Standard path (DIA-NN / closed-search)
        AASequence aa_seq;
        try { aa_seq = AASequence::fromString(pep.sequence); }
        catch (...) { continue; }
        if (aa_seq.size() < 2) continue;

        PeakSpectrum theo;
        tsg.getSpectrum(theo, aa_seq, 1, pep.charge);
        const auto& ion_names    = theo.getStringDataArrays().at(0);
        const auto& frag_charges = theo.getIntegerDataArrays().at(0);

        for (Size ion_i = 0; ion_i < theo.size(); ++ion_i)
        {
          const double mz_mono = theo[ion_i].getMZ();
          const String ion_name = ion_names[ion_i];
          const int    frag_z   = frag_charges[ion_i];

          recordHits(queryAndRecord(mz_mono, ion_name + "[M+0]", pep.rt, pep.im), ion_name + "[M+0]");
          for (int iso = 1; iso <= 3; ++iso)
          {
            const double mz_iso = mz_mono + iso * Constants::C13C12_MASSDIFF_U / frag_z;
            const String iso_label = ion_name + "[M+" + String(iso) + "]";
            recordHits(queryAndRecord(mz_iso, iso_label, pep.rt, pep.im), iso_label);
          }
        }
      }
      else
      {
        // Open-search path: generate ions for each candidate localisation position
        const String delta_tag = formatDeltaTag_(pep.delta_mass);
        const int    pep_len   = static_cast<int>(pep.bare_sequence.size());
        std::set<String> seen_labels;

        for (int cand_pos : pep.localization_candidates)
        {
          std::map<int,double> cand_mods = pep.assigned_mods;
          cand_mods[cand_pos] += pep.delta_mass;

          AASequence cand_seq;
          try { cand_seq = AASequence::fromString(buildOpenSearchSequence_(pep.bare_sequence, cand_mods)); }
          catch (...) { continue; }
          if (cand_seq.size() < 2) continue;

          PeakSpectrum cand_spec;
          tsg.getSpectrum(cand_spec, cand_seq, 1, pep.charge);
          const auto& cand_names   = cand_spec.getStringDataArrays().at(0);
          const auto& cand_charges = cand_spec.getIntegerDataArrays().at(0);

          for (Size ion_i = 0; ion_i < cand_spec.size(); ++ion_i)
          {
            const double mz_mono   = cand_spec[ion_i].getMZ();
            const String base_name = cand_names[ion_i];
            const int    frag_z    = cand_charges[ion_i];
            const bool   is_mod    = ionSpansPosition_(base_name, cand_pos, pep_len);
            const String ion_label = is_mod ? base_name + delta_tag : base_name;

            if (seen_labels.count(ion_label)) continue;
            seen_labels.insert(ion_label);

            recordHits(queryAndRecord(mz_mono, ion_label + "[M+0]", pep.rt, pep.im), ion_label + "[M+0]");
            for (int iso = 1; iso <= 3; ++iso)
            {
              const double mz_iso = mz_mono + iso * Constants::C13C12_MASSDIFF_U / frag_z;
              const String iso_label = ion_label + "[M+" + String(iso) + "]";
              recordHits(queryAndRecord(mz_iso, iso_label, pep.rt, pep.im), iso_label);
            }
          }
        }
      }
    }

    // Collapse raw_claims into semicolon-separated strings
    for (Size i = 0; i < ms2_traces.size(); ++i)
    {
      if (raw_claims[i].empty()) continue;
      String s;
      for (Size k = 0; k < raw_claims[i].size(); ++k)
      {
        if (k > 0) s += ";";
        s += raw_claims[i][k].first + ":" + raw_claims[i][k].second;
      }
      trace_claim_strings[i] = s;
    }
  }

  // -------------------------------------------------------------------------
  // Parameter registration
  // -------------------------------------------------------------------------
  void registerOptionsAndFlags_() override
  {
    // --- Input / output ---
    registerInputFile_("in", "<file>", "", "Input DIA mzML file", true);
    setValidFormats_("in", {"mzML"});

    registerInputFile_("in_ids", "<file>", "",
                       "Peptide identifications at 1% FDR (TSV). "
                       "Accepts DIA-NN report.tsv or MSFragger psm.tsv.", true);
    setValidFormats_("in_ids", {"tsv"});

    registerOutputFile_("out_orphan", "<file>", "",
                        "Output mzML: pseudo spectra retaining only fragment peaks "
                        "NOT explained by any identified peptide.", true);
    setValidFormats_("out_orphan", {"mzML"});

    registerOutputFile_("out_annotated", "<file>", "",
                        "Output mzML: pseudo spectra retaining only fragment peaks "
                        "explained by identified peptides.  Each peak carries a "
                        "StringDataArray 'peptide_ion_claims'.", true);
    setValidFormats_("out_annotated", {"mzML"});

    registerOutputFile_("out_full", "<file>", "",
                        "Output mzML: full pseudo spectra (all peaks, equivalent to diaWeaver output). "
                        "Leave empty to skip (default).", false);
    setValidFormats_("out_full", {"mzML"});

    // --- Ion accounting tolerances ---
    registerDoubleOption_("mz_tolerance", "<ppm>", 20.0,
                          "Fragment ion m/z tolerance for peptide-to-trace mapping (±, ppm)", false);
    setMinFloat_("mz_tolerance", 0.0);

    registerDoubleOption_("rt_tolerance", "<s>", 15.0,
                          "RT tolerance for peptide-to-trace mapping (±, seconds)", false);
    setMinFloat_("rt_tolerance", 0.0);

    registerDoubleOption_("im_tolerance", "<1/K0>", 0.02,
                          "Ion mobility tolerance for peptide-to-trace mapping (±, 1/K0; "
                          "applied only when IM data are present)", false);
    setMinFloat_("im_tolerance", 0.0);

    // --- Flags (mirrored from diaWeaver) ---
    registerFlag_("keep_ms1",
                  "If set, include peak-picked MS1 spectra in the orphan and annotated output files.");

    registerFlag_("aggregate_across_scans",
                  "If set, aggregate signal across neighboring scans using Gaussian weighting "
                  "before peak picking (requires IM data).", false);

    // --- Subsections (identical to diaWeaver) ---
    registerSubsection_("PeakPickerIM",       "Parameters for ion mobility peak picking");
    registerSubsection_("PeakPickerHiRes",    "Parameters for high-resolution peak picking");
    registerSubsection_("FeatureFinderPeptide","Parameters for FeatureFinderPeptide (MS1 precursor detection)");
    registerSubsection_("MassTraceExtractor", "Parameters for MassTraceExtractor (MS2 fragment trace detection)");
    registerSubsection_("ClusterMassTraces",  "Parameters for precursor-fragment clustering");

    registerIntOption_("threads", "<n>", 1, "Total number of threads", false);
    setMinInt_("threads", 1);

    registerIntOption_("threads_outer_loop", "<n>", -1,
                       "Number of threads for the outer (window) loop. "
                       "Remaining threads are used for inner (peak picking) loop. "
                       "Set to -1 to use all threads in the outer loop only.", false);
  }

  Param getSubsectionDefaults_(const String& name) const override
  {
    if (name == "PeakPickerIM")   { PeakPickerIM   p; return p.getDefaults(); }
    if (name == "PeakPickerHiRes"){ PeakPickerHiRes p; return p.getDefaults(); }

    if (name == "FeatureFinderPeptide")
    {
      Param combined;
      Param p_com;
      p_com.setValue("noise_threshold_int", 60.0, "Intensity threshold below which peaks are regarded as noise.");
      p_com.setValue("chrom_peak_snr", 1.0, "Minimum signal-to-noise a mass trace should have.");
      p_com.setValue("chrom_fwhm", 5.0, "Expected chromatographic peak width (in seconds).");
      p_com.setValue("auto_noise_threshold", "true", "Automatically estimate noise threshold from the input map.");
      p_com.setValidStrings("auto_noise_threshold", {"true","false"});
      p_com.setValue("noise_estimation_n_scans", 50, "Number of scans sampled to estimate noise.");
      p_com.setValue("noise_estimation_percentile", 80.0, "Intensity percentile used to define noise level.");
      combined.insert("common:", p_com);
      combined.setSectionDescription("common", "Common parameters for all subsections");

      Param p_mtd = MassTraceDetection().getDefaults();
      p_mtd.setValue("mass_error_ppm", 7.0, "Allowed mass deviation (ppm).");
      p_mtd.setValue("min_trace_length", 5.0, "Minimum expected trace length (seconds).");
      p_mtd.setValue("ion_mobility_tolerance", 0.01, "Allowed ion mobility deviation (1/k0).");
      p_mtd.setValue("reestimate_mt_sd", "false", "Enable dynamic re-estimation of m/z variance.");
      p_mtd.setValue("quant_method", "max_height", "Quantification method for mass traces.");
      p_mtd.setValue("trace_termination_outliers", 2, "Consecutive missing-peak scans before termination.");
      p_mtd.setValue("impute_zeros_missing_scans", "true", "Insert zero-intensity points at missing scans.");
      p_mtd.remove("noise_threshold_int"); p_mtd.remove("chrom_peak_snr");
      p_mtd.remove("auto_noise_threshold"); p_mtd.remove("noise_estimation_n_scans");
      p_mtd.remove("noise_estimation_percentile");
      combined.insert("mtd:", p_mtd);
      combined.setSectionDescription("mtd", "Mass Trace Detection parameters");

      Param p_epd;
      p_epd.setValue("enabled", "true", "Enable splitting of isobaric traces by elution peak detection.");
      p_epd.setValue("width_filtering", "off", "Filter unlikely peak widths.");
      p_epd.setValidStrings("enabled", {"true","false"});
      p_epd.insert("", ElutionPeakDetection().getDefaults());
      p_epd.remove("chrom_peak_snr"); p_epd.remove("chrom_fwhm");
      combined.insert("epd:", p_epd);
      combined.setSectionDescription("epd", "Elution Peak Detection parameters");

      Param p_ffp = FeatureFindingPeptide().getDefaults();
      p_ffp.setValue("local_rt_range", 5.0, "RT range for coeluting mass trace search");
      p_ffp.setValue("local_mz_range", 3.0, "MZ range for isotopic mass trace search");
      p_ffp.setValue("local_im_range", 0.02, "IM range for isotopic mass trace search");
      p_ffp.setValue("charge_lower_bound", 2, "Lowest charge state to consider");
      p_ffp.setValue("charge_upper_bound", 4, "Highest charge state to consider");
      p_ffp.setValue("remove_single_traces", "true", "Remove unassembled single traces.");
      p_ffp.setValue("use_smoothed_intensities", "true", "Use Savitzky-Golay smoothed intensities.");
      p_ffp.setValue("mass_defect_filtering", "true", "Filter hypotheses by peptide mass defect boundaries.");
      p_ffp.setValue("mass_defect_offset", 0.1, "Mass defect tolerance offset.");
      p_ffp.setValue("overlapping_features", "false", "Allow low-confidence hypotheses to reuse traces.");
      p_ffp.setValue("hypothesis_score_quantile", 0.5, "Score quantile for low-confidence classification.");
      p_ffp.setValue("rt_max_lag", 5, "Maximum lag for normalised cross-correlation.");
      p_ffp.setValue("rt_min_pearson_correlation", 0.3, "Minimum Pearson correlation between traces.");
      p_ffp.setValue("rt_peak_overlap_threshold", 0.3, "Minimum FWHM overlap proportion.");
      p_ffp.remove("chrom_fwhm"); p_ffp.remove("report_chromatograms");
      combined.insert("ffp:", p_ffp);
      combined.setSectionDescription("ffp", "FeatureFindingPeptide parameters");
      return combined;
    }

    if (name == "MassTraceExtractor")
    {
      Param combined;
      Param p_com;
      p_com.setValue("noise_threshold_int", 30.0, "Intensity threshold below which peaks are regarded as noise.");
      p_com.setValue("chrom_peak_snr", 1.0, "Minimum signal-to-noise a mass trace should have.");
      p_com.setValue("chrom_fwhm", 3.0, "Expected chromatographic peak width (in seconds).");
      p_com.setValue("auto_noise_threshold", "true", "Automatically estimate noise threshold.");
      p_com.setValidStrings("auto_noise_threshold", {"true","false"});
      p_com.setValue("noise_estimation_n_scans", 50, "Number of scans sampled to estimate noise.");
      p_com.setValue("noise_estimation_percentile", 80.0, "Intensity percentile used to define noise level.");
      combined.insert("common:", p_com);
      combined.setSectionDescription("common", "Common parameters for all subsections");

      Param p_mtd = MassTraceDetection().getDefaults();
      p_mtd.setValue("mass_error_ppm", 7.0, "Allowed mass deviation (ppm).");
      p_mtd.setValue("min_trace_length", 2.0, "Minimum expected trace length (seconds).");
      p_mtd.setValue("ion_mobility_tolerance", 0.01, "Allowed ion mobility deviation (1/k0).");
      p_mtd.setValue("reestimate_mt_sd", "false", "Enable dynamic re-estimation of m/z variance.");
      p_mtd.setValue("quant_method", "max_height", "Quantification method for mass traces.");
      p_mtd.setValue("trace_termination_outliers", 2, "Consecutive missing-peak scans before termination.");
      p_mtd.setValue("impute_zeros_missing_scans", "true", "Insert zero-intensity points at missing scans.");
      p_mtd.remove("noise_threshold_int"); p_mtd.remove("chrom_peak_snr");
      p_mtd.remove("auto_noise_threshold"); p_mtd.remove("noise_estimation_n_scans");
      p_mtd.remove("noise_estimation_percentile");
      combined.insert("mtd:", p_mtd);
      combined.setSectionDescription("mtd", "Mass Trace Detection parameters");

      Param p_epd;
      p_epd.setValue("enabled", "true", "Enable splitting of isobaric traces by elution peak detection.");
      p_epd.setValue("width_filtering", "off", "Filter unlikely peak widths.");
      p_epd.setValidStrings("enabled", {"true","false"});
      p_epd.insert("", ElutionPeakDetection().getDefaults());
      p_epd.remove("chrom_peak_snr"); p_epd.remove("chrom_fwhm");
      combined.insert("epd:", p_epd);
      combined.setSectionDescription("epd", "Elution Peak Detection parameters");
      return combined;
    }

    if (name == "ClusterMassTraces")
    {
      Param p;
      p.setValue("min_pearson_correlation", 0.3, "Minimal Pearson correlation score.");
      p.setValue("max_lag", 1, "Maximal lag (spectra shift).");
      p.setValue("min_nr_ions", 30, "Minimal number of ions per output spectrum.");
      p.setValue("max_rt_apex_difference", 5.0, "Maximal RT apex difference (seconds).");
      p.setValue("im_tolerance", 0.02, "Ion mobility tolerance for precursor-fragment matching.");
      p.setValue("nr_precursors_per_fragment", 50, "Maximum precursors a fragment can be assigned to.");
      p.setValue("rt_tolerance", 2.0, "RT tolerance for trace-point correlation matching.");
      p.setValue("pearson_weight", 1.0, "Weight for the Pearson component in combined score.");
      p.setValue("delta_rt_weight", 1.0, "Weight for the delta RT component in combined score.");
      p.setValue("delta_im_weight", 1.0, "Weight for the delta IM component in combined score.");
      p.setValue("max_nr_ions", 500, "Maximum fragment ions per output spectrum (0 = no limit).");
      p.setValue("assign_unassigned_to_all", "false", "Assign unassigned fragments to all precursors.");
      p.setValidStrings("assign_unassigned_to_all", {"false","true"});
      p.setValue("use_combined_scores", "true", "Rank assignments by combined score if true.");
      p.setValidStrings("use_combined_scores", {"false","true"});
      p.setValue("output_fragment_scores", "false", "Output per-fragment scores as FloatDataArrays.");
      p.setValidStrings("output_fragment_scores", {"false","true"});
      p.setValue("smooth_ms1", "true", "Use EPD-smoothed intensities for MS1 profiles.");
      p.setValidStrings("smooth_ms1", {"false","true"});
      p.setValue("smooth_ms2", "false", "Use EPD-smoothed intensities for MS2 profiles.");
      p.setValidStrings("smooth_ms2", {"false","true"});
      return p;
    }
    return Param();
  }

  // -------------------------------------------------------------------------
  // main_
  // -------------------------------------------------------------------------
  ExitCodes main_(int, const char**) override
  {
    const String in          = getStringOption_("in");
    const String in_ids      = getStringOption_("in_ids");
    const String out_orphan  = getStringOption_("out_orphan");
    const String out_ann     = getStringOption_("out_annotated");
    const String out_full    = getStringOption_("out_full");

    const double mz_tol      = getDoubleOption_("mz_tolerance");
    const double rt_tol      = getDoubleOption_("rt_tolerance");
    const double im_tol      = getDoubleOption_("im_tolerance");
    const bool   keep_ms1    = getFlag_("keep_ms1");
    const bool   aggregate   = getFlag_("aggregate_across_scans");

    const Param ppim_params    = getParam_().copy("PeakPickerIM:", true);
    const Param pphr_params    = getParam_().copy("PeakPickerHiRes:", true);
    const Param ffm_common     = getParam_().copy("FeatureFinderPeptide:common:", true);
    Param ffm_mtd              = getParam_().copy("FeatureFinderPeptide:mtd:", true);
    Param ffm_epd              = getParam_().copy("FeatureFinderPeptide:epd:", true);
    Param ffm_ffp              = getParam_().copy("FeatureFinderPeptide:ffp:", true);
    const Param mte_common     = getParam_().copy("MassTraceExtractor:common:", true);
    Param mte_mtd              = getParam_().copy("MassTraceExtractor:mtd:", true);
    Param mte_epd              = getParam_().copy("MassTraceExtractor:epd:", true);
    const Param cluster_param  = getParam_().copy("ClusterMassTraces:", true);

#ifdef _OPENMP
    const int num_threads        = getIntOption_("threads");
    const int threads_outer_loop = getIntOption_("threads_outer_loop");
#endif

    auto start_time = std::chrono::high_resolution_clock::now();

    // ------------------------------------------------------------------
    // Step 1: Determine DIA windows
    // ------------------------------------------------------------------
    OPENMS_LOG_INFO << "Opening: " << in << std::endl;
    OnDiscMSExperiment on_disc;
    if (!on_disc.openFile(in))
    {
      OPENMS_LOG_ERROR << "Failed to open as indexed mzML: " << in << std::endl;
      return INPUT_FILE_NOT_FOUND;
    }

    DiaWeaver::WindowMap windows;
    DiaWeaver::determineWindows(on_disc, windows);
    const DiaWeaver::IMInfo im_info = DiaWeaver::determineIMInfo(on_disc, windows);

    if (windows.empty())
    {
      OPENMS_LOG_ERROR << "No DIA windows found. Is this a DIA mzML?" << std::endl;
      return INCOMPATIBLE_INPUT_DATA;
    }
    OPENMS_LOG_INFO << "Found " << windows.size() << " DIA windows." << std::endl;

    if (im_info.available)
    {
      OPENMS_LOG_INFO << "Ion mobility data detected. Using PeakPickerIM." << std::endl;
      if (aggregate) OPENMS_LOG_INFO << "Aggregation across scans enabled." << std::endl;
    }
    else
    {
      OPENMS_LOG_INFO << "No ion mobility data detected. Using PeakPickerHiRes." << std::endl;
      if (aggregate) OPENMS_LOG_WARN << "aggregate_across_scans set but no IM data. Aggregation skipped." << std::endl;
    }

    // ------------------------------------------------------------------
    // Step 2: Parse peptide identifications
    // ------------------------------------------------------------------
    String run_id = File::basename(in);
    for (const char* ext : {".mzML",".mzml"})
      if (run_id.hasSuffix(ext)) { run_id = run_id.prefix(run_id.size() - strlen(ext)); break; }
    OPENMS_LOG_INFO << "Run ID: " << run_id << std::endl;

    const std::vector<PeptideEntry> peptides = isOpenSearchPSV_(in_ids)
      ? parseMSFraggerOpenSearchTSV_(in_ids, run_id)
      : parsePeptideTSV_(in_ids, run_id);
    OPENMS_LOG_INFO << "Peptides loaded: " << peptides.size() << std::endl;

    // Theoretical spectrum generator configured once and reused across all windows
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

    // Window vector for indexed OpenMP access
    std::vector<std::pair<DiaWeaver::DIAWindow, std::vector<Size>>> window_vec(
      windows.begin(), windows.end());
    const Size total_windows = window_vec.size();

    // ------------------------------------------------------------------
    // Step 3: Output consumers
    // ------------------------------------------------------------------
    const String output_orphan  = out_orphan;
    const String output_ann     = out_ann;
    const String output_full    = out_full;

    auto makeConsumer = [&](const String& path) -> std::unique_ptr<PlainMSDataWritingConsumer>
    {
      auto c = std::make_unique<PlainMSDataWritingConsumer>(path);
      c->setExpectedSize(0, 0);
      SourceFile sf;
      sf.setNameOfFile(File::basename(in));
      sf.setPathToFile(File::path(in));
      ExperimentalSettings es;
      es.setSourceFiles({sf});
      c->setExperimentalSettings(es);
      return c;
    };

    Size spectra_orphan_written = 0, spectra_ann_written = 0, spectra_full_written = 0;

    {
      auto consumer_orphan = makeConsumer(output_orphan);
      auto consumer_ann    = makeConsumer(output_ann);
      std::unique_ptr<PlainMSDataWritingConsumer> consumer_full;
      if (!output_full.empty()) consumer_full = makeConsumer(output_full);

      Size spectrum_index_orphan = 0, spectrum_index_ann = 0, spectrum_index_full = 0;
      Size processed = 0;

#ifdef _OPENMP
      const int total_threads = num_threads;
      int outer_threads = total_threads;
      int inner_threads = 1;
      if (threads_outer_loop > 0)
      {
        outer_threads = std::min(threads_outer_loop, total_threads);
        inner_threads = std::max(1, total_threads / outer_threads);
        omp_set_nested(1);
        omp_set_dynamic(0);
        OPENMS_LOG_INFO << "Nested parallelism: " << outer_threads << " outer x "
                        << inner_threads << " inner threads." << std::endl;
      }
      else
        OPENMS_LOG_INFO << "Using " << outer_threads << " threads (no nesting)." << std::endl;
      omp_set_num_threads(outer_threads);
#pragma omp parallel for schedule(dynamic, 1) firstprivate(on_disc)
#endif
      for (SignedSize idx = 0; idx < static_cast<SignedSize>(total_windows); ++idx)
      {
        const DiaWeaver::DIAWindow& w       = window_vec[idx].first;
        const std::vector<Size>&    indices = window_vec[idx].second;

        // ---- 3a: Extract + peak-pick MS2 ----
        MSExperiment ms2_exp;
        DiaWeaver::extractSingleMS2Window(on_disc, w, indices, im_info, ms2_exp);

        if (aggregate && im_info.available && !ms2_exp.empty())
        {
          MSExperiment picked;
          picked.resize(ms2_exp.size());
#pragma omp parallel num_threads(inner_threads)
          {
            PeakPickerIM pp; pp.setParameters(ppim_params);
#pragma omp for schedule(dynamic,1)
            for (SignedSize s = 0; s < static_cast<SignedSize>(ms2_exp.size()); ++s)
            {
              if (ms2_exp[s].getIMPeakType() != IMPeakType::IM_CENTROIDED)
              {
                MSSpectrum agg;
                aggregateSpectrum_(ms2_exp, static_cast<Size>(s), pp, agg);
                pp.pickIMTraces(agg);
                picked[s] = std::move(agg);
              }
              else picked[s] = ms2_exp[s];
            }
          }
          ms2_exp = std::move(picked);
        }
        else if (!ms2_exp.empty())
        {
#pragma omp parallel num_threads(inner_threads)
          {
            PeakPickerIM pp_im; PeakPickerHiRes pp_hr;
            if (im_info.available) pp_im.setParameters(ppim_params);
            else pp_hr.setParameters(pphr_params);
#pragma omp for schedule(dynamic,1)
            for (SignedSize s = 0; s < static_cast<SignedSize>(ms2_exp.size()); ++s)
            {
              if (im_info.available)
              {
                if (ms2_exp[s].getIMPeakType() != IMPeakType::IM_CENTROIDED)
                  pp_im.pickIMTraces(ms2_exp[s]);
              }
              else if (ms2_exp[s].getType(true) != SpectrumSettings::SpectrumType::CENTROID)
              {
                MSSpectrum picked; pp_hr.pick(ms2_exp[s], picked); ms2_exp[s] = std::move(picked);
              }
            }
          }
        }

        // ---- 3b: MS2 mass trace extraction ----
        std::vector<MassTrace> ms2_traces;
        if (!ms2_exp.empty())
        {
          Param mtd_c = mte_mtd, epd_c = mte_epd;
          runMassTraceExtractor_(ms2_exp, mte_common, mtd_c, epd_c, ms2_traces);
        }

        // ---- 3c: Ion accounting — build trace_claim_strings ----
        std::vector<String> trace_claim_strings;
        if (!ms2_traces.empty() && !peptides.empty())
        {
          // Filter peptides to this window by precursor m/z
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
          explainTraces_(ms2_traces, peptides, win_pep_indices,
                         im_info.available, mz_tol, rt_tol, im_tol,
                         tsg, trace_claim_strings);
        }
        else
          trace_claim_strings.assign(ms2_traces.size(), String());

        // ---- 3d: Extract + peak-pick MS1 ----
        MSExperiment ms1_exp;
        DiaWeaver::extractSingleMS1Window(on_disc, w, im_info, ms1_exp);

        if (aggregate && im_info.available && !ms1_exp.empty())
        {
          MSExperiment picked;
          picked.resize(ms1_exp.size());
#pragma omp parallel num_threads(inner_threads)
          {
            PeakPickerIM pp; pp.setParameters(ppim_params);
#pragma omp for schedule(dynamic,1)
            for (SignedSize s = 0; s < static_cast<SignedSize>(ms1_exp.size()); ++s)
            {
              if (ms1_exp[s].getIMPeakType() != IMPeakType::IM_CENTROIDED)
              {
                MSSpectrum agg;
                aggregateSpectrum_(ms1_exp, static_cast<Size>(s), pp, agg);
                pp.pickIMTraces(agg);
                picked[s] = std::move(agg);
              }
              else picked[s] = ms1_exp[s];
            }
          }
          ms1_exp = std::move(picked);
        }
        else if (!ms1_exp.empty())
        {
#pragma omp parallel num_threads(inner_threads)
          {
            PeakPickerIM pp_im; PeakPickerHiRes pp_hr;
            if (im_info.available) pp_im.setParameters(ppim_params);
            else pp_hr.setParameters(pphr_params);
#pragma omp for schedule(dynamic,1)
            for (SignedSize s = 0; s < static_cast<SignedSize>(ms1_exp.size()); ++s)
            {
              if (im_info.available)
              {
                if (ms1_exp[s].getIMPeakType() != IMPeakType::IM_CENTROIDED)
                  pp_im.pickIMTraces(ms1_exp[s]);
              }
              else if (ms1_exp[s].getType(true) != SpectrumSettings::SpectrumType::CENTROID)
              {
                MSSpectrum picked; pp_hr.pick(ms1_exp[s], picked); ms1_exp[s] = std::move(picked);
              }
            }
          }
        }

        // ---- 3e: FeatureFinderPeptide on MS1 ----
        if (!ms1_exp.empty() && !ms2_traces.empty())
        {
          FeatureMap ms1_features;
          std::vector<MassTrace> ms1_traces;
          Param mtd_c = ffm_mtd, epd_c = ffm_epd, ffp_c = ffm_ffp;

          if (runFeatureFinderPeptide_(ms1_exp, ffm_common, mtd_c, epd_c, ffp_c, ms1_features, ms1_traces)
              && !ms1_features.empty())
          {
            MSExperiment win_orphan, win_ann, win_full;
            MSExperiment* p_full = consumer_full ? &win_full : nullptr;

            ClusterMassTracesByPrecursor clusterer;
            clusterer.setParameters(cluster_param);
            clusterer.run(ms1_features, ms1_traces, ms2_traces,
                          trace_claim_strings,
                          w.lower_mz, w.upper_mz,
                          win_orphan, win_ann, p_full);

            // Write MS1 if requested (same raw spectra go to both outputs)
            if (keep_ms1)
            {
#pragma omp critical (write_spectra)
              {
                for (auto& sp : ms1_exp)
                {
                  sp.setNativeID("scan=" + String(++spectrum_index_orphan));
                  sp.setType(SpectrumSettings::SpectrumType::CENTROID);
                  consumer_orphan->consumeSpectrum(sp);
                  sp.setNativeID("scan=" + String(++spectrum_index_ann));
                  consumer_ann->consumeSpectrum(sp);
                  if (consumer_full)
                  {
                    sp.setNativeID("scan=" + String(++spectrum_index_full));
                    consumer_full->consumeSpectrum(sp);
                  }
                }
              }
            }

#pragma omp critical (write_spectra)
            {
              for (auto& sp : win_orphan)
              {
                sp.setNativeID("scan=" + String(++spectrum_index_orphan));
                consumer_orphan->consumeSpectrum(sp);
              }
              for (auto& sp : win_ann)
              {
                sp.setNativeID("scan=" + String(++spectrum_index_ann));
                consumer_ann->consumeSpectrum(sp);
              }
              if (consumer_full)
                for (auto& sp : win_full)
                {
                  sp.setNativeID("scan=" + String(++spectrum_index_full));
                  consumer_full->consumeSpectrum(sp);
                }
            }
          }
        }

#ifdef _OPENMP
#pragma omp critical (progress_log)
#endif
        {
          ++processed;
          OPENMS_LOG_INFO << "Processed window " << processed << "/" << total_windows
                          << " (m/z: " << w.lower_mz << "-" << w.upper_mz << ")" << std::endl;
        }
      } // end window loop

#ifdef _OPENMP
      if (threads_outer_loop > 0) omp_set_num_threads(num_threads);
#endif

      spectra_orphan_written = consumer_orphan->getNrSpectraWritten();
      spectra_ann_written    = consumer_ann->getNrSpectraWritten();
      if (consumer_full) spectra_full_written = consumer_full->getNrSpectraWritten();
    } // consumers finalized here

    OPENMS_LOG_INFO << "Written: " << spectra_orphan_written << " orphan, "
                    << spectra_ann_written << " annotated"
                    << (spectra_full_written ? (", " + String(spectra_full_written) + " full") : String(""))
                    << " pseudo spectra." << std::endl;

    // ------------------------------------------------------------------
    // Step 4: Sort, re-ID, compress and rewrite each output
    // ------------------------------------------------------------------
    auto sortAndRewrite = [&](const String& path, Size n_written)
    {
      if (n_written == 0)
      {
        OPENMS_LOG_WARN << "No spectra written to " << path << "." << std::endl;
        return;
      }
      OPENMS_LOG_INFO << "Sorting " << path << " by RT..." << std::endl;
      MSExperiment exp;
      MzMLFile mzml;
      mzml.load(path, exp);
      exp.sortSpectra(false);
      for (Size i = 0; i < exp.size(); ++i)
      {
        exp[i].setNativeID("scan=" + String(i + 1));
        exp[i].getInstrumentSettings().setScanMode(InstrumentSettings::ScanMode::MSNSPECTRUM);
        exp[i].getInstrumentSettings().setPolarity(IonSource::Polarity::POSITIVE);
      }

      MSNumpressCoder::NumpressConfig np_mz, np_int, np_fda;
      np_mz.setCompression("linear"); np_mz.numpressErrorTolerance = 0.0001;
      np_int.setCompression("slof");
      np_fda.setCompression("slof");
      mzml.getOptions().setNumpressConfigurationMassTime(np_mz);
      mzml.getOptions().setNumpressConfigurationIntensity(np_int);
      mzml.getOptions().setNumpressConfigurationFloatDataArray(np_fda);
      mzml.getOptions().setCompression(true);
      mzml.store(path, exp);
    };

    sortAndRewrite(output_orphan, spectra_orphan_written);
    sortAndRewrite(output_ann,    spectra_ann_written);
    if (!output_full.empty()) sortAndRewrite(output_full, spectra_full_written);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto dur = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    auto h = std::chrono::duration_cast<std::chrono::hours>(dur);
    auto m = std::chrono::duration_cast<std::chrono::minutes>(dur % std::chrono::hours(1));
    auto s = dur % std::chrono::minutes(1);
    OPENMS_LOG_INFO << "Done. Total time: " << h.count() << "h " << m.count() << "m " << s.count() << "s" << std::endl;

    return EXECUTION_OK;
  }
};

// ---------------------------------------------------------------------------

int main(int argc, const char** argv)
{
  TOPPDiaWeaverCounter tool;
  return tool.main(argc, argv);
}
