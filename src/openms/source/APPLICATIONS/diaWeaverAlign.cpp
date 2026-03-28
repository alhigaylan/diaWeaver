// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Mohammed Alhigaylan $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/diaWeaverAlign.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <algorithm>
#include <cmath>
#include <numeric>

namespace OpenMS
{

  // =========================================================================
  // Construction / parameter handling
  // =========================================================================

  DiaWeaverAlign::DiaWeaverAlign() :
    DefaultParamHandler("DiaWeaverAlign")
  {
    defaults_.setValue(
      "fragment_mz_tolerance_ppm", 20.0,
      "Fragment m/z tolerance in ppm. Defines the width of each bin on the "
      "logarithmic m/z grid: every bin spans exactly this many ppm, so bin "
      "width in Daltons grows proportionally with m/z. The same tolerance is "
      "used for exact ppm verification during searching.");
    defaults_.setMinFloat("fragment_mz_tolerance_ppm", 0.1);

    defaults_.setValue(
      "max_peaks_per_spectrum", 200,
      "Maximum number of peaks (ranked by descending intensity) to index per "
      "spectrum. Limits index size and suppresses noise peaks. ");
    defaults_.setMinInt("max_peaks_per_spectrum", 1);

    defaults_.setValue(
      "min_matched_peaks", 4,
      "Minimum number of matched fragment peaks required to report a result.");
    defaults_.setMinInt("min_matched_peaks", 1);

    defaultsToParam_();
  }

  void DiaWeaverAlign::updateMembers_()
  {
    fragment_tol_ppm_       = (double)param_.getValue("fragment_mz_tolerance_ppm");
    max_peaks_per_spectrum_ = (Size)(int)param_.getValue("max_peaks_per_spectrum");
    min_matched_peaks_      = (int)param_.getValue("min_matched_peaks");
    // log_step_, min_observed_mz_, and n_bins_ are derived from data in buildIndex.
  }


  // =========================================================================
  // Private helpers — log-space bin grid
  // =========================================================================

  uint32_t DiaWeaverAlign::toBinIdx_(double mz) const
  {
    // Position in log space relative to the grid origin (min_observed_mz_),
    // divided by the bin width (log_step_), gives the bin index.
    return static_cast<uint32_t>(std::log(mz / min_observed_mz_) / log_step_);
  }

  uint16_t DiaWeaverAlign::toMassOffset_(double mz, uint32_t bin_idx) const
  {
    // Fractional position within the bin in log space.
    // The bin occupies [bin_idx * log_step, (bin_idx + 1) * log_step] in
    // log(mz / min_observed_mz_) space.
    const double log_pos = std::log(mz / min_observed_mz_) / log_step_;
    const double frac    = log_pos - static_cast<double>(bin_idx);
    return static_cast<uint16_t>(std::max(0.0, std::min(1.0, frac)) * 65535.0);
  }

  double DiaWeaverAlign::fromBinIdx_(uint32_t bin_idx, uint16_t mass_offset) const
  {
    // Inverse of toBinIdx_ + toMassOffset_: recover approximate m/z from the
    // bin index and fractional log-space offset.
    const double log_pos = static_cast<double>(bin_idx)
                           + static_cast<double>(mass_offset) / 65535.0;
    return min_observed_mz_ * std::exp(log_pos * log_step_);
  }


  // =========================================================================
  // Index building
  // =========================================================================

  void DiaWeaverAlign::buildIndex(const std::vector<String>& mzml_files)
  {
    updateMembers_();

    source_files_.clear();
    spectrum_entries_.clear();
    fragment_entries_.clear();
    bin_offsets_.clear();
    log_step_         = 0.0;
    min_observed_mz_  = 0.0;
    n_bins_           = 0;

    OPENMS_LOG_INFO << "[DiaWeaverAlign] Building fragment index from "
                    << mzml_files.size() << " file(s)." << std::endl;

    // -----------------------------------------------------------------------
    // Phase 1: Load all MS2 spectra.
    //
    // All spectra are collected before sorting so that the final spectrum IDs
    // reflect a globally consistent precursor m/z order across all files,
    // -----------------------------------------------------------------------

    struct RawSpecEntry
    {
      MSSpectrum spectrum;
      double     precursor_mz{0.0};
      double     precursor_im{-1.0};
      int        charge{0};
      String     native_id;
      Size       source_idx{0};
    };

    std::vector<RawSpecEntry> raw_spectra;

    for (Size file_idx = 0; file_idx < mzml_files.size(); ++file_idx)
    {
      source_files_.push_back(mzml_files[file_idx]);

      MSExperiment exp;
      MzMLFile mzml;
      mzml.load(mzml_files[file_idx], exp);

      OPENMS_LOG_INFO << "[DiaWeaverAlign]   " << mzml_files[file_idx]
                      << " — " << exp.size() << " spectra" << std::endl;

      for (auto& spec : exp)
      {
        if (spec.getMSLevel() != 2) continue;
        if (spec.empty()) continue;

        double pmz     = 0.0;
        double pim     = -1.0;
        int    pcharge = 0;

        if (!spec.getPrecursors().empty())
        {
          const auto& prec = spec.getPrecursors()[0];
          pmz     = prec.getMZ();
          pcharge = prec.getCharge();
        }

        // Extract precursor ion mobility from float data arrays if present.
        for (const auto& fda : spec.getFloatDataArrays())
        {
          const String& name = fda.getName();
          if ((name == "mean inverse reduced ion mobility array" ||
               name == "mean drift time array") && !fda.empty())
          {
            pim = static_cast<double>(fda[0]);
            break;
          }
        }

        RawSpecEntry entry;
        entry.precursor_mz  = pmz;
        entry.precursor_im  = pim;
        entry.charge        = pcharge;
        entry.native_id     = spec.getNativeID();
        entry.source_idx    = file_idx;
        entry.spectrum      = std::move(spec);
        raw_spectra.push_back(std::move(entry));
      }
    }

    OPENMS_LOG_INFO << "[DiaWeaverAlign] Loaded " << raw_spectra.size()
                    << " MS2 spectra total." << std::endl;

    if (raw_spectra.empty()) return;

    // -----------------------------------------------------------------------
    // Phase 2: Determine the observed m/z range from all peaks.
    //
    // The log-space bin grid is anchored at min_observed_mz_ and extends to
    // max_observed_mz. No user-supplied range parameter is needed.
    //
    //   log_step  = log(1 + ppm / 1e6)       [bin width in log space]
    //   n_bins    = ceil(log(max / min) / log_step) + 1
    // -----------------------------------------------------------------------

    double max_observed_mz = 0.0;
    min_observed_mz_       = std::numeric_limits<double>::max();

    for (const auto& entry : raw_spectra)
    {
      for (const auto& peak : entry.spectrum)
      {
        const double mz = peak.getMZ();
        if (mz > max_observed_mz) max_observed_mz = mz;
        if (mz < min_observed_mz_) min_observed_mz_ = mz;
      }
    }

    log_step_ = std::log(1.0 + fragment_tol_ppm_ / 1e6);
    n_bins_   = static_cast<uint32_t>(
                  std::ceil(std::log(max_observed_mz / min_observed_mz_) / log_step_)) + 1;

    OPENMS_LOG_INFO << "[DiaWeaverAlign] Observed m/z range: ["
                    << min_observed_mz_ << ", " << max_observed_mz << "] Da" << std::endl;
    OPENMS_LOG_INFO << "[DiaWeaverAlign] Log-space bins: " << n_bins_
                    << "  (each = " << fragment_tol_ppm_ << " ppm)" << std::endl;

    // -----------------------------------------------------------------------
    // Phase 3: Sort spectra by precursor m/z and assign contiguous IDs.
    //
    // entries within each fragment bin will be in
    // spectrum_id order (= precursor m/z order), enabling binary-search
    // restriction to the precursor window during searching.
    // -----------------------------------------------------------------------

    std::sort(raw_spectra.begin(), raw_spectra.end(),
      [](const RawSpecEntry& a, const RawSpecEntry& b) {
        return a.precursor_mz < b.precursor_mz;
      });

    spectrum_entries_.resize(raw_spectra.size());
    for (Size i = 0; i < raw_spectra.size(); ++i)
    {
      spectrum_entries_[i].precursor_mz    = raw_spectra[i].precursor_mz;
      spectrum_entries_[i].precursor_im    = raw_spectra[i].precursor_im;
      spectrum_entries_[i].charge          = raw_spectra[i].charge;
      spectrum_entries_[i].native_id       = raw_spectra[i].native_id;
      spectrum_entries_[i].source_file_idx = raw_spectra[i].source_idx;
    }

    // -----------------------------------------------------------------------
    // Phase 4: Generate fragment entries.
    //
    // For each spectrum (in ID = precursor m/z order), rank peaks by descending
    // intensity, bin each peak on the log-space grid, and record a FragmentEntry.
    // -----------------------------------------------------------------------

    std::vector<std::pair<uint32_t, FragmentEntry>> all_entries;
    all_entries.reserve(raw_spectra.size() *
                        std::min(static_cast<Size>(50), max_peaks_per_spectrum_));

    for (uint32_t spec_id = 0;
         spec_id < static_cast<uint32_t>(raw_spectra.size());
         ++spec_id)
    {
      const MSSpectrum& spec    = raw_spectra[spec_id].spectrum;
      const int         pcharge = raw_spectra[spec_id].charge;
      const uint8_t charge_byte = (pcharge > 0 && pcharge <= 255)
                                  ? static_cast<uint8_t>(pcharge) : 0;

      std::vector<Size> order(spec.size());
      std::iota(order.begin(), order.end(), 0);
      std::sort(order.begin(), order.end(), [&spec](Size a, Size b) {
        return spec[a].getIntensity() > spec[b].getIntensity();
      });

      const Size n_peaks = std::min(spec.size(), max_peaks_per_spectrum_);

      for (Size rank = 0; rank < n_peaks; ++rank)
      {
        const double mz = spec[order[rank]].getMZ();
        if (mz <= 0.0) continue;

        const uint32_t bin_idx = toBinIdx_(mz);
        if (bin_idx >= n_bins_) continue;

        FragmentEntry fe;
        fe.spectrum_id    = spec_id;
        fe.mass_offset    = toMassOffset_(mz, bin_idx);
        fe.intensity_rank = (rank <= 254) ? static_cast<uint8_t>(rank) : 255;
        fe.charge         = charge_byte;

        all_entries.emplace_back(bin_idx, fe);
      }
    }

    OPENMS_LOG_INFO << "[DiaWeaverAlign] Generated " << all_entries.size()
                    << " fragment entries." << std::endl;

    // -----------------------------------------------------------------------
    // Phase 5: Sort by (bin_idx, spectrum_id) and build CSR structure.
    // -----------------------------------------------------------------------

    std::sort(all_entries.begin(), all_entries.end(),
      [](const std::pair<uint32_t, FragmentEntry>& a,
         const std::pair<uint32_t, FragmentEntry>& b) {
        if (a.first != b.first) return a.first < b.first;
        return a.second.spectrum_id < b.second.spectrum_id;
      });

    // Prefix-sum over per-bin counts to build bin_offsets_.
    bin_offsets_.assign(n_bins_ + 1, 0);
    for (const auto& [bin_idx, _] : all_entries)
    {
      ++bin_offsets_[bin_idx + 1];
    }
    for (Size b = 1; b <= n_bins_; ++b)
    {
      bin_offsets_[b] += bin_offsets_[b - 1];
    }

    fragment_entries_.resize(all_entries.size());
    for (Size i = 0; i < all_entries.size(); ++i)
    {
      fragment_entries_[i] = all_entries[i].second;
    }

    OPENMS_LOG_INFO << "[DiaWeaverAlign] Index built: "
                    << spectrum_entries_.size() << " spectra, "
                    << fragment_entries_.size() << " fragment entries, "
                    << n_bins_ << " bins." << std::endl;
  }


  // =========================================================================
  // Searching
  // =========================================================================

  std::vector<DiaWeaverAlign::SearchResult> DiaWeaverAlign::search(
    const MSSpectrum& query,
    double precursor_mz,
    double precursor_mz_tol_ppm,
    double query_im,
    double im_tolerance) const
  {
    if (fragment_entries_.empty() || spectrum_entries_.empty())
    {
      return {};
    }

    // -----------------------------------------------------------------------
    // Step 1: Restrict candidates to the precursor m/z window.
    // -----------------------------------------------------------------------

    const double prec_tol_da = (precursor_mz_tol_ppm / 1e6) * precursor_mz;
    const double prec_lo = precursor_mz - prec_tol_da;
    const double prec_hi = precursor_mz + prec_tol_da;

    const auto spec_lo_it = std::lower_bound(
      spectrum_entries_.begin(), spectrum_entries_.end(), prec_lo,
      [](const SpectrumEntry& e, double v) { return e.precursor_mz < v; });
    const auto spec_hi_it = std::upper_bound(
      spectrum_entries_.begin(), spectrum_entries_.end(), prec_hi,
      [](double v, const SpectrumEntry& e) { return v < e.precursor_mz; });

    if (spec_lo_it == spec_hi_it) return {};

    const uint32_t id_lo = static_cast<uint32_t>(spec_lo_it - spectrum_entries_.begin());
    const uint32_t id_hi = static_cast<uint32_t>(spec_hi_it - spectrum_entries_.begin());

    // -----------------------------------------------------------------------
    // Step 2: Optionally restrict candidates by ion mobility.
    // -----------------------------------------------------------------------

    const uint32_t n_candidates = id_hi - id_lo;
    const bool use_im = (im_tolerance >= 0.0 && query_im >= 0.0);

    std::vector<bool> im_mask(n_candidates, true);
    if (use_im)
    {
      for (uint32_t i = 0; i < n_candidates; ++i)
      {
        const double lib_im = spectrum_entries_[id_lo + i].precursor_im;
        if (lib_im < 0.0 || std::fabs(lib_im - query_im) > im_tolerance)
        {
          im_mask[i] = false;
        }
      }
    }

    // -----------------------------------------------------------------------
    // Step 3: Match query peaks against the fragment index.
    //
    // Each library fragment was assigned to exactly one bin. A query peak looks
    // up its own bin, plus both immediate neighbours. The neighbour check is
    // needed because the query may sit near a bin boundary: a library peak
    // within tolerance on the other side of that boundary lives in the adjacent
    // bin. After the bin lookup, the exact ppm distance is recomputed from the
    // stored mass_offset and only peaks within tol_da are accepted.
    // -----------------------------------------------------------------------

    struct MatchAccum
    {
      uint32_t n_matched{0};
      double   sum_cross{0.0};
      double   sum_lib_sq{0.0};
    };

    std::vector<MatchAccum> accum(n_candidates);

    double query_sq = 0.0;
    for (Size qi = 0; qi < query.size(); ++qi)
    {
      const double qint = query[qi].getIntensity();
      query_sq += qint * qint;
    }

    for (Size qi = 0; qi < query.size(); ++qi)
    {
      const double qmz  = query[qi].getMZ();
      const double qint = query[qi].getIntensity();

      if (qmz <= 0.0 || qmz < min_observed_mz_) continue;

      // Tolerance in Da at this m/z, for exact ppm verification.
      const double tol_da = (fragment_tol_ppm_ / 1e6) * qmz;

      const uint32_t bin_center = toBinIdx_(qmz);
      const uint32_t bin_lo     = (bin_center > 0) ? bin_center - 1 : 0;
      const uint32_t bin_hi     = std::min(bin_center + 1, n_bins_ - 1);

      for (uint32_t b = bin_lo; b <= bin_hi; ++b)
      {
        const uint32_t fe_start = bin_offsets_[b];
        const uint32_t fe_end   = bin_offsets_[b + 1];
        if (fe_start == fe_end) continue;

        const FragmentEntry* b_begin = fragment_entries_.data() + fe_start;
        const FragmentEntry* b_end   = fragment_entries_.data() + fe_end;

        const FragmentEntry* range_lo = std::lower_bound(
          b_begin, b_end, id_lo,
          [](const FragmentEntry& fe, uint32_t v) { return fe.spectrum_id < v; });
        const FragmentEntry* range_hi = std::lower_bound(
          range_lo, b_end, id_hi,
          [](const FragmentEntry& fe, uint32_t v) { return fe.spectrum_id < v; });

        for (const FragmentEntry* fe = range_lo; fe != range_hi; ++fe)
        {
          const uint32_t local_id = fe->spectrum_id - id_lo;
          if (use_im && !im_mask[local_id]) continue;

          const double frag_mz = fromBinIdx_(b, fe->mass_offset);
          if (std::fabs(frag_mz - qmz) > tol_da) continue;

          const double lib_weight = 1.0 / (static_cast<double>(fe->intensity_rank) + 1.0);

          accum[local_id].n_matched++;
          accum[local_id].sum_cross  += qint * lib_weight;
          accum[local_id].sum_lib_sq += lib_weight * lib_weight;
        }
      }
    }

    // -----------------------------------------------------------------------
    // Step 4: Compute cosine scores and collect passing results.
    // -----------------------------------------------------------------------

    std::vector<SearchResult> results;
    results.reserve(16);

    for (uint32_t i = 0; i < n_candidates; ++i)
    {
      if (static_cast<int>(accum[i].n_matched) < min_matched_peaks_) continue;

      const double denom  = std::sqrt(query_sq * accum[i].sum_lib_sq);
      const double cosine = (denom > 0.0) ? accum[i].sum_cross / denom : 0.0;

      SearchResult r;
      r.spectrum_id     = id_lo + i;
      r.n_matched_peaks = accum[i].n_matched;
      r.cosine_score    = cosine;
      results.push_back(r);
    }

    std::sort(results.begin(), results.end(),
      [](const SearchResult& a, const SearchResult& b) {
        return a.cosine_score > b.cosine_score;
      });

    return results;
  }


  // =========================================================================
  // Accessors
  // =========================================================================

  const DiaWeaverAlign::SpectrumEntry& DiaWeaverAlign::getSpectrumEntry(uint32_t id) const
  {
    return spectrum_entries_[id];
  }

  const String& DiaWeaverAlign::getSourceFile(Size idx) const
  {
    return source_files_[idx];
  }

  Size DiaWeaverAlign::getSpectrumCount() const
  {
    return spectrum_entries_.size();
  }

  Size DiaWeaverAlign::getTotalFragmentEntries() const
  {
    return fragment_entries_.size();
  }

  Size DiaWeaverAlign::getBinCount() const
  {
    return n_bins_;
  }

  double DiaWeaverAlign::getBinWidthPpm() const
  {
    return fragment_tol_ppm_;
  }

} // namespace OpenMS
