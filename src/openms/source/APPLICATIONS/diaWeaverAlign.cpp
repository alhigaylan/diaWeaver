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
      "Defines the width of each bin on the logarithmic m/z axis: every bin "
      "spans exactly this many ppm, so bin width in Daltons grows proportionally "
      "with m/z.");
    defaults_.setMinFloat("fragment_mz_tolerance_ppm", 0.1);

    defaults_.setValue(
      "max_peaks_per_spectrum", 200,
      "Maximum number of peaks per spectrum to index, ranked by descending "
      "intensity. Limits index size and suppresses noise peaks.");
    defaults_.setMinInt("max_peaks_per_spectrum", 1);

    defaultsToParam_();
  }

  void DiaWeaverAlign::updateMembers_()
  {
    fragment_tol_ppm_       = (double)param_.getValue("fragment_mz_tolerance_ppm");
    max_peaks_per_spectrum_ = (Size)(int)param_.getValue("max_peaks_per_spectrum");
    // log_step_, min_observed_mz_, and n_bins_ are derived from data in buildIndex.
  }


  // =========================================================================
  // Private helpers — log-space bin axis
  // =========================================================================

  uint32_t DiaWeaverAlign::toBinIdx_(double mz) const
  {
    return static_cast<uint32_t>(std::log(mz / min_observed_mz_) / log_step_);
  }

  uint16_t DiaWeaverAlign::toMassOffset_(double mz, uint32_t bin_idx) const
  {
    const double log_pos = std::log(mz / min_observed_mz_) / log_step_;
    const double frac    = log_pos - static_cast<double>(bin_idx);
    return static_cast<uint16_t>(std::max(0.0, std::min(1.0, frac)) * 65535.0);
  }

  double DiaWeaverAlign::fromBinIdx_(uint32_t bin_idx, uint16_t mass_offset) const
  {
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
    log_step_        = 0.0;
    min_observed_mz_ = 0.0;
    n_bins_          = 0;

    OPENMS_LOG_INFO << "[DiaWeaverAlign] Building fragment index from "
                    << mzml_files.size() << " file(s)." << std::endl;

    // -----------------------------------------------------------------------
    // Phase 1: Load all MS2 spectra from every input file.
    //
    // All spectra are collected first so that RT-based sorting and ID
    // assignment can be done globally across the entire dataset.
    // -----------------------------------------------------------------------

    struct RawSpecEntry
    {
      MSSpectrum spectrum;
      double     retention_time{-1.0};
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
        entry.retention_time = spec.getRT();
        entry.precursor_mz   = pmz;
        entry.precursor_im   = pim;
        entry.charge         = pcharge;
        entry.native_id      = spec.getNativeID();
        entry.source_idx     = file_idx;
        entry.spectrum       = std::move(spec);
        raw_spectra.push_back(std::move(entry));
      }
    }

    OPENMS_LOG_INFO << "[DiaWeaverAlign] Loaded " << raw_spectra.size()
                    << " MS2 spectra total." << std::endl;

    if (raw_spectra.empty()) return;

    // -----------------------------------------------------------------------
    // Phase 2: Scan all peaks to determine the global m/z range.
    //
    // The bin axis is anchored at min_observed_mz_ and spans to
    // max_observed_mz. Both values come from the data; no user parameter
    // is required.
    //
    //   log_step = log(1 + ppm / 1e6)
    //   n_bins   = ceil(log(max_mz / min_mz) / log_step) + 1
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

    OPENMS_LOG_INFO << "[DiaWeaverAlign] m/z range: ["
                    << min_observed_mz_ << ", " << max_observed_mz << "] Da  |  "
                    << n_bins_ << " bins at " << fragment_tol_ppm_ << " ppm each."
                    << std::endl;

    // -----------------------------------------------------------------------
    // Phase 3: Sort spectra by retention time and assign contiguous IDs.
    //
    // spectrum_id 0 is the earliest eluting spectrum across all input files.
    // Within each bin, fragment entries will be ordered by spectrum_id,
    // i.e. by retention time.
    // -----------------------------------------------------------------------

    std::sort(raw_spectra.begin(), raw_spectra.end(),
      [](const RawSpecEntry& a, const RawSpecEntry& b) {
        return a.retention_time < b.retention_time;
      });

    spectrum_entries_.resize(raw_spectra.size());
    for (Size i = 0; i < raw_spectra.size(); ++i)
    {
      spectrum_entries_[i].retention_time  = raw_spectra[i].retention_time;
      spectrum_entries_[i].precursor_mz    = raw_spectra[i].precursor_mz;
      spectrum_entries_[i].precursor_im    = raw_spectra[i].precursor_im;
      spectrum_entries_[i].charge          = raw_spectra[i].charge;
      spectrum_entries_[i].native_id       = raw_spectra[i].native_id;
      spectrum_entries_[i].source_file_idx = raw_spectra[i].source_idx;
    }

    // -----------------------------------------------------------------------
    // Phase 4: Generate fragment entries.
    //
    // For each spectrum (in RT order), take the top max_peaks_per_spectrum
    // peaks by intensity and place each one into its bin as a FragmentEntry.
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

      // Rank peaks by descending intensity.
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
    // Phase 5: Sort by (bin_idx, spectrum_id) and build CSR offsets.
    //
    // After sorting, entries within every bin are in RT order (ascending
    // spectrum_id). bin_offsets_ is built via a prefix sum over per-bin counts.
    // -----------------------------------------------------------------------

    std::sort(all_entries.begin(), all_entries.end(),
      [](const std::pair<uint32_t, FragmentEntry>& a,
         const std::pair<uint32_t, FragmentEntry>& b) {
        if (a.first != b.first) return a.first < b.first;
        return a.second.spectrum_id < b.second.spectrum_id;
      });

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

  double DiaWeaverAlign::getLogStep() const
  {
    return log_step_;
  }

  double DiaWeaverAlign::getMinObservedMz() const
  {
    return min_observed_mz_;
  }

} // namespace OpenMS
