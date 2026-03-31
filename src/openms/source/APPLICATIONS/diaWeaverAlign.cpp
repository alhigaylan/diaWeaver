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

  DiaWeaverAlign::DiaWeaverAlign() :
    DefaultParamHandler("DiaWeaverAlign")
  {
    defaults_.setValue("lower_mz",   400.0, "Lower m/z bound of the fragment index.");
    defaults_.setValue("upper_mz",  2000.0, "Upper m/z bound of the fragment index.");
    defaults_.setValue("bin_width",    0.1, "Bin width in Daltons.");
    defaults_.setMinFloat("lower_mz",   1.0);
    defaults_.setMinFloat("upper_mz",   1.0);
    defaults_.setMinFloat("bin_width", 1e-3);
    defaultsToParam_();
  }

  void DiaWeaverAlign::updateMembers_()
  {
    lower_mz_  = (double)param_.getValue("lower_mz");
    upper_mz_  = (double)param_.getValue("upper_mz");
    bin_width_ = (double)param_.getValue("bin_width");
    n_bins_    = static_cast<uint32_t>(std::ceil((upper_mz_ - lower_mz_) / bin_width_));
  }

  uint32_t DiaWeaverAlign::toBinIdx_(double mz) const
  {
    return static_cast<uint32_t>((mz - lower_mz_) / bin_width_);
  }

  uint16_t DiaWeaverAlign::toMassOffset_(double mz, uint32_t bin_idx) const
  {
    const double bin_start = lower_mz_ + bin_idx * bin_width_;
    const double frac      = (mz - bin_start) / bin_width_;
    return static_cast<uint16_t>(std::max(0.0, std::min(1.0, frac)) * 65535.0);
  }

  void DiaWeaverAlign::buildIndex(const std::vector<String>& mzml_files)
  {
    updateMembers_();

    source_files_.clear();
    spectrum_entries_.clear();
    fragment_entries_.clear();
    bin_offsets_.clear();

    OPENMS_LOG_INFO << "[DiaWeaverAlign] Building fragment index  "
                    << "range=[" << lower_mz_ << ", " << upper_mz_ << "] Da  "
                    << "bin_width=" << bin_width_ << " Da  "
                    << "n_bins=" << n_bins_ << std::endl;

    // -----------------------------------------------------------------------
    // Phase 1: Load all MS2 spectra.
    // -----------------------------------------------------------------------

    struct RawSpecEntry
    {
      MSSpectrum spectrum;
      double     retention_time{-1.0};
      int        charge{0};
      String     native_id;
      Size       source_idx{0};
    };

    std::vector<RawSpecEntry> raw_spectra;

    for (Size file_idx = 0; file_idx < mzml_files.size(); ++file_idx)
    {
      source_files_.push_back(mzml_files[file_idx]);

      MSExperiment exp;
      MzMLFile().load(mzml_files[file_idx], exp);

      for (auto& spec : exp)
      {
        if (spec.getMSLevel() != 2 || spec.empty()) continue;

        int pcharge = 0;
        if (!spec.getPrecursors().empty())
          pcharge = spec.getPrecursors()[0].getCharge();

        RawSpecEntry entry;
        entry.retention_time = spec.getRT();
        entry.charge         = pcharge;
        entry.native_id      = spec.getNativeID();
        entry.source_idx     = file_idx;
        entry.spectrum       = std::move(spec);
        raw_spectra.push_back(std::move(entry));
      }
    }

    OPENMS_LOG_INFO << "[DiaWeaverAlign] Loaded " << raw_spectra.size()
                    << " MS2 spectra." << std::endl;

    if (raw_spectra.empty()) return;

    // -----------------------------------------------------------------------
    // Phase 2: Sort by retention time and assign contiguous spectrum IDs.
    // -----------------------------------------------------------------------

    std::sort(raw_spectra.begin(), raw_spectra.end(),
      [](const RawSpecEntry& a, const RawSpecEntry& b) {
        return a.retention_time < b.retention_time;
      });

    spectrum_entries_.resize(raw_spectra.size());
    for (Size i = 0; i < raw_spectra.size(); ++i)
    {
      spectrum_entries_[i].retention_time  = raw_spectra[i].retention_time;
      spectrum_entries_[i].native_id       = raw_spectra[i].native_id;
      spectrum_entries_[i].source_file_idx = raw_spectra[i].source_idx;
    }

    // -----------------------------------------------------------------------
    // Phase 3: Generate fragment entries.
    //
    // For each spectrum, iterate all peaks and place each one into its bin.
    // Peaks outside [lower_mz_, upper_mz_] are skipped.
    // Peaks within a spectrum are ranked by descending intensity so the
    // intensity_rank field reflects how prominent each peak is.
    // -----------------------------------------------------------------------

    std::vector<std::pair<uint32_t, FragmentEntry>> all_entries;
    all_entries.reserve(raw_spectra.size() * 50);

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

      for (Size rank = 0; rank < spec.size(); ++rank)
      {
        const double mz = spec[order[rank]].getMZ();
        if (mz < lower_mz_ || mz >= upper_mz_) continue;

        const uint32_t bin_idx = toBinIdx_(mz);

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
    // Phase 4: Sort by (bin_idx, spectrum_id) and build CSR offsets.
    //
    // Entries within each bin are in RT order (ascending spectrum_id).
    // -----------------------------------------------------------------------

    std::sort(all_entries.begin(), all_entries.end(),
      [](const std::pair<uint32_t, FragmentEntry>& a,
         const std::pair<uint32_t, FragmentEntry>& b) {
        if (a.first != b.first) return a.first < b.first;
        return a.second.spectrum_id < b.second.spectrum_id;
      });

    bin_offsets_.assign(n_bins_ + 1, 0);
    for (const auto& [bin_idx, _] : all_entries)
      ++bin_offsets_[bin_idx + 1];
    for (Size b = 1; b <= n_bins_; ++b)
      bin_offsets_[b] += bin_offsets_[b - 1];

    fragment_entries_.resize(all_entries.size());
    for (Size i = 0; i < all_entries.size(); ++i)
      fragment_entries_[i] = all_entries[i].second;

    OPENMS_LOG_INFO << "[DiaWeaverAlign] Index complete: "
                    << spectrum_entries_.size() << " spectra, "
                    << fragment_entries_.size() << " fragment entries." << std::endl;
  }

  // -------------------------------------------------------------------------
  // Accessors
  // -------------------------------------------------------------------------

  const DiaWeaverAlign::SpectrumEntry& DiaWeaverAlign::getSpectrumEntry(uint32_t id) const
  {
    return spectrum_entries_[id];
  }

  const String& DiaWeaverAlign::getSourceFile(Size idx) const
  {
    return source_files_[idx];
  }

  std::pair<const DiaWeaverAlign::FragmentEntry*, const DiaWeaverAlign::FragmentEntry*>
  DiaWeaverAlign::getBinEntries(uint32_t bin_idx) const
  {
    const FragmentEntry* base = fragment_entries_.data();
    return { base + bin_offsets_[bin_idx], base + bin_offsets_[bin_idx + 1] };
  }

  Size   DiaWeaverAlign::getSpectrumCount()        const { return spectrum_entries_.size(); }
  Size   DiaWeaverAlign::getTotalFragmentEntries() const { return fragment_entries_.size(); }
  Size   DiaWeaverAlign::getBinCount()             const { return n_bins_; }
  double DiaWeaverAlign::getBinWidth()             const { return bin_width_; }
  double DiaWeaverAlign::getLowerMz()              const { return lower_mz_; }
  double DiaWeaverAlign::getUpperMz()              const { return upper_mz_; }

} // namespace OpenMS
