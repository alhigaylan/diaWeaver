// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Mohammed Alhigaylan $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/diaWeaverAlign.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <algorithm>
#include <cmath>

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
    precursor_im_detected_ = false;

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
      double     precursor_mz{-1.0};
      double     drift_time{-1.0};
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

        if (raw_spectra.empty())  // first MS2 spectrum sets the IM expectation
        {
          precursor_im_detected_ = (spec.getDriftTime() != IMTypes::DRIFTTIME_NOT_SET);
          OPENMS_LOG_INFO << "[DiaWeaverAlign] Ion mobility: "
                          << (precursor_im_detected_ ? "detected" : "not detected") << std::endl;
        }
        else if (precursor_im_detected_ && spec.getDriftTime() == IMTypes::DRIFTTIME_NOT_SET)
        {
          throw Exception::MissingInformation(
            __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Ion mobility was detected in earlier spectra but spectrum '" +
            spec.getNativeID() + "' in '" + mzml_files[file_idx] + "' has no drift time.");
        }

        RawSpecEntry entry;
        entry.retention_time = spec.getRT();
        if (precursor_im_detected_)
          entry.drift_time = spec.getDriftTime();
        entry.charge         = spec.getPrecursors()[0].getCharge();
        entry.precursor_mz   = spec.getPrecursors()[0].getMZ();
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
    // Phase 2: Assign contiguous spectrum IDs in load order
    //          (file order, then within-file order).
    //
    // Each spectrum_id is an index into spectrum_entries_, giving O(1) metadata
    // lookup from a FragmentEntry without needing a (file_idx, per-file-id) pair.
    // -----------------------------------------------------------------------

    spectrum_entries_.resize(raw_spectra.size());
    for (Size i = 0; i < raw_spectra.size(); ++i)
    {
      spectrum_entries_[i].retention_time   = raw_spectra[i].retention_time;
      spectrum_entries_[i].precursor_mz     = raw_spectra[i].precursor_mz;
      spectrum_entries_[i].drift_time       = raw_spectra[i].drift_time;
      spectrum_entries_[i].precursor_charge = raw_spectra[i].charge;
      spectrum_entries_[i].native_id        = raw_spectra[i].native_id;
      spectrum_entries_[i].source_file_idx  = raw_spectra[i].source_idx;
    }

    // -----------------------------------------------------------------------
    // Phase 3: Generate fragment entries.
    //
    // We scan linearly: when a peak crosses into a new
    // bin we emit the best (highest intensity) entry for the completed bin and
    // start tracking the new one.
    // -----------------------------------------------------------------------

    std::vector<std::pair<uint32_t, FragmentEntry>> all_entries;
    all_entries.reserve(raw_spectra.size() * 50);

    for (uint32_t spec_id = 0;
         spec_id < static_cast<uint32_t>(raw_spectra.size());
         ++spec_id)
    {
      const MSSpectrum& spec    = raw_spectra[spec_id].spectrum;
      const uint8_t charge_byte = static_cast<uint8_t>(raw_spectra[spec_id].charge);

      uint32_t     current_bin = std::numeric_limits<uint32_t>::max();
      float        best_intensity = 0.0f;
      FragmentEntry best_fe{};

      auto emit = [&]() {
        if (current_bin != std::numeric_limits<uint32_t>::max())
          all_entries.emplace_back(current_bin, best_fe);
      };

      for (Size i = 0; i < spec.size(); ++i)
      {
        const double mz        = spec[i].getMZ();
        const float  intensity = spec[i].getIntensity();

        if (mz < lower_mz_ || mz >= upper_mz_) continue;

        const uint32_t bin_idx = toBinIdx_(mz);

        if (bin_idx != current_bin)
        {
          emit();
          current_bin       = bin_idx;
          best_intensity    = intensity;
          best_fe.spectrum_id = spec_id;
          best_fe.mass_offset = toMassOffset_(mz, bin_idx);
          best_fe.charge      = charge_byte;
          best_fe.reserved    = 0;
        }
        else if (intensity > best_intensity)
        {
          best_intensity      = intensity;
          best_fe.mass_offset = toMassOffset_(mz, bin_idx);
        }
      }

      emit(); // flush the last bin
    }

    OPENMS_LOG_INFO << "[DiaWeaverAlign] Generated " << all_entries.size()
                    << " fragment entries." << std::endl;

    // -----------------------------------------------------------------------
    // Phase 4: Sort by (bin_idx, spectrum_id) and build CSR offsets.
    //
    // Entries within each bin are in load order (ascending spectrum_id).
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
  bool   DiaWeaverAlign::hasPrecursorIM()          const { return precursor_im_detected_; }

} // namespace OpenMS
