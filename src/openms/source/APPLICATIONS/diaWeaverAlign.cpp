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
    defaults_.setValue("lower_mz",      400.0, "Lower m/z bound of the fragment index.");
    defaults_.setValue("upper_mz",     2000.0, "Upper m/z bound of the fragment index.");
    defaults_.setValue("bin_width",       0.1, "m/z bin width in Daltons.");
    defaults_.setValue("lower_im",       0.60, "Lower ion mobility bound (Vs/cm^2).");
    defaults_.setValue("upper_im",       1.70, "Upper ion mobility bound (Vs/cm^2).");
    defaults_.setValue("bin_width_im",   0.01, "Ion mobility bin width (Vs/cm^2).");
    defaults_.setMinFloat("lower_mz",     1.0);
    defaults_.setMinFloat("upper_mz",     1.0);
    defaults_.setMinFloat("bin_width",   1e-3);
    defaults_.setMinFloat("lower_im",     0.0);
    defaults_.setMinFloat("upper_im",     0.0);
    defaults_.setMinFloat("bin_width_im", 1e-4);
    defaultsToParam_();
  }

  void DiaWeaverAlign::updateMembers_()
  {
    lower_mz_      = (double)param_.getValue("lower_mz");
    upper_mz_      = (double)param_.getValue("upper_mz");
    bin_width_     = (double)param_.getValue("bin_width");
    n_bins_        = static_cast<uint32_t>(std::ceil((upper_mz_ - lower_mz_) / bin_width_));
    lower_im_      = (double)param_.getValue("lower_im");
    upper_im_      = (double)param_.getValue("upper_im");
    bin_width_im_  = (double)param_.getValue("bin_width_im");
    n_im_bins_     = static_cast<uint32_t>(std::ceil((upper_im_ - lower_im_) / bin_width_im_));
  }

  uint32_t DiaWeaverAlign::toBinIdx_(double mz) const
  {
    return static_cast<uint32_t>((mz - lower_mz_) / bin_width_);
  }

  uint32_t DiaWeaverAlign::toBinIdx_im_(double im) const
  {
    return static_cast<uint32_t>((im - lower_im_) / bin_width_im_);
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
    im_detected_ = false;

    OPENMS_LOG_INFO << "[DiaWeaverAlign] Building fragment index  "
                    << "mz=[" << lower_mz_ << ", " << upper_mz_ << "] Da  "
                    << "bin_width=" << bin_width_ << " Da  n_mz_bins=" << n_bins_ << "  "
                    << "im=[" << lower_im_ << ", " << upper_im_ << "] Vs/cm^2  "
                    << "bin_width_im=" << bin_width_im_ << "  n_im_bins=" << n_im_bins_ << std::endl;

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
          im_detected_ = (spec.getDriftTime() != IMTypes::DRIFTTIME_NOT_SET);
          if (im_detected_ != spec.containsIMData())
          {
            throw Exception::MissingInformation(
              __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Spectrum '" + spec.getNativeID() + "' in '" + mzml_files[file_idx] +
              "': precursor drift time and fragment-level IM data are inconsistent.");
          }
          OPENMS_LOG_INFO << "[DiaWeaverAlign] Ion mobility: "
                          << (im_detected_ ? "detected" : "not detected") << std::endl;
        }
        else if (im_detected_ &&
                 (spec.getDriftTime() == IMTypes::DRIFTTIME_NOT_SET || !spec.containsIMData()))
        {
          throw Exception::MissingInformation(
            __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Ion mobility was detected in earlier spectra but spectrum '" +
            spec.getNativeID() + "' in '" + mzml_files[file_idx] + "' has incomplete IM data.");
        }

        RawSpecEntry entry;
        entry.retention_time = spec.getRT();
        if (im_detected_)
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
    // Peaks are m/z-sorted, so mz_bin values are already sorted in ascending order
    // The only unordered dimension (when IM is present) is IM within the same
    // mz_bin.  We exploit this with a flat array of n_im_bins_ slots that
    // accumulates the highest intensity peak per IM bin for the current mz_bin.
    // Only occupied slots are flushed and reset, giving O(n_peaks) total work.
    //
    // When IM is absent the index is 1D and a simple streaming scan suffices.
    // -----------------------------------------------------------------------

    std::vector<std::pair<uint32_t, FragmentEntry>> all_entries;
    /// TODO: diaWeaver by default outputs 500 peaks per spectrum maximum. Is there a better way?
    all_entries.reserve(raw_spectra.size() * 500);

    if (im_detected_)
    {
      // Per-IM-bin accumulators — allocated once, reused across all spectra.
      std::vector<float>         im_best_intensity(n_im_bins_, -1.0f);
      std::vector<FragmentEntry> im_best_fe(n_im_bins_);
      std::vector<uint32_t>      occupied_im_bins;
      occupied_im_bins.reserve(n_im_bins_);

      for (uint32_t spec_id = 0;
           spec_id < static_cast<uint32_t>(raw_spectra.size());
           ++spec_id)
      {
        const MSSpectrum& spec    = raw_spectra[spec_id].spectrum;
        const uint8_t charge_byte = static_cast<uint8_t>(raw_spectra[spec_id].charge);
        const auto& im_array      = spec.getFloatDataArrays()[spec.getIMData().first];

        uint32_t current_mz_bin = std::numeric_limits<uint32_t>::max();

        // Flush occupied IM slots for the current mz_bin and reset them.
        auto flush_mz_bin = [&]()
        {
          for (uint32_t im_b : occupied_im_bins)
          {
            all_entries.emplace_back(current_mz_bin * n_im_bins_ + im_b, im_best_fe[im_b]);
            im_best_intensity[im_b] = -1.0f;
          }
          occupied_im_bins.clear();
        };

        for (Size i = 0; i < spec.size(); ++i)
        {
          const double mz        = spec[i].getMZ();
          const float  intensity = spec[i].getIntensity();
          const float  im        = im_array[i];

          if (mz < lower_mz_ || mz >= upper_mz_) continue;
          if (im < lower_im_ || im >= upper_im_) continue;

          const uint32_t mz_bin = toBinIdx_(mz);
          const uint32_t im_bin = toBinIdx_im_(im);

          if (mz_bin != current_mz_bin)
          {
            if (current_mz_bin != std::numeric_limits<uint32_t>::max())
              flush_mz_bin();
            current_mz_bin = mz_bin;
          }

          if (im_best_intensity[im_bin] < 0.0f)  // first peak in this IM slot
          {
            occupied_im_bins.push_back(im_bin);
            im_best_intensity[im_bin] = intensity;
            im_best_fe[im_bin] = {spec_id, toMassOffset_(mz, mz_bin), charge_byte, 0};
          }
          else if (intensity > im_best_intensity[im_bin])
          {
            im_best_intensity[im_bin]              = intensity;
            im_best_fe[im_bin].mass_offset = toMassOffset_(mz, mz_bin);
          }
        }

        if (current_mz_bin != std::numeric_limits<uint32_t>::max())
          flush_mz_bin();
      }
    }
    else
    {
      // 1D index: peaks are m/z-sorted so bin transitions are already in order.
      for (uint32_t spec_id = 0;
           spec_id < static_cast<uint32_t>(raw_spectra.size());
           ++spec_id)
      {
        const MSSpectrum& spec    = raw_spectra[spec_id].spectrum;
        const uint8_t charge_byte = static_cast<uint8_t>(raw_spectra[spec_id].charge);

        uint32_t      current_bin    = std::numeric_limits<uint32_t>::max();
        float         best_intensity = 0.0f;
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
            current_bin         = bin_idx;
            best_intensity      = intensity;
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

        emit();
      }
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

    const Size n_flat_bins = im_detected_
      ? static_cast<Size>(n_bins_) * n_im_bins_
      : n_bins_;

    bin_offsets_.assign(n_flat_bins + 1, 0);
    for (const auto& [flat_idx, _] : all_entries)
      ++bin_offsets_[flat_idx + 1];
    for (Size b = 1; b <= n_flat_bins; ++b)
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
  DiaWeaverAlign::getBinEntries(uint32_t mz_bin, uint32_t im_bin) const
  {
    const uint32_t flat_idx = im_detected_
      ? mz_bin * n_im_bins_ + im_bin
      : mz_bin;
    const FragmentEntry* base = fragment_entries_.data();
    return { base + bin_offsets_[flat_idx], base + bin_offsets_[flat_idx + 1] };
  }

  Size   DiaWeaverAlign::getSpectrumCount()        const { return spectrum_entries_.size(); }
  Size   DiaWeaverAlign::getTotalFragmentEntries() const { return fragment_entries_.size(); }
  Size   DiaWeaverAlign::getMZBinCount()           const { return n_bins_; }
  Size   DiaWeaverAlign::getIMBinCount()           const { return n_im_bins_; }
  double DiaWeaverAlign::getBinWidth()             const { return bin_width_; }
  double DiaWeaverAlign::getBinWidthIM()           const { return bin_width_im_; }
  double DiaWeaverAlign::getLowerMz()              const { return lower_mz_; }
  double DiaWeaverAlign::getUpperMz()              const { return upper_mz_; }
  double DiaWeaverAlign::getLowerIM()              const { return lower_im_; }
  double DiaWeaverAlign::getUpperIM()              const { return upper_im_; }
  bool   DiaWeaverAlign::hasIM()                   const { return im_detected_; }

} // namespace OpenMS
