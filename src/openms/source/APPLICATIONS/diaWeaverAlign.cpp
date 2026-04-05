// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Mohammed Alhigaylan $
// --------------------------------------------------------------------------

#include "OpenMS/CONCEPT/LogStream.h"
#include <OpenMS/APPLICATIONS/diaWeaverAlign.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <unordered_map>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS
{

  DiaWeaverAlign::DiaWeaverAlign() :
    DefaultParamHandler("DiaWeaverAlign")
  {
    defaults_.setValue("lower_mz",                400.0, "Lower m/z bound of the fragment index.");
    defaults_.setValue("upper_mz",               2000.0, "Upper m/z bound of the fragment index.");
    defaults_.setValue("bin_width",                 0.1, "Fragment m/z bin width in Daltons.");
    defaults_.setValue("lower_im",                 0.60, "Lower ion mobility bound (Vs/cm^2).");
    defaults_.setValue("upper_im",                 1.70, "Upper ion mobility bound (Vs/cm^2).");
    defaults_.setValue("bin_width_im",             0.01, "Ion mobility bin width (Vs/cm^2).");
    defaults_.setValue("lower_precursor_mz",      400.0, "Lower m/z bound of the precursor index.");
    defaults_.setValue("upper_precursor_mz",     1200.0, "Upper m/z bound of the precursor index.");
    defaults_.setValue("precursor_bin_width",      0.02, "Precursor m/z bin width in Daltons.");
    defaults_.setMinFloat("lower_mz",               1.0);
    defaults_.setMinFloat("upper_mz",               1.0);
    defaults_.setMinFloat("bin_width",             1e-3);
    defaults_.setMinFloat("lower_im",               0.0);
    defaults_.setMinFloat("upper_im",               0.0);
    defaults_.setMinFloat("bin_width_im",          1e-4);
    defaults_.setMinFloat("lower_precursor_mz",     1.0);
    defaults_.setMinFloat("upper_precursor_mz",     1.0);
    defaults_.setMinFloat("precursor_bin_width",   1e-4);
    defaultsToParam_();
  }

  DiaWeaverAlign::~DiaWeaverAlign()
  {
    // Close any streaming FILE handles that may be open (e.g. after an exception).
    if (stream_frags_fp_) { std::fclose(stream_frags_fp_); stream_frags_fp_ = nullptr; }
    if (stream_meta_fp_)  { std::fclose(stream_meta_fp_);  stream_meta_fp_  = nullptr; }
  }

  void DiaWeaverAlign::updateMembers_()
  {
    lower_mz_             = (double)param_.getValue("lower_mz");
    upper_mz_             = (double)param_.getValue("upper_mz");
    bin_width_            = (double)param_.getValue("bin_width");
    n_bins_               = static_cast<uint32_t>(std::ceil((upper_mz_ - lower_mz_) / bin_width_));
    lower_im_             = (double)param_.getValue("lower_im");
    upper_im_             = (double)param_.getValue("upper_im");
    bin_width_im_         = (double)param_.getValue("bin_width_im");
    n_im_bins_            = static_cast<uint32_t>(std::ceil((upper_im_ - lower_im_) / bin_width_im_));
    lower_precursor_mz_   = (double)param_.getValue("lower_precursor_mz");
    upper_precursor_mz_   = (double)param_.getValue("upper_precursor_mz");
    precursor_bin_width_  = (double)param_.getValue("precursor_bin_width");
    n_precursor_bins_     = static_cast<uint32_t>(std::ceil((upper_precursor_mz_ - lower_precursor_mz_) / precursor_bin_width_));
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

  uint32_t DiaWeaverAlign::toPrecursorBinIdx_(double mz) const
  {
    return static_cast<uint32_t>((mz - lower_precursor_mz_) / precursor_bin_width_);
  }

  void DiaWeaverAlign::buildPrecursorIndex()
  {
    updateMembers_();

    precursor_entries_.clear();
    precursor_offsets_.clear();

    if (spectrum_entries_.empty()) return;

    const bool use_im  = hasIM();
    const Size n_flat  = use_im
      ? static_cast<Size>(n_precursor_bins_) * n_im_bins_
      : n_precursor_bins_;

    OPENMS_LOG_INFO << "[DiaWeaverAlign] Building precursor index  "
                    << "mz=[" << lower_precursor_mz_ << ", " << upper_precursor_mz_ << "] Da  "
                    << "precursor_bin_width=" << precursor_bin_width_ << " Da  "
                    << "n_precursor_bins=" << n_precursor_bins_ << std::endl;

    // Build (flat_idx, spectrum_id) pairs from spectrum metadata.
    std::vector<std::pair<uint32_t, uint32_t>> entries;
    entries.reserve(spectrum_entries_.size());

    Size n_out_of_range = 0;
    for (Size i = 0; i < spectrum_entries_.size(); ++i)
    {
      const SpectrumEntry& se = spectrum_entries_[i];
      if (se.precursor_mz < lower_precursor_mz_ || se.precursor_mz >= upper_precursor_mz_)
      {
        ++n_out_of_range;
        OPENMS_LOG_DEBUG << "[DiaWeaverAlign::buildPrecursorIndex] spectrum_id=" << i
                         << " skipped: prec_mz=" << se.precursor_mz
                         << " outside [" << lower_precursor_mz_ << ", "
                         << upper_precursor_mz_ << ")" << std::endl;
        continue;
      }

      const uint32_t mz_bin = toPrecursorBinIdx_(se.precursor_mz);

      uint32_t flat_idx;
      if (use_im)
      {
        if (se.drift_time < lower_im_ || se.drift_time >= upper_im_)
        {
          ++n_out_of_range;
          OPENMS_LOG_DEBUG << "[DiaWeaverAlign::buildPrecursorIndex] spectrum_id=" << i
                           << " skipped: drift_time=" << se.drift_time
                           << " outside [" << lower_im_ << ", " << upper_im_ << ")" << std::endl;
          continue;
        }
        const uint32_t im_bin = toBinIdx_im_(se.drift_time);
        flat_idx = mz_bin * n_im_bins_ + im_bin;
        OPENMS_LOG_DEBUG << "[DiaWeaverAlign::buildPrecursorIndex] spectrum_id=" << i
                         << "  prec_mz=" << se.precursor_mz
                         << "  -> mz_bin=" << mz_bin
                         << "  drift_time=" << se.drift_time
                         << "  -> im_bin=" << im_bin
                         << "  flat_idx=" << flat_idx << std::endl;
      }
      else
      {
        flat_idx = mz_bin;
        OPENMS_LOG_DEBUG << "[DiaWeaverAlign::buildPrecursorIndex] spectrum_id=" << i
                         << "  prec_mz=" << se.precursor_mz
                         << "  -> mz_bin=" << mz_bin << std::endl;
      }

      entries.emplace_back(flat_idx, static_cast<uint32_t>(i));
    }

    OPENMS_LOG_DEBUG << "[DiaWeaverAlign::buildPrecursorIndex] indexed=" << entries.size()
                     << "  skipped(out-of-range)=" << n_out_of_range << std::endl;

    std::sort(entries.begin(), entries.end());

    precursor_offsets_.assign(n_flat + 1, 0);
    for (const auto& [flat_idx, _] : entries)
      ++precursor_offsets_[flat_idx + 1];
    for (Size b = 1; b <= n_flat; ++b)
      precursor_offsets_[b] += precursor_offsets_[b - 1];

    precursor_entries_.resize(entries.size());
    for (Size i = 0; i < entries.size(); ++i)
      precursor_entries_[i] = entries[i].second;

    OPENMS_LOG_INFO << "[DiaWeaverAlign] Precursor index complete: "
                    << precursor_entries_.size() << " entries." << std::endl;
  }

  // ---------------------------------------------------------------------------
  // Private helper: populate RawSpecEntry_ from a single MSExperiment.
  // Accepts pseudo-MS2 spectra regardless of MS level (diaWeaver emits them
  // as MS level 1). A spectrum is accepted when it is non-empty and carries
  // at least one precursor.
  // ---------------------------------------------------------------------------
  bool DiaWeaverAlign::collectSpectraFromExperiment_(
    const MSExperiment& exp,
    Size file_idx,
    const String& file_label,
    bool& im_detected,
    bool im_initialised,
    std::vector<RawSpecEntry_>& raw_spectra)
  {
    for (const MSSpectrum& spec_ref : exp)
    {
      // Pseudo-MS2 spectra from diaWeaver are MS level 2. Skip anything else.
      if (spec_ref.getMSLevel() != 2 || spec_ref.empty()) continue;

      MSSpectrum spec = spec_ref; // need a mutable copy for std::move below

      if (!im_initialised)
      {
        im_detected = (spec.getDriftTime() != IMTypes::DRIFTTIME_NOT_SET);
        if (im_detected != spec.containsIMData())
        {
          throw Exception::MissingInformation(
            __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Spectrum '" + spec.getNativeID() + "' in '" + file_label +
            "': precursor drift time and fragment-level IM data are inconsistent.");
        }
        OPENMS_LOG_INFO << "[DiaWeaverAlign] Ion mobility: "
                        << (im_detected ? "detected" : "not detected") << std::endl;
        im_initialised = true;
      }
      else if (im_detected &&
               (spec.getDriftTime() == IMTypes::DRIFTTIME_NOT_SET || !spec.containsIMData()))
      {
        throw Exception::MissingInformation(
          __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Ion mobility was detected in earlier spectra but spectrum '" +
          spec.getNativeID() + "' in '" + file_label + "' has incomplete IM data.");
      }

      DiaWeaverAlign::RawSpecEntry_ entry;
      entry.retention_time = spec.getRT();
      if (im_detected)
        entry.drift_time = spec.getDriftTime();
      entry.charge       = spec.getPrecursors()[0].getCharge();
      entry.precursor_mz = spec.getPrecursors()[0].getMZ();
      entry.native_id    = spec.getNativeID();
      entry.source_idx   = file_idx;
      entry.spectrum     = std::move(spec);
      raw_spectra.push_back(std::move(entry));
    }
    return im_initialised;
  }

  void DiaWeaverAlign::buildIndex(const std::vector<String>& mzml_files)
  {
    updateMembers_();

    source_files_.clear();
    spectrum_entries_.clear();
    fragment_entries_.clear();
    bin_offsets_.clear();

    bool im_detected     = false;
    bool im_initialised  = false;

    OPENMS_LOG_INFO << "[DiaWeaverAlign] Building fragment index  "
                    << "mz=[" << lower_mz_ << ", " << upper_mz_ << "] Da  "
                    << "bin_width=" << bin_width_ << " Da  n_mz_bins=" << n_bins_ << "  "
                    << "im=[" << lower_im_ << ", " << upper_im_ << "] Vs/cm^2  "
                    << "bin_width_im=" << bin_width_im_ << "  n_im_bins=" << n_im_bins_ << std::endl;

    // -----------------------------------------------------------------------
    // Phase 1: Load all pseudo-MS2 spectra.
    // -----------------------------------------------------------------------

    std::vector<RawSpecEntry_> raw_spectra;

    for (Size file_idx = 0; file_idx < mzml_files.size(); ++file_idx)
    {
      source_files_.push_back(mzml_files[file_idx]);

      MSExperiment exp;
      MzMLFile().load(mzml_files[file_idx], exp);

      im_initialised = collectSpectraFromExperiment_(
        exp, file_idx, mzml_files[file_idx], im_detected, im_initialised, raw_spectra);
    }

    OPENMS_LOG_INFO << "[DiaWeaverAlign] Loaded " << raw_spectra.size()
                    << " pseudo-MS2 spectra." << std::endl;

    if (raw_spectra.empty()) return;

    buildIndexCore_(raw_spectra);
  }

  void DiaWeaverAlign::buildIndex(const std::vector<MSExperiment>& experiments,
                                   const std::vector<String>& source_names)
  {
    updateMembers_();

    source_files_.clear();
    spectrum_entries_.clear();
    fragment_entries_.clear();
    bin_offsets_.clear();

    bool im_detected    = false;
    bool im_initialised = false;

    OPENMS_LOG_INFO << "[DiaWeaverAlign] Building fragment index (in-memory)  "
                    << "mz=[" << lower_mz_ << ", " << upper_mz_ << "] Da  "
                    << "bin_width=" << bin_width_ << " Da  n_mz_bins=" << n_bins_ << "  "
                    << "im=[" << lower_im_ << ", " << upper_im_ << "] Vs/cm^2  "
                    << "bin_width_im=" << bin_width_im_ << "  n_im_bins=" << n_im_bins_ << std::endl;

    std::vector<RawSpecEntry_> raw_spectra;

    for (Size exp_idx = 0; exp_idx < experiments.size(); ++exp_idx)
    {
      const String& label = (exp_idx < source_names.size())
        ? source_names[exp_idx]
        : String("experiment[") + String(exp_idx) + "]";
      source_files_.push_back(label);

      im_initialised = collectSpectraFromExperiment_(
        experiments[exp_idx], exp_idx, label, im_detected, im_initialised, raw_spectra);
    }

    OPENMS_LOG_INFO << "[DiaWeaverAlign] Loaded " << raw_spectra.size()
                    << " pseudo-MS2 spectra from " << experiments.size() << " experiment(s)." << std::endl;

    if (raw_spectra.empty()) return;

    buildIndexCore_(raw_spectra);
  }

  // -------------------------------------------------------------------------
  // binSingleSpectrum_: bin fragment peaks of one spectrum.
  // Used by appendSpectrumToStream; buildIndexCore_ has its own optimised loop.
  // -------------------------------------------------------------------------
  void DiaWeaverAlign::binSingleSpectrum_(
    const MSSpectrum&                                spec,
    uint32_t                                         spec_id,
    uint8_t                                          charge_byte,
    bool                                             im_detected,
    std::vector<std::pair<uint32_t, FragmentEntry>>& out_entries) const
  {
    if (im_detected)
    {
      const auto& im_array = spec.getFloatDataArrays()[spec.getIMData().first];

      std::vector<float>         im_best_intensity(n_im_bins_, -1.0f);
      std::vector<FragmentEntry> im_best_fe(n_im_bins_);
      std::vector<uint32_t>      occupied_im_bins;
      occupied_im_bins.reserve(32);

      uint32_t current_mz_bin = std::numeric_limits<uint32_t>::max();

      auto flush_mz_bin = [&]()
      {
        for (uint32_t im_b : occupied_im_bins)
        {
          out_entries.emplace_back(current_mz_bin * n_im_bins_ + im_b, im_best_fe[im_b]);
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

        if (im_best_intensity[im_bin] < 0.0f)
        {
          occupied_im_bins.push_back(im_bin);
          im_best_intensity[im_bin] = intensity;
          im_best_fe[im_bin] = {spec_id, toMassOffset_(mz, mz_bin), charge_byte, 0};
        }
        else if (intensity > im_best_intensity[im_bin])
        {
          im_best_intensity[im_bin]     = intensity;
          im_best_fe[im_bin].mass_offset = toMassOffset_(mz, mz_bin);
        }
      }

      if (current_mz_bin != std::numeric_limits<uint32_t>::max())
        flush_mz_bin();
    }
    else
    {
      uint32_t      current_bin    = std::numeric_limits<uint32_t>::max();
      float         best_intensity = 0.0f;
      FragmentEntry best_fe{};

      auto emit = [&]() {
        if (current_bin != std::numeric_limits<uint32_t>::max())
          out_entries.emplace_back(current_bin, best_fe);
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

  // -------------------------------------------------------------------------
  // buildCSR_: sort entries and build bin_offsets_ / fragment_entries_.
  // Called by buildIndexCore_ (phase 4) and finalizeFromStream.
  // -------------------------------------------------------------------------
  void DiaWeaverAlign::buildCSR_(
    std::vector<std::pair<uint32_t, FragmentEntry>>& all_entries,
    bool                                              im_detected)
  {
    OPENMS_LOG_DEBUG << "[DiaWeaverAlign::buildCSR_] sorting " << all_entries.size()
                     << " entries, im=" << im_detected << std::endl;

    std::sort(all_entries.begin(), all_entries.end(),
      [](const std::pair<uint32_t, FragmentEntry>& a,
         const std::pair<uint32_t, FragmentEntry>& b) {
        if (a.first != b.first) return a.first < b.first;
        return a.second.spectrum_id < b.second.spectrum_id;
      });

    const Size n_flat_bins = im_detected
      ? static_cast<Size>(n_bins_) * n_im_bins_
      : n_bins_;

    bin_offsets_.assign(n_flat_bins + 1, 0);
    for (const auto& [flat_idx, _] : all_entries)
      ++bin_offsets_[flat_idx + 1];
    for (Size b = 1; b <= n_flat_bins; ++b)
      bin_offsets_[b] += bin_offsets_[b - 1];

    // Count occupied bins for debug summary
    Size occupied = 0;
    for (Size b = 0; b < n_flat_bins; ++b)
      if (bin_offsets_[b + 1] > bin_offsets_[b]) ++occupied;

    OPENMS_LOG_DEBUG << "[DiaWeaverAlign::buildCSR_] CSR complete: "
                     << all_entries.size() << " entries across "
                     << occupied << "/" << n_flat_bins << " occupied bins" << std::endl;

    fragment_entries_.resize(all_entries.size());
    for (Size i = 0; i < all_entries.size(); ++i)
      fragment_entries_[i] = all_entries[i].second;
  }

  // -----------------------------------------------------------------------
  // Phase 2-4: Shared index construction.
  // Assigns spectrum IDs, bins fragments, builds CSR.
  // Called by both buildIndex overloads after they populate raw_spectra.
  // -----------------------------------------------------------------------
  void DiaWeaverAlign::buildIndexCore_(std::vector<RawSpecEntry_>& raw_spectra)
  {
    // Determine IM presence from the loaded data. hasIM() cannot be used here
    // because bin_offsets_ is still empty at this point in construction.
    const bool im_detected =
      !raw_spectra.empty() &&
      raw_spectra.front().spectrum.getDriftTime() != IMTypes::DRIFTTIME_NOT_SET;

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

    if (im_detected)
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
    // Phase 4: Sort and build CSR.
    // -----------------------------------------------------------------------
    buildCSR_(all_entries, im_detected);

    OPENMS_LOG_INFO << "[DiaWeaverAlign] Index complete: "
                    << spectrum_entries_.size() << " spectra, "
                    << fragment_entries_.size() << " fragment entries." << std::endl;
  }

  // -------------------------------------------------------------------------
  // Streaming index construction
  // -------------------------------------------------------------------------

  void DiaWeaverAlign::openStream(const String& frags_path, const String& meta_path)
  {
    updateMembers_();
    OPENMS_LOG_DEBUG << "[DiaWeaverAlign::openStream] frags=" << frags_path
                     << "  meta=" << meta_path << std::endl;
    OPENMS_LOG_DEBUG << "[DiaWeaverAlign::openStream] fragment axis: mz=["
                     << lower_mz_ << ", " << upper_mz_ << "] Da  bin_width=" << bin_width_
                     << "  n_mz_bins=" << n_bins_ << std::endl;
    OPENMS_LOG_DEBUG << "[DiaWeaverAlign::openStream] IM axis: im=["
                     << lower_im_ << ", " << upper_im_ << "]  bin_width_im=" << bin_width_im_
                     << "  n_im_bins=" << n_im_bins_ << std::endl;

    // Reset index state
    source_files_.clear();
    spectrum_entries_.clear();
    fragment_entries_.clear();
    bin_offsets_.clear();
    precursor_entries_.clear();
    precursor_offsets_.clear();

    // Reset streaming counters
    stream_next_spec_id_   = 0;
    stream_im_detected_    = false;
    stream_im_initialised_ = false;

    // Open raw fragment file
    stream_frags_fp_ = std::fopen(frags_path.c_str(), "wb");
    if (!stream_frags_fp_)
      throw Exception::FileNotFound(
        __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, frags_path);

    // Write frags header: magic (4) + n_spectra placeholder (4) + has_im placeholder (4)
    const char frags_magic[4] = {'D', 'W', 'A', 'F'};
    uint32_t   placeholder    = 0;
    std::fwrite(frags_magic, 1, 4, stream_frags_fp_);
    std::fwrite(&placeholder, 4, 1, stream_frags_fp_);  // n_spectra — filled by closeStream
    std::fwrite(&placeholder, 4, 1, stream_frags_fp_);  // has_im    — filled by closeStream

    // Open raw metadata file
    stream_meta_fp_ = std::fopen(meta_path.c_str(), "wb");
    if (!stream_meta_fp_)
    {
      std::fclose(stream_frags_fp_);
      stream_frags_fp_ = nullptr;
      throw Exception::FileNotFound(
        __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, meta_path);
    }

    // Write meta header: magic (4) + n_spectra placeholder (4)
    const char meta_magic[4] = {'D', 'W', 'A', 'M'};
    std::fwrite(meta_magic, 1, 4, stream_meta_fp_);
    std::fwrite(&placeholder, 4, 1, stream_meta_fp_);   // n_spectra — filled by closeStream
  }

  void DiaWeaverAlign::appendSpectrumToStream(const MSSpectrum& spec, Size source_file_idx)
  {
    if (spec.getMSLevel() != 2 || spec.empty() || spec.getPrecursors().empty()) return;

    if (!stream_frags_fp_ || !stream_meta_fp_)
      throw Exception::MissingInformation(
        __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "appendSpectrumToStream() called without a prior openStream().");

    // Detect IM on first spectrum seen
    if (!stream_im_initialised_)
    {
      stream_im_detected_    = (spec.getDriftTime() != IMTypes::DRIFTTIME_NOT_SET);
      stream_im_initialised_ = true;
      OPENMS_LOG_INFO << "[DiaWeaverAlign] Streaming IM: "
                      << (stream_im_detected_ ? "detected" : "not detected") << std::endl;
    }

    const uint32_t spec_id     = stream_next_spec_id_++;
    const uint8_t  charge_byte = static_cast<uint8_t>(spec.getPrecursors()[0].getCharge());

    // Bin fragment peaks and write records: uint32_t flat_idx (4) + FragmentEntry (8) = 12 bytes each
    std::vector<std::pair<uint32_t, FragmentEntry>> entries;
    binSingleSpectrum_(spec, spec_id, charge_byte, stream_im_detected_, entries);

    OPENMS_LOG_DEBUG << "[DiaWeaverAlign::appendSpectrumToStream] spec_id=" << spec_id
                     << "  native_id='" << spec.getNativeID() << "'"
                     << "  RT=" << spec.getRT()
                     << "  prec_mz=" << spec.getPrecursors()[0].getMZ()
                     << (stream_im_detected_
                           ? ("  drift_t=" + String(spec.getDriftTime()))
                           : "")
                     << "  n_input_peaks=" << spec.size()
                     << "  n_indexed_entries=" << entries.size() << std::endl;

    for (const auto& [flat_idx, fe] : entries)
    {
      std::fwrite(&flat_idx, 4, 1, stream_frags_fp_);
      std::fwrite(&fe,       8, 1, stream_frags_fp_);
    }

    // Write metadata record to meta file
    const double   rt      = spec.getRT();
    const double   prec_mz = spec.getPrecursors()[0].getMZ();
    const double   dt      = stream_im_detected_ ? spec.getDriftTime() : -1.0;
    const int32_t  charge  = static_cast<int32_t>(spec.getPrecursors()[0].getCharge());
    const uint32_t src     = static_cast<uint32_t>(source_file_idx);

    std::fwrite(&rt,      8, 1, stream_meta_fp_);
    std::fwrite(&prec_mz, 8, 1, stream_meta_fp_);
    std::fwrite(&dt,      8, 1, stream_meta_fp_);
    std::fwrite(&charge,  4, 1, stream_meta_fp_);
    std::fwrite(&src,     4, 1, stream_meta_fp_);

    const String&  native_id = spec.getNativeID();
    const uint32_t id_len    = static_cast<uint32_t>(native_id.size());
    std::fwrite(&id_len, 4, 1, stream_meta_fp_);
    if (id_len > 0)
      std::fwrite(native_id.c_str(), 1, id_len, stream_meta_fp_);
  }

  void DiaWeaverAlign::closeStream()
  {
    if (!stream_frags_fp_ && !stream_meta_fp_) return;

    const uint32_t n_spectra = stream_next_spec_id_;
    const uint32_t has_im    = stream_im_detected_ ? 1u : 0u;

    if (stream_frags_fp_)
    {
      // Seek back to header offsets 4..11 and write actual n_spectra and has_im
      std::fseek(stream_frags_fp_, 4, SEEK_SET);
      std::fwrite(&n_spectra, 4, 1, stream_frags_fp_);
      std::fwrite(&has_im,    4, 1, stream_frags_fp_);
      std::fclose(stream_frags_fp_);
      stream_frags_fp_ = nullptr;
    }

    if (stream_meta_fp_)
    {
      // Seek back to header offset 4 and write actual n_spectra
      std::fseek(stream_meta_fp_, 4, SEEK_SET);
      std::fwrite(&n_spectra, 4, 1, stream_meta_fp_);
      std::fclose(stream_meta_fp_);
      stream_meta_fp_ = nullptr;
    }

    OPENMS_LOG_DEBUG << "[DiaWeaverAlign::closeStream] closed  n_spectra=" << n_spectra
                     << "  has_im=" << has_im << std::endl;

    stream_next_spec_id_   = 0;
    stream_im_initialised_ = false;
    stream_im_detected_    = false;
  }

  void DiaWeaverAlign::finalizeFromStream(
    const String&              frags_path,
    const String&              meta_path,
    const std::vector<String>& source_file_names)
  {
    source_files_ = source_file_names;

    // -----------------------------------------------------------------------
    // Read raw fragment entries file
    // -----------------------------------------------------------------------
    FILE* fp = std::fopen(frags_path.c_str(), "rb");
    if (!fp)
      throw Exception::FileNotFound(
        __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, frags_path);

    char     magic_f[4];
    uint32_t n_spectra_f, has_im_flag;
    std::fread(magic_f, 1, 4, fp);
    if (magic_f[0] != 'D' || magic_f[1] != 'W' || magic_f[2] != 'A' || magic_f[3] != 'F')
    {
      std::fclose(fp);
      throw Exception::ParseError(
        __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        frags_path, "Bad magic number; expected DWAF.");
    }
    std::fread(&n_spectra_f, 4, 1, fp);
    std::fread(&has_im_flag, 4, 1, fp);

    const bool im_detected = (has_im_flag != 0);

    // Derive record count from file size (each record is 12 bytes: 4+8)
    std::fseek(fp, 0, SEEK_END);
    const long file_size = std::ftell(fp);
    std::fseek(fp, 12, SEEK_SET);   // rewind to first record

    const uint64_t n_records = (file_size > 12)
      ? static_cast<uint64_t>(file_size - 12) / 12
      : 0;

    OPENMS_LOG_INFO << "[DiaWeaverAlign] Finalizing: reading " << n_records
                    << " fragment records, im=" << im_detected
                    << ", from " << frags_path << std::endl;

    std::vector<std::pair<uint32_t, FragmentEntry>> all_entries(n_records);
    for (uint64_t i = 0; i < n_records; ++i)
    {
      uint32_t     flat_idx;
      FragmentEntry fe;
      std::fread(&flat_idx, 4, 1, fp);
      std::fread(&fe,       8, 1, fp);
      all_entries[i] = {flat_idx, fe};
    }
    std::fclose(fp);

    // -----------------------------------------------------------------------
    // Read spectrum metadata file
    // -----------------------------------------------------------------------
    FILE* fpm = std::fopen(meta_path.c_str(), "rb");
    if (!fpm)
      throw Exception::FileNotFound(
        __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, meta_path);

    char     magic_m[4];
    uint32_t n_spectra_m;
    std::fread(magic_m, 1, 4, fpm);
    if (magic_m[0] != 'D' || magic_m[1] != 'W' || magic_m[2] != 'A' || magic_m[3] != 'M')
    {
      std::fclose(fpm);
      throw Exception::ParseError(
        __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        meta_path, "Bad magic number; expected DWAM.");
    }
    std::fread(&n_spectra_m, 4, 1, fpm);

    OPENMS_LOG_INFO << "[DiaWeaverAlign] Reading " << n_spectra_m
                    << " spectrum metadata entries." << std::endl;

    spectrum_entries_.resize(n_spectra_m);
    for (uint32_t i = 0; i < n_spectra_m; ++i)
    {
      SpectrumEntry& se = spectrum_entries_[i];
      int32_t  charge;
      uint32_t src_idx, id_len;

      std::fread(&se.retention_time, 8, 1, fpm);
      std::fread(&se.precursor_mz,   8, 1, fpm);
      std::fread(&se.drift_time,     8, 1, fpm);
      std::fread(&charge,            4, 1, fpm);
      std::fread(&src_idx,           4, 1, fpm);
      std::fread(&id_len,            4, 1, fpm);

      se.precursor_charge = static_cast<int>(charge);
      se.source_file_idx  = static_cast<Size>(src_idx);

      if (id_len > 0)
      {
        std::string buf(id_len, '\0');
        std::fread(&buf[0], 1, id_len, fpm);
        se.native_id = buf;
      }

      OPENMS_LOG_DEBUG << "[DiaWeaverAlign::finalizeFromStream] spectrum_id=" << i
                       << "  RT=" << se.retention_time
                       << "  prec_mz=" << se.precursor_mz
                       << "  drift_time=" << se.drift_time
                       << "  charge=" << se.precursor_charge
                       << "  src_file=" << se.source_file_idx
                       << "  native_id='" << se.native_id << "'" << std::endl;
    }
    std::fclose(fpm);

    // -----------------------------------------------------------------------
    // Sort entries and build CSR fragment index
    // -----------------------------------------------------------------------
    OPENMS_LOG_INFO << "[DiaWeaverAlign] Building CSR from "
                    << all_entries.size() << " entries..." << std::endl;

    buildCSR_(all_entries, im_detected);

    OPENMS_LOG_INFO << "[DiaWeaverAlign] Finalization complete: "
                    << spectrum_entries_.size() << " spectra, "
                    << fragment_entries_.size() << " fragment entries." << std::endl;
  }

  // -------------------------------------------------------------------------
  // Index serialization / deserialization
  //
  // File format (.dwaindex):
  //   magic[4]  "DWAI"
  //   doubles/uint32_ts for all bin-axis scalars + has_im flag
  //   source_files_   (count + length-prefixed strings)
  //   spectrum_entries_  (count + per-entry fixed+variable data)
  //   fragment_entries_  (count uint64 + raw array)
  //   bin_offsets_       (count uint64 + raw array)
  //   precursor_entries_ (count uint64 + raw array)
  //   precursor_offsets_ (count uint64 + raw array)
  // -------------------------------------------------------------------------

  void DiaWeaverAlign::saveIndex(const String& path) const
  {
    FILE* fp = std::fopen(path.c_str(), "wb");
    if (!fp)
      throw Exception::FileNotFound(
        __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path);

    // Magic
    const char magic[4] = {'D', 'W', 'A', 'I'};
    std::fwrite(magic, 1, 4, fp);

    // Bin-axis scalars
    std::fwrite(&lower_mz_,            8, 1, fp);
    std::fwrite(&upper_mz_,            8, 1, fp);
    std::fwrite(&bin_width_,           8, 1, fp);
    std::fwrite(&n_bins_,              4, 1, fp);
    std::fwrite(&lower_im_,            8, 1, fp);
    std::fwrite(&upper_im_,            8, 1, fp);
    std::fwrite(&bin_width_im_,        8, 1, fp);
    std::fwrite(&n_im_bins_,           4, 1, fp);
    std::fwrite(&lower_precursor_mz_,  8, 1, fp);
    std::fwrite(&upper_precursor_mz_,  8, 1, fp);
    std::fwrite(&precursor_bin_width_, 8, 1, fp);
    std::fwrite(&n_precursor_bins_,    4, 1, fp);
    const uint32_t has_im = hasIM() ? 1u : 0u;
    std::fwrite(&has_im, 4, 1, fp);

    // DIA window bounds — makes the file self-describing about which precursor
    // isolation window this index was built for (independent of fragment m/z range)
    std::fwrite(&window_lower_mz_, 8, 1, fp);
    std::fwrite(&window_upper_mz_, 8, 1, fp);
    std::fwrite(&window_lower_im_, 8, 1, fp);
    std::fwrite(&window_upper_im_, 8, 1, fp);

    // Source files
    const uint32_t n_src = static_cast<uint32_t>(source_files_.size());
    std::fwrite(&n_src, 4, 1, fp);
    for (const String& sf : source_files_)
    {
      const uint32_t len = static_cast<uint32_t>(sf.size());
      std::fwrite(&len,        4,   1,   fp);
      std::fwrite(sf.c_str(), 1,   len, fp);
    }

    // Spectrum entries
    const uint32_t n_spec = static_cast<uint32_t>(spectrum_entries_.size());
    std::fwrite(&n_spec, 4, 1, fp);
    for (const SpectrumEntry& se : spectrum_entries_)
    {
      const int32_t  charge  = static_cast<int32_t>(se.precursor_charge);
      const uint32_t src_idx = static_cast<uint32_t>(se.source_file_idx);
      const uint32_t id_len  = static_cast<uint32_t>(se.native_id.size());
      std::fwrite(&se.retention_time, 8, 1, fp);
      std::fwrite(&se.precursor_mz,   8, 1, fp);
      std::fwrite(&se.drift_time,     8, 1, fp);
      std::fwrite(&charge,            4, 1, fp);
      std::fwrite(&src_idx,           4, 1, fp);
      std::fwrite(&id_len,            4, 1, fp);
      if (id_len > 0)
        std::fwrite(se.native_id.c_str(), 1, id_len, fp);
    }

    // Fragment entries (8 bytes each)
    const uint64_t n_frag = static_cast<uint64_t>(fragment_entries_.size());
    std::fwrite(&n_frag, 8, 1, fp);
    if (n_frag > 0)
      std::fwrite(fragment_entries_.data(), 8, n_frag, fp);

    // Bin offsets (4 bytes each)
    const uint64_t n_bin_off = static_cast<uint64_t>(bin_offsets_.size());
    std::fwrite(&n_bin_off, 8, 1, fp);
    if (n_bin_off > 0)
      std::fwrite(bin_offsets_.data(), 4, n_bin_off, fp);

    // Precursor entries (4 bytes each)
    const uint64_t n_prec = static_cast<uint64_t>(precursor_entries_.size());
    std::fwrite(&n_prec, 8, 1, fp);
    if (n_prec > 0)
      std::fwrite(precursor_entries_.data(), 4, n_prec, fp);

    // Precursor offsets (4 bytes each)
    const uint64_t n_prec_off = static_cast<uint64_t>(precursor_offsets_.size());
    std::fwrite(&n_prec_off, 8, 1, fp);
    if (n_prec_off > 0)
      std::fwrite(precursor_offsets_.data(), 4, n_prec_off, fp);

    std::fclose(fp);
    OPENMS_LOG_INFO << "[DiaWeaverAlign] Index saved to " << path << std::endl;
  }

  void DiaWeaverAlign::loadIndex(const String& path)
  {
    FILE* fp = std::fopen(path.c_str(), "rb");
    if (!fp)
      throw Exception::FileNotFound(
        __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path);

    // Magic
    char magic[4];
    std::fread(magic, 1, 4, fp);
    if (magic[0] != 'D' || magic[1] != 'W' || magic[2] != 'A' || magic[3] != 'I')
    {
      std::fclose(fp);
      throw Exception::ParseError(
        __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        path, "Bad magic number; expected DWAI.");
    }

    // Bin-axis scalars — restore directly (bypasses updateMembers_)
    std::fread(&lower_mz_,            8, 1, fp);
    std::fread(&upper_mz_,            8, 1, fp);
    std::fread(&bin_width_,           8, 1, fp);
    std::fread(&n_bins_,              4, 1, fp);
    std::fread(&lower_im_,            8, 1, fp);
    std::fread(&upper_im_,            8, 1, fp);
    std::fread(&bin_width_im_,        8, 1, fp);
    std::fread(&n_im_bins_,           4, 1, fp);
    std::fread(&lower_precursor_mz_,  8, 1, fp);
    std::fread(&upper_precursor_mz_,  8, 1, fp);
    std::fread(&precursor_bin_width_, 8, 1, fp);
    std::fread(&n_precursor_bins_,    4, 1, fp);
    uint32_t has_im_flag;
    std::fread(&has_im_flag, 4, 1, fp);

    // DIA window bounds
    std::fread(&window_lower_mz_, 8, 1, fp);
    std::fread(&window_upper_mz_, 8, 1, fp);
    std::fread(&window_lower_im_, 8, 1, fp);
    std::fread(&window_upper_im_, 8, 1, fp);

    // Source files
    uint32_t n_src;
    std::fread(&n_src, 4, 1, fp);
    source_files_.resize(n_src);
    for (uint32_t i = 0; i < n_src; ++i)
    {
      uint32_t len;
      std::fread(&len, 4, 1, fp);
      if (len > 0)
      {
        std::string buf(len, '\0');
        std::fread(&buf[0], 1, len, fp);
        source_files_[i] = buf;
      }
    }

    // Spectrum entries
    uint32_t n_spec;
    std::fread(&n_spec, 4, 1, fp);
    spectrum_entries_.resize(n_spec);
    for (uint32_t i = 0; i < n_spec; ++i)
    {
      SpectrumEntry& se = spectrum_entries_[i];
      int32_t  charge;
      uint32_t src_idx, id_len;
      std::fread(&se.retention_time, 8, 1, fp);
      std::fread(&se.precursor_mz,   8, 1, fp);
      std::fread(&se.drift_time,     8, 1, fp);
      std::fread(&charge,            4, 1, fp);
      std::fread(&src_idx,           4, 1, fp);
      std::fread(&id_len,            4, 1, fp);
      se.precursor_charge = static_cast<int>(charge);
      se.source_file_idx  = static_cast<Size>(src_idx);
      if (id_len > 0)
      {
        std::string buf(id_len, '\0');
        std::fread(&buf[0], 1, id_len, fp);
        se.native_id = buf;
      }
    }

    // Fragment entries
    uint64_t n_frag;
    std::fread(&n_frag, 8, 1, fp);
    fragment_entries_.resize(n_frag);
    if (n_frag > 0)
      std::fread(fragment_entries_.data(), 8, n_frag, fp);

    // Bin offsets
    uint64_t n_bin_off;
    std::fread(&n_bin_off, 8, 1, fp);
    bin_offsets_.resize(n_bin_off);
    if (n_bin_off > 0)
      std::fread(bin_offsets_.data(), 4, n_bin_off, fp);

    // Precursor entries
    uint64_t n_prec;
    std::fread(&n_prec, 8, 1, fp);
    precursor_entries_.resize(n_prec);
    if (n_prec > 0)
      std::fread(precursor_entries_.data(), 4, n_prec, fp);

    // Precursor offsets
    uint64_t n_prec_off;
    std::fread(&n_prec_off, 8, 1, fp);
    precursor_offsets_.resize(n_prec_off);
    if (n_prec_off > 0)
      std::fread(precursor_offsets_.data(), 4, n_prec_off, fp);

    std::fclose(fp);
    OPENMS_LOG_INFO << "[DiaWeaverAlign] Index loaded from " << path
                    << "  spectra=" << spectrum_entries_.size()
                    << "  frags=" << fragment_entries_.size() << std::endl;
  }

  // -------------------------------------------------------------------------
  // DIA window bounds
  // -------------------------------------------------------------------------

  void DiaWeaverAlign::setWindowBounds(
    double lower_mz, double upper_mz,
    double lower_im, double upper_im)
  {
    window_lower_mz_ = lower_mz;
    window_upper_mz_ = upper_mz;
    window_lower_im_ = lower_im;
    window_upper_im_ = upper_im;
  }

  double DiaWeaverAlign::getWindowLowerMz() const { return window_lower_mz_; }
  double DiaWeaverAlign::getWindowUpperMz() const { return window_upper_mz_; }
  double DiaWeaverAlign::getWindowLowerIM() const { return window_lower_im_; }
  double DiaWeaverAlign::getWindowUpperIM() const { return window_upper_im_; }

  // -------------------------------------------------------------------------
  // Spectrum matching
  // -------------------------------------------------------------------------

  std::vector<DiaWeaverAlign::MatchResult>
  DiaWeaverAlign::matchSpectrum(const MSSpectrum& query) const
  {
    // -----------------------------------------------------------------------
    // Step 1: Collect raw spectrum_id hits.
    //
    // For every query peak that falls within the index bounds, look up its
    // (mz_bin, im_bin) cell and append every spectrum_id found there.
    // The same spectrum_id is appended once per query peak that hits it —
    // that repetition is intentional: it allows us to count matched peaks per spectrum.
    // -----------------------------------------------------------------------
    
    const bool query_has_im = query.containsIMData();

    if (query_has_im != hasIM())
    {
      throw Exception::MissingInformation(
        __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Ion mobility mismatch: the fragment index was built " +
        String(hasIM() ? "with" : "without") +
        " IM data but the query spectrum '" + query.getNativeID() +
        "' " + (query_has_im ? "has" : "has no") + " per-peak ion mobility.");
    }

    std::vector<uint32_t> hits;

    Size im_array_idx = 0;
    if (query_has_im)
      im_array_idx = query.getIMData().first;

    for (Size i = 0; i < query.size(); ++i)
    {
      const double mz = query[i].getMZ();
      if (mz < lower_mz_ || mz >= upper_mz_) continue;

      const uint32_t mz_bin = toBinIdx_(mz);

      if (query_has_im)
      {
        const float im = query.getFloatDataArrays()[im_array_idx][i];
        if (im < lower_im_ || im >= upper_im_) continue;

        const uint32_t im_bin = toBinIdx_im_(im);
        auto [begin, end] = getBinEntries(mz_bin, im_bin);
        for (const FragmentEntry* fe = begin; fe != end; ++fe)
          hits.push_back(fe->spectrum_id);
      }
      else
      {
        auto [begin, end] = getBinEntries(mz_bin);
        for (const FragmentEntry* fe = begin; fe != end; ++fe)
          hits.push_back(fe->spectrum_id);
      }
    }

    // -----------------------------------------------------------------------
    // Step 2: Sort hits so identical spectrum_ids are adjacent.
    // -----------------------------------------------------------------------

    std::sort(hits.begin(), hits.end());

    // -----------------------------------------------------------------------
    // Step 3: Single linear scan — count runs of the same spectrum_id.
    // -----------------------------------------------------------------------

    std::vector<MatchResult> results;

    Size i = 0;
    while (i < hits.size())
    {
      const uint32_t spec_id = hits[i];
      uint32_t count = 0;
      while (i < hits.size() && hits[i] == spec_id)
      {
        ++count;
        ++i;
      }
      results.push_back({spec_id, count});
      OPENMS_LOG_DEBUG << "[DiaWeaverAlign::matchSpectrum] query='" << query.getNativeID()
                       << "'  matched spectrum_id=" << spec_id
                       << "  matched_peaks=" << count << std::endl;
    }

    OPENMS_LOG_DEBUG << "[DiaWeaverAlign::matchSpectrum] query='" << query.getNativeID()
                     << "'  n_query_peaks=" << query.size()
                     << "  total_raw_hits=" << hits.size()
                     << "  distinct_spectra_hit=" << results.size() << std::endl;

    return results;
  }

  std::vector<DiaWeaverAlign::ScoreTrace>
  DiaWeaverAlign::matchExperiment(const MSExperiment& experiment, uint32_t min_matched_peaks) const
  {
    // Collect pointers to query spectra in load order (ascending RT).
    //
    // Peak-picked DIA MS2 data produced by diaWeaver's PeakPickerIM pipeline
    // carries MS level 1 (the DIA extraction step sets it so that downstream
    // mass-trace tools see "MS1" data). This is not clean but it works for now.
    // every peak-picked MS2 spectrum has a
    // precursor (isolation window)
    std::vector<const MSSpectrum*> query_spectra;
    for (const auto& spec : experiment)
    {
      if (!spec.empty() && !spec.getPrecursors().empty())
        query_spectra.push_back(&spec);
    }

    if (query_spectra.empty()) return {};

    OPENMS_LOG_DEBUG << "[DiaWeaverAlign::matchExperiment] "
                     << query_spectra.size() << " query spectra  "
                     << "min_matched_peaks=" << min_matched_peaks
                     << "  index_has_im=" << hasIM() << std::endl;

    // Validate IM consistency once before the parallel region.
    // matchSpectrum() would also catch this per call, but exceptions thrown
    // inside an OpenMP parallel region do not propagate reliably to the caller.
    const bool query_has_im = query_spectra.front()->containsIMData();
    if (query_has_im != hasIM())
    {
      throw Exception::MissingInformation(
        __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Ion mobility mismatch: the fragment index was built " +
        String(hasIM() ? "with" : "without") +
        " IM data but the query experiment spectra " +
        (query_has_im ? "have" : "have no") + " per-peak ion mobility.");
    }

    // -----------------------------------------------------------------------
    // First pass: parallel fragment scoring.
    //
    // Each thread accumulates into its own private map to avoid contention.
    // LocalEntry also records the query spectrum index alongside rt and score
    // so the apex query spectrum can be identified during merge without a
    // separate RT-to-index lookup.
    // -----------------------------------------------------------------------

    struct LocalEntry
    {
      std::vector<double>   rts;
      std::vector<uint32_t> scores;
      std::vector<int>      query_indices; ///< index into query_spectra parallel to rts/scores
    };
    using LocalTraceMap = std::unordered_map<uint32_t, LocalEntry>;
    std::vector<LocalTraceMap> thread_traces;

    auto process_query = [&](int q, LocalTraceMap& local_map)
    {
      const double rt = query_spectra[q]->getRT();
      for (const MatchResult& r : matchSpectrum(*query_spectra[q]))
      {
        if (r.matched_peaks < min_matched_peaks) continue;
        auto& entry = local_map[r.spectrum_id];
        entry.rts.push_back(rt);
        entry.scores.push_back(r.matched_peaks);
        entry.query_indices.push_back(q);
      }
    };

#ifdef _OPENMP
    #pragma omp parallel
    {
      #pragma omp single
      thread_traces.resize(static_cast<Size>(omp_get_num_threads()));

      auto& local = thread_traces[omp_get_thread_num()];

      #pragma omp for schedule(static)
      for (int q = 0; q < static_cast<int>(query_spectra.size()); ++q)
        process_query(q, local);
    }
#else
    thread_traces.resize(1);
    for (int q = 0; q < static_cast<int>(query_spectra.size()); ++q)
      process_query(q, thread_traces[0]);
#endif

    // -----------------------------------------------------------------------
    // Merge per-thread traces and track apex (best score) for each spectrum_id.
    //
    // Threads processed contiguous RT chunks in ascending order, so appending
    // in thread-index order preserves RT ordering within each trace.
    // -----------------------------------------------------------------------

    struct ApexInfo
    {
      uint32_t score{0};
      double   rt{-1.0};
      int      query_idx{-1};
    };

    std::unordered_map<uint32_t, ScoreTrace> merged;
    std::unordered_map<uint32_t, ApexInfo>   apex_info;

    for (Size t = 0; t < thread_traces.size(); ++t)
    {
      for (auto& [sid, local_entry] : thread_traces[t])
      {
        ScoreTrace& trace = merged[sid];
        trace.spectrum_id = sid;
        trace.rts.insert(trace.rts.end(),
                         local_entry.rts.begin(), local_entry.rts.end());
        trace.scores.insert(trace.scores.end(),
                            local_entry.scores.begin(), local_entry.scores.end());

        ApexInfo& best = apex_info[sid];
        for (Size j = 0; j < local_entry.scores.size(); ++j)
        {
          if (local_entry.scores[j] > best.score)
          {
            best.score     = local_entry.scores[j];
            best.rt        = local_entry.rts[j];
            best.query_idx = local_entry.query_indices[j];
          }
        }
      }
    }

    std::vector<ScoreTrace> result;
    result.reserve(merged.size());
    for (auto& [sid, trace] : merged)
      result.push_back(std::move(trace));

    std::sort(result.begin(), result.end(),
      [](const ScoreTrace& a, const ScoreTrace& b)
      { return a.spectrum_id < b.spectrum_id; });

    OPENMS_LOG_DEBUG << "[DiaWeaverAlign::matchExperiment] merged "
                     << result.size() << " spectrum_ids with at least one above-threshold query hit"
                     << std::endl;

    // Populate apex fields and group spectrum_ids by their apex query spectrum index.
    std::vector<std::vector<uint32_t>> apex_groups(query_spectra.size());
    for (ScoreTrace& trace : result)
    {
      const ApexInfo& best = apex_info[trace.spectrum_id];
      trace.apex_score = best.score;
      trace.apex_rt    = best.rt;
      apex_groups[best.query_idx].push_back(trace.spectrum_id);

      OPENMS_LOG_DEBUG << "[DiaWeaverAlign::matchExperiment] spectrum_id=" << trace.spectrum_id
                       << "  apex_rt=" << trace.apex_rt
                       << "  apex_score=" << trace.apex_score
                       << "  n_query_hits=" << trace.rts.size() << std::endl;
    }

    // -----------------------------------------------------------------------
    // Second pass: collect apex fingerprints.
    //
    // For each query spectrum that is the apex for at least one spectrum_id,
    // re-scan its peaks and record (flat_bin_idx, experimental_intensity) for
    // every fragment index hit that belongs to a target spectrum_id.
    //
    // is_target uses uint8_t (not bool) to avoid std::vector<bool> bit-packing,
    // which would prevent safe per-element writes from the inner loop.
    // spec_to_result_idx maps spectrum_id -> index in result for direct access.
    // Both arrays are O(N_index) and reset via target_list to avoid O(N_index) scans.
    // -----------------------------------------------------------------------

    std::vector<int>     spec_to_result_idx(spectrum_entries_.size(), -1);
    std::vector<uint8_t> is_target(spectrum_entries_.size(), 0);
    std::vector<uint32_t> target_list;

    for (int r = 0; r < static_cast<int>(result.size()); ++r)
      spec_to_result_idx[result[r].spectrum_id] = r;

    for (int q = 0; q < static_cast<int>(query_spectra.size()); ++q)
    {
      if (apex_groups[q].empty()) continue;

      for (uint32_t sid : apex_groups[q])
      {
        is_target[sid] = 1;
        target_list.push_back(sid);
      }

      const MSSpectrum& apex_spec = *query_spectra[q];
      Size im_array_idx = 0;
      if (query_has_im)
        im_array_idx = apex_spec.getIMData().first;

      for (Size i = 0; i < apex_spec.size(); ++i)
      {
        const double mz = apex_spec[i].getMZ();
        if (mz < lower_mz_ || mz >= upper_mz_) continue;
        const float    intensity = apex_spec[i].getIntensity();
        const uint32_t mz_bin   = toBinIdx_(mz);

        uint32_t im_bin = 0;
        if (query_has_im)
        {
          const float im = apex_spec.getFloatDataArrays()[im_array_idx][i];
          if (im < lower_im_ || im >= upper_im_) continue;
          im_bin = toBinIdx_im_(im);
        }

        const uint32_t flat_idx = query_has_im
          ? mz_bin * n_im_bins_ + im_bin
          : mz_bin;

        auto [begin, end] = getBinEntries(mz_bin, im_bin);
        for (const FragmentEntry* fe = begin; fe != end; ++fe)
        {
          if (is_target[fe->spectrum_id])
            result[spec_to_result_idx[fe->spectrum_id]].apex_fingerprint.push_back({flat_idx, intensity});
        }
      }

      for (uint32_t sid : target_list) is_target[sid] = 0;
      target_list.clear();
    }

    for (const ScoreTrace& trace : result)
    {
      OPENMS_LOG_DEBUG << "[DiaWeaverAlign::matchExperiment] apex_fingerprint  spectrum_id="
                       << trace.spectrum_id
                       << "  n_apex_fragment_hits=" << trace.apex_fingerprint.size()
                       << std::endl;
    }

    return result;
  }

  std::vector<DiaWeaverAlign::FeatureGroup>
  DiaWeaverAlign::groupFeatures(const std::vector<ScoreTrace>& traces,
                                 double   rt_tolerance,
                                 uint32_t min_fragment_overlap) const
  {
    if (traces.empty() || precursor_offsets_.empty()) return {};

    // -----------------------------------------------------------------------
    // Pre-computation
    //
    // spec_to_trace_idx: direct-address lookup from spectrum_id to index in traces.
    // fingerprint_bins:  per-trace sorted, deduplicated flat_bin_idx vector for
    //                    O(F_A + F_B) two-pointer overlap counting.
    // -----------------------------------------------------------------------

    std::vector<int> spec_to_trace_idx(spectrum_entries_.size(), -1);
    for (int i = 0; i < static_cast<int>(traces.size()); ++i)
      spec_to_trace_idx[traces[i].spectrum_id] = i;

    std::vector<std::vector<uint32_t>> fingerprint_bins(traces.size());
    for (Size i = 0; i < traces.size(); ++i)
    {
      auto& bins = fingerprint_bins[i];
      bins.reserve(traces[i].apex_fingerprint.size());
      for (const ApexFragment& af : traces[i].apex_fingerprint)
        bins.push_back(af.flat_bin_idx);
      std::sort(bins.begin(), bins.end());
      bins.erase(std::unique(bins.begin(), bins.end()), bins.end());
    }

    // Two-pointer overlap count on sorted flat_bin_idx vectors.
    auto count_overlap = [](const std::vector<uint32_t>& a,
                            const std::vector<uint32_t>& b) -> uint32_t
    {
      uint32_t count = 0;
      Size i = 0, j = 0;
      while (i < a.size() && j < b.size())
      {
        if      (a[i] < b[j]) ++i;
        else if (b[j] < a[i]) ++j;
        else { ++count; ++i; ++j; }
      }
      return count;
    };

    // -----------------------------------------------------------------------
    // Union-Find with path-halving and union by rank.
    // -----------------------------------------------------------------------

    std::vector<int> parent(traces.size());
    std::vector<int> uf_rank(traces.size(), 0);
    std::iota(parent.begin(), parent.end(), 0);

    auto find = [&](int x) -> int
    {
      while (parent[x] != x)
      {
        parent[x] = parent[parent[x]]; // path halving
        x = parent[x];
      }
      return x;
    };

    auto unite = [&](int x, int y)
    {
      x = find(x); y = find(y);
      if (x == y) return;
      if (uf_rank[x] < uf_rank[y]) std::swap(x, y);
      parent[y] = x;
      if (uf_rank[x] == uf_rank[y]) ++uf_rank[x];
    };

    // -----------------------------------------------------------------------
    // Pair validation.
    //
    // For each trace i, query the ±1 precursor bin neighborhood to find
    // candidate spectrum_ids. Pairs are only processed when tj > ti to avoid
    // double-counting. Three criteria must all pass: different source file,
    // apex RT within tolerance, and sufficient fragment overlap.
    // -----------------------------------------------------------------------

    const bool use_im = hasIM();

    for (Size ti = 0; ti < traces.size(); ++ti)
    {
      const SpectrumEntry& se_i = getSpectrumEntry(traces[ti].spectrum_id);
      if (se_i.precursor_mz < lower_precursor_mz_ || se_i.precursor_mz >= upper_precursor_mz_) continue;

      const uint32_t mz_bin_i = toPrecursorBinIdx_(se_i.precursor_mz);
      const uint32_t im_bin_i = use_im ? toBinIdx_im_(se_i.drift_time) : 0;

      const int mz_lo = std::max(0, static_cast<int>(mz_bin_i) - 1);
      const int mz_hi = std::min(static_cast<int>(n_precursor_bins_) - 1, static_cast<int>(mz_bin_i) + 1);
      const int im_lo = use_im ? std::max(0, static_cast<int>(im_bin_i) - 1) : 0;
      const int im_hi = use_im ? std::min(static_cast<int>(n_im_bins_) - 1, static_cast<int>(im_bin_i) + 1) : 0;

      for (int mb = mz_lo; mb <= mz_hi; ++mb)
      {
        for (int ib = im_lo; ib <= im_hi; ++ib)
        {
          auto [begin, end] = getPrecursorBinEntries(
            static_cast<uint32_t>(mb), static_cast<uint32_t>(ib));

          for (const uint32_t* sid_ptr = begin; sid_ptr != end; ++sid_ptr)
          {
            const int tj = spec_to_trace_idx[*sid_ptr];
            if (tj <= static_cast<int>(ti)) continue; // process each pair once only

            const SpectrumEntry& se_j = getSpectrumEntry(*sid_ptr);

            if (se_i.source_file_idx == se_j.source_file_idx)
            {
              OPENMS_LOG_DEBUG << "[DiaWeaverAlign::groupFeatures] REJECT pair ("
                               << traces[ti].spectrum_id << ", " << *sid_ptr
                               << "): same source_file=" << se_i.source_file_idx << std::endl;
              continue;
            }

            const double rt_diff = std::abs(traces[ti].apex_rt - traces[tj].apex_rt);
            if (rt_diff > rt_tolerance)
            {
              OPENMS_LOG_DEBUG << "[DiaWeaverAlign::groupFeatures] REJECT pair ("
                               << traces[ti].spectrum_id << ", " << *sid_ptr
                               << "): |RT_diff|=" << rt_diff
                               << " > rt_tolerance=" << rt_tolerance << std::endl;
              continue;
            }

            const uint32_t overlap = count_overlap(fingerprint_bins[ti], fingerprint_bins[tj]);
            if (overlap < min_fragment_overlap)
            {
              OPENMS_LOG_DEBUG << "[DiaWeaverAlign::groupFeatures] REJECT pair ("
                               << traces[ti].spectrum_id << ", " << *sid_ptr
                               << "): fragment_overlap=" << overlap
                               << " < min=" << min_fragment_overlap << std::endl;
              continue;
            }

            OPENMS_LOG_DEBUG << "[DiaWeaverAlign::groupFeatures] MERGE pair ("
                             << traces[ti].spectrum_id << ", " << *sid_ptr
                             << "): |RT_diff|=" << rt_diff
                             << "  fragment_overlap=" << overlap
                             << "  src_files=(" << se_i.source_file_idx
                             << ", " << se_j.source_file_idx << ")" << std::endl;

            unite(static_cast<int>(ti), tj);
          }
        }
      }
    }

    // -----------------------------------------------------------------------
    // Collect Union-Find components into FeatureGroups.
    //
    // Singletons (no valid pair was found) are discarded.
    // For each group, accumulate per-fragment (count, max_intensity) to
    // identify shared fragments (count >= 2) and their best intensity.
    // -----------------------------------------------------------------------

    std::unordered_map<int, std::vector<int>> components;
    for (Size i = 0; i < traces.size(); ++i)
      components[find(static_cast<int>(i))].push_back(static_cast<int>(i));

    std::vector<FeatureGroup> result;
    result.reserve(components.size());

    for (auto& [root, members] : components)
    {
      if (members.size() < 2) continue;

      FeatureGroup group;
      group.spectrum_ids.reserve(members.size());

      double rt_sum = 0.0;
      std::unordered_map<uint32_t, std::pair<uint32_t, float>> frag_acc; // flat_idx -> (count, max_intensity)

      for (int mi : members)
      {
        group.spectrum_ids.push_back(traces[mi].spectrum_id);
        rt_sum += traces[mi].apex_rt;

        for (const ApexFragment& af : traces[mi].apex_fingerprint)
        {
          auto& acc = frag_acc[af.flat_bin_idx];
          ++acc.first;
          if (af.experimental_intensity > acc.second)
            acc.second = af.experimental_intensity;
        }
      }

      group.mean_apex_rt = rt_sum / static_cast<double>(members.size());

      for (const auto& [flat_idx, acc] : frag_acc)
      {
        if (acc.first >= 2)
          group.shared_fragments.push_back({flat_idx, acc.second});
      }

      std::sort(group.shared_fragments.begin(), group.shared_fragments.end(),
        [](const SharedFragment& a, const SharedFragment& b)
        { return a.flat_bin_idx < b.flat_bin_idx; });

      OPENMS_LOG_DEBUG << "[DiaWeaverAlign::groupFeatures] group: n_members="
                       << group.spectrum_ids.size()
                       << "  mean_apex_rt=" << group.mean_apex_rt
                       << "  n_shared_frags=" << group.shared_fragments.size()
                       << "  spectrum_ids=[";
      for (Size m = 0; m < group.spectrum_ids.size(); ++m)
      {
        OPENMS_LOG_DEBUG << group.spectrum_ids[m];
        if (m + 1 < group.spectrum_ids.size()) OPENMS_LOG_DEBUG << ",";
      }
      OPENMS_LOG_DEBUG << "]" << std::endl;

      result.push_back(std::move(group));
    }

    OPENMS_LOG_DEBUG << "[DiaWeaverAlign::groupFeatures] total groups (>=2 members): "
                     << result.size() << std::endl;

    std::sort(result.begin(), result.end(),
      [](const FeatureGroup& a, const FeatureGroup& b)
      { return a.mean_apex_rt < b.mean_apex_rt; });

    return result;
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
    const uint32_t flat_idx = hasIM()
      ? mz_bin * n_im_bins_ + im_bin
      : mz_bin;
    const FragmentEntry* base = fragment_entries_.data();
    return { base + bin_offsets_[flat_idx], base + bin_offsets_[flat_idx + 1] };
  }

  std::pair<const uint32_t*, const uint32_t*>
  DiaWeaverAlign::getPrecursorBinEntries(uint32_t mz_bin, uint32_t im_bin) const
  {
    const uint32_t flat_idx = hasIM()
      ? mz_bin * n_im_bins_ + im_bin
      : mz_bin;
    const uint32_t* base = precursor_entries_.data();
    return { base + precursor_offsets_[flat_idx], base + precursor_offsets_[flat_idx + 1] };
  }

  Size   DiaWeaverAlign::getPrecursorBinCount()  const { return n_precursor_bins_; }
  double DiaWeaverAlign::getPrecursorBinWidth()  const { return precursor_bin_width_; }
  double DiaWeaverAlign::getLowerPrecursorMz()   const { return lower_precursor_mz_; }
  double DiaWeaverAlign::getUpperPrecursorMz()   const { return upper_precursor_mz_; }

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
  bool   DiaWeaverAlign::hasIM() const
  {
    return !bin_offsets_.empty() &&
           bin_offsets_.size() == static_cast<Size>(n_bins_) * n_im_bins_ + 1;
  }

} // namespace OpenMS
