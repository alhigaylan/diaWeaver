// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Mohammed Alhigaylan $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OpenMSConfig.h>

#include <cstdint>
#include <vector>

namespace OpenMS
{
  /**
    @brief Builds a fragment index from diaWeaver pseudo MS2 spectra.

    Iterates over a set of pseudo MS2 mzML files (output of diaWeaver) and maps
    every fragment peak onto a common logarithmic m/z bin axis. The axis spans
    the minimum and maximum fragment m/z observed across the entire dataset; no
    user-supplied range is required.

    **Logarithmic bin axis**

    Because ppm tolerance scales with m/z, a fixed-Da bin width would be too
    coarse at high m/z or unnecessarily fine at low m/z. Each bin therefore
    spans a constant number of ppm:

    @code
      log_step = log(1 + ppm / 1e6)          // bin width in log space
      bin_idx  = floor(log(mz / min_mz) / log_step)
      n_bins   = ceil(log(max_mz / min_mz) / log_step) + 1
    @endcode

    Each fragment peak is assigned to exactly one bin. The @c mass_offset field
    in FragmentEntry stores its fractional position within that bin (in log m/z
    space), so the exact m/z can be recovered at any time.

    **Index storage — CSR (Compressed Sparse Row) format**

    @code
      bin_offsets_[b]     → first index in fragment_entries_ for bin b
      bin_offsets_[b + 1] → one past the last index for bin b
    @endcode

    Within each bin, entries are ordered by @c spectrum_id, which corresponds
    to retention time order because spectra are sorted by RT before IDs are
    assigned.

    **Spectrum ordering**

    All spectra across all input files are sorted by retention time before being
    assigned contiguous 32-bit IDs. This means @c spectrum_id is a global RT
    rank across the entire dataset.
  */
  class OPENMS_DLLAPI DiaWeaverAlign : public DefaultParamHandler
  {
  public:

    // -----------------------------------------------------------------------
    // Data structures
    // -----------------------------------------------------------------------

    /**
      @brief A single entry in the fragment index (exactly 8 bytes).

      Each peak in a pseudo MS2 spectrum produces one FragmentEntry, placed in
      the bin corresponding to its m/z. The @c mass_offset encodes the peak's
      fractional log-space position within the bin, enabling m/z recovery:

        frag_mz = min_observed_mz * exp((bin_idx + mass_offset / 65535.0) * log_step)
    */
#pragma pack(push, 1)
    struct FragmentEntry
    {
      uint32_t spectrum_id;    ///< RT-ordered index of the parent spectrum
      uint16_t mass_offset;    ///< Fractional log-space position within bin, scaled to [0, 65535]
      uint8_t  intensity_rank; ///< 0-based rank by descending intensity within parent spectrum
      uint8_t  charge;         ///< Precursor charge state (0 = unknown)
    };
#pragma pack(pop)
    static_assert(sizeof(FragmentEntry) == 8, "FragmentEntry must be exactly 8 bytes");

    /**
      @brief Metadata for one indexed pseudo MS2 spectrum.

      Stored in retention time ascending order; @c spectrum_id is the array index.
    */
    struct SpectrumEntry
    {
      double retention_time{-1.0};  ///< Retention time in seconds
      double precursor_mz{0.0};     ///< Precursor m/z
      double precursor_im{-1.0};    ///< Precursor ion mobility (-1 = not available)
      int    charge{0};             ///< Precursor charge (0 = unknown)
      String native_id;             ///< Original native ID from source mzML
      Size   source_file_idx{0};    ///< Index into the source file list
    };


    // -----------------------------------------------------------------------
    // Construction / parameter handling
    // -----------------------------------------------------------------------

    DiaWeaverAlign();

    ~DiaWeaverAlign() override = default;


    // -----------------------------------------------------------------------
    // Index building
    // -----------------------------------------------------------------------

    /**
      @brief Build the fragment index from pseudo MS2 mzML files.

      - Loads all MS2 spectra from the given files.
      - Scans all peaks to determine the global m/z range and derive the
        logarithmic bin axis.
      - Sorts spectra by retention time and assigns contiguous 32-bit IDs.
      - For each spectrum, indexes the top @c max_peaks_per_spectrum peaks
        (by intensity) as FragmentEntry records in the CSR structure.

      @param mzml_files  Paths to pseudo MS2 mzML files produced by diaWeaver.
      @pre All files must be centroided mzML files containing MS2 spectra.
    */
    void buildIndex(const std::vector<String>& mzml_files);


    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    /// Returns the metadata for the given spectrum ID.
    const SpectrumEntry& getSpectrumEntry(uint32_t id) const;

    /// Returns the source mzML file path for the given file index.
    const String& getSourceFile(Size idx) const;

    /// Returns the total number of indexed spectra.
    Size getSpectrumCount() const;

    /// Returns the total number of fragment entries across all bins.
    Size getTotalFragmentEntries() const;

    /// Returns the number of bins in the fragment index.
    Size getBinCount() const;

    /// Returns the ppm value that defines the width of each bin.
    double getBinWidthPpm() const;

    /// Returns the log-space step size (= log(1 + ppm/1e6)).
    double getLogStep() const;

    /// Returns the minimum observed fragment m/z (anchor of the bin axis).
    double getMinObservedMz() const;


  protected:

    void updateMembers_() override;


  private:

    // -----------------------------------------------------------------------
    // Index storage (CSR format)
    // -----------------------------------------------------------------------

    /// All fragment entries, sorted by (bin_idx, spectrum_id).
    std::vector<FragmentEntry> fragment_entries_;

    /// bin_offsets_[b] = first index in fragment_entries_ for bin b.
    /// Size is n_bins_ + 1; bin_offsets_[n_bins_] == fragment_entries_.size().
    std::vector<uint32_t> bin_offsets_;

    /// Spectrum metadata sorted by retention time; index == spectrum_id.
    std::vector<SpectrumEntry> spectrum_entries_;

    /// Source mzML file paths, indexed by SpectrumEntry::source_file_idx.
    std::vector<String> source_files_;

    // -----------------------------------------------------------------------
    // Bin axis parameters — derived from data during buildIndex
    // -----------------------------------------------------------------------

    double   log_step_{0.0};         ///< Bin width in log space = log(1 + ppm/1e6)
    double   min_observed_mz_{0.0};  ///< Anchor (lowest observed fragment m/z)
    uint32_t n_bins_{0};             ///< Total number of bins

    // -----------------------------------------------------------------------
    // Cached parameter values
    // -----------------------------------------------------------------------

    double fragment_tol_ppm_{20.0};
    Size   max_peaks_per_spectrum_{200};

    // -----------------------------------------------------------------------
    // Private helpers
    // -----------------------------------------------------------------------

    /// Map a fragment m/z to its bin index.
    uint32_t toBinIdx_(double mz) const;

    /// Encode the fractional log-space position within a bin as a uint16.
    uint16_t toMassOffset_(double mz, uint32_t bin_idx) const;

    /// Recover the approximate fragment m/z from a bin index and mass_offset.
    double fromBinIdx_(uint32_t bin_idx, uint16_t mass_offset) const;
  };

} // namespace OpenMS
