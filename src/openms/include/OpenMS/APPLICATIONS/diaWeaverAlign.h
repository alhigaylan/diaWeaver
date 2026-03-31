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

    Reads MS2 spectra from a set of mzML files and maps every fragment peak
    onto a fixed-width Da bin axis defined by a user-specified m/z range and
    bin width. The index is stored in CSR format for efficient column access.

    **Bin axis**

      bin_idx(mz) = floor((mz - lower_mz) / bin_width)
      n_bins      = ceil((upper_mz - lower_mz) / bin_width)

    Peaks outside [lower_mz, upper_mz] are ignored.

    **CSR storage**

      bin_offsets_[b]     → first entry in fragment_entries_ for bin b
      bin_offsets_[b + 1] → one past the last entry for bin b

    Within each bin, entries are in retention time order (spectrum_id order).
  */
  class OPENMS_DLLAPI DiaWeaverAlign : public DefaultParamHandler
  {
  public:

    /**
      @brief One entry in the fragment index (exactly 8 bytes).

      Each indexed peak produces one FragmentEntry stored in the bin that
      contains its m/z.
    */
#pragma pack(push, 1)
    struct FragmentEntry
    {
      uint32_t spectrum_id; ///< RT-ordered index of the parent spectrum
      uint16_t mass_offset; ///< Sub-bin Da offset scaled to [0, 65535]:
                            ///  offset_da = (mass_offset / 65535.0) * bin_width
      uint8_t  charge;      ///< Precursor charge (0 = unknown)
      uint8_t  reserved{0}; ///< Reserved for future use
    };
#pragma pack(pop)
    static_assert(sizeof(FragmentEntry) == 8, "FragmentEntry must be exactly 8 bytes");

    /**
      @brief Metadata for one indexed pseudo MS2 spectrum.

      Stored in retention time ascending order; @c spectrum_id is the index.
    */
    struct SpectrumEntry
    {
      double retention_time{-1.0};  ///< Retention time in seconds
      String native_id;             ///< Native ID from the source mzML
      Size   source_file_idx{0};    ///< Index into the source file list
    };


    DiaWeaverAlign();

    ~DiaWeaverAlign() override = default;

    /**
      @brief Build the fragment index from a list of pseudo MS2 mzML files.

      Loads all MS2 spectra, sorts them by retention time, assigns contiguous
      32-bit IDs, and populates the CSR fragment index. Peaks outside
      [lower_mz, upper_mz] are skipped.

      @param mzml_files  Paths to centroided pseudo MS2 mzML files.
    */
    void buildIndex(const std::vector<String>& mzml_files);

    /// Returns the metadata for the given spectrum ID.
    const SpectrumEntry& getSpectrumEntry(uint32_t id) const;

    /// Returns the source mzML file path for the given file index.
    const String& getSourceFile(Size idx) const;

    /// Returns the fragment entries for a given bin (contiguous span in fragment_entries_).
    /// The returned pointers are valid until the next call to buildIndex().
    std::pair<const FragmentEntry*, const FragmentEntry*> getBinEntries(uint32_t bin_idx) const;

    Size   getSpectrumCount()        const;
    Size   getTotalFragmentEntries() const;
    Size   getBinCount()             const;
    double getBinWidth()             const;
    double getLowerMz()              const;
    double getUpperMz()              const;


  protected:

    void updateMembers_() override;


  private:

    std::vector<FragmentEntry> fragment_entries_;  ///< All entries, sorted by (bin, spectrum_id)
    std::vector<uint32_t>      bin_offsets_;        ///< CSR offsets, size = n_bins_ + 1
    std::vector<SpectrumEntry> spectrum_entries_;   ///< Spectrum metadata, RT-sorted
    std::vector<String>        source_files_;       ///< Source file paths

    // Derived from parameters in updateMembers_
    double   lower_mz_{400.0};
    double   upper_mz_{2000.0};
    double   bin_width_{0.1};
    uint32_t n_bins_{0};

    uint32_t toBinIdx_(double mz) const;
    uint16_t toMassOffset_(double mz, uint32_t bin_idx) const;
  };

} // namespace OpenMS
