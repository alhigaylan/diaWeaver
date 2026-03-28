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
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/OpenMSConfig.h>

#include <cstdint>
#include <vector>

namespace OpenMS
{
  /**
    @brief Fragment index for spectral library alignment of diaWeaver pseudo-spectra.

    Builds a fragment index from a set of pseudo MS2 spectra (output of diaWeaver),
    mirroring MSFragger's fragment index design. Instead of theoretical peptide
    fragments, the index is built from observed peaks in pseudo MS2 spectra.

    **Logarithmic bin grid**

    Because ppm tolerance scales linearly with m/z, the bin width in Daltons
    grows with m/z. A fixed-Da bin width would be either too coarse at high m/z
    or unnecessarily fine at low m/z. This class uses a logarithmic bin grid
    instead, where every bin spans exactly @c fragment_mz_tolerance_ppm ppm:

    @code
      log_step  = log(1 + ppm / 1e6)
      bin_idx   = floor( log(mz / min_mz) / log_step )
      n_bins    = ceil( log(max_mz / min_mz) / log_step ) + 1
    @endcode

    The minimum and maximum m/z boundaries are determined automatically by
    scanning every peak in the input data; no user-supplied range is required.

    Each fragment peak is assigned to exactly one bin. The @c mass_offset field
    in FragmentEntry stores the fractional position of the peak within its bin
    (in log m/z space), so the exact m/z can always be recovered for verification.

    **Why adjacent bins are checked during search**

    When searching, a query peak is looked up in the bin it falls in. However,
    if the query sits near a bin boundary, a library peak that is within
    @c ppm tolerance of the query may have been assigned to the neighbouring bin.
    For example, with 20 ppm bins and edges at 999.98 / 1000.00 / 1000.02 Da:
    a query at 1000.001 Da and a library peak at 999.983 Da are 18 ppm apart
    (within tolerance) but land in different bins. Checking the query's bin and
    both immediate neighbours guarantees no match is missed. After the bin
    lookup, the exact ppm distance is verified using the recovered m/z so that
    only true matches within tolerance are accepted.

    **Index storage — CSR format**

    - @c fragment_entries_: flat array of 8-byte FragmentEntry structs, sorted
      first by bin index then by spectrum_id (and hence by precursor m/z).
    - @c bin_offsets_: array of size n_bins + 1 where @c bin_offsets_[b] is the
      index of the first entry for bin @c b in @c fragment_entries_.

    Each FragmentEntry is exactly 8 bytes (matching MSFragger's entry size):
    - spectrum_id (4 bytes): index into the spectrum array sorted by precursor m/z
    - mass_offset (2 bytes): fractional position within the bin in log space,
      scaled to [0, 65535]; allows recovery of the approximate fragment m/z as
      @c min_mz * exp((bin_idx + offset/65535) * log_step)
    - intensity_rank (1 byte): 0-based rank by descending intensity in the
      parent spectrum (0 = most intense)
    - charge (1 byte): precursor charge state (0 = unknown)

    **Searching**

    Candidates are first restricted to a precursor m/z window via binary search
    on the sorted spectrum array (analogous to MSFragger's precursor-mass filter).
    An optional ion mobility tolerance further narrows the candidate set. For each
    query peak the three relevant bins are scanned; within each bin a binary search
    restricts to the candidate spectrum ID range, and exact ppm distances are
    verified before accumulating the match score.

    @ingroup Analysis
  */
  class OPENMS_DLLAPI DiaWeaverAlign : public DefaultParamHandler
  {
  public:

    // -----------------------------------------------------------------------
    // Data structures
    // -----------------------------------------------------------------------

    /**
      @brief A single entry in the fragment index (exactly 8 bytes).

      Mirrors MSFragger's 8-byte fragment entry. The mass_offset field encodes
      the fractional position of the peak within its bin in log m/z space,
      allowing approximate m/z recovery:

        frag_mz ≈ min_observed_mz * exp((bin_idx + mass_offset / 65535.0) * log_step)
    */
#pragma pack(push, 1)
    struct FragmentEntry
    {
      uint32_t spectrum_id;    ///< Index of the parent spectrum (precursor m/z order)
      uint16_t mass_offset;    ///< Fractional log-space position within bin, scaled to [0, 65535]
      uint8_t  intensity_rank; ///< 0-based intensity rank within parent spectrum (0 = most intense)
      uint8_t  charge;         ///< Precursor charge state (0 = unknown)
    };
#pragma pack(pop)
    static_assert(sizeof(FragmentEntry) == 8, "FragmentEntry must be exactly 8 bytes");

    /**
      @brief Metadata for one indexed pseudo MS2 spectrum.

      Stored in precursor m/z ascending order; spectrum_id is the array index.
    */
    struct SpectrumEntry
    {
      double precursor_mz{0.0};    ///< Precursor m/z
      double precursor_im{-1.0};   ///< Precursor ion mobility (-1 = not available)
      int    charge{0};            ///< Precursor charge (0 = unknown)
      String native_id;            ///< Original native ID from source mzML
      Size   source_file_idx{0};   ///< Index into the source file list
    };

    /**
      @brief Result of a single spectrum search.
    */
    struct SearchResult
    {
      uint32_t spectrum_id{0};     ///< Index into the spectrum entry array
      uint32_t n_matched_peaks{0}; ///< Number of fragment peaks matched within tolerance
      double   cosine_score{0.0};  ///< Rank-weighted cosine similarity (see note in search())
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

      Scans all peaks to determine the observed m/z range, derives the log-space
      bin grid from that range and the ppm tolerance, sorts spectra by precursor
      m/z, and populates the CSR fragment index. Only the top
      @c max_peaks_per_spectrum peaks (by intensity) per spectrum are indexed.

      @param mzml_files  Paths to pseudo MS2 mzML files produced by diaWeaver.
      @pre All files must be centroided mzML files containing MS2 spectra.
    */
    void buildIndex(const std::vector<String>& mzml_files);


    // -----------------------------------------------------------------------
    // Searching
    // -----------------------------------------------------------------------

    /**
      @brief Search the index for spectra matching a query spectrum.

      1. Binary-search the precursor m/z-sorted spectrum array to obtain the
         candidate spectrum ID range within @c precursor_mz_tol_ppm.
      2. Optionally narrow candidates by ion mobility tolerance.
      3. For each query peak, check its bin and the two adjacent bins. Within
         each bin, binary-search for the candidate ID range and verify the exact
         ppm distance.
      4. Return results with >= @c min_matched_peaks matches, sorted by
         @c cosine_score descending.

      @note  The cosine score uses 1/(intensity_rank + 1) as the library peak
             weight because actual intensities are not stored in the 8-byte
             FragmentEntry. This is a rank-weighted approximation.

      @param query                 Centroided query spectrum.
      @param precursor_mz          Precursor m/z of the query spectrum.
      @param precursor_mz_tol_ppm  Precursor m/z tolerance in ppm.
      @param query_im              Precursor ion mobility of the query (-1 = not available).
      @param im_tolerance          Ion mobility tolerance (-1 = skip IM filtering).
      @return Matches sorted by cosine_score descending; empty if none pass the threshold.
    */
    std::vector<SearchResult> search(
      const MSSpectrum& query,
      double precursor_mz,
      double precursor_mz_tol_ppm,
      double query_im = -1.0,
      double im_tolerance = -1.0) const;


    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    /// Returns the metadata entry for the given spectrum ID.
    const SpectrumEntry& getSpectrumEntry(uint32_t id) const;

    /// Returns the source mzML file path for the given file index.
    const String& getSourceFile(Size idx) const;

    /// Returns the total number of indexed spectra.
    Size getSpectrumCount() const;

    /// Returns the total number of fragment entries across all bins.
    Size getTotalFragmentEntries() const;

    /// Returns the number of bins in the fragment index.
    Size getBinCount() const;

    /// Returns the ppm tolerance used to build the log-space bin grid.
    double getBinWidthPpm() const;


  protected:

    void updateMembers_() override;


  private:

    // -----------------------------------------------------------------------
    // Index storage (CSR format)
    // -----------------------------------------------------------------------

    /// Flat array of all fragment entries, sorted by (bin_idx, spectrum_id).
    std::vector<FragmentEntry> fragment_entries_;

    /// bin_offsets_[b] = first index in fragment_entries_ for bin b.
    /// Size is n_bins_ + 1; bin_offsets_[n_bins_] == fragment_entries_.size().
    std::vector<uint32_t> bin_offsets_;

    /// Spectrum metadata (precursor m/z ascending); index == spectrum_id.
    std::vector<SpectrumEntry> spectrum_entries_;

    /// Source mzML file paths, indexed by SpectrumEntry::source_file_idx.
    std::vector<String> source_files_;

    // -----------------------------------------------------------------------
    // Log-space bin grid parameters — set during buildIndex from the data
    // -----------------------------------------------------------------------

    double   log_step_{0.0};          ///< log(1 + ppm/1e6); the bin width in log space
    double   min_observed_mz_{0.0};   ///< Base of the log grid (smallest observed peak m/z)
    uint32_t n_bins_{0};              ///< Total number of bins

    // -----------------------------------------------------------------------
    // Cached parameter values (updated by updateMembers_)
    // -----------------------------------------------------------------------

    double fragment_tol_ppm_{20.0};
    Size   max_peaks_per_spectrum_{200};
    int    min_matched_peaks_{4};

    // -----------------------------------------------------------------------
    // Private helpers
    // -----------------------------------------------------------------------

    /// Map a fragment m/z to a bin index on the log grid.
    uint32_t toBinIdx_(double mz) const;

    /// Encode the fractional log-space position within a bin as a uint16.
    uint16_t toMassOffset_(double mz, uint32_t bin_idx) const;

    /// Recover the approximate fragment m/z from a bin index and mass_offset.
    double fromBinIdx_(uint32_t bin_idx, uint16_t mass_offset) const;
  };

} // namespace OpenMS
