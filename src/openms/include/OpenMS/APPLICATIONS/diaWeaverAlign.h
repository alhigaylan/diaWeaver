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
  class MSSpectrum;   // full definition in MSSpectrum.h, included by diaWeaverAlign.cpp
  class MSExperiment; // full definition in MSExperiment.h, included by diaWeaverAlign.cpp

  /**
    @brief Builds a fragment index from diaWeaver pseudo MS2 spectra.

    Reads MS2 spectra from a set of mzML files and maps every fragment peak
    onto a 2D bin grid (m/z × ion mobility) or a 1D m/z bin axis when no IM
    data is present. The index is stored in flat Compressed Sparse Row (CSR) format
    for O(1) bin access. In CSR, a one-dimensional offset array encodes where each bin's
    entries begin and end inside a contiguous entry array, avoiding per-bin dynamic allocation.

    **Bin axes**

      mz_bin(mz) = floor((mz  - lower_mz)  / bin_width)
      im_bin(im) = floor((im  - lower_im)  / bin_width_im)
      n_mz_bins  = ceil((upper_mz - lower_mz) / bin_width)
      n_im_bins  = ceil((upper_im - lower_im) / bin_width_im)

    Peaks outside [lower_mz, upper_mz] or [lower_im, upper_im] are ignored.

    **Flat 2D CSR storage**

      flat_idx                   = mz_bin * n_im_bins + im_bin
      bin_offsets_[flat_idx]     → first entry in fragment_entries_ for that cell
      bin_offsets_[flat_idx + 1] → one past the last entry

    When IM data is absent the index degenerates to 1D (flat_idx = mz_bin).
    Within each cell, entries are in load order (file order, then within-file order).
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
      double retention_time{-1.0};   ///< Retention time in seconds
      double precursor_mz{-1.0};     ///< Precursor m/z (-1 if unknown)
      double drift_time{};            ///< Ion mobility drift time. Only valid when hasIM() is true.
      String native_id;              ///< Native ID from the source mzML
      Size   source_file_idx{0};     ///< Index into the source file list
      int    precursor_charge{0};    ///< Precursor charge (0 = unknown)
    };

    /**
      @brief Score of one pseudo-MS2 spectrum against an experimental query spectrum.

      @c matched_peaks is the number of query peaks whose (m/z, IM) bin contained
      at least one fragment entry for this spectrum in the index.
    */
    struct MatchResult
    {
      uint32_t spectrum_id{0};     ///< Index into spectrum_entries_
      uint32_t matched_peaks{0};   ///< Number of query peaks that hit this spectrum's bins
    };

    /**
      @brief One fragment hit recorded at the apex query spectrum for a spectrum_id.

      Stores the flat 2D bin index and the experimental query peak intensity that
      produced the hit. Used for fragment-level similarity comparison between
      spectrum_ids that share the same precursor identity bucket.
    */
    struct ApexFragment
    {
      uint32_t flat_bin_idx{0};          ///< Flat 2D fragment bin (mz_bin * n_im_bins + im_bin)
      float    experimental_intensity{0}; ///< Intensity of the query peak that produced this hit
    };

    /**
      @brief Per-spectrum-id score trace across an MSExperiment of query spectra.

      @c rts and @c scores are parallel vectors in ascending RT order.
      Only query spectra whose matched_peaks count met the threshold passed to
      matchExperiment() contribute an entry. @c apex_rt and @c apex_score identify
      the single highest-scoring query RT. @c apex_fingerprint holds the fragment
      hits recorded at that apex query spectrum, used for isobaric disambiguation.
    */
    struct ScoreTrace
    {
      uint32_t             spectrum_id{0}; ///< Index into spectrum_entries_
      std::vector<double>   rts;           ///< Query spectrum RTs at which score >= threshold, ascending
      std::vector<uint32_t> scores;        ///< matched_peaks count at each RT, parallel to rts
      double   apex_rt{-1.0};             ///< Query RT at which matched_peaks was highest
      uint32_t apex_score{0};             ///< Highest matched_peaks count across all query RTs
      std::vector<ApexFragment> apex_fingerprint; ///< Fragment hits at the apex query spectrum
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

    /**
      @brief Score an experimental peak-picked MS2 spectrum against the fragment index.

      For each query peak, all spectrum_ids stored in the matching (m/z, IM) bin are
      collected into a flat list. That list is sorted and scanned once to count how many
      query peaks mapped to each spectrum_id — the matched_peaks score.

      If the index was built with IM data the query spectrum must also carry per-peak IM
      (containsIMData()); a mismatch throws Exception::MissingInformation. If the index
      is 1D, IM on the query is ignored.

      @param query  Centroided experimental MS2 spectrum.
      @return       One MatchResult per pseudo-MS2 spectrum that shared at least one bin
                    with the query, in ascending spectrum_id order.
    */
    std::vector<MatchResult> matchSpectrum(const MSSpectrum& query) const;

    /**
      @brief Score every MS2 spectrum in @p experiment against the fragment index.

      Iterates over all MS2 spectra in load order (ascending RT), calling
      matchSpectrum() for each. Scores that meet @p min_matched_peaks are
      accumulated into per-spectrum_id traces. Query spectra are distributed
      across available OpenMP threads; each thread accumulates into a private
      map and results are merged after all threads complete. Because threads
      process contiguous RT chunks, merging in thread-index order preserves
      RT ordering within each trace.

      IM consistency is validated before the parallel loop: if the index was
      built with IM the experiment must carry per-peak IM on every MS2 spectrum,
      and vice versa.

      @param experiment         Peak-picked experimental MSExperiment.
      @param min_matched_peaks  Minimum matched_peaks count to include in a trace (default 2).
      @return                   One ScoreTrace per spectrum_id that accumulated at least one
                                above-threshold observation, sorted by spectrum_id.
    */
    std::vector<ScoreTrace> matchExperiment(const MSExperiment& experiment,
                                             uint32_t min_matched_peaks = 2) const;

    /**
      @brief Build the precursor index from the loaded spectrum metadata.

      Maps every spectrum_id onto a 2D bin (precursor m/z × ion mobility) or a
      1D m/z bin when no IM data is present, using the same CSR layout as the
      fragment index. Must be called after buildIndex().

      Spectrum_ids whose precursor_mz falls outside [lower_mz, upper_mz] are skipped.
      The precursor bin width is controlled by the @c precursor_bin_width parameter
      (default 0.02 Da), which is independent of the fragment @c bin_width.
    */
    void buildPrecursorIndex();

    /**
      @brief Returns the spectrum_ids stored in a given precursor (m/z bin, IM bin) cell.

      When hasIM() is false, @p im_bin must be 0 (the index is 1D).
      The returned pointers are valid until the next call to buildIndex() or buildPrecursorIndex().
    */
    std::pair<const uint32_t*, const uint32_t*>
      getPrecursorBinEntries(uint32_t mz_bin, uint32_t im_bin = 0) const;

    Size   getPrecursorBinCount()  const;
    double getPrecursorBinWidth()  const;
    double getLowerPrecursorMz()   const;
    double getUpperPrecursorMz()   const;

    /// Returns the metadata for the given spectrum ID.
    const SpectrumEntry& getSpectrumEntry(uint32_t id) const;

    /// Returns the source mzML file path for the given file index.
    const String& getSourceFile(Size idx) const;

    /**
      @brief Returns the fragment entries for a given (m/z bin, IM bin) cell.

      When hasIM() is false, @p im_bin must be 0 (the index is 1D).
      The returned pointers are valid until the next call to buildIndex().
    */
    std::pair<const FragmentEntry*, const FragmentEntry*>
      getBinEntries(uint32_t mz_bin, uint32_t im_bin = 0) const;

    Size   getSpectrumCount()        const;
    Size   getTotalFragmentEntries() const;
    Size   getMZBinCount()           const;
    Size   getIMBinCount()           const;
    double getBinWidth()             const;
    double getBinWidthIM()           const;
    double getLowerMz()              const;
    double getUpperMz()              const;
    double getLowerIM()              const;
    double getUpperIM()              const;

    /// Returns true if the index was built with ion mobility data (derived from the index structure, not a stored flag).
    bool   hasIM()                   const;


  protected:

    void updateMembers_() override;


  private:

    std::vector<FragmentEntry> fragment_entries_;   ///< All entries, sorted by (flat_idx, spectrum_id)
    std::vector<uint32_t>      bin_offsets_;         ///< CSR offsets, size = n_bins_*n_im_bins_ + 1
    std::vector<uint32_t>      precursor_entries_;   ///< spectrum_ids sorted by precursor flat_idx
    std::vector<uint32_t>      precursor_offsets_;   ///< CSR offsets, size = n_precursor_bins_*n_im_bins_ + 1
    std::vector<SpectrumEntry> spectrum_entries_;    ///< Spectrum metadata, in load order
    std::vector<String>        source_files_;        ///< Source file paths

    // Derived from parameters in updateMembers_
    double   lower_mz_{400.0};
    double   upper_mz_{2000.0};
    double   bin_width_{0.1};
    uint32_t n_bins_{0};
    double   lower_im_{0.60};
    double   upper_im_{1.70};
    double   bin_width_im_{0.01};
    uint32_t n_im_bins_{0};
    double   lower_precursor_mz_{400.0};
    double   upper_precursor_mz_{1200.0};
    double   precursor_bin_width_{0.02};
    uint32_t n_precursor_bins_{0};

    uint32_t toBinIdx_(double mz) const;
    uint32_t toBinIdx_im_(double im) const;
    uint16_t toMassOffset_(double mz, uint32_t bin_idx) const;
    uint32_t toPrecursorBinIdx_(double mz) const;
  };

} // namespace OpenMS
