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
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/OpenMSConfig.h>

#include <cstdint>
#include <cstdio>
#include <vector>

namespace OpenMS
{

  /**
    @brief Builds a fragment index from diaWeaver pseudo MS2 spectra.

    Reads MS2 spectra from a set of mzML files and maps every fragment peak
    onto a 2D bin grid (m/z × ion mobility) or a 1D m/z bin axis when no IM
    data is present. The index is stored in flat Compressed Sparse Row (CSR) format
    for O(1) bin access. In CSR, a one-dimensional offset array encodes where each bin's
    entries begin and end inside a contiguous entry array, avoiding per-bin dynamic allocation.

    **Bin axes**

      mz_bin(mz) = floor(log(mz / lower_mz) / log(1 + ppm/1e6))
      im_bin(im) = floor((im  - lower_im)  / bin_width_im)
      n_mz_bins  = floor(log(upper_mz / lower_mz) / log(1 + ppm/1e6)) + 1
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
      @brief One entry in the fragment index (exactly 9 bytes).

      Each indexed peak produces one FragmentEntry stored in the bin that
      contains its m/z.

      log2_intensity encodes log2(peak_intensity) as a fixed-point uint16_t
      scaled by 2048.  Decodes as: log2_val = log2_intensity / 2048.0f.
      Covers intensities up to 2^32 (≈ 4e9) with ~0.0005 log2 resolution.
    */
#pragma pack(push, 1)
    struct FragmentEntry
    {
      uint32_t spectrum_id;       ///< RT-ordered index of the parent spectrum
      uint16_t mass_offset;       ///< Fractional offset within the bin, scaled to [0, 65535]:
                                  ///  frac = mass_offset / 65535.0  (0 = bin start, 1 = bin end)
      uint8_t  charge;            ///< Precursor charge (0 = unknown)
      uint16_t log2_intensity{0}; ///< log2(peak intensity) * 2048, clamped to [0, 65535]
    };
#pragma pack(pop)
    static_assert(sizeof(FragmentEntry) == 9, "FragmentEntry must be exactly 9 bytes");

    /**
      @brief Metadata for one indexed pseudo MS2 spectrum.

      Stored in retention time ascending order; @c spectrum_id is the index.
    */
    struct SpectrumEntry
    {
      double   retention_time{-1.0};      ///< Retention time in seconds
      double   precursor_mz{-1.0};        ///< Precursor m/z (-1 if unknown)
      double   drift_time{};              ///< Ion mobility drift time. Only valid when hasIM() is true.
      String   native_id;                 ///< Native ID from the source mzML
      Size     source_file_idx{0};        ///< Index into the source file list
      int      precursor_charge{0};       ///< Precursor charge (0 = unknown)
      uint32_t n_indexed_fragments{0};    ///< Number of fragment bins this spectrum occupies in the index
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
      uint32_t flat_bin_idx{0};            ///< Flat 2D fragment bin (mz_bin * n_im_bins + im_bin)
      float    experimental_intensity{0};  ///< Intensity of the query peak that produced this hit
      float    index_log2_intensity{0.0f}; ///< log2(index peak intensity) decoded from FragmentEntry::log2_intensity
    };

    /**
      @brief One shared fragment entry in a FeatureGroup.

      Represents a fragment bin seen in at least two grouped spectrum_ids, carrying
      the maximum experimental intensity observed across all group members at that bin.
    */
    struct SharedFragment
    {
      uint32_t flat_bin_idx{0};   ///< Flat 2D fragment bin (mz_bin * n_im_bins + im_bin)
      float    max_intensity{0};  ///< Max experimental intensity across all group members
    };

    /**
      @brief A group of spectrum_ids aligned to the same compound across source files.

      Members share the same precursor bin, have close apex RTs in the query experiment,
      and have sufficient fragment overlap. @c shared_fragments lists every fragment bin
      seen in at least two members, with the maximum intensity across the group.
    */
    struct FeatureGroup
    {
      std::vector<uint32_t>       spectrum_ids;          ///< Grouped spectrum_ids, one per source file ideally
      std::vector<SharedFragment> shared_fragments;      ///< Fragment bins shared by >=2 members within this group, sorted by flat_bin_idx
      std::vector<uint32_t>       unique_fragment_bins;  ///< Fragment bins present in this sub-group's union fingerprint but absent from every other sub-group produced by the same split. Empty when sub-clustering was not triggered.
      double                      mean_apex_rt{-1.0};    ///< Mean apex RT across all members (seconds)
      double                      min_internal_overlap{1.0}; ///< Minimum pairwise Overlap Coefficient (|A∩B|/min(|A|,|B|)) within this group (1.0 when sub-clustering was not triggered).
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

    ~DiaWeaverAlign() override;

    /**
      @brief Build the fragment index from a list of pseudo MS2 mzML files.

      Loads all MS2 spectra, sorts them by retention time, assigns contiguous
      32-bit IDs, and populates the CSR fragment index. Peaks outside
      [lower_mz, upper_mz] are skipped.

      @param mzml_files  Paths to centroided pseudo MS2 mzML files.
    */
    void buildIndex(const std::vector<String>& mzml_files);

    /**
      @brief Build the fragment index from pre-loaded MSExperiment objects.

      Same as the file-path overload but accepts experiments already in memory,
      avoiding a second round of disk I/O when the caller has already loaded the
      data (e.g. after window-based partitioning). Each element of @p experiments
      corresponds to one logical source file; @p source_names provides the display
      name stored in the source-file list (may be empty strings).

      Spectra without precursors or with empty peak lists are skipped.

      @param experiments   Pre-loaded pseudo-MS2 experiments, one per source file.
      @param source_names  Display names parallel to @p experiments (may be shorter;
                           missing entries default to an empty string).
    */
    void buildIndex(const std::vector<MSExperiment>& experiments,
                    const std::vector<String>& source_names = {});

    /**
      @brief Open a streaming session for incremental index construction.

      Initialises the m/z and IM axis parameters, opens binary raw files for
      writing, and resets all streaming counters. Must be called before any
      calls to appendSpectrumToStream().

      @param frags_path  Path to the raw fragment entries file to create/overwrite.
      @param meta_path   Path to the raw spectrum metadata file to create/overwrite.
    */
    void openStream(const String& frags_path, const String& meta_path);

    /**
      @brief Append one pseudo-MS2 spectrum to the open streaming session.

      Bins every fragment peak of @p spec and writes the resulting (flat_idx,
      FragmentEntry) pairs to the raw fragment file. Spectrum metadata is written
      to the raw meta file. Must be called after openStream().

      Spectra with MS level != 2, empty peak lists, or missing precursors are
      silently skipped.

      @param spec             Pseudo-MS2 spectrum to index.
      @param source_file_idx  Index of the source file in the final source-file list.
    */
    void appendSpectrumToStream(const MSSpectrum& spec, Size source_file_idx);

    /**
      @brief Finalise the open streaming session.

      Seeks back to the raw file headers and writes the actual spectrum count and
      IM flag, then closes both FILE* handles. After this call the raw files are
      complete and ready for finalizeFromStream().
    */
    void closeStream();

    /**
      @brief Build the in-memory fragment index from raw streaming files.

      Reads the raw fragment and metadata files produced by openStream /
      appendSpectrumToStream / closeStream, sorts the entries, and constructs
      the full CSR index and spectrum_entries_. The precursor index is NOT built
      here; call buildPrecursorIndex() afterwards.

      @param frags_path        Path to the raw fragment entries file (.dwaf).
      @param meta_path         Path to the raw spectrum metadata file (.dwam).
      @param source_file_names Source file display names (one per logical source file).
    */
    void finalizeFromStream(const String& frags_path,
                            const String& meta_path,
                            const std::vector<String>& source_file_names);

    /**
      @brief Serialize the complete index (fragment + precursor + metadata) to a binary file.

      Must be called after buildIndex() (or finalizeFromStream() + buildPrecursorIndex()).
      The resulting file can be reloaded with loadIndex() without re-building.

      @param path  Output file path (conventionally .dwaindex).
    */
    /**
      @brief Serialize the complete index (fragment + precursor + metadata) to a binary file.

      Stores all bin parameters, the DIA window bounds set via setWindowBounds(),
      source files, spectrum metadata, and all CSR arrays. The resulting file is
      fully self-describing and can be reloaded with loadIndex().

      @param path  Output file path (conventionally .dwaindex).
    */
    void saveIndex(const String& path) const;

    /**
      @brief Deserialize a previously saved index from a binary file.

      Restores all private data members including DIA window bounds, making the
      object immediately usable for matchSpectrum(), matchExperiment(), and
      groupFeatures() without calling buildIndex() or buildPrecursorIndex().

      @param path  Input file path produced by saveIndex().
    */
    void loadIndex(const String& path);

    /**
      @brief Store the DIA isolation window this index was built for.

      Must be called before saveIndex() so the window bounds are embedded in the
      binary file. The values are retrieved with getWindowLowerMz() etc. after
      loadIndex().
    */
    void setWindowBounds(double lower_mz, double upper_mz,
                         double lower_im = -1.0, double upper_im = -1.0);

    double getWindowLowerMz() const;
    double getWindowUpperMz() const;
    double getWindowLowerIM() const;
    double getWindowUpperIM() const;

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

      Spectrum_ids whose precursor_mz falls outside [lower_precursor_mz, upper_precursor_mz] are skipped.
      The precursor bin width uses logarithmic ppm binning, controlled by the
      @c precursor_ppm_tolerance parameter (default 20 ppm), independent of fragment binning.
    */
    void buildPrecursorIndex();

    /**
      @brief Group aligned spectrum_ids into feature groups across source files.

      Operates on the output of matchExperiment(). For each ScoreTrace, queries
      the precursor index in a ±1 bin neighborhood to find candidate spectrum_ids
      with a similar precursor identity. Candidate pairs are validated against three
      criteria: different source files, apex RT within @p rt_tolerance, and at least
      @p min_fragment_overlap shared flat_bin_idx values in their apex fingerprints.
      Valid pairs are merged with Union-Find; the final connected components become
      FeatureGroups.

      For each FeatureGroup, shared_fragments lists every flat_bin_idx seen in >=2
      members, carrying the maximum experimental intensity across the group.

      Must be called after buildPrecursorIndex().

      @param traces                  Output of matchExperiment().
      @param rt_tolerance            Maximum apex RT difference in seconds (default 10.0).
      @param min_fragment_overlap    Minimum shared flat_bin_idx count for a valid pair (default 5).
      @param im_tolerance            Maximum drift time difference (Vs/cm^2) between two
                                     spectrum-ids to be considered the same precursor identity
                                     (default 0.01). Ignored when the index has no IM data.
      @param precursor_ppm_tolerance Maximum precursor m/z difference in ppm (default 20.0).
      @return                        FeatureGroups with >=2 members, sorted by mean_apex_rt.
    */
    std::vector<FeatureGroup> groupFeatures(const std::vector<ScoreTrace>& traces,
                                             double   rt_tolerance            = 10.0,
                                             uint32_t min_fragment_overlap    = 5,
                                             double   im_tolerance            = 0.01,
                                             double   precursor_ppm_tolerance = 20.0,
                                             uint32_t isotope_error_tol       = 3,
                                             double   min_overlap_similarity  = 0.7) const;

    /**
      @brief Returns the spectrum_ids stored in a given precursor (m/z bin, IM bin) cell.

      When hasIM() is false, @p im_bin must be 0 (the index is 1D).
      The returned pointers are valid until the next call to buildIndex() or buildPrecursorIndex().
    */
    std::pair<const uint32_t*, const uint32_t*>
      getPrecursorBinEntries(uint32_t mz_bin, uint32_t im_bin = 0) const;

    Size   getPrecursorBinCount()  const;
    double getPrecursorPPMTolerance() const;
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
    double getFragmentPPMTolerance() const;
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
    double   fragment_ppm_tolerance_{20.0};
    double   log_bin_ratio_{0.0};  ///< Pre-computed log(1 + fragment_ppm_tolerance_/1e6)
    uint32_t n_bins_{0};
    double   lower_im_{0.60};
    double   upper_im_{1.70};
    double   bin_width_im_{0.01};
    uint32_t n_im_bins_{0};
    double   lower_precursor_mz_{400.0};
    double   upper_precursor_mz_{1200.0};
    double   precursor_ppm_tolerance_{20.0};
    double   log_precursor_bin_ratio_{0.0};  ///< Pre-computed log(1 + precursor_ppm_tolerance_/1e6)
    uint32_t n_precursor_bins_{0};

    /// Internal per-spectrum record used during index construction.
    struct RawSpecEntry_
    {
      MSSpectrum spectrum;
      double     retention_time{-1.0};
      double     precursor_mz{-1.0};
      double     drift_time{-1.0};
      int        charge{0};
      String     native_id;
      Size       source_idx{0};
    };

    /// Shared phases 2-4 of index construction (spectrum-ID assignment,
    /// fragment binning, CSR build). Called by both buildIndex overloads.
    void buildIndexCore_(std::vector<RawSpecEntry_>& raw_spectra);

    /// Append spectra from one MSExperiment to @p raw_spectra.
    /// Returns the updated @p im_initialised flag.
    static bool collectSpectraFromExperiment_(
      const MSExperiment& exp,
      Size                file_idx,
      const String&       file_label,
      bool&               im_detected,
      bool                im_initialised,
      std::vector<RawSpecEntry_>& raw_spectra);

    uint32_t toBinIdx_(double mz) const;
    uint32_t toBinIdx_im_(double im) const;
    uint16_t toMassOffset_(double mz, uint32_t bin_idx) const;
    uint32_t toPrecursorBinIdx_(double mz) const;

    /**
      @brief Bin one spectrum's fragment peaks into (flat_idx, FragmentEntry) pairs.

      Applies intensity-max-per-(mz_bin, im_bin) deduplication. Results are
      appended to @p out_entries. Used by appendSpectrumToStream(); buildIndexCore_()
      uses its own optimised batch loop.

      @param spec         Centroided pseudo-MS2 spectrum.
      @param spec_id      Spectrum ID to store in each FragmentEntry.
      @param charge_byte  Precursor charge truncated to uint8_t.
      @param im_detected  True if the index is being built with ion mobility data.
      @param out_entries  Output vector; entries are appended (not cleared).
    */
    void binSingleSpectrum_(const MSSpectrum&                              spec,
                            uint32_t                                       spec_id,
                            uint8_t                                        charge_byte,
                            bool                                           im_detected,
                            std::vector<std::pair<uint32_t, FragmentEntry>>& out_entries) const;

    /**
      @brief Sort @p all_entries and build the CSR fragment index arrays.

      Sorts by (flat_idx, spectrum_id) then populates fragment_entries_ and
      bin_offsets_. The @p im_detected flag determines whether a 2D or 1D flat
      index is used for the offset array size.
    */
    void buildCSR_(std::vector<std::pair<uint32_t, FragmentEntry>>& all_entries,
                   bool                                              im_detected);

    /**
      @brief Split one FeatureGroup into sub-groups using complete-linkage
             re-clustering on the Overlap Coefficient of apex fingerprints.

      For each pair of members the pairwise positional tolerances are checked
      first (apex RT within @p rt_tolerance, drift time within @p im_tolerance
      when IM data is present, isotope-corrected precursor m/z within
      @p precursor_ppm_tolerance). Pairs that violate any tolerance are
      assigned a similarity of 0 and can never be merged. For pairs that pass
      the positional checks the Overlap Coefficient (|A∩B| / min(|A|,|B|)) is
      computed over their sorted flat_bin_idx sets; this metric is insensitive
      to the size asymmetry that arises when one member's fingerprint is
      inflated by being from the same LC-MS run as the query. Complete-linkage
      agglomerative clustering merges the pair with the highest similarity; a
      merge is accepted only when the similarity exceeds @p min_overlap_similarity.
      For each resulting sub-group the unique_fragment_bins (bins absent from
      every other sub-group's union fingerprint) and min_internal_overlap are
      populated. Singleton sub-clusters are discarded. Returns an empty vector
      when all members become singletons.
    */
    static std::vector<FeatureGroup> splitByOverlap_(
      const FeatureGroup&                         group,
      const std::vector<ScoreTrace>&              traces,
      const std::vector<std::vector<uint32_t>>&   fingerprint_bins,
      const std::vector<int>&                     spec_to_trace_idx,
      const std::vector<SpectrumEntry>&            spectrum_entries,
      double                                       min_overlap_similarity,
      double                                       rt_tolerance,
      double                                       im_tolerance,
      double                                       precursor_ppm_tolerance,
      uint32_t                                     isotope_error_tol,
      bool                                         use_im);

    // --- DIA window bounds (stored in .dwaindex for self-describing files) ---
    double window_lower_mz_{-1.0};
    double window_upper_mz_{-1.0};
    double window_lower_im_{-1.0};
    double window_upper_im_{-1.0};

    // --- Streaming state ---
    FILE*    stream_frags_fp_{nullptr};      ///< Raw fragment entries file (write, open between openStream/closeStream)
    FILE*    stream_meta_fp_{nullptr};       ///< Raw spectrum metadata file (write, open between openStream/closeStream)
    uint32_t stream_next_spec_id_{0};        ///< Running spectrum_id counter
    bool     stream_im_detected_{false};     ///< IM detected in first streamed spectrum
    bool     stream_im_initialised_{false};  ///< True once the first spectrum has been processed
  };

} // namespace OpenMS
