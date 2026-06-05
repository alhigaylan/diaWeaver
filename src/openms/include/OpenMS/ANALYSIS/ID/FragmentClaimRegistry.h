// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: OpenMS Developers $
// $Authors: OpenMS Developers $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <cstdint>
#include <unordered_map>
#include <vector>

namespace OpenMS
{

  /**
   * @brief Tracks which peptide has claimed ownership of each MassTrace fragment.
   *
   * Fragment traces are globally identified by a (window_id, trace_id) pair encoded
   * as a single int64_t key: `(window_id << 32) | static_cast<uint32_t>(trace_id)`.
   * This guarantees uniqueness across all DIA windows regardless of the number of
   * traces per window (supports up to 2^32 traces per window and 2^32 windows).
   *
   * During iterative re-scoring in diaWeaverPeptide, PSMs are sorted by score
   * (descending). For each PSM the registry attempts to claim the associated
   * fragment trace keys. A trace can only be claimed by one peptide — the first
   * (highest-scoring) one to attempt it wins.
   */
  class OPENMS_DLLAPI FragmentClaimRegistry
  {
  public:
    using TraceKey = int64_t;

    struct ClaimRecord
    {
      String   peptide_seq;
      double   score;
      Size     spectrum_idx;
    };

    /// Encode a (window_id, local trace_id) pair into a globally unique 64-bit key.
    static TraceKey makeKey(Int window_id, Int trace_id) noexcept
    {
      return (static_cast<int64_t>(window_id) << 32) |
             static_cast<int64_t>(static_cast<uint32_t>(trace_id));
    }

    /// Attempt to claim all keys in @p keys for @p peptide_seq.
    /// Only unclaimed keys are claimed; already-claimed ones are silently skipped.
    /// Returns the number of keys that were successfully claimed.
    Size tryClaim(const std::vector<TraceKey>& keys,
                  const String& peptide_seq,
                  double score,
                  Size spectrum_idx);

    /// True if @p key has been claimed by any peptide.
    bool isClaimed(TraceKey key) const;

    /// Return the ClaimRecord for a key, or nullptr if unclaimed.
    const ClaimRecord* getClaimRecord(TraceKey key) const;

    /**
     * @brief Build a boolean exclusion mask for a spectrum.
     *
     * For each peak in @p spec, the mask entry is `true` if the corresponding
     * (window_id, trace_id) key has been claimed by a peptide OTHER than
     * @p query_peptide_seq.  Peaks whose key is unclaimed, or claimed by the
     * query peptide itself, are NOT excluded (mask entry = `false`).
     *
     * Requires both "fragment_trace_id" and "fragment_window_id" IntegerDataArrays.
     * Returns an empty vector if either array is absent.
     *
     * @param spec               Pseudo spectrum carrying both IntegerDataArrays
     * @param query_peptide_seq  Sequence of the peptide being re-scored
     */
    std::vector<bool> buildExclusionMask(const MSSpectrum& spec,
                                         const String& query_peptide_seq) const;

    /// Merge all claims from @p other into this registry.
    /// If a key is already claimed in this registry, the existing record is kept (first claimer wins).
    void merge(const FragmentClaimRegistry& other)
    {
      for (const auto& [key, record] : other.claims_)
        claims_.try_emplace(key, record);
    }

    Size claimedCount() const noexcept { return claims_.size(); }

    void clear() noexcept { claims_.clear(); }

  private:
    std::unordered_map<TraceKey, ClaimRecord> claims_;
  };

}  // namespace OpenMS
