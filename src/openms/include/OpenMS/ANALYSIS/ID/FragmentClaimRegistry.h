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

#include <unordered_map>
#include <vector>

namespace OpenMS
{

  /**
   * @brief Tracks which peptide has claimed ownership of each MassTrace fragment.
   *
   * During iterative re-scoring in diaWeaverPeptide, PSMs are sorted by score
   * (descending). For each PSM the registry attempts to claim the associated
   * fragment trace IDs. A trace can only be claimed by one peptide — the first
   * (highest-scoring) one to attempt it wins.
   *
   * Once claiming is complete, buildExclusionMask() produces a per-peak boolean
   * vector for a given spectrum, indicating which peaks have already been claimed
   * by a *different* peptide. Those peaks are excluded when re-scoring the spectrum.
   */
  class OPENMS_DLLAPI FragmentClaimRegistry
  {
  public:
    struct ClaimRecord
    {
      String   peptide_seq;
      double   score;
      Size     spectrum_idx;
    };

    /// Attempt to claim all trace IDs listed in @p trace_ids for @p peptide_seq.
    /// Only unclaimed trace IDs are claimed; already-claimed ones are silently skipped.
    /// Returns the number of IDs that were successfully claimed (not previously taken).
    Size tryClaim(const std::vector<Int>& trace_ids,
                  const String& peptide_seq,
                  double score,
                  Size spectrum_idx);

    /// True if @p trace_id has been claimed by any peptide.
    bool isClaimed(Int trace_id) const;

    /// Return the ClaimRecord for a trace, or nullptr if unclaimed.
    const ClaimRecord* getClaimRecord(Int trace_id) const;

    /**
     * @brief Build a boolean exclusion mask for a spectrum.
     *
     * For each peak in @p spec, the mask entry is `true` if the corresponding
     * trace_id has been claimed by a peptide OTHER than @p query_peptide_seq.
     * Peaks whose trace_id is unclaimed, or claimed by the query peptide itself,
     * are NOT excluded (mask entry = `false`).
     *
     * The returned vector has the same length as the number of peaks in @p spec.
     * Returns an empty vector if the spectrum has no "fragment_trace_id" array.
     *
     * @param spec            Pseudo spectrum (must carry "fragment_trace_id" IntegerDataArray)
     * @param query_peptide_seq  Sequence of the peptide being re-scored
     */
    std::vector<bool> buildExclusionMask(const MSSpectrum& spec,
                                         const String& query_peptide_seq) const;

    Size claimedCount() const noexcept { return claims_.size(); }

    void clear() noexcept { claims_.clear(); }

  private:
    std::unordered_map<Int, ClaimRecord> claims_;
  };

}  // namespace OpenMS
