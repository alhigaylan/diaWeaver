// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: OpenMS Developers $
// $Authors: OpenMS Developers $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/FragmentClaimRegistry.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

namespace OpenMS
{

  Size FragmentClaimRegistry::tryClaim(const std::vector<Int>& trace_ids,
                                        const String& peptide_seq,
                                        double score,
                                        Size spectrum_idx)
  {
    Size claimed = 0;
    for (Int tid : trace_ids)
    {
      if (claims_.emplace(tid, ClaimRecord{peptide_seq, score, spectrum_idx}).second)
        ++claimed;
    }
    return claimed;
  }

  bool FragmentClaimRegistry::isClaimed(Int trace_id) const
  {
    return claims_.find(trace_id) != claims_.end();
  }

  const FragmentClaimRegistry::ClaimRecord*
  FragmentClaimRegistry::getClaimRecord(Int trace_id) const
  {
    auto it = claims_.find(trace_id);
    if (it == claims_.end()) return nullptr;
    return &it->second;
  }

  std::vector<bool> FragmentClaimRegistry::buildExclusionMask(const MSSpectrum& spec,
                                                               const String& query_peptide_seq) const
  {
    const MSSpectrum::IntegerDataArray* trace_id_array = nullptr;
    for (const auto& arr : spec.getIntegerDataArrays())
    {
      if (arr.getName() == "fragment_trace_id")
      {
        trace_id_array = &arr;
        break;
      }
    }

    if (trace_id_array == nullptr) return {};

    std::vector<bool> mask(trace_id_array->size(), false);
    for (Size i = 0; i < trace_id_array->size(); ++i)
    {
      Int tid = (*trace_id_array)[i];
      auto it = claims_.find(tid);
      if (it != claims_.end() && it->second.peptide_seq != query_peptide_seq)
        mask[i] = true;
    }
    return mask;
  }

}  // namespace OpenMS
