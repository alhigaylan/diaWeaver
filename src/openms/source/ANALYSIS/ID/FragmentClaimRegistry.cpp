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

  Size FragmentClaimRegistry::tryClaim(const std::vector<TraceKey>& keys,
                                        const String& peptide_seq,
                                        double score,
                                        Size spectrum_idx)
  {
    Size claimed = 0;
    for (TraceKey key : keys)
    {
      if (claims_.emplace(key, ClaimRecord{peptide_seq, score, spectrum_idx}).second)
        ++claimed;
    }
    return claimed;
  }

  bool FragmentClaimRegistry::isClaimed(TraceKey key) const
  {
    return claims_.find(key) != claims_.end();
  }

  const FragmentClaimRegistry::ClaimRecord*
  FragmentClaimRegistry::getClaimRecord(TraceKey key) const
  {
    auto it = claims_.find(key);
    if (it == claims_.end()) return nullptr;
    return &it->second;
  }

  std::vector<bool> FragmentClaimRegistry::buildExclusionMask(const MSSpectrum& spec,
                                                               const String& query_peptide_seq) const
  {
    const MSSpectrum::IntegerDataArray* trace_arr  = nullptr;
    const MSSpectrum::IntegerDataArray* window_arr = nullptr;
    for (const auto& arr : spec.getIntegerDataArrays())
    {
      if      (arr.getName() == "fragment_trace_id")  trace_arr  = &arr;
      else if (arr.getName() == "fragment_window_id") window_arr = &arr;
    }

    if (trace_arr == nullptr || window_arr == nullptr) return {};

    const Size n = trace_arr->size();
    std::vector<bool> mask(n, false);
    for (Size i = 0; i < n; ++i)
    {
      TraceKey key = makeKey((*window_arr)[i], (*trace_arr)[i]);
      auto it = claims_.find(key);
      if (it != claims_.end() && it->second.peptide_seq != query_peptide_seq)
        mask[i] = true;
    }
    return mask;
  }

}  // namespace OpenMS
