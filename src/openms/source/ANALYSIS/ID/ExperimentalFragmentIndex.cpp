// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: OpenMS Developers $
// $Authors: OpenMS Developers $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/ExperimentalFragmentIndex.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

namespace OpenMS
{

  const ExperimentalFragmentIndex::SpectrumIdxList ExperimentalFragmentIndex::empty_;

  void ExperimentalFragmentIndex::build(const MSExperiment& pseudo_spectra)
  {
    index_.clear();

    for (Size spec_idx = 0; spec_idx < pseudo_spectra.size(); ++spec_idx)
    {
      const MSSpectrum& spec = pseudo_spectra[spec_idx];
      const auto& int_arrays = spec.getIntegerDataArrays();

      // Locate the fragment_trace_id array
      const MSSpectrum::IntegerDataArray* trace_id_array = nullptr;
      for (const auto& arr : int_arrays)
      {
        if (arr.getName() == "fragment_trace_id")
        {
          trace_id_array = &arr;
          break;
        }
      }

      if (trace_id_array == nullptr) continue;

      for (Int trace_id : *trace_id_array)
      {
        index_[trace_id].push_back(spec_idx);
      }
    }
  }

  const ExperimentalFragmentIndex::SpectrumIdxList&
  ExperimentalFragmentIndex::getSpectraForTrace(Int trace_id) const
  {
    auto it = index_.find(trace_id);
    if (it == index_.end()) return empty_;
    return it->second;
  }

}  // namespace OpenMS
