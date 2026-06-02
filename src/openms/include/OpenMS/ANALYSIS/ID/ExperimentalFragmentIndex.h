// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: OpenMS Developers $
// $Authors: OpenMS Developers $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <unordered_map>
#include <vector>

namespace OpenMS
{

  /**
   * @brief Maps MassTrace-derived fragment IDs to spectrum indices within a DIA window.
   *
   * Built from the pseudo-spectra MSExperiment produced by ClusterMassTracesByPrecursor,
   * using the "fragment_trace_id" IntegerDataArray emitted per spectrum.
   *
   * Because DIA windows are m/z-disjoint, trace IDs are local integers (positions in the
   * window's ms2_traces vector) and never collide across windows. No global encoding is needed.
   *
   * Usage:
   * @code
   *   ExperimentalFragmentIndex efi;
   *   efi.build(pseudo_spectra);
   *   const auto& spectra = efi.getSpectraForTrace(42);  // {spec_idx, ...}
   * @endcode
   */
  class OPENMS_DLLAPI ExperimentalFragmentIndex
  {
  public:
    /// spectrum_idx values for all pseudo spectra containing a given trace_id
    using SpectrumIdxList = std::vector<Size>;

    /// Build the index from pseudo spectra that carry "fragment_trace_id" IntegerDataArrays.
    void build(const MSExperiment& pseudo_spectra);

    /// Return all spectrum indices that contain fragment with the given trace_id.
    /// Returns an empty range if the trace_id is unknown.
    const SpectrumIdxList& getSpectraForTrace(Int trace_id) const;

    /// Total number of unique trace IDs registered.
    Size size() const noexcept { return index_.size(); }

    void clear() noexcept { index_.clear(); }

  private:
    std::unordered_map<Int, SpectrumIdxList> index_;
    static const SpectrumIdxList empty_;
  };

}  // namespace OpenMS
