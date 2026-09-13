// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>

#include <algorithm>
#include <atomic>
#include <thread>
#include <vector>

namespace OpenMS
{
  /**
    @brief Runs @p body(state, i) for i in [0, n), across up to @p num_threads
           std::thread workers, joined before this function returns.

    Each worker constructs its own thread-local state once (via @p make_state)
    and then repeatedly claims the next unclaimed index from a shared atomic
    counter, approximating OpenMP's `schedule(dynamic, 1)` load balancing.

    Unlike a nested `#pragma omp parallel` region, the synchronization here
    (std::thread::join()) is a plain-C++-guaranteed happens-before edge, not
    dependent on the OpenMP runtime's own (nested) team-join implementation.
    Use this instead of nested OpenMP when a caller needs a reliable "wait
    until this specific batch of workers is done" guarantee from within code
    that may itself already be running inside an OpenMP parallel region.

    @param n Number of work items
    @param num_threads Worker thread budget (clamped to [1, n]; 1 runs inline, no threads spawned)
    @param make_state Nullary callable returning one worker's thread-local state
    @param body Callable taking (state&, Size i) for one work item
  */
  template <typename StateFactory, typename Body>
  void parallelForWithState(Size n, Size num_threads, StateFactory&& make_state, Body&& body)
  {
    if (n == 0) return;

    const Size n_threads = std::max<Size>(1, std::min<Size>(num_threads, n));

    if (n_threads == 1)
    {
      // no threads spawned; caller's own thread does the work
      auto state = make_state();
      for (Size i = 0; i < n; ++i) body(state, i);
      return;
    }

    std::atomic<Size> next{0};
    auto run_worker = [&]()
    {
      auto state = make_state();
      Size i;
      while ((i = next.fetch_add(1, std::memory_order_relaxed)) < n)
      {
        body(state, i);
      }
    };

    std::vector<std::thread> workers;
    workers.reserve(n_threads);
    for (Size t = 0; t < n_threads; ++t)
    {
      workers.emplace_back(run_worker);
    }
    for (auto& w : workers) w.join();
  }

  /// Same as parallelForWithState(), for callers that need no per-thread state.
  template <typename Body>
  void parallelFor(Size n, Size num_threads, Body&& body)
  {
    parallelForWithState(n, num_threads, [](){ return 0; },
                          [&body](int&, Size i) { body(i); });
  }

} // namespace OpenMS
