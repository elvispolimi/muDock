#pragma once

#include <cassert>
#include <mudock/alpaka_implementation/queue_alpaka.hpp>
#include <type_traits>
#include <utility>

namespace mudock {
  template<class F, class... Args>
  inline void queue_alpaka::invoke_kernel(const index3D gridDim, const index3D blockDim, Args&&... args) {
    assert(gridDim.size_x() > 0 && blockDim.size_x() > 0);
    assert(gridDim.size_y() == 1 && blockDim.size_y() == 1 &&
           "queue_alpaka currently supports only 1D launches");
    assert(gridDim.size_z() == 1 && blockDim.size_z() == 1 &&
           "queue_alpaka currently supports only 1D launches");

    const auto grid  = alpaka::Vec<dim, idx>{static_cast<idx>(gridDim.size_x())};
    const auto block = alpaka::Vec<dim, idx>{static_cast<idx>(block_threads(blockDim.size_x()))};
    const auto elems = alpaka::Vec<dim, idx>{idx{1}};

    const auto work_div = alpaka::WorkDivMembers<dim, idx>{grid, block, elems};
    auto task           = alpaka::createTaskKernel<acc>(work_div,
                                              F{},
                                              static_cast<std::decay_t<Args>>(std::forward<Args>(args))...);
    alpaka::enqueue(native_queue(), task);
  }

  template<class F, class... Args>
  inline void queue_alpaka::invoke_kernel(const int gridDim, const int blockDim, Args&&... args) {
    assert(gridDim >= 0 && blockDim >= 0);
    invoke_kernel<F>(index3D{gridDim, 1, 1}, index3D{blockDim, 1, 1}, std::forward<Args>(args)...);
  }

  inline constexpr int queue_alpaka::block_threads(const int requested_threads) {
#if defined(MUDOCK_ALPAKA_BACKEND_SERIAL) || defined(MUDOCK_ALPAKA_BACKEND_TBB) || \
    defined(MUDOCK_ALPAKA_BACKEND_OMP2)
    (void) requested_threads;
    return 1;
#else
    return requested_threads;
#endif
  }
} // namespace mudock
