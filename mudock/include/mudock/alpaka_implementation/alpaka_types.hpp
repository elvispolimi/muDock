#pragma once

#include <alpaka/alpaka.hpp>

#include <cstddef>

namespace mudock::alpaka_backend {
  using dim = alpaka::DimInt<1u>;
  using idx = std::size_t;

#if defined(MUDOCK_ALPAKA_BACKEND_SERIAL)
  using acc = alpaka::AccCpuSerial<dim, idx>;
#elif defined(MUDOCK_ALPAKA_BACKEND_THREADS)
  using acc = alpaka::AccCpuThreads<dim, idx>;
#elif defined(MUDOCK_ALPAKA_BACKEND_TBB)
  using acc = alpaka::AccCpuTbbBlocks<dim, idx>;
#elif defined(MUDOCK_ALPAKA_BACKEND_OMP2)
  using acc = alpaka::AccCpuOmp2Blocks<dim, idx>;
#elif defined(MUDOCK_ALPAKA_BACKEND_CUDA)
  using acc = alpaka::AccGpuCudaRt<dim, idx>;
#elif defined(MUDOCK_ALPAKA_BACKEND_HIP)
  using acc = alpaka::AccGpuHipRt<dim, idx>;
#elif defined(MUDOCK_ALPAKA_BACKEND_SYCL)
  using acc = alpaka::AccGpuSyclIntel<dim, idx>;
#else
  #error "MUDOCK_ALPAKA_BACKEND_* compile definition is required when MUDOCK_USE_ALPAKA is enabled"
#endif

  using dev_acc = alpaka::Dev<acc>;
  using queue_acc = alpaka::Queue<acc, alpaka::Blocking>;
} // namespace mudock::alpaka_backend
