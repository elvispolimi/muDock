#pragma once

#include <hip/hip_runtime.h>
#include <mudock/compute/batch_multiple.hpp>
#include <mudock/log.hpp>
#include <stdexcept>

// TODO check -Wterminate
#define MUDOCK_CHECK(call)                                                            \
  {                                                                                   \
    const hipError_t err = call;                                                      \
    if (err != hipSuccess) {                                                          \
      mudock::error(hipGetErrorString(err), " in ", __FILE__, " at line ", __LINE__); \
      throw std::runtime_error("HIP call failed, see log for details");               \
    }                                                                                 \
  }

#define MUDOCK_CHECK_KERNELCALL()                                                     \
  {                                                                                   \
    const hipError_t err = hipGetLastError();                                         \
    if (err != hipSuccess) {                                                          \
      mudock::error(hipGetErrorString(err), " in ", __FILE__, " at line ", __LINE__); \
      throw std::runtime_error("HIP kernel call failed, see log for details");        \
    }                                                                                 \
  }

namespace mudock {
  template<auto Kernel>
  inline batch_multiple get_kernel_batch_multiple_hip(const int device_id,
                                                      const int block_size,
                                                      const size_t dynamic_shared_mem = 0,
                                                      const char* kernel_label = "unknown_kernel") {
    MUDOCK_CHECK(hipSetDevice(device_id));
    hipDeviceProp_t props;
    MUDOCK_CHECK(hipGetDeviceProperties(&props, device_id));

    int num_blocks_per_sm = 0;
    MUDOCK_CHECK(
        hipOccupancyMaxActiveBlocksPerMultiprocessor(&num_blocks_per_sm, Kernel, block_size, dynamic_shared_mem));
    (void) kernel_label;
    return {num_blocks_per_sm, props.multiProcessorCount};
  }

#if defined(__HIP_PLATFORM_NVCC__) || defined(__NVCC__)
  #define BITLANE_MASK 0xFFFFFFFF
#elif defined(__HIP_PLATFORM_AMD__)
  #define BITLANE_MASK 0xFFFFFFFFFFFFFFFF
#endif
#if defined(__HIP_PLATFORM_AMD__)
  #define SHFL_DOWN(mask, v, delta, width) __shfl_down((v), (delta), (width))
  #define SHFL(mask, v, src, width)        __shfl((v), (src), (width))
#else
  #define SHFL_DOWN(mask, v, delta, width) __shfl_down_sync((mask), (v), (delta), (width))
  #define SHFL(mask, v, src, width)        __shfl_sync((mask), (v), (src), (width))
#endif
} // namespace mudock
