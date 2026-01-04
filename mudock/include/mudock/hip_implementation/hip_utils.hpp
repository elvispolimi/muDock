#pragma once

#include <hip/hip_runtime.h>
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
#if defined(__HIP_PLATFORM_NVCC__) || defined(__NVCC__)
  #define BITLANE_MASK 0xFFFFFFFF
#elif defined(__HIP_PLATFORM_AMD__)
  #define BITLANE_MASK 0xFFFFFFFFFFFFFFFF
#endif
} // namespace mudock
