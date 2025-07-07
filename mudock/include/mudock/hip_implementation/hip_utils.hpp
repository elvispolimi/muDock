#pragma once

namespace mudock {
#if defined(__HIP_PLATFORM_NVCC__) || defined(__NVCC__)
  #define BITLANE_MASK 0xFFFFFFFF
#elif defined(__HIP_PLATFORM_AMD__)
  #define BITLANE_MASK 0xFFFFFFFFFFFFFFFF
#endif
} // namespace mudock
