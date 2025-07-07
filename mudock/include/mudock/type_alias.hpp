#pragma once

#include <type_traits>
namespace mudock {

  using fp_type = float;
#if defined(MUDOCK_USE_HIP)
  static_assert(std::is_same_v<fp_type, float>, "HIP supports only float data type");
#endif
#if defined(MUDOCK_USE_CUDA)
  static_assert(std::is_same_v<fp_type, float>, "CUDA supports only float data type");
#endif
#if defined(MUDOCK_USE_XSIMD)
  static_assert(std::is_same_v<fp_type, float>, "XSIMD does not supports double data type");
#endif
#if defined(MUDOCK_USE_GH)
  static_assert(std::is_same_v<fp_type, float>, "GH does not supports double data type");
#endif
} // namespace mudock
