#pragma once

#include <hip/hip_runtime.h>
#include <mudock/hip_implementation/hip_check_error_macro.hpp>
#include <mudock/hip_implementation/hip_wrapper.hpp>

namespace mudock {
  struct hipTexture_wrapper {
    hipTextureObject_t tex;
    ~hipTexture_wrapper() noexcept(false) {
      hipResourceDesc resDesc;
      MUDOCK_CHECK(hipGetTextureObjectResourceDesc(&resDesc, tex));
      MUDOCK_CHECK(hipDestroyTextureObject(tex));
      if (resDesc.resType == hipResourceTypeArray) {
        hipArray_t mem = resDesc.res.array.array;
        MUDOCK_CHECK(hipFreeArray(mem));
      }
    };
    __host__ __device__ const hipTextureObject_t& operator()() const { return tex; };
    hipTextureObject_t& operator()() { return tex; };
  };
} // namespace mudock
