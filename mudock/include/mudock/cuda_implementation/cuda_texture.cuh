#pragma once

#include <cuda_runtime.h>
#include <mudock/cuda_implementation/cuda_check_error_macro.cuh>
#include <mudock/cuda_implementation/cuda_wrapper.cuh>

namespace mudock {
  struct cudaTexture_wrapper {
    cudaTextureObject_t tex;
    ~cudaTexture_wrapper() noexcept(false) {
      cudaResourceDesc resDesc;
      MUDOCK_CHECK(cudaGetTextureObjectResourceDesc(&resDesc, tex));
      MUDOCK_CHECK(cudaDestroyTextureObject(tex));
      if (resDesc.resType == cudaResourceTypeArray) {
        cudaArray_t mem = resDesc.res.array.array;
        MUDOCK_CHECK(cudaFreeArray(mem));
      }
    };
    __host__ __device__ const cudaTextureObject_t& operator()() const { return tex; };
    cudaTextureObject_t& operator()() { return tex; };
  };
} // namespace mudock
