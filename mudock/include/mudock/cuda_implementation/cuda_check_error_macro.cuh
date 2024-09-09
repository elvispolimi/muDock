#pragma once

#include <cuda_runtime.h>
#include <mudock/log.hpp>
#include <stdexcept>

// TODO check -Wterminate
#define MUDOCK_CHECK(call)                                                             \
  {                                                                                    \
    const cudaError_t err = call;                                                      \
    if (err != cudaSuccess) {                                                          \
      mudock::error(cudaGetErrorString(err), " in ", __FILE__, " at line ", __LINE__); \
      throw std::runtime_error("CUDA call failed, see log for details");               \
    }                                                                                  \
  }

#define MUDOCK_CHECK_KERNELCALL()                                                      \
  {                                                                                    \
    const cudaError_t err = cudaGetLastError();                                        \
    if (err != cudaSuccess) {                                                          \
      mudock::error(cudaGetErrorString(err), " in ", __FILE__, " at line ", __LINE__); \
      throw std::runtime_error("CUDA call failed, see log for details");               \
    }                                                                                  \
  }
