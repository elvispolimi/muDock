#pragma once

#include <cuda_runtime.h>
#include <mudock/compute/batch_multiple.hpp>
#include <mudock/log.hpp>
#include <stdexcept>

// TODO check -Wterminate
// TODO add check if DEBUG mode is enable
// TODO add better checks if function invoked after driver teardown
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

namespace mudock {
  template<auto Kernel>
  inline batch_multiple get_kernel_batch_multiple_cuda(const int device_id,
                                                       const int block_size,
                                                       const size_t dynamic_shared_mem = 0,
                                                       const char* kernel_label = "unknown_kernel") {
    MUDOCK_CHECK(cudaSetDevice(device_id));
    cudaDeviceProp props;
    MUDOCK_CHECK(cudaGetDeviceProperties(&props, device_id));

    int num_blocks_per_sm = 0;
    MUDOCK_CHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&num_blocks_per_sm,
                                                               Kernel,
                                                               block_size,
                                                               dynamic_shared_mem));
    (void) kernel_label;
    return {num_blocks_per_sm, props.multiProcessorCount};
  }
} // namespace mudock
