#pragma once

#include <cuda_runtime.h>
#include <mudock/cuda_implementation/cuda_utils.cuh>
#include <mudock/type_alias.hpp>

namespace mudock {
  struct cuda_texture_devices {
    fp_type* tex_dev       = nullptr;
    const int dev_id       = 0;

    cuda_texture_devices(const int dev_id_, const int xyz, const int num_tex, const fp_type* src): dev_id(dev_id_) {
      MUDOCK_CHECK(cudaSetDevice(dev_id));
      const int num_elements = xyz * num_tex;
      MUDOCK_CHECK(cudaMalloc(&tex_dev, num_elements * sizeof(fp_type)));
      MUDOCK_CHECK(cudaMemcpy(tex_dev, src, num_elements * sizeof(fp_type), cudaMemcpyHostToDevice));
      MUDOCK_CHECK(cudaDeviceSynchronize());
    };

    cuda_texture_devices(const cuda_texture_devices&)            = delete;
    cuda_texture_devices& operator=(const cuda_texture_devices&) = delete;
    cuda_texture_devices(cuda_texture_devices&& other)            = delete;
    cuda_texture_devices& operator=(cuda_texture_devices&& other) = delete;

    ~cuda_texture_devices() noexcept {
      if (tex_dev) {
        cudaSetDevice(dev_id);
        cudaFree(tex_dev);
        tex_dev = nullptr;
      }
    };
  };
} // namespace mudock
