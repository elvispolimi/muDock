#pragma once

#include <cuda_runtime.h>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/cuda_implementation/cuda_check_error_macro.cuh>
#include <mudock/cuda_implementation/cuda_texture.cuh>
#include <mudock/cuda_implementation/cuda_wrapper.cuh>
#include <mudock/grid.hpp>

namespace mudock {
  struct cudaStream_wrapper {
    cudaStream_wrapper(const int id) {
      MUDOCK_CHECK(cudaSetDevice(static_cast<int>(id)));
      MUDOCK_CHECK(cudaStreamCreate(&stream););
    };
    ~cudaStream_wrapper() noexcept(false) { MUDOCK_CHECK(cudaStreamDestroy(stream)) };
    cudaStream_t& operator()() { return stream; };

  private:
    cudaStream_t stream;
  };

  struct device {
    // Device ID
    const std::size_t id;
    // Grid Maps
    const point<fp_type, 3> center_maps;

    const autodock_protein& adt_protein;

    device(const std::size_t gpu_id, const autodock_protein& adt_protein);

    cudaStream_wrapper get_stream() const;
    const auto* get_tex_dev_pointer() const { return atom_tex.dev_pointer(); };

  private:
    cudaStream_wrapper stream;
    cuda_wrapper<std::vector, cudaTexture_wrapper> atom_tex;
  };
} // namespace mudock
