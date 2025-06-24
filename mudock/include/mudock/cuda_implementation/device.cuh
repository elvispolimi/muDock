#pragma once

#include <cuda_runtime.h>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/cuda_implementation/cuda_check_error_macro.cuh>
#include <mudock/cuda_implementation/cuda_texture.cuh>
#include <mudock/cuda_implementation/cuda_wrapper.cuh>
#include <mudock/grid.hpp>

namespace mudock {
  struct device {
    // Device ID
    const std::size_t id;
    // Grid Maps
    const point<fp_type, 3> center_maps;
    struct cudaStream_wrapper {
      cudaStream_t stream;
      cudaStream_wrapper(const cudaStream_t stream): stream(stream) {};
      ~cudaStream_wrapper() noexcept(false) { MUDOCK_CHECK(cudaStreamDestroy(stream)) };
      cudaStream_t& operator()() { return stream; };
    } stream;

    cuda_wrapper<std::vector, cudaTexture_wrapper> atom_tex;
    const autodock_protein& adt_protein;

    device(const std::size_t gpu_id, const autodock_protein& adt_protein);

    cudaStream_t get_stream() const;
  };
} // namespace mudock
