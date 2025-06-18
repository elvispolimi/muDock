#pragma once

#include <cuda_runtime.h>
#include <memory>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/cuda_implementation/cuda_wrapper.cuh>
#include <mudock/grid.hpp>

namespace mudock {

  struct device {
    // Device ID
    const std::size_t id;
    // Grid Maps
    const point<fp_type, 3> center_maps;
    const cudaStream_t stream;
    cudaTextureObject_t electro_tex, desolv_tex;
    cuda_wrapper<std::vector, cudaTextureObject_t> atom_texs;
    const autodock_protein& adt_protein;

    device(const std::size_t gpu_id, const autodock_protein& adt_protein);

    cudaStream_t get_stream() const;
  };
} // namespace mudock
