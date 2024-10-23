#pragma once

#include <cuda_runtime.h>
#include <memory>
#include <mudock/cuda_implementation/cuda_wrapper.cuh>
#include <mudock/grid.hpp>

namespace mudock {

  struct device {
    // Device ID
    const std::size_t id;
    // Grid Maps
    const point3D center_maps;
    const cudaStream_t stream;
    cudaTextureObject_t electro_tex, desolv_tex;
    cuda_wrapper<std::vector, cudaTextureObject_t> atom_texs;

    device(const std::size_t gpu_id,
           std::shared_ptr<const grid_atom_mapper>& grid_atom_maps,
           std::shared_ptr<const grid_map>& electro_map,
           std::shared_ptr<const grid_map>& desolv_map);

    cudaStream_t get_stream() const;
  };
} // namespace mudock
