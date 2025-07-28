#pragma once

#include <hip/hip_runtime.h>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/hip_implementation/hip_check_error_macro.hpp>
#include <mudock/hip_implementation/hip_wrapper.hpp>

namespace mudock {
  struct hipTexture_wrapper {
    hip_wrapper<std::vector, fp_type> tex;
    ~hipTexture_wrapper() = default;
    const hip_wrapper<std::vector, fp_type>& operator()() const { return tex; };
    hip_wrapper<std::vector, fp_type>& operator()() { return tex; };
    hipTexture_wrapper(hipStream_t& stream, const autodock_protein& adt_protein): tex(stream) {
      const fp_type* grid_map = adt_protein.get_maps_pointer();
      const int num_elements  = adt_protein.get_map_flat_size() * num_autodock_grids();
      tex.alloc(num_elements);

      std::memcpy(tex(), grid_map, num_elements * sizeof(fp_type));

      tex.copy_host2device();
    }
  };
} // namespace mudock
