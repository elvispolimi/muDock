#pragma once

#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/grid/space_grid.hpp>
#include <mudock/sycl_implementation/sycl_wrapper.hpp>
#include <mudock/type_alias.hpp>
// #include <sycl/image.hpp>
#include <sycl/sycl.hpp>
#include <vector>

// TODO oneAPI images
namespace mudock {
  struct syclTexture_wrapper {
    // sycl::image<3> tex;
    sycl_wrapper<std::vector, fp_type> tex;
    ~syclTexture_wrapper() = default;
    const sycl_wrapper<std::vector, fp_type>& operator()() const { return tex; };
    sycl_wrapper<std::vector, fp_type>& operator()() { return tex; };
    syclTexture_wrapper(sycl::queue& queue, const autodock_protein& adt_protein): tex(queue) {
      const fp_type* grid_map = adt_protein.get_maps_pointer();
      const int num_elements  = adt_protein.get_map_flat_size() * num_autodock_grids();
      tex.alloc(num_elements);

      std::memcpy(tex.host_pointer(), grid_map, num_elements * sizeof(fp_type));

      tex.copy_host2device();
    }
  };
} // namespace mudock
