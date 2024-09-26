#pragma once

#include <mudock/batch.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/grid.hpp>
#include <mudock/knobs.hpp>

namespace mudock {

  class virtual_screen_sycl {
    // the configuration of the GA algorithm
    knobs configuration;

    // Grid Maps
    const index3D index_maps;
    const point3D center_maps;

  public:
    virtual_screen_sycl(const knobs k,
                        std::shared_ptr<const grid_atom_mapper>& grid_atom_maps,
                        std::shared_ptr<const grid_map>& electro_map,
                        std::shared_ptr<const grid_map>& desolv_map);

    void operator()(batch& incoming_batch);
  };
} // namespace mudock
