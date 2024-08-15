#pragma once

#include <memory>
#include <mudock/compute.hpp>
#include <mudock/grid.hpp>
#include <mudock/molecule.hpp>
#include <string_view>

namespace mudock {

  // this function will configure and create (if needed) cpp workers to the threadpool
  void manage_cpp(std::string_view configuration,
                  threadpool& pool,
                  std::shared_ptr<const dynamic_molecule> protein,
                  std::shared_ptr<const grid_atom_mapper> grid_atom_maps,
                  std::shared_ptr<const grid_map> electro_map,
                  std::shared_ptr<const grid_map> desolv_map,
                  std::shared_ptr<safe_stack<static_molecule>> input_molecules,
                  std::shared_ptr<safe_stack<static_molecule>> output_molecules);
} // namespace mudock
