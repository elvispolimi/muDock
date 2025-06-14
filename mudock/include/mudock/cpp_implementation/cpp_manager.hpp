#pragma once

#include <memory>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/compute.hpp>
#include <mudock/grid.hpp>
#include <mudock/grid/grid_map.hpp>
#include <mudock/knobs.hpp>
#include <mudock/molecule.hpp>

namespace mudock {

  // this function will configure and create (if needed) cpp workers to the threadpool
  void manage_cpp(const std::vector<std::string>& configurations,
                  threadpool& pool,
                  const autodock_protein& adt_protein,
                  const knobs knobs,
                  std::shared_ptr<safe_stack<static_molecule>>& input_molecules,
                  std::shared_ptr<safe_stack<static_molecule>>& output_molecules);
} // namespace mudock
