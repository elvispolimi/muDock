#pragma once

#include <memory>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/compute.hpp>
#include <mudock/grid.hpp>
#include <mudock/knobs.hpp>
#include <mudock/molecule.hpp>

namespace mudock {

  // this function will configure and create (if needed) cuda workers to the threadpool
  void manage_cuda(const std::vector<std::string>& configurations,
                   threadpool& pool,
                   const knobs knobs,
                   const autodock_protein& adt_protein,
                   std::shared_ptr<safe_stack<autodock_ligand>>& input_molecules,
                   std::shared_ptr<safe_stack<static_molecule>>& output_molecules);
} // namespace mudock
