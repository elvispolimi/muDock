#pragma once

#include <memory>
#include <mudock/compute.hpp>
#include <mudock/knobs.hpp>
#include <mudock/molecule.hpp>
#include <string_view>

namespace mudock {

  // this function will configure and create (if needed) cuda workers to the threadpool
  void manage_cuda(std::string_view configuration,
                   threadpool& pool,
                   std::shared_ptr<dynamic_molecule> protein,
                   const knobs knobs,
                   std::shared_ptr<safe_stack<static_molecule>> input_molecules,
                   std::shared_ptr<safe_stack<static_molecule>> output_molecules);
} // namespace mudock
