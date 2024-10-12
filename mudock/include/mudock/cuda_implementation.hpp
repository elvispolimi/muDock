#pragma once

#include <mudock/knobs.hpp>
#ifdef MUDOCK_USE_POLY || MUDOCK_USE_POLY
  #include <mudock/cuda_implementation/cuda_batch_sizer.cuh>
  #include <mudock/cuda_implementation/cuda_check_error_macro.cuh>
  #include <mudock/cuda_implementation/cuda_manager.hpp>
  #include <mudock/cuda_implementation/cuda_object.cuh>
  #include <mudock/cuda_implementation/cuda_wrapper.cuh>
  #include <mudock/cuda_implementation/virtual_screen.cuh>
#else
  #include <mudock/log.hpp>
namespace mudock {
  inline void manage_cuda(std::string_view,
                          threadpool&,
                          const knobs,
                          [[maybe_unused]] std::shared_ptr<const grid_atom_mapper> grid_atom_maps,
                          [[maybe_unused]] std::shared_ptr<const grid_map> electro_map,
                          [[maybe_unused]] std::shared_ptr<const grid_map> desolv_map,
                          std::shared_ptr<safe_stack<static_molecule> >,
                          std::shared_ptr<safe_stack<static_molecule> >) {
    warning("The CUDA implementation is disabled");
  }
} // namespace mudock
#endif
