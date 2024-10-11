#pragma once

#include <mudock/knobs.hpp>
#ifdef MUDOCK_USE_OMP
  #include <mudock/omp_implementation/omp_batch_sizer.hpp>
  // #include <mudock/cuda_implementation/cuda_check_error_macro.cuh>
  #include <mudock/omp_implementation/omp_manager.hpp>
  // #include <mudock/cuda_implementation/cuda_object.cuh>
  // #include <mudock/cuda_implementation/cuda_wrapper.cuh>
  #include <mudock/omp_implementation/virtual_screen.hpp>
#else
  #include <mudock/log.hpp>
namespace mudock {
  inline void manage_omp(std::string_view,
                         threadpool&,
                         const knobs,
                         [[maybe_unused]] std::shared_ptr<const grid_atom_mapper> grid_atom_maps,
                         [[maybe_unused]] std::shared_ptr<const grid_map> electro_map,
                         [[maybe_unused]] std::shared_ptr<const grid_map> desolv_map,
                         std::shared_ptr<safe_stack<static_molecule> >,
                         std::shared_ptr<safe_stack<static_molecule> >) {
    warning("The OpenMP implementation is disabled");
  }
} // namespace mudock
#endif
