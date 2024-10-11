#include <mudock/compute.hpp>
#include <mudock/omp_implementation/omp_batch_sizer.hpp>
#include <mudock/omp_implementation/omp_manager.hpp>
#include <mudock/omp_implementation/omp_worker.hpp>
#include <stdexcept>

namespace mudock {

  static constexpr auto omp_token = std::string_view{"OMP"};

  void manage_omp(std::string_view configuration,
                  threadpool& pool,
                  const knobs knobs,
                  std::shared_ptr<const grid_atom_mapper>& grid_atom_maps,
                  std::shared_ptr<const grid_map>& electro_map,
                  std::shared_ptr<const grid_map>& desolv_map,
                  std::shared_ptr<safe_stack<static_molecule>>& input_molecules,
                  std::shared_ptr<safe_stack<static_molecule>>& output_molecules) {
    // single out the OpenMP description
    const auto begin_omp_description = configuration.find(omp_token);
    const auto end_omp_description   = configuration.find(";", begin_omp_description);
    configuration                    = configuration.substr(begin_omp_description + omp_token.size(),
                                         end_omp_description - begin_omp_description);

    // parse the OpenMP description (if any)
    if (!configuration.empty()) {
      // parse the ID of the gpus that we target
      const auto device_ids = parse_ids(configuration);

      // now we need to allocate reorder buffers for all the OpenMP wrappers. In theory we can use a single
      // reorder buffer for all of them, but it can become a bottleneck. In the current implementation we
      // go for this solution, but we need to investigate better approaches
      auto rob = std::make_shared<reorder_buffer>(&compute_batch_size);

      // add the workers that we found parsing the configuration
      // we spawn workers for each device
      // TODO check what to do with CPUs and GPUs
      for (const auto id: device_ids) {
        pool.add_worker<mudock::omp_worker>(knobs,
                                             grid_atom_maps,
                                             electro_map,
                                             desolv_map,
                                             input_molecules,
                                             output_molecules,
                                             rob,
                                             id);
        pool.add_worker<mudock::omp_worker>(knobs,
                                             grid_atom_maps,
                                             electro_map,
                                             desolv_map,
                                             input_molecules,
                                             output_molecules,
                                             rob,
                                             id);
      }
    }
  }
} // namespace mudock
