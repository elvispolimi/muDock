#include <mudock/compute.hpp>
#include <mudock/cuda_implementation/cuda_batch_sizer.cuh>
#include <mudock/cuda_implementation/cuda_manager.hpp>
#include <mudock/cuda_implementation/cuda_worker.hpp>
#include <stdexcept>

namespace mudock {

  static constexpr auto cuda_token = std::string_view{"CUDA"};

  void manage_cuda(std::string_view configuration,
                   threadpool& pool,
                   const knobs knobs,
                   std::shared_ptr<const grid_atom_mapper>& grid_atom_maps,
                   std::shared_ptr<const grid_map>& electro_map,
                   std::shared_ptr<const grid_map>& desolv_map,
                   std::shared_ptr<safe_stack<static_molecule>>& input_molecules,
                   std::shared_ptr<safe_stack<static_molecule>>& output_molecules) {
    // single out the CUDA description
    const auto begin_cuda_description = configuration.find(cuda_token);
    const auto end_cuda_description   = configuration.find(";", begin_cuda_description);
    configuration                     = configuration.substr(begin_cuda_description + cuda_token.size(),
                                         end_cuda_description - begin_cuda_description);

    // parse the CUDA description (if any)
    if (!configuration.empty()) {
      // the description should start with a colon
      if (configuration.front() != ':') [[unlikely]] {
        throw std::runtime_error(std::string{"CUDA description should start with ':' ("} +
                                 std::string{configuration} + std::string{")"});
      }
      configuration = configuration.substr(1);

      // make sure that the device is the CPU
      const auto colon_index = configuration.find(':');
      const auto device_name = configuration.substr(0, colon_index);
      if (device_name != gpu_token) [[unlikely]] {
        throw std::runtime_error(std::string{"Unsupported device '"} + std::string{device_name} +
                                 std::string{"' for the CUDA implementation"});
      }
      configuration = configuration.substr(colon_index);

      // the core counts description should start with a colon
      if (configuration.front() != ':') [[unlikely]] {
        throw std::runtime_error(std::string{"GPU count description should start with ':' ("} +
                                 std::string{configuration} + std::string{")"});
      }
      configuration = configuration.substr(1);

      // parse the ID of the gpus that we target
      const auto gpu_ids = parse_ids(configuration);

      // now we need to allocate reorder buffers for all the CUDA wrappers. In theory we can use a single
      // reorder buffer for all of them, but it can become a bottleneck. In the current implementation we
      // go for this solution, but we need to investigate better approaches
      auto rob = std::make_shared<reorder_buffer>(&compute_batch_size);

      // add the workers that we found parsing the configuration
      for (const auto id: gpu_ids) {
        // we spawn two workers for each GPU to implement the double buffer
        pool.add_worker<mudock::cuda_worker>(knobs,
                                             grid_atom_maps,
                                             electro_map,
                                             desolv_map,
                                             input_molecules,
                                             output_molecules,
                                             rob,
                                             id);
        pool.add_worker<mudock::cuda_worker>(knobs,
                                             grid_atom_maps,
                                             electro_map,
                                             desolv_map,
                                             input_molecules,
                                             output_molecules,
                                             rob,
                                             id);
        pool.add_worker<mudock::cuda_worker>(knobs,
                                             grid_atom_maps,
                                             electro_map,
                                             desolv_map,
                                             input_molecules,
                                             output_molecules,
                                             rob,
                                             id);
        pool.add_worker<mudock::cuda_worker>(knobs,
                                             grid_atom_maps,
                                             electro_map,
                                             desolv_map,
                                             input_molecules,
                                             output_molecules,
                                             rob,
                                             id);
        pool.add_worker<mudock::cuda_worker>(knobs,
                                             grid_atom_maps,
                                             electro_map,
                                             desolv_map,
                                             input_molecules,
                                             output_molecules,
                                             rob,
                                             id);
        pool.add_worker<mudock::cuda_worker>(knobs,
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
