#include <mudock/hip_implementation/hip_batch_sizer.hpp>
#include <mudock/hip_implementation/hip_manager.hpp>
#include <mudock/hip_implementation/hip_worker.hpp>
#include <stdexcept>

namespace mudock {

  static constexpr auto hip_token = std::string_view{"HIP"};

  void manage_hip(std::string_view configuration,
                  threadpool& pool,
                  const knobs knobs,
                  std::shared_ptr<const grid_atom_mapper>& grid_atom_maps,
                  std::shared_ptr<const grid_map>& electro_map,
                  std::shared_ptr<const grid_map>& desolv_map,
                  std::shared_ptr<safe_stack<static_molecule>>& input_molecules,
                  std::shared_ptr<safe_stack<static_molecule>>& output_molecules) {
    // single out the HIP description
    const auto begin_hip_description = configuration.find(hip_token);
    const auto end_hip_description   = configuration.find(";", begin_hip_description);
    configuration                    = configuration.substr(begin_hip_description + hip_token.size(),
                                         end_hip_description - begin_hip_description);

    // parse the HIP description (if any)
    if (!configuration.empty()) {
      // the description should start with a colon
      if (configuration.front() != ':') [[unlikely]] {
        throw std::runtime_error(std::string{"HIP description should start with ':' ("} +
                                 std::string{configuration} + std::string{")"});
      }
      configuration = configuration.substr(1);

      // make sure that the device is the GPU
      // TODO extend to other devices
      const auto colon_index = configuration.find(':');
      const auto device_name = configuration.substr(0, colon_index);
      if (device_name != gpu_token) [[unlikely]] {
        throw std::runtime_error(std::string{"Unsupported device '"} + std::string{device_name} +
                                 std::string{"' for the HIP implementation"});
      }
      configuration = configuration.substr(colon_index);

      // the core counts description should start with a colon
      if (configuration.front() != ':') [[unlikely]] {
        throw std::runtime_error(std::string{"Device count description should start with ':' ("} +
                                 std::string{configuration} + std::string{")"});
      }
      configuration = configuration.substr(1);

      // parse the ID of the gpus that we target
      const auto gpu_ids = parse_ids(configuration);

      // now we need to allocate reorder buffers for all the HIP wrappers. In theory we can use a single
      // reorder buffer for all of them, but it can become a bottleneck. In the current implementation we
      // go for this solution, but we need to investigate better approaches
      auto rob = std::make_shared<reorder_buffer>(&compute_batch_size);

      // add the workers that we found parsing the configuration
      for (const auto id: gpu_ids) {
        // we spawn two workers for each GPU to implement the double buffer
        pool.add_worker<mudock::hip_worker>(knobs,
                                            grid_atom_maps,
                                            electro_map,
                                            desolv_map,
                                            input_molecules,
                                            output_molecules,
                                            rob,
                                            id);
        pool.add_worker<mudock::hip_worker>(knobs,
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
