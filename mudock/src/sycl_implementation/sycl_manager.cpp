#include <mudock/sycl_implementation/sycl_batch_sizer.hpp>
#include <mudock/sycl_implementation/sycl_manager.hpp>
#include <mudock/sycl_implementation/sycl_worker.hpp>
#include <stdexcept>
#include <sycl/sycl.hpp>

namespace mudock {

  static constexpr auto sycl_token = std::string_view{"SYCL"};

  void manage_sycl(std::string_view configuration,
                   threadpool& pool,
                   const knobs knobs,
                   std::shared_ptr<const grid_atom_mapper>& grid_atom_maps,
                   std::shared_ptr<const grid_map>& electro_map,
                   std::shared_ptr<const grid_map>& desolv_map,
                   std::shared_ptr<safe_stack<static_molecule>>& input_molecules,
                   std::shared_ptr<safe_stack<static_molecule>>& output_molecules) {
    // single out the SYCL description
    const auto begin_sycl_description = configuration.find(sycl_token);
    const auto end_sycl_description   = configuration.find(";", begin_sycl_description);
    configuration                     = configuration.substr(begin_sycl_description + sycl_token.size(),
                                         end_sycl_description - begin_sycl_description);

    // parse the SYCL description (if any)
    if (!configuration.empty()) {
      // the description should start with a colon
      if (configuration.front() != ':') [[unlikely]] {
        throw std::runtime_error(std::string{"SYCL description should start with ':' ("} +
                                 std::string{configuration} + std::string{")"});
      }
      configuration = configuration.substr(1);

      // make sure that the device is the CPU
      const auto colon_index = configuration.find(':');
      const auto device_name = configuration.substr(0, colon_index);
      if (device_name != gpu_token || device_name != cpu_token) [[unlikely]] {
        throw std::runtime_error(std::string{"Unsupported device '"} + std::string{device_name} +
                                 std::string{"' for the SYCL implementation"});
      }
      configuration = configuration.substr(colon_index);

      // the core counts description should start with a colon
      if (configuration.front() != ':') [[unlikely]] {
        throw std::runtime_error(std::string{"Device count description should start with ':' ("} +
                                 std::string{configuration} + std::string{")"});
      }
      configuration = configuration.substr(1);

      // parse the ID of the devices that we target
      const auto device_ids = parse_ids(configuration);

      // Select the correct device
      // Get a list of available devices
      std::vector<sycl::device> devices = sycl::device::get_devices();
      // Vector to hold devices of the requested type
      std::vector<sycl::device> filtered_devices;
      // Filter devices based on the input type (CPU or device)
      for (const auto& dev: devices) {
        if (device_name == cpu_token && dev.is_cpu()) {
          filtered_devices.push_back(dev);
        } else if (device_name == gpu_token && dev.is_gpu()) {
          filtered_devices.push_back(dev);
        }
      }
      // Check if the requested device index is valid
      if (filtered_devices.empty()) {
        throw std::runtime_error(std::string{"No devices of type "} + std::string{device_name} +
                                 std::string{" found."});
      } else {
        if (std::all_of(device_ids.begin(), device_ids.end(), [&](const int num) {
              return num >= 0 && num <= static_cast<int>(filtered_devices.size());
            }))
          throw std::runtime_error(std::string{"Invalid SYCL device numbers."});
      }

      // now we need to allocate reorder buffers for all the SYCL wrappers. In theory we can use a single
      // reorder buffer for all of them, but it can become a bottleneck. In the current implementation we
      // go for this solution, but we need to investigate better approaches
      auto rob = std::make_shared<reorder_buffer>(&compute_batch_size);

      // add the workers that we found parsing the configuration
      for (const auto id: device_ids) {
        // we spawn two workers for each GPU to implement the double buffer
        pool.add_worker<mudock::sycl_worker>(knobs,
                                             grid_atom_maps,
                                             electro_map,
                                             desolv_map,
                                             input_molecules,
                                             output_molecules,
                                             rob,
                                             id,
                                             filtered_devices[id]);
        pool.add_worker<mudock::sycl_worker>(knobs,
                                             grid_atom_maps,
                                             electro_map,
                                             desolv_map,
                                             input_molecules,
                                             output_molecules,
                                             rob,
                                             id,
                                             filtered_devices[id]);
      }
    }
  }
} // namespace mudock
