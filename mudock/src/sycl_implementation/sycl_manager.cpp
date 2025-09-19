#include <mudock/sycl_implementation/sycl_batch_sizer.hpp>
#include <mudock/sycl_implementation/sycl_manager.hpp>
#include <mudock/sycl_implementation/sycl_worker.hpp>
#include <stdexcept>
#include <sycl/sycl.hpp>

namespace mudock {

  static constexpr auto sycl_token = std::string_view{"SYCL"};

  void manage_sycl(const std::vector<std::string>& configurations,
                   threadpool& pool,
                   const knobs knobs,
                   const autodock_protein& adt_protein,
                   std::shared_ptr<safe_stack<static_molecule>>& input_molecules,
                   std::shared_ptr<safe_stack<static_molecule>>& output_molecules) {
    // single out the SYCL description
    const auto it =
        std::find_if(configurations.begin(), configurations.end(), [](const std::string_view& str) {
          return str.find(sycl_token) != std::string::npos; // Check if the target is a substring
        });

    // parse the SYCL description (if any)
    if (it != configurations.end()) {
      auto configuration = *it;

      configuration = configuration.substr(sycl_token.size());

      // the description should start with a colon
      if (configuration.front() != ':') [[unlikely]] {
        throw std::runtime_error(std::string{"SYCL description should start with ':' ("} +
                                 std::string{configuration} + std::string{")"});
      }
      configuration = configuration.substr(1);

      // make sure that the device is supported
      const auto colon_index = configuration.find(':');
      const auto device_name = configuration.substr(0, colon_index);
      if (device_name != gpu_token && device_name != cpu_token) [[unlikely]] {
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
      std::vector<sycl::device> gpu_devices;
      // TODO it can happen?
      std::optional<sycl::device> cpu_device;
      // Filter devices based on the input type (CPU or device)
      for (const auto& dev: devices) {
        if (device_name == cpu_token && dev.is_cpu()) {
          cpu_device.emplace(dev);
        } else if (device_name == gpu_token && dev.is_gpu()) {
          gpu_devices.push_back(dev);
        }
      }
      // Check if the requested device index is valid
      sycl::device dev_ook;
      if (device_name == gpu_token) {
        if (gpu_devices.empty()) {
          throw std::runtime_error(std::string{"No devices of type "} + std::string{device_name} +
                                   std::string{" found."});
        } else {
          if (!std::all_of(device_ids.begin(), device_ids.end(), [&](const int num) {
                return num >= 0 && num <= static_cast<int>(gpu_devices.size());
              }))
            throw std::runtime_error(std::string{"Invalid SYCL device numbers."});
        }
        dev_ook = gpu_devices[0];
      } else if (device_name == cpu_token) {
        if (!cpu_device.has_value())
          throw std::runtime_error(std::string{"No devices of type "} + std::string{device_name} +
                                   std::string{" found."});
        dev_ook = cpu_device.value();
      }

      // now we need to allocate reorder buffers for all the SYCL wrappers. In theory we can use a single
      // reorder buffer for all of them, but it can become a bottleneck. In the current implementation we
      // go for this solution, but we need to investigate better approaches
      auto rob = std::make_shared<reorder_buffer>(
          [dev_ook](const int num_atoms) { return compute_batch_size(dev_ook, num_atoms); });

      // add the workers that we found parsing the configuration
      for (const auto id: device_ids) {
        if (device_name == gpu_token) {
          const auto dev = std::make_shared<device>(gpu_devices[id], adt_protein);

          // we spawn two workers for each GPU to implement the double buffer
          pool.add_worker<mudock::sycl_worker>(knobs, input_molecules, output_molecules, rob, dev);

          pool.add_worker<mudock::sycl_worker>(knobs, input_molecules, output_molecules, rob, dev);
        } else if (device_name == cpu_token) {
          const auto dev = std::make_shared<device>(cpu_device.value(), adt_protein);

          pool.add_worker<mudock::sycl_worker>(knobs, input_molecules, output_molecules, rob, dev);
        }
      }
    }
  }
} // namespace mudock
