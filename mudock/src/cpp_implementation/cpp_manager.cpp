#include <mudock/compute.hpp>
#include <mudock/cpp_implementation/cpp_manager.hpp>
#include <mudock/cpp_implementation/cpp_worker.hpp>
#include <stdexcept>

namespace mudock {

  static constexpr auto cpp_token = std::string_view{"CPP"};

  void manage_cpp(const std::vector<std::string>& configurations,
                  threadpool& pool,
                  std::shared_ptr<const grid_atom_mapper>& grid_atom_maps,
                  std::shared_ptr<const grid_map>& electro_map,
                  std::shared_ptr<const grid_map>& desolv_map,
                  const knobs knobs,
                  std::shared_ptr<safe_stack<static_molecule>>& input_molecules,
                  std::shared_ptr<safe_stack<static_molecule>>& output_molecules) {
    // single out the CPP description
    const auto it =
        std::find_if(configurations.begin(), configurations.end(), [](const std::string_view& str) {
          return str.find(cpp_token) != std::string::npos; // Check if the target is a substring
        });

    // parse the CPP description (if any)
    if (it != configurations.end()) {
      auto configuration = *it;

      configuration = configuration.substr(cpp_token.size());

      // the description should start with a colon
      if (configuration.front() != ':') [[unlikely]] {
        throw std::runtime_error(std::string{"CPP description should start with ':' ("} +
                                 std::string{configuration} + std::string{")"});
      }
      configuration = configuration.substr(1);

      // make sure that the device is the CPU
      const auto colon_index = configuration.find(':');
      const auto device_name = configuration.substr(0, colon_index);
      if (device_name != cpu_token) [[unlikely]] {
        throw std::runtime_error(std::string{"Unsupported device '"} + std::string{device_name} +
                                 std::string{"' for the CPP implementation"});
      }
      configuration = configuration.substr(colon_index);

      // the core counts description should start with a colon
      if (configuration.front() != ':') [[unlikely]] {
        throw std::runtime_error(std::string{"Core count description should start with ':' ("} +
                                 std::string{configuration} + std::string{")"});
      }
      configuration = configuration.substr(1);

      // add the workers that we found parsing the configuration
      for (const auto id: parse_ids(configuration)) {
        pool.add_worker<mudock::cpp_worker>(knobs,
                                            grid_atom_maps,
                                            electro_map,
                                            desolv_map,
                                            input_molecules,
                                            output_molecules,
                                            id);
      }
    }
  }
} // namespace mudock
