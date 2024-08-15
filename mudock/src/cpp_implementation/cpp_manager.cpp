#include <mudock/compute.hpp>
#include <mudock/cpp_implementation/cpp_manager.hpp>
#include <mudock/cpp_implementation/cpp_worker.hpp>
#include <stdexcept>

namespace mudock {

  static constexpr auto cpp_token = std::string_view{"CPP"};

  void manage_cpp(std::string_view configuration,
                  threadpool& pool,
                  std::shared_ptr<const dynamic_molecule> protein,
                  std::shared_ptr<const grid_atom_mapper> grid_atom_maps,
                  std::shared_ptr<const grid_map> electro_map,
                  std::shared_ptr<const grid_map> desolv_map,
                  std::shared_ptr<safe_stack<static_molecule>> input_molecules,
                  std::shared_ptr<safe_stack<static_molecule>> output_molecules) {
    // single out the CPP description
    const auto begin_cpp_description = configuration.find(cpp_token);
    const auto end_cpp_description   = configuration.find(";", begin_cpp_description);
    configuration                    = configuration.substr(begin_cpp_description + cpp_token.size(),
                                         end_cpp_description - begin_cpp_description);

    // parse the CPP description (if any)
    if (!configuration.empty()) {
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
        pool.add_worker<mudock::cpp_worker>(protein,
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
