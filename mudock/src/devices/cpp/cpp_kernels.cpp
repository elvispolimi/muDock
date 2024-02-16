#include <mudock/devices/cpp/cpp_kernels.hpp>

namespace mudock {
  template<>
  void dock_and_score<language_types::CPP>(molecules_scratchpad<language_types::CPP>& mol_scratchpad,
                                           const device_context<language_types::CPP>& dev_context) {
    // TODO Add kernel
    for (auto& mol: mol_scratchpad.get_ligands()) { 
      std::string score = "N/A";
      mol->properties.assign(property_type::SCORE, score); }
  };
} // namespace mudock
