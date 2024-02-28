#include <mudock/devices/cpp/kernel.hpp>

namespace mudock {
  template<>
  void dock_and_score<language_types::CPP>(language_scratchpad_impl<language_types::CPP>& mol_scratchpad,
                                           const device_context<language_types::CPP>& dev_context) {
    // TODO Add kernel
    for (auto& mol: mol_scratchpad.get_ligands()) { 
      std::string score = "N/A";
      mol->properties.assign(property_type::SCORE, score); }
  };
} // namespace mudock
