#include <mudock/devices/cpp/kernel.hpp>

namespace mudock {
  template<>
  void dock_and_score<kernel_type::CPP>(kernel_scratchpad_impl<kernel_type::CPP>& mol_scratchpad,
                                        const device_context<kernel_type::CPP>& dev_context) {
    // TODO Add kernel
    for (auto& mol: mol_scratchpad.get_ligands()) {
      std::string score = "N/A";
      mol->properties.assign(property_type::SCORE, score);
    }
  };
} // namespace mudock
