#include <array>
#include <mudock/devices/kernel_types.hpp>
#include <stdexcept>

namespace mudock {

  static std::array<std::string_view, 1> language_kind_names = {{
      "CPP",
  }};

  kernel_type parse_kernel(const std::string& language_kind) {
    for (std::size_t i{0}; i < language_kind_names.size(); ++i) {
      if (language_kind == language_kind_names[i]) {
        return static_cast<kernel_type>(i);
      }
    }
    throw std::runtime_error{"Invalid kernel  " + language_kind + " !"};
  }
} // namespace mudock
