#include <array>
#include <mudock/devices/device_types.hpp>
#include <stdexcept>

namespace mudock {

  static std::array<std::string_view, 2> device_kind_names = {{
      "CPU",
      "GPU",
  }};

  device_type parse_device(const std::string& device_kind) {
    for (std::size_t i{0}; i < device_kind_names.size(); ++i) {
      if (device_kind == device_kind_names[i]) {
        return static_cast<device_type>(i);
      }
    }
    throw std::runtime_error{"Invalid device  " + device_kind + " !"};
  }
} // namespace mudock
