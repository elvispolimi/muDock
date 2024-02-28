#pragma once

#include <string>

namespace mudock {
  // Enumerates the devices that are support by muDock
  enum class device_type { CPU = 0, GPU };

  // From a string you get the device type, if available
  // Otherwise it throws a runtime exception
  device_type parse_device(const std::string& device_kind);

} // namespace mudock
