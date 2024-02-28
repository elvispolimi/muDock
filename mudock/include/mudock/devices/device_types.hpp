#pragma once

#include <string>

namespace mudock {
  // Enumerates the devices that are support by muDock
  enum class device_types { CPU = 0, GPU };

  // From a string you get the device type, if available
  // Otherwise it throws a runtime exception 
  device_types parse_device(std::string);

} // namespace mudock
