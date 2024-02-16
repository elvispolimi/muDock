#pragma once

#include <string>

namespace mudock {
  enum class device_types { CPU = 0, GPU };

  device_types parse_device(std::string);

} // namespace mudock
