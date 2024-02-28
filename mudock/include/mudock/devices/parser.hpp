#pragma once

#include <cstdint>
#include <mudock/devices/device_types.hpp>
#include <mudock/devices/kernel_types.hpp>
#include <mudock/type_alias.hpp>
#include <string_view>
#include <vector>

namespace mudock {

  // This structure contains the devices's configurations obtained from command line:
  //  - the device type, for example CPU
  //  - the kernel type, for example CPP
  //  - the ids of the device
  struct device_conf {
    const device_type d_t;
    const kernel_type l_t;
    const std::vector<std::size_t> ids;
  };

  // Get the device configurations defined by the user from command line
  device_conf parse_conf(const std::string& conf);
} // namespace mudock
