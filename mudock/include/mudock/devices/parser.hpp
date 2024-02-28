#pragma once

#include <mudock/devices/device_types.hpp>
#include <mudock/devices/language_types.hpp>
#include <mudock/type_alias.hpp>
#include <string_view>
#include <vector>

namespace mudock {

  // This structure contains the devices's configurations obtained from command line:
  //  - the device type, for example CPU
  //  - the language types, for example CPP
  //  - the ids of the device
  struct device_conf {
    const device_types d_t;
    const language_types l_t;
    const std::vector<index_type> ids;
  };

  // Get the device configurations defined by the user from command line
  device_conf parse_conf(const std::string);
} // namespace mudock
