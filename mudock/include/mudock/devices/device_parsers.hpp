#pragma once

#include <mudock/devices/device_types.hpp>
#include <mudock/devices/language_types.hpp>
#include <mudock/type_alias.hpp>
#include <string_view>
#include <vector>

namespace mudock {

  struct device_conf {
    const device_types d_t;
    const language_types l_t;
    const std::vector<index_type> ids;
  };

  device_conf parse_conf(const std::string);
} // namespace mudock
