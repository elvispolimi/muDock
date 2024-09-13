#pragma once

#include <cstdint>
#include <string_view>
#include <vector>

namespace mudock {

  std::vector<int> parse_ids(std::string_view description);

}
