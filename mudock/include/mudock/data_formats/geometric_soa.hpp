#pragma once

#include <mudock/soa_template/soa_blocks.hpp>
#include <mudock/soa_template/soa_layout.hpp>

namespace mudock {

  GENERATE_SOA_LAYOUT(GeometryFragments,
                      SOA_COLUMN(int, values))

}  // namespace mudock
