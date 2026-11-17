#pragma once

#include <mudock/soa_template/soa_layout.hpp>

namespace mudock {

  GENERATE_SOA_LAYOUT(BatchInfo,
                      SOA_COLUMN(int, num_atoms),
                      SOA_COLUMN(int, num_rotamers))

}