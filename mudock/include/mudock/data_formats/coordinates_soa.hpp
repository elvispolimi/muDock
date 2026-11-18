#pragma once

#include <mudock/soa_template/soa_layout.hpp>

namespace mudock {

  GENERATE_SOA_LAYOUT(Coordinates,
                      SOA_COLUMN(fp_type, x),
                      SOA_COLUMN(fp_type, y),
                      SOA_COLUMN(fp_type, z))

}