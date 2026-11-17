#pragma once

#include <mudock/soa_template/soa_blocks.hpp>
#include <mudock/soa_template/soa_layout.hpp>

namespace mudock {

  GENERATE_SOA_LAYOUT(GeometryFragments,
                      SOA_COLUMN(int, values))

  GENERATE_SOA_BLOCKS(GeometryData,
                      SOA_BLOCK(fragments, GeometryFragments),
                      SOA_BLOCK(fragments_start, GeometryFragments),
                      SOA_BLOCK(fragments_start_atoms, GeometryFragments),
                      SOA_BLOCK(fragments_stop_atoms, GeometryFragments),
                      SOA_BLOCK(fragments_indices_start, GeometryFragments))

}  // namespace mudock
