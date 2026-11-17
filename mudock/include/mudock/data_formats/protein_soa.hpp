#pragma once

#include <mudock/soa_template/soa_layout.hpp>

namespace mudock {

  GENERATE_SOA_LAYOUT(ProteinBounds,
                      SOA_COLUMN(fp_type, minimum),
                      SOA_COLUMN(fp_type, maximum),
                      SOA_COLUMN(fp_type, center))

  GENERATE_SOA_LAYOUT(ProteinDimensions,
                      SOA_COLUMN(int, size_x),
                      SOA_COLUMN(int, size_xy),
                      SOA_COLUMN(int, size_xyz))

  GENERATE_SOA_LAYOUT(ProteinMaps,
                      SOA_COLUMN(fp_type, maps))

  GENERATE_SOA_BLOCKS(ProteinData,
                      SOA_BLOCK(protein_bounds, ProteinBounds),
                      SOA_BLOCK(protein_dimensions, ProteinDimensions),
                      SOA_BLOCK(protein_maps, ProteinMaps))

}
