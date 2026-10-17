#pragma once

#include <mudock/soa_template/soa_blocks.hpp>
#include <mudock/soa_template/soa_layout.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  GENERATE_SOA_LAYOUT(AdtAtomData,
                      SOA_COLUMN(fp_type, vols),
                      SOA_COLUMN(fp_type, solpars),
                      SOA_COLUMN(fp_type, charges),
                      SOA_COLUMN(int, map_offsets))

  GENERATE_SOA_LAYOUT(AdtNumNonbondData,
                      SOA_COLUMN(int, num_nonbond))

  GENERATE_SOA_LAYOUT(AdtNonbondData,
                      SOA_COLUMN(int, a1),
                      SOA_COLUMN(int, a2),
                      SOA_COLUMN(fp_type, cA),
                      SOA_COLUMN(fp_type, cB),
                      SOA_COLUMN(int, xB))

  GENERATE_SOA_BLOCKS(AdtScoreData,
                      SOA_BLOCK(atom_data, AdtAtomData),
                      SOA_BLOCK(num_nonbond_data, AdtNumNonbondData),
                      SOA_BLOCK(nonbond_data, AdtNonbondData))

}  // namespace mudock