#pragma once

#include <array>
#include <mudock/chem/drugscore/drugscore_types.hpp>
#include <mudock/type_alias.hpp>

//===------------------------------------------------------------------------------------------------------
// WARNING: This file has been automatically generated from chem/drugscore_scores.json
//===------------------------------------------------------------------------------------------------------

namespace mudock {

  static constexpr fp_type maxdist() { return 6.0f; }
  static constexpr fp_type mindist() { return 1.1f; }
  static constexpr unsigned int num_buckets() { return 50; }

  // Indexed by : DRUGSCORE_MODEL[ligand][protein][distance] = score
  extern const fp_type DRUGSCORE_MODEL[num_valued_ds_types()][num_valued_ds_types()][num_buckets()];

} // namespace mudock
