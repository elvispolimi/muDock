#pragma once

#include <mudock/chem/drugscore/drugscore_model_data.hpp>

namespace mudock {

  static constexpr const fp_type *DRUGSCORE_MODEL_FLAT() { return DRUGSCORE_MODEL[0][0]; }

  static constexpr unsigned int bucketize(fp_type distance) {
    return static_cast<unsigned int>(distance * 10 + 0.5f);
  }

  static constexpr unsigned int max_bucket() { return bucketize(maxdist()); }
  static constexpr unsigned int min_bucket() { return bucketize(mindist()); }

  static constexpr unsigned int outer_stride() { return num_valued_ds_types() * num_buckets(); }
  static constexpr unsigned int inner_stride() { return num_buckets(); }

} // namespace mudock
