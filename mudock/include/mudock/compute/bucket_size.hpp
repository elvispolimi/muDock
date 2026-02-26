#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <mudock/log.hpp>
#include <stdexcept>

namespace mudock {
  template<class get_multiple_t>
  inline int resolve_bucket_size(const char* backend_name,
                                 const int atoms,
                                 const size_t max_bucket_size,
                                 get_multiple_t&& get_default_multiple) {
#ifdef MUDOCK_ADT_BUCKET_OVERRIDE
    const int capped = std::min<int>(MUDOCK_ADT_BUCKET_OVERRIDE, max_bucket_size);
    mudock::info(backend_name,
                 " Bucket size for ",
                 atoms,
                 " atoms override -> ",
                 MUDOCK_ADT_BUCKET_OVERRIDE,
                 ", capped -> ",
                 capped);
    return capped;
#else
    const int bucket_multiple = get_default_multiple();
    if (bucket_multiple <= 0)
      throw std::runtime_error(
          "Compilation error: there is a bucket of atoms number which it is not handled.");

    int bucket_size = max_bucket_size / bucket_multiple;
    bucket_size     = bucket_size == 0 ? max_bucket_size : bucket_size * bucket_multiple;

  #ifdef MUDOCK_ADT_BUCKET_MULTIPLE_OVERRIDE
    const int base_bucket_size = bucket_size;
    const double override_multiplier = static_cast<double>(MUDOCK_ADT_BUCKET_MULTIPLE_OVERRIDE);
    if (override_multiplier <= 0.0) {
      throw std::runtime_error(
          "Compilation error: MUDOCK_ADT_BUCKET_MULTIPLE_OVERRIDE must be > 0.");
    }
    const int scaled_bucket_size =
        static_cast<int>(std::llround(static_cast<double>(bucket_size) * override_multiplier));
    bucket_size = std::max<int>(1, scaled_bucket_size);
    mudock::info(backend_name,
                 " Bucket size for ",
                 atoms,
                 " atoms base bucket size ",
                 base_bucket_size,
                 " scaled by override ",
                 override_multiplier,
                 " -> ",
                 bucket_size);
  #else
    mudock::info(backend_name,
                 " Bucket size for ",
                 atoms,
                 " atoms ",
                 bucket_multiple,
                 " bucket multiple, ",
                 max_bucket_size,
                 " max bucket size -> ",
                 bucket_size);
  #endif
    return bucket_size;
#endif
  }
} // namespace mudock
