#pragma once

#include <algorithm>
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
  #ifdef MUDOCK_ADT_BUCKET_MULTIPLE_OVERRIDE
    const int bucket_multiple = MUDOCK_ADT_BUCKET_MULTIPLE_OVERRIDE;
  #else
    const int bucket_multiple = get_default_multiple();
    if (bucket_multiple == 0)
      throw std::runtime_error(
          "Compilation error: there is a bucket of atoms number which it is not handled.");
  #endif

    int bucket_size = max_bucket_size / bucket_multiple;
    bucket_size     = bucket_size == 0 ? max_bucket_size : bucket_size * bucket_multiple;
    mudock::info(backend_name,
                 " Bucket size for ",
                 atoms,
                 " atoms ",
                 bucket_multiple,
                 " bucket multiple, ",
                 max_bucket_size,
                 " max bucket size -> ",
                 bucket_size);
    return bucket_size;
#endif
  }
} // namespace mudock
