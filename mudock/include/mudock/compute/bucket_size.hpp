#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <mudock/compute/batch_multiple.hpp>
#include <mudock/log.hpp>
#include <stdexcept>

namespace mudock {
  template<class get_multiple_t>
  inline int resolve_bucket_size(const char* backend_name,
                                 const int atoms,
                                 const size_t max_bucket_size,
                                 const size_t mem_per_ligand_bytes,
                                 get_multiple_t&& get_default_multiple) {
#ifdef MUDOCK_ADT_BUCKET_OVERRIDE
    (void) get_default_multiple;
    static_assert(MUDOCK_ADT_BUCKET_OVERRIDE > 0,
                  "MUDOCK_ADT_BUCKET_OVERRIDE must be > 0.");
    const int capped                 = std::min<int>(MUDOCK_ADT_BUCKET_OVERRIDE, max_bucket_size);
    const size_t estimated_mem_bytes = static_cast<size_t>(capped) * mem_per_ligand_bytes;
    const double estimated_mem_mib   = static_cast<double>(estimated_mem_bytes) / (1024.0 * 1024.0);
    mudock::info(backend_name,
                 " Bucket size for ",
                 atoms,
                 " atoms override -> ",
                 MUDOCK_ADT_BUCKET_OVERRIDE,
                 ", capped -> ",
                 capped,
                 ", estimated device memory -> ",
                 estimated_mem_bytes,
                 " B (",
                 estimated_mem_mib,
                 " MiB)");
    return capped;
#else
    const auto base_multiple_info = normalize_batch_multiple(get_default_multiple());
    const int base_multiple       = base_multiple_info.total_multiple();
    if (base_multiple <= 0)
      throw std::runtime_error(
          "Compilation error: there is a bucket of atoms number which it is not handled.");

    double effective_multiplier = 1.0;
    int bucket_size             = 1;

  #ifdef MUDOCK_ADT_BUCKET_MULTIPLE_OVERRIDE
    static_assert(MUDOCK_ADT_BUCKET_MULTIPLE_OVERRIDE > 0,
                  "MUDOCK_ADT_BUCKET_MULTIPLE_OVERRIDE must be > 0.");
    const double override_multiplier = static_cast<double>(MUDOCK_ADT_BUCKET_MULTIPLE_OVERRIDE);
    effective_multiplier            = override_multiplier;
    bucket_size = static_cast<int>(std::llround(static_cast<double>(base_multiple) * effective_multiplier));
    if (bucket_size <= 0)
      bucket_size = 1;

    const size_t estimated_mem_bytes = static_cast<size_t>(bucket_size) * mem_per_ligand_bytes;
    const double estimated_mem_mib   = static_cast<double>(estimated_mem_bytes) / (1024.0 * 1024.0);
    mudock::info(backend_name,
                 " Bucket size for ",
                 atoms,
                 " atoms base multiple ",
                 base_multiple,
                 " with override multiplier ",
                 effective_multiplier,
                 " -> ",
                 bucket_size,
                 ", estimated device memory -> ",
                 estimated_mem_bytes,
                 " B (",
                 estimated_mem_mib,
                 " MiB)");
  #else
  #ifdef MUDOCK_ADT_BUCKET_POLICY_MEMORY_ONLY
    bucket_size = std::max<int>(1, static_cast<int>(max_bucket_size));
    effective_multiplier =
        static_cast<double>(bucket_size) / static_cast<double>(base_multiple);
  #elif defined(MUDOCK_ADT_BUCKET_POLICY_ALIGNED_PER_SM)
    const int per_sm_multiple = std::max(1, base_multiple_info.active_blocks_per_sm);
    bucket_size = static_cast<int>((max_bucket_size / static_cast<size_t>(per_sm_multiple)) *
                                   static_cast<size_t>(per_sm_multiple));
    if (bucket_size <= 0) {
      // If budget is smaller than one alignment unit, keep minimum legal batch.
      bucket_size = 1;
    }
    effective_multiplier =
        static_cast<double>(bucket_size) / static_cast<double>(base_multiple);
  #else
    bucket_size = std::max<int>(1, std::min(base_multiple, static_cast<int>(max_bucket_size)));
    effective_multiplier =
        static_cast<double>(bucket_size) / static_cast<double>(base_multiple);
  #endif

    const size_t estimated_mem_bytes = static_cast<size_t>(bucket_size) * mem_per_ligand_bytes;
    const double estimated_mem_mib   = static_cast<double>(estimated_mem_bytes) / (1024.0 * 1024.0);
    mudock::info(backend_name,
                 " Bucket size for ",
                 atoms,
                 " atoms base multiple ",
                 base_multiple,
                 " (active_blocks_per_sm=",
                 base_multiple_info.active_blocks_per_sm,
                 ", num_sms=",
                 base_multiple_info.num_sms,
                 ")",
  #ifdef MUDOCK_ADT_BUCKET_POLICY_MEMORY_ONLY
                 ", memory-only multiplier ",
  #elif defined(MUDOCK_ADT_BUCKET_POLICY_ALIGNED_PER_SM)
                 ", aligned-per-sm multiplier ",
  #else
                 ", aligned-per-gpu multiplier ",
  #endif
                 effective_multiplier,
                 ", max bucket size ",
                 max_bucket_size,
                 " -> ",
                 bucket_size,
                 ", estimated device memory -> ",
                 estimated_mem_bytes,
                 " B (",
                 estimated_mem_mib,
                 " MiB)");
  #endif
    return bucket_size;
#endif
  }
} // namespace mudock
