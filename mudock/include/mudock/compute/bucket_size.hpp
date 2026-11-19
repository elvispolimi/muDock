#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <mudock/compute/batch_multiple.hpp>
#include <mudock/log.hpp>
#include <stdexcept>
#include <utility>

namespace mudock {
  template<class get_multiple_t>
  inline int resolve_stage_bucket_size(const char* stage_name,
                                       const int atoms,
                                       const size_t max_bucket_size,
                                       const size_t mem_per_ligand_bytes,
                                       const bool honors_stage_bucket_policy,
                                       get_multiple_t&& get_default_multiple) {
    (void) honors_stage_bucket_policy;
#ifdef MUDOCK_STAGE_BUCKET_OVERRIDE
    (void) get_default_multiple;
#endif
#ifdef MUDOCK_STAGE_BUCKET_OVERRIDE
    static_assert(MUDOCK_STAGE_BUCKET_OVERRIDE > 0,
                  "MUDOCK_STAGE_BUCKET_OVERRIDE must be > 0.");
    constexpr int stage_bucket_override = MUDOCK_STAGE_BUCKET_OVERRIDE;
    const size_t capped_sz           = std::min(static_cast<size_t>(stage_bucket_override), max_bucket_size);
    const int capped                 = static_cast<int>(std::min(capped_sz, static_cast<size_t>(std::numeric_limits<int>::max())));
    const size_t estimated_mem_bytes = static_cast<size_t>(capped) * mem_per_ligand_bytes;
    const double estimated_mem_mib   = static_cast<double>(estimated_mem_bytes) / (1024.0 * 1024.0);
    mudock::stage_bucket_trace(stage_name,
                 " stage bucket for ",
                 atoms,
                 " atoms override -> ",
                 stage_bucket_override,
                 ", capped -> ",
                 capped,
                 ", estimated worker batch memory -> ",
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
    const int max_bucket = std::max(
        1,
        static_cast<int>(std::min(max_bucket_size,
                                  static_cast<size_t>(std::numeric_limits<int>::max()))));

  #ifdef MUDOCK_STAGE_BUCKET_MULTIPLE_OVERRIDE
    // An explicit multiplier selects the requested scaling, but the resulting
    // bucket remains capped by the configured worker memory budget.
    static_assert(MUDOCK_STAGE_BUCKET_MULTIPLE_OVERRIDE > 0,
                  "MUDOCK_STAGE_BUCKET_MULTIPLE_OVERRIDE must be > 0.");
    const double override_multiplier = static_cast<double>(MUDOCK_STAGE_BUCKET_MULTIPLE_OVERRIDE);
    effective_multiplier            = override_multiplier;
    const double requested_bucket = static_cast<double>(base_multiple) * effective_multiplier;
    if (!std::isfinite(requested_bucket) || requested_bucket <= 0.0) {
      throw std::runtime_error(
          "MUDOCK_STAGE_BUCKET_MULTIPLE_OVERRIDE must produce a finite positive bucket.");
    }
    const bool memory_cap_applied = requested_bucket >= static_cast<double>(max_bucket);
    bucket_size = memory_cap_applied ? max_bucket : static_cast<int>(std::llround(requested_bucket));
    if (bucket_size <= 0)
      bucket_size = 1;

    const size_t estimated_mem_bytes = static_cast<size_t>(bucket_size) * mem_per_ligand_bytes;
    const double estimated_mem_mib   = static_cast<double>(estimated_mem_bytes) / (1024.0 * 1024.0);
    mudock::stage_bucket_trace(stage_name,
                 " stage bucket for ",
                 atoms,
                 " atoms base multiple ",
                 base_multiple,
                 " with override multiplier ",
                 effective_multiplier,
                 ", requested bucket ",
                 requested_bucket,
                 ", selected bucket ",
                 bucket_size,
                 ", memory cap applied=",
                 memory_cap_applied,
                 ", estimated worker batch memory -> ",
                 estimated_mem_bytes,
                 " B (",
                 estimated_mem_mib,
                 " MiB)");
  #else
  #if defined(MUDOCK_STAGE_BUCKET_POLICY_MAX_UTILIZATION)
    if (honors_stage_bucket_policy) {
      bucket_size = max_bucket;
      effective_multiplier =
          static_cast<double>(bucket_size) / static_cast<double>(base_multiple);
    } else {
      bucket_size = std::min(base_multiple, max_bucket);
      effective_multiplier =
          static_cast<double>(bucket_size) / static_cast<double>(base_multiple);
    }
  #elif defined(MUDOCK_STAGE_BUCKET_POLICY_SM_ALIGNED)
    const int per_sm_multiple = std::max(1, base_multiple_info.active_blocks_per_sm);
    const int sm_aligned_size = (max_bucket / per_sm_multiple) * per_sm_multiple;
    bucket_size = sm_aligned_size;
    if (bucket_size <= 0) {
      bucket_size = 1;
    }
    effective_multiplier =
        static_cast<double>(bucket_size) / static_cast<double>(base_multiple);
  #elif defined(MUDOCK_STAGE_BUCKET_POLICY_DEVICE_ALIGNED)
    // Without an explicit multiplier, prefer the largest whole-GPU multiple
    // that fits the memory assigned to this worker. If no whole-GPU multiple
    // fits, fall back to the largest SM-aligned multiple before using an
    // unaligned memory-safe bucket.
    const int gpu_multiple = std::max(1, base_multiple);
    const int aligned_size = (max_bucket / gpu_multiple) * gpu_multiple;
    const int sm_multiple = std::max(1, base_multiple_info.active_blocks_per_sm);
    const int sm_aligned_size = (max_bucket / sm_multiple) * sm_multiple;
    const int selected_size = aligned_size > 0 ? aligned_size
                                               : (sm_aligned_size > 0 ? sm_aligned_size : max_bucket);
    bucket_size = std::max(1, selected_size);
    if (bucket_size <= 0) {
      // If budget is smaller than one alignment unit, keep minimum legal batch.
      bucket_size = 1;
    }
    effective_multiplier =
        static_cast<double>(bucket_size) / static_cast<double>(base_multiple);
  #else
    bucket_size = std::min(base_multiple, max_bucket);
    effective_multiplier =
        static_cast<double>(bucket_size) / static_cast<double>(base_multiple);
  #endif

    const size_t estimated_mem_bytes = static_cast<size_t>(bucket_size) * mem_per_ligand_bytes;
    const double estimated_mem_mib   = static_cast<double>(estimated_mem_bytes) / (1024.0 * 1024.0);
    const char* alignment_label = "device-aligned";
  #if defined(MUDOCK_STAGE_BUCKET_POLICY_MAX_UTILIZATION)
    alignment_label = honors_stage_bucket_policy ? "max-utilization" : "device-aligned fallback";
  #elif defined(MUDOCK_STAGE_BUCKET_POLICY_SM_ALIGNED)
    alignment_label = "sm-aligned";
  #elif defined(MUDOCK_STAGE_BUCKET_POLICY_DEVICE_ALIGNED)
    alignment_label = aligned_size > 0
                          ? "device-aligned"
                          : (sm_aligned_size > 0 ? "sm-aligned fallback" : "max-size fallback");
  #endif
    mudock::stage_bucket_trace(stage_name,
                 " stage bucket for ",
                 atoms,
                 " atoms base multiple ",
                 base_multiple,
                 " (active_blocks_per_sm=",
                 base_multiple_info.active_blocks_per_sm,
                 ", num_sms=",
                 base_multiple_info.num_sms,
                 ")",
                 ", ",
                 alignment_label,
                 " multiplier ",
                 effective_multiplier,
                 ", max bucket size ",
                 max_bucket_size,
                 " -> ",
                 bucket_size,
                 ", estimated worker batch memory -> ",
                 estimated_mem_bytes,
                 " B (",
                 estimated_mem_mib,
                 " MiB)");
  #endif
    return bucket_size;
#endif
  }
} // namespace mudock
