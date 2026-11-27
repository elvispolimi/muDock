#pragma once

#include <algorithm>
#include <cstddef>
#include <limits>
#include <mudock/log.hpp>

namespace mudock {
  inline int resolve_stage_bucket_size(const char* stage_name,
                                      const int atoms,
                                      const std::size_t memory_limited_size,
                                      const std::size_t memory_per_ligand) {
    const auto memory_cap = std::max(
        std::size_t{1},
        std::min(memory_limited_size, static_cast<std::size_t>(std::numeric_limits<int>::max())));
    std::size_t selected_size = memory_cap;
#ifdef MUDOCK_STAGE_BUCKET_OVERRIDE
    static_assert(MUDOCK_STAGE_BUCKET_OVERRIDE > 0,
                  "MUDOCK_STAGE_BUCKET_OVERRIDE must be greater than zero.");
    selected_size = std::min(static_cast<std::size_t>(MUDOCK_STAGE_BUCKET_OVERRIDE), memory_cap);
#endif
    const auto estimated_memory = selected_size * memory_per_ligand;
    mudock::stage_bucket_trace(stage_name,
                               " batch size for ",
                               atoms,
                               " atoms: ",
                               selected_size,
                               " ligands (memory cap ",
                               memory_cap,
                               ", estimated batch memory ",
                               estimated_memory,
                               " B)");
    return static_cast<int>(selected_size);
  }
} // namespace mudock
