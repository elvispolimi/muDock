#pragma once

#include <cstdint>
#include <mpi.h>
#include <mudock/format/supported_format.hpp>
#include <mudock/mpi_implementation/byte_range.hpp>
#include <mudock/mpi_implementation/stream_partition.hpp>
#include <string>
#include <vector>

namespace mudock {

inline byte_range broadcast_range(const std::vector<byte_range>& ranges, int rank_id, int num_ranks,
                                  MPI_Comm comm = MPI_COMM_WORLD) {
  if (num_ranks <= 0 || rank_id < 0 || rank_id >= num_ranks) return {};

  std::vector<std::uint64_t> offsets(static_cast<std::size_t>(2 * num_ranks), 0);
  if (rank_id == 0) {
    for (int i = 0; i < num_ranks; ++i) {
      offsets[static_cast<std::size_t>(2 * i)] = ranges[static_cast<std::size_t>(i)].begin;
      offsets[static_cast<std::size_t>(2 * i + 1)] = ranges[static_cast<std::size_t>(i)].end;
    }
  }

  MPI_Bcast(offsets.data(), 2 * num_ranks, MPI_UINT64_T, 0, comm);
  return {offsets[static_cast<std::size_t>(2 * rank_id)],
          offsets[static_cast<std::size_t>(2 * rank_id + 1)]};
}

template<supported_format format>
inline byte_range distribute_aligned_ranges(const std::string& path, int rank_id, int num_ranks,
                                            MPI_Comm comm = MPI_COMM_WORLD) {
  const auto ranges =
      rank_id == 0 ? compute_aligned_ranges<format>(path, num_ranks) : std::vector<byte_range>{};
  return broadcast_range(ranges, rank_id, num_ranks, comm);
}

template<supported_format format>
inline byte_range mpi_splitter_bcast(const std::string& path, int rank_id, int num_ranks,
                                     MPI_Comm comm = MPI_COMM_WORLD) {
  return distribute_aligned_ranges<format>(path, rank_id, num_ranks, comm);
}

}  // namespace mudock
