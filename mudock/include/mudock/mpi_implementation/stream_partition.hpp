#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <istream>
#include <mudock/format/supported_format.hpp>
#include <mudock/mpi_implementation/byte_range.hpp>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace mudock {

namespace detail {

inline std::size_t find_previous_marker(std::istream& file, std::size_t pos,
                                        const std::string_view marker) {
  if (pos == 0) return 0;

  constexpr std::size_t search_block_size = 2000;
  std::string buffer(search_block_size, '\0');

  while (pos > 0) {
    const std::size_t read_size = std::min(search_block_size, pos);
    const std::size_t start = pos - read_size;

    file.clear();
    file.seekg(static_cast<std::streamoff>(start), std::ios::beg);
    file.read(buffer.data(), static_cast<std::streamsize>(read_size));

    const std::size_t got = static_cast<std::size_t>(file.gcount());
    if (got == 0) break;

    const std::string_view view(buffer.data(), got);
    const std::size_t cut = view.rfind(marker);
    if (cut != std::string_view::npos) {
      return start + cut;
    }

    pos = start;
  }

  return 0;
}

}  // namespace detail

template<supported_format format>
inline std::vector<byte_range> compute_aligned_ranges(const std::string& path, int num_ranks) {
  if (num_ranks <= 0) return {};

  std::vector<byte_range> ranges(static_cast<std::size_t>(num_ranks));
  const std::uint64_t file_size = static_cast<std::uint64_t>(std::filesystem::file_size(path));
  if (file_size == 0) return ranges;

  if (num_ranks == 1) {
    ranges.front() = {0, file_size};
    return ranges;
  }

  const std::uint64_t chunk_size = file_size / static_cast<std::uint64_t>(num_ranks);
  std::ifstream file(path, std::ios::binary);
  if (!file) {
    throw std::runtime_error("Failed to open file: " + path);
  }

  for (int rank_id = 0; rank_id < num_ranks; ++rank_id) {
    const std::uint64_t begin = static_cast<std::uint64_t>(rank_id) * chunk_size;
    const std::uint64_t end =
        rank_id == num_ranks - 1 ? file_size : begin + chunk_size;
    ranges[static_cast<std::size_t>(rank_id)] = {begin, end};
  }

  for (int rank_id = 1; rank_id < num_ranks; ++rank_id) {
    ranges[static_cast<std::size_t>(rank_id)].begin = static_cast<std::uint64_t>(
        detail::find_previous_marker(
            file,
            static_cast<std::size_t>(ranges[static_cast<std::size_t>(rank_id)].begin),
            type_of_format<format>::MOLECULE_TOKEN));
  }

  for (int rank_id = 0; rank_id < num_ranks - 1; ++rank_id) {
    ranges[static_cast<std::size_t>(rank_id)].end = static_cast<std::uint64_t>(
        detail::find_previous_marker(
            file,
            static_cast<std::size_t>(ranges[static_cast<std::size_t>(rank_id)].end),
            type_of_format<format>::MOLECULE_TOKEN));
  }

  for (auto& range : ranges) {
    if (range.end < range.begin) range.end = range.begin;
    if (range.begin > file_size) range.begin = file_size;
    if (range.end > file_size) range.end = file_size;
  }

  return ranges;
}

}  // namespace mudock
