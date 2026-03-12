#pragma once

#include <algorithm>
#include <cstdint>
#include <cstddef>
#include <fstream>
#include <filesystem>
#include <istream>
#include <mpi.h>
#include <mudock/format/adt_mol2.hpp>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mudock {

// helper: offset of previous marker token
inline std::size_t find_previous_marker(std::istream& file, std::size_t pos, const std::string_view marker) {
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

// Compute all ranges on rank 0; broadcast to everyone; return range for this rank.
inline std::pair<std::uint64_t, std::uint64_t>
mpi_splitter_bcast(const std::string& path, int rank_id, int num_ranks, MPI_Comm comm = MPI_COMM_WORLD)
{
    if (num_ranks <= 0 || rank_id < 0 || rank_id >= num_ranks)
        return {0, 0};

    // flat buffer: [start0,end0,start1,end1,...]
    std::vector<std::uint64_t> offsets;
    offsets.resize(static_cast<std::size_t>(2 * num_ranks), 0);

    if (rank_id == 0) {
        const std::uint64_t file_size = static_cast<std::uint64_t>(std::filesystem::file_size(path));
        if (file_size == 0) {

        } else if (num_ranks == 1) {
            offsets[0] = 0;
            offsets[1] = file_size;
        } else {
            const std::uint64_t chunk_size = file_size / static_cast<std::uint64_t>(num_ranks);

            std::ifstream file(path, std::ios::binary);
            if (!file)
                throw std::runtime_error("Failed to open file: " + path);

            // first pass: naive chunks
            for (int r = 0; r < num_ranks; ++r) {
                std::uint64_t start_pos = static_cast<std::uint64_t>(r) * chunk_size;
                std::uint64_t end_pos   = (r == num_ranks - 1)
                                        ? file_size
                                        : start_pos + chunk_size;

                offsets[2 * r + 0] = start_pos;
                offsets[2 * r + 1] = end_pos;
            }

            // second pass: align to previous marker
            // rank 0 start stays 0
            for (int r = 1; r < num_ranks; ++r) {
                offsets[2 * r + 0] =
                    static_cast<std::uint64_t>(find_previous_marker(file, 
                        static_cast<std::size_t>(offsets[2 * r + 0]), adt_mol2_tokens::MOLECULE_TOKEN));
            }
            // last rank end stays file_size
            for (int r = 0; r < num_ranks - 1; ++r) {
                offsets[2 * r + 1] =
                    static_cast<std::uint64_t>(find_previous_marker(file, 
                        static_cast<std::size_t>(offsets[2 * r + 1]), adt_mol2_tokens::MOLECULE_TOKEN));
            }

            // optional sanity: enforce non-decreasing and not inverted
            for (int r = 0; r < num_ranks; ++r) {
                auto& s = offsets[2 * r + 0];
                auto& e = offsets[2 * r + 1];
                if (e < s) e = s;
                if (e > file_size) e = file_size;
                if (s > file_size) s = file_size;
            }
        }
    }

    // Broadcast flat offsets to all ranks
    MPI_Bcast(offsets.data(), 2 * num_ranks, MPI_UINT64_T, 0, comm);

    return { offsets[2 * rank_id + 0], offsets[2 * rank_id + 1] };
}

} // namespace mudock
