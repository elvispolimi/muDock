#include <mudock/format/pdb.hpp>
#include <stdexcept>
#include <string_view>

namespace mudock {

  std::string_view::size_type pdb::next_molecule_start_index(std::string_view text) const {
    // throw std::runtime_error("PDB Splitter not implemented yet!");
    const std::string_view start_token{PDB_START_TOKEN};
    const auto index_first_token = text.find(start_token);
    return index_first_token != std::string_view::npos
               ? text.find(PDB_START_TOKEN, index_first_token + start_token.size())
               : std::string_view::npos;
  }

} // namespace mudock
