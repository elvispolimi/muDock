#include <mudock/format/pdb.hpp>
#include <stdexcept>
#include <string_view>

namespace mudock {

  std::string_view::size_type pdb::next_molecule_start_index(std::string_view) const {
    throw std::runtime_error("PDB Splitter not implemented yet!");
  }

} // namespace mudock
