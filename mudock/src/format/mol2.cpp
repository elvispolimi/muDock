#include <mudock/format/mol2.hpp>

namespace mudock {

  std::string_view::size_type mol2::next_molecule_start_index(std::string_view text) const {
    return text.find("@<TRIPOS>MOLECULE");
  }

} // namespace mudock
