#include <mudock/format/pdbqt.hpp>
#include <string_view>

namespace mudock {

  std::string_view::size_type pdbqt::next_molecule_start_index(std::string_view text) const {
    const std::string_view start_token{PDBQT_START_TOKEN};
    const auto index_first_token = text.find(start_token);
    return index_first_token != std::string_view::npos
               ? text.find(PDBQT_START_TOKEN, index_first_token + start_token.size())
               : std::string_view::npos;
  }

} // namespace mudock
