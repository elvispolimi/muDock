#include <mudock/format/mol2x.hpp>
#include <string_view>

namespace mudock {

  std::string_view::size_type mol2x::next_molecule_start_index(std::string_view text) const {
    static constexpr auto molecule_token = std::string_view{"@<TRIPOS>MOLECULE"};
    const auto index_first_token         = text.find(molecule_token);
    return index_first_token != std::string_view::npos
               ? text.find(molecule_token, index_first_token + molecule_token.size())
               : std::string_view::npos;
  }

} // namespace mudock
