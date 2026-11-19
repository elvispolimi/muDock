#include <mudock/format/pdb.hpp>
#include <stdexcept>
#include <string_view>

namespace mudock {

  std::string_view::size_type pdb::next_molecule_start_index(std::string_view text) const {
    constexpr std::string_view model_token = "MODEL ";
    const auto first_model = text.find(model_token);
    if (first_model != std::string_view::npos) {
      const auto second_model = text.find(model_token, first_model + model_token.size());
      if (second_model != std::string_view::npos) return second_model;
    }
    const std::string_view start_token{PDB_START_TOKEN};
    const auto index_first_token = text.find(start_token);
    return index_first_token != std::string_view::npos
               ? text.find(PDB_START_TOKEN, index_first_token + start_token.size())
               : std::string_view::npos;
  }

} // namespace mudock
