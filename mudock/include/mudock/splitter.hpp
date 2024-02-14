#pragma once

#include <mudock/format.hpp>
#include <string>
#include <vector>

namespace mudock {

  template<class format_type>
    requires can_split<format_type>
  class splitter {
    std::string input_text;
    format_type format_splitter;

  public:
    template<class... args_type>
    splitter(args_type&&... args): format_splitter(args...) {}

    std::vector<std::string> operator()(std::string new_text);
    std::string flush() {
      auto&& text_remainder = std::move(input_text);
      input_text.clear();
      return text_remainder;
    }
  };

  //===------------------------------------------------------------------------------------------------------
  // Out-of-class method definitions
  //===------------------------------------------------------------------------------------------------------

  template<class format_type>
    requires can_split<format_type>
  std::vector<std::string> splitter<format_type>::operator()(std::string new_text) {
    std::vector<std::string> result;

    // append the new text to the old one
    input_text += new_text;

    // find all the molecule in the text
    auto found_molecule_descriptions = std::vector<std::string>{};
    auto parsing_index               = std::string::size_type{0};
    while (parsing_index != std::string::npos) {
      const auto index_next_molecule = format_splitter.next_molecule_start_index(new_text);
      if (index_next_molecule != std::string::npos) { // we found a description
        const auto description_size = index_next_molecule - parsing_index;
        result.emplace_back(input_text.substr(parsing_index, description_size));
      } else { // there is no description
        input_text = input_text.substr(parsing_index);
      }
      parsing_index = index_next_molecule;
    }

    return result;
  }

} // namespace mudock
