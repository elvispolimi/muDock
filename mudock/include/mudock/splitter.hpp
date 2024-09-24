#pragma once

#include <mudock/format.hpp>
#include <string>
#include <string_view>
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

    [[nodiscard]] std::vector<std::string> operator()(std::string new_text);
    [[nodiscard]] std::string flush() {
      auto text_remainder = std::move(input_text);
      input_text.clear();
      return text_remainder;
    }
  };

  //===------------------------------------------------------------------------------------------------------
  // Out-of-class method definitions
  //===------------------------------------------------------------------------------------------------------

  template<class format_type>
  requires can_split<format_type>
  [[nodiscard]] std::vector<std::string> splitter<format_type>::operator()(std::string new_text) {
    std::vector<std::string> result;

    // append the new text to the old one
    input_text += new_text;

    // find all the molecule in the text
    auto input_view = std::string_view{input_text};
    while (!input_view.empty()) {
      const auto index_next_molecule = format_splitter.next_molecule_start_index(input_view);
      if (index_next_molecule != std::string::npos) { // we found a description
        result.emplace_back(input_view.substr(std::size_t{0}, index_next_molecule));
        input_view = input_view.substr(index_next_molecule);
      } else { // there is no description
        input_text = std::string{input_view};
        input_view = std::string_view{};
      }
    }

    return result;
  }

} // namespace mudock
