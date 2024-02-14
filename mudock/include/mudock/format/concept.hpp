#pragma once

#include <concepts>
#include <mudock/molecule.hpp>
#include <string_view>

namespace mudock {

  template<class T>
  concept can_split = requires(T t) {
    { t.next_molecule_start_index(std::string_view()) } -> std::convertible_to<std::string_view::size_type>;
  };

  template<class T, class molecule_type>
  concept can_parse = requires(T t) {
    { t.parse(molecule_type(), std::string_view()) };
  } && is_molecule<molecule_type>;

} // namespace mudock
