#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <stdexcept>
#include <string_view>

namespace mudock {
  //TODO add CMake conf file
  constexpr bool is_gh_enabled =
#ifdef MUDOCK_USE_GH
      true;
#else
      false;
#endif

  constexpr bool is_xsimd_enabled =
#ifdef MUDOCK_USE_XSIMD
      true;
#else
      false;
#endif

  enum class cpu_vectorization { AUTO = 0, GH, XSIMD };

  struct cpu_vect_description {
    cpu_vectorization value;
    std::string_view name;
    bool is_enabled;
  };

  static constexpr auto num_vectorization_type() { return 3; }

  static constexpr std::array<cpu_vect_description, num_vectorization_type()> VECT_DESCRIPTION = {
      {{cpu_vectorization::AUTO, "CPP", true},
       {cpu_vectorization::GH, "GH", is_gh_enabled},
       {cpu_vectorization::XSIMD, "XSIMD", is_xsimd_enabled}}};

  inline constexpr const cpu_vect_description& get_description(const cpu_vectorization e) {
    assert(VECT_DESCRIPTION[static_cast<int>(e)].value == e);
    return VECT_DESCRIPTION[static_cast<int>(e)];
  }

  inline cpu_vectorization parse_vect_type(const std::string_view vect) {
    const auto element_it = std::find_if(std::begin(VECT_DESCRIPTION),
                                         std::end(VECT_DESCRIPTION),
                                         [&vect](const auto& e) { return e.name == vect; });
    if (element_it != std::end(VECT_DESCRIPTION))
      return element_it->value;
    else
      throw std::runtime_error("Missing vectorization type");
  }
} // namespace mudock
