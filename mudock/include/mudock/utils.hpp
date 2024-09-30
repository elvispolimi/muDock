#pragma once

#include <type_traits>

template<auto Start, auto End, auto Inc, class F>
constexpr void constexpr_for(F&& f) {
  if constexpr (Start < End) {
    f(std::integral_constant<decltype(Start), Start>());
    constexpr_for<Start + Inc, End, Inc>(f);
  }
}

constexpr auto is_debug() {
  // TODO this should work only with CMAKE
#ifdef DEBUG_MODE
  return true;
#else
  return false;
#endif
}