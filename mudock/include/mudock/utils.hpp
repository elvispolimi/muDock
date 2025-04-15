#pragma once

#include <cassert>
#include <iostream>
#include <type_traits>

template<auto Start, auto End, auto Inc, class F>
constexpr void constexpr_for(F&& f) {
  if constexpr (Start < End) {
    f(std::integral_constant<decltype(Start), Start>());
    constexpr_for<Start + Inc, End, Inc>(f);
  }
}

template<auto Start, auto End, auto Inc, class T, class F>
constexpr void constexpr_switch(F&& f, T value) {
  if constexpr (Start < End) {
    if (static_cast<T>(Start) == value)
      f(std::integral_constant<decltype(Start), Start>());
    else
      constexpr_switch<Start + Inc, End, Inc>(f, value);
  }
}

// utility function that reads the whole content of a stream
template<class stream_type>
inline auto read_from_stream(stream_type&& in) {
  assert(in.good());
  return std::string{std::istreambuf_iterator<std::string::value_type>{in},
                     std::istreambuf_iterator<std::string::value_type>{}};
}
constexpr auto is_debug() {
  // TODO this should work only with CMAKE
#ifdef DEBUG_MODE
  return false;
#else
  return false;
#endif
}
