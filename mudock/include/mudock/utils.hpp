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

// TODO fix me: there is no real need to have them separated, it could be template for the type of comparison
template<auto Start, auto End, auto Inc, class T, class V, class F>
constexpr void constexpr_switch_bucket(F&& f, T value, V* values) {
  if constexpr (Start < End) {
    if (value <= static_cast<T>(values[Start]))
      f(std::integral_constant<decltype(Start), Start>());
    else
      constexpr_switch_bucket<Start + Inc, End, Inc>(f, value, values);
  }
}

// Unroll control for kernels
#define MUDOCK_STRINGIFY_INNER(x) #x
#define MUDOCK_STRINGIFY(x) MUDOCK_STRINGIFY_INNER(x)
#ifndef MUDOCK_UNROLL_FACTOR
  #define MUDOCK_UNROLL_FACTOR 8
#endif
#define MUDOCK_ATOM_LOOP_UNROLL_FACTOR(MAX_ATOMS, STEP) (((MAX_ATOMS) + (STEP)-1) / (STEP))
#ifdef MUDOCK_DISABLE_UNROLL
  #define MUDOCK_PRAGMA_UNROLL(factor) _Pragma("unroll 1")
#else
  #define MUDOCK_PRAGMA_UNROLL(factor) _Pragma(MUDOCK_STRINGIFY(unroll factor))
#endif

// utility function that reads the whole content of a stream
template<class stream_type>
inline auto read_from_stream(stream_type&& in) {
  assert(in.good());
  return std::string{std::istreambuf_iterator<std::string::value_type>{in},
                     std::istreambuf_iterator<std::string::value_type>{}};
}
constexpr auto is_debug() {
#ifdef MUDOCK_TEST
  return true;
#else
  return false;
#endif
}

template<class T>
constexpr T big_bound() {
  if constexpr (std::is_same_v<T, float>)
    return T(1e30f); // safe for adds/mults in many cases
  if constexpr (std::is_same_v<T, double>)
    return T(1e300);
  return std::numeric_limits<T>::max() / T(4); // fallback
}
