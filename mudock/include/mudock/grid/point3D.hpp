#pragma once

#include <cmath>
#include <mudock/grid/pi.hpp>
#include <mudock/type_alias.hpp>
#include <type_traits>
#include <utility>

namespace mudock {

  struct point3D {
    fp_type x = fp_type{0.0};
    fp_type y = fp_type{0.0};
    fp_type z = fp_type{0.0};
  };

  template<class point_type>
  constexpr std::remove_reference_t<point_type> add(point_type&& a, point_type&& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
  }

  template<class point_type, class... point_types>
  constexpr std::remove_reference_t<point_type> add(point_type&& a, point_types&&... bs) {
    return add(std::forward<point_type>(a), std::forward<point_type>(bs...));
  }

  template<class point_type>
  constexpr std::remove_reference_t<point_type> difference(point_type&& a, point_type&& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
  }

  template<class point_type, class... point_types>
  constexpr std::remove_reference_t<point_type> difference(point_type&& a, point_types&&... bs) {
    return difference(std::forward<point_type>(a), std::forward<point_type>(bs...));
  }

  template<class point_type>
  constexpr std::remove_reference_t<point_type> scalar_division(point_type&& a, const fp_type value) {
    return {a.x / value, a.y / value, a.z / value};
  }

  template<class point_type>
  constexpr std::remove_reference_t<point_type> scale(point_type&& a, const fp_type value) {
    return {a.x * value, a.y * value, a.z * value};
  }

  template<class point_type>
  constexpr std::remove_reference_t<point_type> square(point_type&& point) {
    return {point.x * point.x, point.y * point.y, point.z * point.z};
  }

  template<class point_type>
  constexpr fp_type sum_components(point_type&& point) {
    return point.x + point.y + point.z;
  }

  template<class point_type>
  constexpr std::remove_reference_t<point_type> normalize(point_type&& a, point_type&& b) {
    point3D diff            = difference(point3D{a.x, a.y, a.z}, point3D{b.x, b.y, b.z});
    fp_type square_distance = sum_components(square(diff));
    if (square_distance == fp_type{0}) {
      // TODO @Davide we should log something
      // error("Attempt to divide by zero was just prevented.");
      square_distance = std::numeric_limits<fp_type>::epsilon();
    }
    // TODO ask @Davide about this
    fp_type inv_rd = fp_type{1} / sqrtf(square_distance);
    return scale(diff, inv_rd);
  }

  template<class point_type>
  constexpr std::remove_reference_t<point_type> normalize(point_type&& diff) {
    fp_type square_distance = sum_components(square(diff));
    if (square_distance == fp_type{0}) {
      // TODO @Davide we should log something
      // error("Attempt to divide by zero was just prevented.");
      square_distance = std::numeric_limits<fp_type>::epsilon();
    }
    // TODO ask @Davide about this
    fp_type inv_rd = fp_type{1} / sqrtf(square_distance);
    return scale(diff, inv_rd);
  }

  template<class point_type>
  constexpr fp_type diff_components(point_type&& point) {
    return point.x - point.y - point.z;
  }

  template<class point_type>
  constexpr fp_type distance2(point_type&& a, point_type&& b) {
    return sum_components(square(difference(std::forward<point_type>(a), std::forward<point_type>(b))));
  }

  template<class point_type>
  constexpr fp_type distance(point_type&& a, point_type&& b) {
    return sqrtf(
        sum_components(square(difference(std::forward<point_type>(a), std::forward<point_type>(b)))));
  }

  template<class point_type>
  constexpr fp_type inner_product(point_type&& a, point_type&& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
  }

  template<class point_type>
  constexpr std::remove_reference_t<point_type> product(point_type&& a, point_type&& b) {
    return {a.x * b.x, a.y * b.y, a.z * b.z};
  }

  // this function consider the given point as the end of a vector that starts in the origin
  template<class point_type>
  constexpr fp_type magnitude(point_type&& a) {
    return std::sqrt(sum_components(square(std::forward<point_type>(a))));
  }

  template<class point_type>
  constexpr fp_type angle(point_type&& origin, point_type&& a, point_type&& b) {
    auto v1 = difference(std::forward<point_type>(a), std::forward<point_type>(origin));
    auto v2 = difference(std::forward<point_type>(b), std::forward<point_type>(origin));
    return std::acos(inner_product(v1, v2) / (magnitude(v1) * magnitude(v2)));
  }

} // namespace mudock
