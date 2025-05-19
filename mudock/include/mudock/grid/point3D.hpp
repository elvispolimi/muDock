#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <mudock/grid/pi.hpp>
#include <mudock/type_alias.hpp>
#include <ranges>
#include <type_traits>
#include <utility>

namespace mudock {
  template<typename T, std::size_t n>
  class point {
    static_assert(n > 0, "A multi dimensional point must at least a dimension");

  protected:
    // the storage of sizes and coefficients for dealing with md indexes
    std::array<T, n> _sizes;

  public:
    static constexpr auto num_dimensions = n;

    ~point()                  = default;
    point(point&& other)      = default;
    point(const point& other) = default;
    point& operator=(point&& other) {
      _sizes = std::move(other._sizes);
      return *this;
    };
    point& operator=(const point& other) {
      _sizes = other._sizes;
      return *this;
    };

    template<typename... I>
    point(I... sizes)
      requires(std::conjunction_v<std::is_same<T, I>...> && sizeof...(sizes) > 0 && sizeof...(sizes) <= n)
    {
      const auto size_list = std::initializer_list{static_cast<std::size_t>(sizes)...};
      std::copy(std::cbegin(size_list), std::cend(size_list), std::begin(_sizes));
    }

    point() {};

    template<std::size_t index>
    [[nodiscard]] std::size_t size() const {
      return _sizes[index];
    }

    bool operator<(const point<T, n>& other) const {
      return [&]<std::size_t... I>(const std::array<T, n>& a,
                                   const std::array<T, n>& b,
                                   std::index_sequence<I...>) {
        return ((a[I] < b[I]) || ...);
      }(_sizes, other._sizes, std::make_index_sequence<n>{});
    }
    bool sum_components() const {
      return [&]<std::size_t... I>(const std::array<T, n>& a, std::index_sequence<I...>) {
        return (a[I] + ...);
      }(_sizes, std::make_index_sequence<n>{});
    }
    bool operator>(const point<T, n>& other) const { return other < this; }
    T distance(const point<T, n>& other) const {
      return std::sqrt((*this - other).square().sum_components());
    }
    point<T, n> operator-(const point<T, n>& other) const {
      return [&]<std::size_t... I>(const std::array<T, n>& a,
                                   const std::array<T, n>& b,
                                   std::index_sequence<I...>) {
        return point<T, n>{(a[I] - b[I])...};
      }(_sizes, other._sizes, std::make_index_sequence<n>{});
    }
    point<T, n> operator+(const point<T, n>& other) const {
      return [&]<std::size_t... I>(const std::array<T, n>& a,
                                   const std::array<T, n>& b,
                                   std::index_sequence<I...>) {
        return point<T, n>{(a[I] + b[I])...};
      }(_sizes, other._sizes, std::make_index_sequence<n>{});
    }
    point<T, n> square() const {
      return [&]<std::size_t... I>(const std::array<T, n>& a, std::index_sequence<I...>) {
        return point<T, n>{(a[I] * a[I])...};
      }(_sizes, std::make_index_sequence<n>{});
    }

    point<T, n> truncate() {
      return [&]<std::size_t... I>(const std::array<T, n>& a, std::index_sequence<I...>) {
        return point<T, n>{std::trunc(a[I])...};
      }(_sizes, std::make_index_sequence<n>{});
    }

    point<T, n> operator*(const T scale) const {
      return [&]<std::size_t... I>(const std::array<T, n>& a, std::index_sequence<I...>) {
        return point<T, n>{(a[I] * scale)...};
      }(_sizes, std::make_index_sequence<n>{});
    }
    point<T, n> operator*(const point<T, n>& other) const {
      return [&]<std::size_t... I>(const std::array<T, n>& a,
                                   const std::array<T, n>& b,
                                   std::index_sequence<I...>) {
        return point<T, n>{(a[I] * b[I])...};
      }(_sizes, other._sizes, std::make_index_sequence<n>{});
    }
    template<typename A>
    point<T, n> product(const A other) const {
      return *this * other;
    }
  };

  struct point3D: public point<fp_type, 3> {
    fp_type& x = _sizes[0];
    fp_type& y = _sizes[1];
    fp_type& z = _sizes[2];

    point3D(const fp_type _x, const fp_type _y, const fp_type _z): point<fp_type, 3>(_x, _y, _z) {};
    point3D(const fp_type v): point<fp_type, 3>(v, v, v) {};
    point3D(): point<fp_type, 3>() {};

    ~point3D()                    = default;
    point3D(point3D&& other)      = default;
    point3D(const point3D& other) = default;
    point3D& operator=(point3D&& other) {
      point::operator=(other);
      return *this;
    }
    point3D& operator=(const point3D& other) {
      point::operator=(other);
      return *this;
    }

    inline std::array<fp_type, 3> get_array() const { return {x, y, z}; }
  };
  // FIX ME make them deprected
  template<class point_type>
  [[nodiscard]] constexpr std::remove_reference_t<point_type> add(point_type&& a, point_type&& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
  }

  template<class point_type>
  [[nodiscard]] constexpr std::remove_reference_t<point_type> divide(point_type&& a, point_type&& b) {
    return {a.x / b.x, a.y / b.y, a.z / b.z};
  }

  template<class point_type>
  [[nodiscard]] constexpr std::remove_reference_t<point_type> difference(point_type&& a, point_type&& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
  }

  template<class point_type, class... point_types>
  [[nodiscard]] constexpr std::remove_reference_t<point_type> difference(point_type&& a,
                                                                         point_types&&... bs) {
    return difference(std::forward<point_type>(a), std::forward<point_type>(bs...));
  }

  template<class point_type>
  [[nodiscard]] constexpr std::remove_reference_t<point_type> scalar_division(point_type&& a,
                                                                              const fp_type value) {
    return {a.x / value, a.y / value, a.z / value};
  }

  template<class point_type>
  [[nodiscard]] constexpr std::remove_reference_t<point_type> scale(point_type&& a, const fp_type value) {
    return {a.x * value, a.y * value, a.z * value};
  }

  template<class point_type>
  [[nodiscard]] constexpr std::remove_reference_t<point_type> square(point_type&& point) {
    return {point.x * point.x, point.y * point.y, point.z * point.z};
  }

  template<class point_type>
  [[nodiscard]] constexpr fp_type sum_components(point_type&& point) {
    return point.x + point.y + point.z;
  }

  template<class point_type>
  [[nodiscard]] constexpr std::remove_reference_t<point_type> normalize(point_type&& diff) {
    fp_type square_distance = sum_components(square(diff));
    if (square_distance == fp_type{0}) {
      // TODO @Davide we should log something
      // error("Attempt to divide by zero was just prevented.");
      square_distance = std::numeric_limits<fp_type>::epsilon();
    }
    // TODO ask @Davide about this
    fp_type inv_rd = fp_type{1} / std::sqrt(square_distance);
    return scale(diff, inv_rd);
  }

  // TODO fix the constexpr here
  template<class point_type>
  [[nodiscard]] std::remove_reference_t<point_type> normalize(point_type&& a, point_type&& b) {
    point3D diff = difference(point3D{a.x, a.y, a.z}, point3D{b.x, b.y, b.z});
    return normalize(diff);
  }

  template<class point_type>
  [[nodiscard]] constexpr std::remove_reference_t<point_type> cross_product(point_type&& a, point_type&& b) {
    return point3D{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
  }

  template<class point_type>
  [[nodiscard]] constexpr fp_type distance2(point_type&& a, point_type&& b) {
    return sum_components(square(difference(std::forward<point_type>(a), std::forward<point_type>(b))));
  }

  template<class point_type>
  [[nodiscard]] constexpr fp_type distance(point_type&& a, point_type&& b) {
    return std::sqrt(
        sum_components(square(difference(std::forward<point_type>(a), std::forward<point_type>(b)))));
  }

  template<class point_type>
  [[nodiscard]] constexpr fp_type inner_product(point_type&& a, point_type&& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
  }

  template<class point_type>
  [[nodiscard]] constexpr std::remove_reference_t<point_type> product(point_type&& a, point_type&& b) {
    return {a.x * b.x, a.y * b.y, a.z * b.z};
  }

  template<class point_type, class... point_types>
  [[nodiscard]] constexpr std::remove_reference_t<point_type> product(point_type&& a, point_types&&... bs) {
    return product(std::forward<point_type>(a), std::forward<point_type>(bs...));
  }

  // this function consider the given point as the end of a vector that starts in the origin
  template<class point_type>
  [[nodiscard]] constexpr fp_type magnitude(point_type&& a) {
    return std::sqrt(sum_components(square(std::forward<point_type>(a))));
  }

  template<class point_type>
  [[nodiscard]] constexpr fp_type angle(point_type&& origin, point_type&& a, point_type&& b) {
    auto v1 = difference(std::forward<point_type>(a), std::forward<point_type>(origin));
    auto v2 = difference(std::forward<point_type>(b), std::forward<point_type>(origin));
    return std::acos(inner_product(v1, v2) / (magnitude(v1) * magnitude(v2)));
  }

} // namespace mudock
