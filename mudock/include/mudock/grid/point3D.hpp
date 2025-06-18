#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <mudock/grid/pi.hpp>
#include <mudock/type_alias.hpp>
#include <type_traits>
#include <utility>

namespace mudock {
  template<typename T, std::size_t n>
  class point {
    static_assert(n > 0, "A multi dimensional point must at least a dimension");

  protected:
    std::array<T, n> _components;

  public:
    static constexpr auto num_dimensions = n;

    ~point()                  = default;
    point(point&& other)      = default;
    point(const point& other) = default;
    point& operator=(point&& other) {
      _components = std::move(other._components);
      return *this;
    };
    point& operator=(const point& other) {
      _components = other._components;
      return *this;
    };
    auto begin() { return _components.begin(); }
    auto end() { return _components.end(); }

    auto* data() const { return _components.data(); }

    template<typename... I,
             typename =
                 std::enable_if_t<(sizeof...(I) > 0) && (sizeof...(I) == n) && (std::is_same_v<T, I> && ...)>>
    point(I... components)
    // requires((std::is_same_v<T, I> && ...) && sizeof...(components) > 0 && sizeof...(components) <= n)
    {
      const auto size_list = std::initializer_list<T>{static_cast<T>(components)...};
      std::copy(std::cbegin(size_list), std::cend(size_list), std::begin(_components));
    }

    template<typename = std::enable_if_t<(n > 0)>>
    point(const T component): _components() {
      _components.fill(component);
    }

    point(): _components() {};

    template<std::size_t index>
    [[nodiscard]] T& component() {
      return _components[index];
    }
    template<std::size_t index>
    [[nodiscard]] T component() const {
      return _components[index];
    }
    [[nodiscard]] const std::array<T, n>& get_component() const { return _components; }
    [[nodiscard]] const T* get_component_p() const { return _components.data(); }
    [[nodiscard]] inline T& x() { return component<0>(); }
    [[nodiscard]] inline T& y()
      requires(n > 1)
    {
      return component<1>();
    }
    [[nodiscard]] T& z()
      requires(n > 2)
    {
      return component<2>();
    }
    [[nodiscard]] inline T x() const { return component<0>(); }
    [[nodiscard]] inline T y() const
      requires(n > 1)
    {
      return component<1>();
    }
    [[nodiscard]] inline T z() const
      requires(n > 2)
    {
      return component<2>();
    }

    [[nodiscard]] bool operator<(const point<T, n>& other) const {
      return [&]<std::size_t... I>(const std::array<T, n>& a,
                                   const std::array<T, n>& b,
                                   std::index_sequence<I...>) {
        return ((a[I] < b[I]) || ...);
      }(_components, other._components, std::make_index_sequence<n>{});
    }
    // Note cross product make sense only for 3 and 7 dimensional space
    [[nodiscard]] point<T, n> cross(const point<T, n>& other) const
      requires(n == 3)
    {
      const auto a_x = this->x();
      const auto a_y = this->y();
      const auto a_z = this->z();
      const auto b_x = other.x();
      const auto b_y = other.y();
      const auto b_z = other.z();
      return point<T, n>{a_y * b_z - a_z * b_y, a_z * b_x - a_x * b_z, a_x * b_y - a_y * b_x};
    }
    [[nodiscard]] T sum_components() const {
      return [&]<std::size_t... I>(const std::array<T, n>& a, std::index_sequence<I...>) {
        return (a[I] + ...);
      }(_components, std::make_index_sequence<n>{});
    }
    [[nodiscard]] point<T, n> normalize() const {
      const fp_type square_distance =
          std::max(this->square().sum_components(), std::numeric_limits<fp_type>::epsilon());
      const fp_type inv_rd = fp_type{1} / std::sqrt(square_distance);
      return this->product(point<T, n>{inv_rd, inv_rd, inv_rd});
    }
    [[nodiscard]] point<T, n> normalize(const point<T, n>& other) {
      return normalize(this->difference(other));
    }
    [[nodiscard]] bool operator>(const point<T, n>& other) const { return other < this; }
    T distance(const point<T, n>& other) const {
      return std::sqrt((*this - other).square().sum_components());
    }
    [[nodiscard]] point<T, n> operator-(const point<T, n>& other) const {
      return [&]<std::size_t... I>(const std::array<T, n>& a,
                                   const std::array<T, n>& b,
                                   std::index_sequence<I...>) {
        return point<T, n>{(a[I] - b[I])...};
      }(_components, other._components, std::make_index_sequence<n>{});
    }

    [[nodiscard]] point<T, n> operator/(const point<T, n>& other) const {
      return [&]<std::size_t... I>(const std::array<T, n>& a,
                                   const std::array<T, n>& b,
                                   std::index_sequence<I...>) {
        return point<T, n>{(a[I] / b[I])...};
      }(_components, other._components, std::make_index_sequence<n>{});
    }
    [[nodiscard]] point<T, n> divide(const point<T, n>& other) const { return *this / other; }
    [[nodiscard]] point<T, n> add(const point<T, n>& other) const { return *this + other; }
    [[nodiscard]] point<T, n> difference(const point<T, n>& other) const { return *this - other; }
    [[nodiscard]] point<T, n> operator+(const point<T, n>& other) const {
      return [&]<std::size_t... I>(const std::array<T, n>& a,
                                   const std::array<T, n>& b,
                                   std::index_sequence<I...>) {
        return point<T, n>{(a[I] + b[I])...};
      }(_components, other._components, std::make_index_sequence<n>{});
    }
    [[nodiscard]] point<T, n> square() const {
      return [&]<std::size_t... I>(const std::array<T, n>& a, std::index_sequence<I...>) {
        return point<T, n>{(a[I] * a[I])...};
      }(_components, std::make_index_sequence<n>{});
    }

    [[nodiscard]] point<T, n> truncate() {
      return [&]<std::size_t... I>(const std::array<T, n>& a, std::index_sequence<I...>) {
        return point<T, n>{std::trunc(a[I])...};
      }(_components, std::make_index_sequence<n>{});
    }

    [[nodiscard]] point<T, n> operator*(const T scale) const {
      return [&]<std::size_t... I>(const std::array<T, n>& a, std::index_sequence<I...>) {
        return point<T, n>{(a[I] * scale)...};
      }(_components, std::make_index_sequence<n>{});
    }
    [[nodiscard]] point<T, n> operator*(const point<T, n>& other) const {
      return [&]<std::size_t... I>(const std::array<T, n>& a,
                                   const std::array<T, n>& b,
                                   std::index_sequence<I...>) {
        return point<T, n>{(a[I] * b[I])...};
      }(_components, other._components, std::make_index_sequence<n>{});
    }
    [[nodiscard]] T inner_product(const point<T, n>& other) const {
      return [&]<std::size_t... I>(const std::array<T, n>& a,
                                   const std::array<T, n>& b,
                                   std::index_sequence<I...>) {
        return ((a[I] * b[I]) + ...);
      }(_components, other._components, std::make_index_sequence<n>{});
    }
    template<typename A>
    [[nodiscard]] point<T, n> product(const A other) const {
      return *this * other;
    }

    [[nodiscard]] T distance2(const point<T, n>& other) const {
      return this->difference(other).square().sum_components();
    }

    [[nodiscard]] T angle(const point<T, n>& a, const point<T, n>& b) const {
      const auto v1 = a.difference(*this);
      const auto v2 = b.difference(*this);
      return std::acos(v1.inner_product(v2) / (v1.magnitude() * v2.magnitude()));
    }

    // this function consider the given point as the end of a vector that starts in the origin
    [[nodiscard]] T magnitude() const { return std::sqrt(this->square().sum_components()); }

    template<typename F>
    [[nodiscard]] point<T, n> apply(F&& f) const {
      // use a lambda + index_sequence to expand components
      return [&]<std::size_t... I>(std::index_sequence<I...>) {
        return point{T(f(_components[I]))...};
      }(std::make_index_sequence<n>{});
    }
  };

  using point3D = point<fp_type, 3>;
  // struct point3D: public point<fp_type, 3> {
  //   fp_type& x;
  //   fp_type& y;
  //   fp_type& z;
  //
  //   point3D(const fp_type _x, const fp_type _y, const fp_type _z)
  //       : point<fp_type, 3>(_x, _y, _z),
  //         x(this->_components[0]),
  //         y(this->_components[1]),
  //         z(this->_components[2]) {};
  //   point3D(const fp_type v)
  //       : point<fp_type, 3>(v, v, v),
  //         x(this->_components[0]),
  //         y(this->_components[1]),
  //         z(this->_components[2]) {};
  //   point3D()
  //       : point<fp_type, 3>(), x(this->_components[0]), y(this->_components[1]), z(this->_components[2]) {};
  //
  //   ~point3D() = default;
  //   point3D(point3D&& other)
  //       : point<fp_type, 3>(other),
  //         x(this->_components[0]),
  //         y(this->_components[1]),
  //         z(this->_components[2]) {};
  //   point3D(const point3D& other)
  //       : point<fp_type, 3>(other),
  //         x(this->_components[0]),
  //         y(this->_components[1]),
  //         z(this->_components[2]) {};
  //
  //   point3D& operator=(point3D&& other) {
  //     point::operator=(other);
  //     return *this;
  //   }
  //   point3D& operator=(const point3D& other) {
  //     point::operator=(other);
  //     return *this;
  //   }
  //
  //   inline std::array<fp_type, 3> get_array() const { return {x, y, z}; }
  // };
  // FIX ME make them deprected
  // template<class point_type>
  // [[nodiscard]] constexpr std::remove_reference_t<point_type> add(const point_type& a, const point_type& b) {
  //   return {a.x() + b.x(), a.y() + b.y(), a.z() + b.z()};
  // }

  // template<class point_type>
  // [[nodiscard]] constexpr std::remove_reference_t<point_type> divide(point_type&& a, point_type&& b) {
  //   return {a.x() / b.x(), a.y() / b.y(), a.z() / b.z()};
  // }
  //
  // template<class point_type>
  // [[nodiscard]] constexpr std::remove_reference_t<point_type> difference(point_type&& a, point_type&& b) {
  //   return {a.x() - b.x(), a.y() - b.y(), a.z() - b.z()};
  // }
  //
  // template<class point_type, class... point_types>
  // [[nodiscard]] constexpr std::remove_reference_t<point_type> difference(point_type&& a,
  //                                                                        point_types&&... bs) {
  //   return difference(std::forward<point_type>(a), std::forward<point_type>(bs...));
  // }
  //
  // template<class point_type>
  // [[nodiscard]] constexpr std::remove_reference_t<point_type> scalar_division(point_type&& a,
  //                                                                             const fp_type value) {
  //   return {a.x / value, a.y / value, a.z / value};
  // }
  //
  // template<class point_type>
  // [[nodiscard]] constexpr std::remove_reference_t<point_type> scale(point_type&& a, const fp_type value) {
  //   return {a.x() * value, a.y() * value, a.z() * value};
  // }
  //
  // template<class point_type>
  // [[nodiscard]] constexpr std::remove_reference_t<point_type> square(point_type&& point) {
  //   return {point.x() * point.x(), point.y() * point.y(), point.z() * point.z()};
  // }
  //
  // template<class point_type>
  // [[nodiscard]] constexpr fp_type sum_components(point_type&& point) {
  //   return point.x() + point.y() + point.z();
  // }
  //
  // template<class point_type>
  // [[nodiscard]] constexpr std::remove_reference_t<point_type> normalize(point_type&& diff) {
  //   fp_type square_distance = sum_components(square(diff));
  //   if (square_distance == fp_type{0}) {
  //     // TODO @Davide we should log something
  //     // error("Attempt to divide by zero was just prevented.");
  //     square_distance = std::numeric_limits<fp_type>::epsilon();
  //   }
  //   // TODO ask @Davide about this
  //   fp_type inv_rd = fp_type{1} / std::sqrt(square_distance);
  //   return scale(diff, inv_rd);
  // }
  //
  // // TODO fix the constexpr here
  // template<class point_type>
  // [[nodiscard]] std::remove_reference_t<point_type> normalize(point_type&& a, point_type&& b) {
  //   point3D diff = difference(point3D{a.x, a.y, a.z}, point3D{b.x, b.y, b.z});
  //   return normalize(diff);
  // }
  //
  // template<class point_type>
  // [[nodiscard]] constexpr std::remove_reference_t<point_type> cross_product(point_type&& a, point_type&& b) {
  //   return point3D{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
  // }
  //
  // template<class point_type>
  // [[nodiscard]] constexpr fp_type distance2(point_type&& a, point_type&& b) {
  //   return sum_components(square(difference(std::forward<point_type>(a), std::forward<point_type>(b))));
  // }
  //
  // template<class point_type>
  // [[nodiscard]] constexpr fp_type distance(point_type&& a, point_type&& b) {
  //   return std::sqrt(
  //       sum_components(square(difference(std::forward<point_type>(a), std::forward<point_type>(b)))));
  // }
  //
  // template<class point_type>
  // [[nodiscard]] constexpr fp_type inner_product(point_type&& a, point_type&& b) {
  //   return a.x() * b.x() + a.y() * b.y() + a.z() * b.z();
  // }
  //
  // template<class point_type>
  // [[nodiscard]] constexpr std::remove_reference_t<point_type> product(point_type&& a, point_type&& b) {
  //   return {a.x() * b.x(), a.y() * b.y(), a.z() * b.z()};
  // }
  //
  // template<class point_type, class... point_types>
  // [[nodiscard]] constexpr std::remove_reference_t<point_type> product(point_type&& a, point_types&&... bs) {
  //   return product(std::forward<point_type>(a), std::forward<point_type>(bs...));
  // }
  //
  // // this function consider the given point as the end of a vector that starts in the origin
  // template<class point_type>
  // [[nodiscard]] constexpr fp_type magnitude(point_type&& a) {
  //   return std::sqrt(sum_components(square(std::forward<point_type>(a))));
  // }
  //
  // template<class point_type>
  // [[nodiscard]] constexpr fp_type angle(point_type&& origin, point_type&& a, point_type&& b) {
  //   auto v1 = difference(std::forward<point_type>(a), std::forward<point_type>(origin));
  //   auto v2 = difference(std::forward<point_type>(b), std::forward<point_type>(origin));
  //   return std::acos(inner_product(v1, v2) / (magnitude(v1) * magnitude(v2)));
  // }

} // namespace mudock
