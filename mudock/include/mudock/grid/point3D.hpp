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

} // namespace mudock
