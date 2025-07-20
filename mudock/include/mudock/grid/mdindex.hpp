#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <concepts>
#include <initializer_list>
#include <mudock/grid/point3D.hpp>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <type_traits>

namespace mudock {
  template<typename T>
  concept Numeric = std::integral<std::remove_cvref_t<T>>;
  /**
   * This class converts a multidimensional index to a flat index and viceversa. The spinning direction is
   * the lefmost index, i.e. we use the following formula to compute the multindex
   *   flat_index = i0 + i1*d0 + i2*d0*d1 + i3*d0*d1*d2 + ... + in*d0*d1*d2*...*d(n-1)
   */
  template<std::size_t n>
  class md_index {
    static_assert(n > 0, "A multi dimensional index must at least a dimension");

    // the storage of sizes and coefficients for dealing with md indexes
    std::array<std::size_t, n> _sizes;
    std::array<std::size_t, n> _coefs;

  public:
    static constexpr auto num_dimensions = n;

    // by default we assume that we need a single element
    md_index() {
      _sizes.fill(1);
      _coefs.fill(1);
    }
    template<typename T>
    md_index(const point<T, n>& point) {
      const auto& size_list = point.get_component();
      std::copy(std::cbegin(size_list), std::cend(size_list), std::begin(_sizes));
      _coefs[0] = std::size_t{1};
      for (std::size_t i = 1; i < num_dimensions; ++i) { _coefs[i] = _coefs[i - 1] * _sizes[i - 1]; }
    }
    ~md_index()                                = default;
    md_index(md_index&& other)                 = default;
    md_index(const md_index& other)            = default;
    md_index& operator=(md_index&& other)      = default;
    md_index& operator=(const md_index& other) = default;

    // otherwise we initialize the index as requested by the user
    template<Numeric... T>
    md_index(T&&... sizes) {
      static_assert(sizeof...(sizes) == n, "Mismatch between sizes and dimension numbers");
      const auto size_list = std::initializer_list<std::size_t>{static_cast<std::size_t>(sizes)...};
      std::copy(std::cbegin(size_list), std::cend(size_list), std::begin(_sizes));
      _coefs[0] = std::size_t{1};
      for (std::size_t i = 1; i < num_dimensions; ++i) { _coefs[i] = _coefs[i - 1] * _sizes[i - 1]; }
    }

    // function to perform the index conversion
    template<Numeric... T>
    [[nodiscard]] std::size_t to1D(T&&... indexes) const {
      static_assert(sizeof...(indexes) == n, "Mismatch between indexes and dimension numbers");
      assert(is_inside(indexes...));
      const auto index_list  = std::initializer_list<std::size_t>{static_cast<std::size_t>(indexes)...};
      auto index_it          = std::begin(index_list);
      const auto begin_coefs = std::begin(_coefs);
      return std::accumulate(
          begin_coefs + 1,
          std::end(_coefs),
          *index_it,
          [&index_it](const auto sum, const auto coef) { return sum + coef * (*++index_it); });
    }
    [[nodiscard]] auto toND(const std::size_t index) const {
      assert(index < flat_size());
      auto result = std::array<std::size_t, n>{};
      result[0]   = index % _sizes[0];
      for (std::size_t i = 0; i < n; ++i) { result[i] = (index / _coefs[i]) % _sizes[i]; }
      return result;
    }

    [[nodiscard]] bool operator==(const md_index<n>& other) const {
      return [&]<std::size_t... I>(const std::array<std::size_t, n>& a,
                                   const std::array<std::size_t, n>& b,
                                   std::index_sequence<I...>) {
        return ((a[I] == b[I]) && ...);
      }(_sizes, other._sizes, std::make_index_sequence<n>{});
    }

    // utility functions to get sizes in the index
    template<std::size_t index>
    [[nodiscard]] std::size_t size() const {
      return _sizes[index];
    }
    [[nodiscard]] std::size_t size_x() const { return _sizes[0]; }

    template<typename = std::enable_if<(n > 1)>>
    [[nodiscard]] std::size_t size_y() const {
      return _sizes[1];
    }
    template<typename = std::enable_if<(n > 2)>>
    [[nodiscard]] std::size_t size_z() const {
      return _sizes[2];
    }

    template<typename = std::enable_if<(n > 2)>>
    [[nodiscard]] std::size_t size_xy() const {
      return _coefs[2];
    }

    [[nodiscard]] std::size_t flat_size() const { return _sizes[n - 1] * _coefs[n - 1]; }

    // utility function to perform boundaries check
    template<Numeric... T>
    [[nodiscard]] auto is_inside(T&&... indexes) const {
      static_assert(sizeof...(indexes) == n, "Mismatch between indexes and dimension numbers");
      const auto index_list = std::initializer_list<std::size_t>{static_cast<std::size_t>(indexes)...};
      return std::none_of(
          std::begin(index_list),
          std::end(index_list),
          [size_it = std::cbegin(_sizes)](const auto index) mutable { return index >= *size_it++; });
    }
  };

  //===------------------------------------------------------------------------------------------------------
  // Helper class to convert a 1D index to 2D
  //===------------------------------------------------------------------------------------------------------

  struct index {
    virtual int to1D(const int index_x, const int index_y, const int index_z) const = 0;
    virtual int to1D(const int index_x, const int index_y) const                    = 0;
    virtual int get_dim() const                                                     = 0;

    virtual ~index() = default;
  };

  class index2D: public index {
    int n_x = 0;
    int n_y = 0;

  public:
    inline index2D(const int size_x = 0, const int size_y = 0): n_x(size_x), n_y(size_y) {}
    [[nodiscard]] inline int to1D(const int index_x, const int index_y) const override {
      assert(index_x < n_x);
      assert(index_y < n_y);
      return index_y * n_x + index_x;
    }
    [[nodiscard]] inline int to1D([[maybe_unused]] const int index_x,
                                  [[maybe_unused]] const int index_y,
                                  [[maybe_unused]] const int index_z) const override {
      throw std::logic_error("This method is not applicable for 2D indexing.");
    }
    [[nodiscard]] inline auto to2D(const int index1D) const {
      assert(index1D < n_x * n_y);
      return std::make_tuple(index1D % n_x, index1D / n_x);
    }
    [[nodiscard]] inline auto size_x() const { return n_x; }
    [[nodiscard]] inline auto size_y() const { return n_y; }

    [[nodiscard]] inline int get_dim() const override { return n_x * n_y; }
  };

  //===------------------------------------------------------------------------------------------------------
  // Helper class to convert a 1D index to 3D
  //===------------------------------------------------------------------------------------------------------

  class index3D: public index {
    int n_x  = 0;
    int n_y  = 0;
    int n_z  = 0;
    int n_xy = 0;

  public:
    inline index3D(const int size_x = 0, const int size_y = 0, const int size_z = 0)
        : n_x(size_x), n_y(size_y), n_z(size_z), n_xy(size_x * size_y) {}
    [[nodiscard]] inline int to1D([[maybe_unused]] const int index_x,
                                  [[maybe_unused]] const int index_y) const override {
      throw std::logic_error("This method is not applicable for 3D indexing.");
    }

    bool operator==(const index3D& d) const { return n_x == d.n_x && n_y == d.n_y && n_z == d.n_z; }

    [[nodiscard]] inline int to1D(const int index_x, const int index_y, const int index_z) const override {
      assert(index_x < n_x);
      assert(index_y < n_y);
      assert(index_z < n_z);
      return n_xy * index_z + index_y * n_x + index_x;
    }
    [[nodiscard]] inline auto to3D(const int index1D) const {
      assert(index1D < n_x * n_y * n_z);
      return std::make_tuple(index1D % n_x, (index1D % n_xy) / n_y, index1D / n_xy);
    }
    [[nodiscard]] inline auto size_x() const { return n_x; }
    [[nodiscard]] inline auto size_y() const { return n_y; }
    [[nodiscard]] inline auto size_z() const { return n_z; }
    [[nodiscard]] inline auto size_xy() const { return n_xy; }

    [[nodiscard]] inline int get_dim() const override { return n_xy * n_z; }
  };

  template<class T>
  concept is_index =
      (std::same_as<std::remove_cvref_t<T>, index2D> || std::same_as<std::remove_cvref_t<T>, index3D>);

} // namespace mudock
