#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <initializer_list>
#include <numeric>
#include <stdexcept>
#include <tuple>

namespace mudock {

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

    // otherwise we initialize the index as requested by the user
    template<class... T>
    md_index(T&&... sizes) {
      static_assert(sizeof...(sizes) == n, "Mismatch between sizes and dimension numbers");
      const auto size_list = std::initializer_list{static_cast<std::size_t>(sizes)...};
      std::copy(std::cbegin(size_list), std::cend(size_list), std::begin(_sizes));
      _coefs[0] = std::size_t{1};
      for (std::size_t i = 1; i < num_dimensions; ++i) { _coefs[i] = _coefs[i - 1] * _sizes[i - 1]; }
    }

    // function to perform the index conversion
    template<class... T>
    [[nodiscard]] std::size_t to1D(T&&... indexes) const {
      static_assert(sizeof...(indexes) == n, "Mismatch between indexes and dimension numbers");
      assert(is_inside(indexes...));
      const auto index_list  = std::initializer_list{static_cast<std::size_t>(indexes)...};
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

    // utility functions to get sizes in the index
    template<std::size_t index>
    [[nodiscard]] std::size_t size() const {
      return _sizes[index];
    }
    [[nodiscard]] std::size_t flat_size() const { return _sizes[0] * _coefs[n - 1]; }

    // utility function to perform boundaries check
    template<class... T>
    [[nodiscard]] auto is_inside(T&&... indexes) const {
      static_assert(sizeof...(indexes) == n, "Mismatch between indexes and dimension numbers");
      const auto index_list = std::initializer_list{static_cast<std::size_t>(indexes)...};
      return std::none_of(std::begin(index_list),
                          std::end(index_list),
                          [coef_it = std::begin(_coefs)](const auto index) { return index >= *coef_it++; });
    }
  };

  //===------------------------------------------------------------------------------------------------------
  // Helper class to convert a 1D index to 2D
  //===------------------------------------------------------------------------------------------------------

  struct index {
    virtual std::size_t
        to1D(const std::size_t index_x, const std::size_t index_y, const std::size_t index_z) const = 0;
    virtual std::size_t to1D(const std::size_t index_x, const std::size_t index_y) const            = 0;
    virtual std::size_t get_dim() const                                                             = 0;

    virtual ~index() = default;
  };

  class index2D: public index {
    std::size_t n_x = 0;
    std::size_t n_y = 0;

  public:
    inline index2D(const std::size_t size_x = 0, const std::size_t size_y = 0): n_x(size_x), n_y(size_y) {}
    [[nodiscard]] inline std::size_t to1D(const std::size_t index_x,
                                          const std::size_t index_y) const override {
      assert(index_x < n_x);
      assert(index_y < n_y);
      return index_y * n_x + index_x;
    }
    [[nodiscard]] inline std::size_t to1D([[maybe_unused]] const std::size_t index_x,
                                          [[maybe_unused]] const std::size_t index_y,
                                          [[maybe_unused]] const std::size_t index_z) const override {
      throw std::logic_error("This method is not applicable for 2D indexing.");
    }
    [[nodiscard]] inline auto to2D(const std::size_t index1D) const {
      assert(index1D < n_x * n_y);
      return std::make_tuple(index1D % n_x, index1D / n_x);
    }
    [[nodiscard]] inline auto size_x() const { return n_x; }
    [[nodiscard]] inline auto size_y() const { return n_y; }

    [[nodiscard]] inline std::size_t get_dim() const override { return n_x * n_y; }
  };

  //===------------------------------------------------------------------------------------------------------
  // Helper class to convert a 1D index to 3D
  //===------------------------------------------------------------------------------------------------------

  class index3D: public index {
    std::size_t n_x  = 0;
    std::size_t n_y  = 0;
    std::size_t n_z  = 0;
    std::size_t n_xy = 0;

  public:
    inline index3D(const std::size_t size_x = 0, const std::size_t size_y = 0, const std::size_t size_z = 0)
        : n_x(size_x), n_y(size_y), n_z(size_z), n_xy(size_x * size_y) {}
    [[nodiscard]] inline std::size_t to1D([[maybe_unused]] const std::size_t index_x,
                                          [[maybe_unused]] const std::size_t index_y) const override {
      throw std::logic_error("This method is not applicable for 3D indexing.");
    }
    [[nodiscard]] inline std::size_t
        to1D(const std::size_t index_x, const std::size_t index_y, const std::size_t index_z) const override {
      assert(index_x < n_x);
      assert(index_y < n_y);
      assert(index_z < n_z);
      return n_xy * index_z + index_y * n_x + index_x;
    }
    [[nodiscard]] inline auto to3D(const std::size_t index1D) const {
      assert(index1D < n_x * n_y * n_z);
      return std::make_tuple(index1D % n_x, (index1D % n_xy) / n_y, index1D / n_xy);
    }
    [[nodiscard]] inline auto size_x() const { return n_x; }
    [[nodiscard]] inline auto size_y() const { return n_y; }
    [[nodiscard]] inline auto size_z() const { return n_z; }

    [[nodiscard]] inline std::size_t get_dim() const override { return n_xy * n_z; }
  };

  template<class T>
  concept is_index =
      (std::same_as<std::remove_cvref_t<T>, index2D> || std::same_as<std::remove_cvref_t<T>, index3D>);

} // namespace mudock
