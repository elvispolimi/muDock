#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <tuple>

namespace mudock {
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
    inline index2D(const std::size_t size_x, const std::size_t size_y): n_x(size_x), n_y(size_y) {}
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
    inline index3D(const std::size_t size_x, const std::size_t size_y, const std::size_t size_z)
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
