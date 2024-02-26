#pragma once

#include <cassert>
#include <cstdint>
#include <tuple>

namespace mudock {

  //===------------------------------------------------------------------------------------------------------
  // Helper class to convert a 1D index to 2D
  //===------------------------------------------------------------------------------------------------------

  class index2D {
    std::size_t n_x = 0;
    std::size_t n_y = 0;

  public:
    inline index2D(const std::size_t size_x, const std::size_t size_y): n_x(size_x), n_y(size_y) {}
    [[nodiscard]] inline auto to1D(const std::size_t index_x, const std::size_t index_y) const {
      assert(index_x < n_x);
      assert(index_y < n_y);
      return index_y * n_x + index_x;
    }
    [[nodiscard]] inline auto to2D(const std::size_t index1D) const {
      return std::make_tuple(index1D % n_x, index1D / n_x);
    }
    [[nodiscard]] inline auto size_x() const { return n_x; }
    [[nodiscard]] inline auto size_y() const { return n_y; }
  };

  //===------------------------------------------------------------------------------------------------------
  // Helper class to convert a 1D index to 3D
  //===------------------------------------------------------------------------------------------------------

  class index3D {
    std::size_t n_x  = 0;
    std::size_t n_y  = 0;
    std::size_t n_z  = 0;
    std::size_t n_xy = 0;

  public:
    inline index3D(const std::size_t size_x, const std::size_t size_y, const std::size_t size_z)
        : n_x(size_x), n_y(size_y), n_z(size_z), n_xy(size_x * size_y) {}
    [[nodiscard]] inline auto
        to1D(const std::size_t index_x, const std::size_t index_y, const std::size_t index_z) const {
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
  };

} // namespace mudock
