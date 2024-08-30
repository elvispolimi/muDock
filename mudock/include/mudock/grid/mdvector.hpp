#pragma once

#include <cstdint>
#include <mudock/grid/mdindex.hpp>
#include <vector>

namespace mudock {

  /**
   * This class represents a multidimensional vector that adapt a std::vector with more dimensions
   */
  template<class T, std::size_t n>
  class md_vector: public md_index<n> {
  protected:
    std::vector<T> _data;

  public:
    // data constructor
    md_vector(): md_index<n>() {}
    template<class... Y>
    md_vector(Y&&... sizes): md_index<n>(sizes...), _data(md_index<n>::flat_size()) {}

    // operator to retrieve data
    template<class... Y>
    [[nodiscard]] auto& get(Y&&... indexes) {
      return _data[md_index<n>::to1D(indexes...)];
    }
    template<class... Y>
    [[nodiscard]] const auto& get(Y&&... indexes) const {
      return _data[md_index<n>::to1D(indexes...)];
    }
  };

} // namespace mudock
