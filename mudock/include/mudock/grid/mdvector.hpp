#pragma once

#include <cstdint>
#include <mudock/grid/mdindex.hpp>
#include <vector>

namespace mudock {

  template<class T, std::size_t n>
  class md_vector: public md_index<n> {
    std::vector<T> data;

  public:
    // data constructor
    md_vector(): md_index<n>() {}
    template<class... Y>
    md_vector(Y&&... sizes): md_index<n>(sizes...), data(md_index<n>::flat_size()) {}

    // operator to retrieve data
    template<class... Y>
    auto& get(Y&&... indexes) {
      return data[md_index<n>::to1D(indexes...)];
    }
    template<class... Y>
    const auto& get(Y&&... indexes) const {
      return data[md_index<n>::to1D(indexes...)];
    }
  };

} // namespace mudock
