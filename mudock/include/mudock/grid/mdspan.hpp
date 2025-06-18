#pragma once

#include <mudock/grid/mdindex.hpp>
#include <vector>

namespace mudock {
  template<typename T>
  concept owns_data = requires(T t, std::size_t i) {
    { T(i) };                 // constructible from size
    { t[i] };                 // indexable
  } && !std::is_pointer_v<T>; // not a raw pointer

  /**
   * This class represents a multidimensional vector that adapt a std::vector with more dimensions
   */
  template<class T, std::size_t n>
  class md_span: public md_index<n> {
  protected:
    T* _data;

  public:
    md_span()                                = default;
    ~md_span()                               = default;
    md_span(md_span&& other)                 = default;
    md_span(const md_span& other)            = default;
    md_span& operator=(md_span&& other)      = default;
    md_span& operator=(const md_span& other) = default;

    template<class... Y>
    md_span(T* data, Y&&... sizes): md_index<n>(sizes...), _data(data) {}
    // operator to retrieve data
    template<class... Y>
    [[nodiscard]] auto& get(Y&&... indexes) {
      return _data[md_index<n>::to1D(indexes...)];
    }
    template<class... Y>
    [[nodiscard]] const auto& get(Y&&... indexes) const {
      return _data[md_index<n>::to1D(indexes...)];
    }
    [[nodiscard]] T* data() const { return _data; }
  };

  template<class T, std::size_t n>
    requires owns_data<T>
  class md_container: public md_index<n> {
  protected:
    T _data;

  public:
    md_container() = default;

    ~md_container()                                    = default;
    md_container(md_container&& other)                 = default;
    md_container(const md_container& other)            = default;
    md_container& operator=(md_container&& other)      = default;
    md_container& operator=(const md_container& other) = default;

    template<class... Y>
    md_container(Y&&... sizes): md_index<n>(sizes...), _data(md_index<n>::flat_size()) {}
    // operator to retrieve data
    template<class... Y>
    [[nodiscard]] auto& get(Y&&... indexes) {
      return _data[md_index<n>::to1D(indexes...)];
    }
    template<class... Y>
    [[nodiscard]] const auto& get(Y&&... indexes) const {
      return _data[md_index<n>::to1D(indexes...)];
    }

    template<std::size_t a, std::size_t b>
    [[nodiscard]] auto get_slice(const md_index<a> begin, const md_index<b> size) {
      return md_span<typename T::value_type, b>(_data.data() + begin.flat_size(), size);
    }
    template<std::size_t a, std::size_t b>
    [[nodiscard]] auto get_slice(const md_index<a> begin, const md_index<b> size) const {
      return md_span<const typename T::value_type, b>(_data.data() + begin.flat_size(), size);
    }
    [[nodiscard]] auto data() const { return _data.data(); }
  };

  template<class T, std::size_t n>
  using md_vector = md_container<std::vector<T>, n>;
} // namespace mudock
