#pragma once

#include <concepts>
#include <cstring>
#include <memory>
#include <mudock/chem/elements.hpp>
#include <mudock/compute/object.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/type_alias.hpp>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <vector>

namespace mudock {

  enum class buffer_data_type {
    X_COORDS,
    Y_COORDS,
    Z_COORDS,
    ELEMENTS,
    CHROMOSOMES,
    SCORES,
    X_SCRATCH,
    Y_SCRATCH,
    Z_SCRATCH,
    NUM_ATOMS,
    NUM_ROTAMERS,
    PROT_MIN,
    PROT_MAX,
    PROT_CENTER,
    PROT_SIZE_X,
    PROT_SIZE_XY,
    PROT_SIZE_XYZ,
    PROT_GRID_MAPS
  };

  using buffer_type_list = std::tuple<fp_type, int, chromosome>;

  template<typename T, typename... Ts>
  constexpr bool is_in_tuple_v = (std::same_as<T, Ts> || ...);

  template<typename T>
  concept buffer_type_allowed =
      []<typename... Ts>(std::tuple<Ts...>*) { return is_in_tuple_v<T, Ts...>; }((buffer_type_list*) nullptr);

  template<buffer_data_type buff_t, class T>
    requires buffer_type_allowed<T>
  struct buffer_type_traits_impl {
    using type = T;
  };

  template<buffer_data_type buff_t>
  struct buffer_type_traits {};
  template<>
  struct buffer_type_traits<buffer_data_type::X_COORDS> {
    using type = buffer_type_traits_impl<buffer_data_type::X_COORDS, fp_type>::type;
  };
  template<>
  struct buffer_type_traits<buffer_data_type::Y_COORDS> {
    using type = buffer_type_traits_impl<buffer_data_type::Y_COORDS, fp_type>::type;
  };
  template<>
  struct buffer_type_traits<buffer_data_type::Z_COORDS> {
    using type = buffer_type_traits_impl<buffer_data_type::Z_COORDS, fp_type>::type;
  };
  template<>
  struct buffer_type_traits<buffer_data_type::ELEMENTS> {
    using type = buffer_type_traits_impl<buffer_data_type::ELEMENTS, int>::type;
  };
  template<>
  struct buffer_type_traits<buffer_data_type::CHROMOSOMES> {
    using type = buffer_type_traits_impl<buffer_data_type::CHROMOSOMES, chromosome>::type;
  };
  template<>
  struct buffer_type_traits<buffer_data_type::SCORES> {
    using type = buffer_type_traits_impl<buffer_data_type::SCORES, fp_type>::type;
  };
  template<>
  struct buffer_type_traits<buffer_data_type::X_SCRATCH> {
    using type = buffer_type_traits_impl<buffer_data_type::X_SCRATCH, fp_type>::type;
  };
  template<>
  struct buffer_type_traits<buffer_data_type::Y_SCRATCH> {
    using type = buffer_type_traits_impl<buffer_data_type::Y_SCRATCH, fp_type>::type;
  };
  template<>
  struct buffer_type_traits<buffer_data_type::Z_SCRATCH> {
    using type = buffer_type_traits_impl<buffer_data_type::Z_SCRATCH, fp_type>::type;
  };
  template<>
  struct buffer_type_traits<buffer_data_type::NUM_ATOMS> {
    using type = buffer_type_traits_impl<buffer_data_type::NUM_ATOMS, int>::type;
  };
  template<>
  struct buffer_type_traits<buffer_data_type::NUM_ROTAMERS> {
    using type = buffer_type_traits_impl<buffer_data_type::NUM_ROTAMERS, int>::type;
  };
  template<>
  struct buffer_type_traits<buffer_data_type::PROT_MAX> {
    using type = buffer_type_traits_impl<buffer_data_type::PROT_MAX, fp_type>::type;
  };
  template<>
  struct buffer_type_traits<buffer_data_type::PROT_MIN> {
    using type = buffer_type_traits_impl<buffer_data_type::PROT_MIN, fp_type>::type;
  };
  template<>
  struct buffer_type_traits<buffer_data_type::PROT_CENTER> {
    using type = buffer_type_traits_impl<buffer_data_type::PROT_CENTER, fp_type>::type;
  };
  template<>
  struct buffer_type_traits<buffer_data_type::PROT_SIZE_X> {
    using type = buffer_type_traits_impl<buffer_data_type::PROT_SIZE_X, int>::type;
  };
  template<>
  struct buffer_type_traits<buffer_data_type::PROT_SIZE_XY> {
    using type = buffer_type_traits_impl<buffer_data_type::PROT_SIZE_XY, int>::type;
  };
  template<>
  struct buffer_type_traits<buffer_data_type::PROT_SIZE_XYZ> {
    using type = buffer_type_traits_impl<buffer_data_type::PROT_SIZE_XYZ, int>::type;
  };
  template<>
  struct buffer_type_traits<buffer_data_type::PROT_GRID_MAPS> {
    using type = buffer_type_traits_impl<buffer_data_type::PROT_SIZE_XYZ, fp_type>::type;
  };

  template<template<class...> class container_type, typename T, class queue_t, class... args>
  struct buffer_impl {
  private:
    std::unique_ptr<object<T>> obj;
    container_type<T, args...> host;
    std::shared_ptr<queue_t> q;
    bool valid;

  public:
    buffer_impl(std::shared_ptr<queue_t> _queue, const int dim = 0): valid(false) {
      q = _queue;
      if (q->obj_required()) {
        obj = std::make_unique<object<T>>(_queue);
      }
      if (!dim)
        this->alloc(dim);
    }
    buffer_impl(const buffer_impl& other) = delete;

    buffer_impl& operator=(const buffer_impl& other) = delete;

    // Move
    buffer_impl(buffer_impl&&) = default;

    buffer_impl& operator=(buffer_impl&&) = default;

    bool is_valid() { return valid; }
    void set_valid() { valid = true; }
    void set_not_valid() { valid = false; }

    ~buffer_impl() = default;
    [[nodiscard]] inline auto host_pointer() const { return host.data(); }
    [[nodiscard]] inline auto host_pointer() { return host.data(); }
    [[nodiscard]] inline auto num_elements() const { return host.size(); };

    auto* operator()() { return host.data(); }

    inline void alloc(const std::size_t num_elements) {
      host.resize(num_elements);
      if (obj)
        (*obj).alloc(num_elements);
    };
    inline void alloc(const std::size_t num_elements, const int value) {
      host.resize(num_elements, value);
      if (obj) {
        (*obj).alloc(num_elements);
        (*obj).set_to_value(value);
      }
    };

    inline void copy_host2device(const std::size_t copy_size = 0) {
      valid = true;
      if (obj) {
        (*obj).alloc(copy_size ? copy_size : host.size());
        (*obj).copy_host2device(host.data(), copy_size);
      }
    };
    inline void copy_device2host() {
      valid = false;
      if (obj)
        (*obj).copy_device2host(host.data());
    };
    inline void copy_device2device(const buffer_impl<container_type, T, queue_t, args...>& other,
                                   const int n = 0) {
      valid = true;
      if (obj && other.obj) {
        (*obj).alloc(other.num_elements());
        (*obj).copy_device2device(*(other.obj), n);
        host.resize(other.num_elements());
      } else if (!obj && !other.obj) {
        host = other.host;
      } else {
        throw std::runtime_error("Device to device copy with one of the two buffers without object");
      }
    }

    [[nodiscard]] inline auto dev_pointer() {
      if (obj)
        return (*obj).dev_pointer();
      else
        return host_pointer();
    }

    [[nodiscard]] T** dev_pointer_ref() {
      if (obj)
        return (*obj).dev_pointer_ref();
      else
        throw std::runtime_error("Reference to a buffer host pointer not yet implemented");
    };
  };

  template<typename T, typename queue_t>
  using buffer_vector = buffer_impl<std::vector, T, queue_t>;

} // namespace mudock
