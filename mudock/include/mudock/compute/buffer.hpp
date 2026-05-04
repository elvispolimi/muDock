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

#ifdef MUDOCK_USE_ALPAKA
  #include <alpaka/alpaka.hpp>
  #include <cstdint>
  #include <mudock/alpaka_implementation/queue_alpaka.hpp>
  #include <optional>
#endif

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

  template<typename T, typename tuple_t>
  struct is_in_tuple: std::false_type {};

  template<typename T, typename... Ts>
  struct is_in_tuple<T, std::tuple<Ts...>>: std::bool_constant<(std::same_as<T, Ts> || ...)> {};

  template<typename T>
  concept buffer_type_allowed = is_in_tuple<T, buffer_type_list>::value;

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
    std::unique_ptr<object<T, queue_t>> obj;
    container_type<T, args...> host;
    std::shared_ptr<queue_t> q;
    bool valid;

  public:
    buffer_impl(std::shared_ptr<queue_t> _queue, const std::size_t num_elements = 0): valid(false) {
      q = _queue;
      if (q->obj_required()) {
        obj = std::make_unique<object<T, queue_t>>(_queue);
      }
      if (num_elements)
        this->alloc(num_elements);
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

#ifdef MUDOCK_USE_ALPAKA
  template<template<class...> class container_type, typename T, class... args>
  struct buffer_impl<container_type, T, queue_alpaka, args...> {
  private:
    using dim         = queue_alpaka::dim;
    using idx         = queue_alpaka::idx;
    using dev_acc     = queue_alpaka::dev_acc;
    using buffer_type = alpaka::Buf<dev_acc, T, dim, idx>;
    using byte_type   = std::uint8_t;

    container_type<T, args...> host;
    std::shared_ptr<queue_alpaka> q;
    std::optional<buffer_type> buffer;
    T* ptr                 = nullptr;
    std::size_t size       = 0;
    std::size_t alloc_size = 0;
    bool valid             = false;

    static auto extent(const std::size_t num_elements) {
      return alpaka::Vec<dim, idx>{static_cast<idx>(num_elements)};
    }

    static auto host_device() { return alpaka::getDevByIdx(alpaka::PlatformCpu{}, 0); }

    void alloc_device(const std::size_t num_elements) {
      if (num_elements > alloc_size) {
        q->synchronize();
        buffer.emplace(alpaka::allocBuf<T, idx>(q->native_device(), extent(num_elements)));
        ptr        = alpaka::getPtrNative(*buffer);
        alloc_size = num_elements;
      }
      size = num_elements;
    }

  public:
    buffer_impl(std::shared_ptr<queue_alpaka> _queue, const std::size_t num_elements = 0): q(_queue) {
      if (num_elements)
        this->alloc(num_elements);
    }
    buffer_impl(const buffer_impl& other) = delete;

    buffer_impl& operator=(const buffer_impl& other) = delete;

    buffer_impl(buffer_impl&&) = default;

    buffer_impl& operator=(buffer_impl&&) = default;

    bool is_valid() { return valid; }
    void set_valid() { valid = true; }
    void set_not_valid() { valid = false; }

    ~buffer_impl() {
      if (buffer) {
        q->synchronize();
      }
    }
    [[nodiscard]] inline auto host_pointer() const { return host.data(); }
    [[nodiscard]] inline auto host_pointer() { return host.data(); }
    [[nodiscard]] inline auto num_elements() const { return host.size(); };

    auto* operator()() { return host.data(); }

    inline void alloc(const std::size_t num_elements) {
      host.resize(num_elements);
      alloc_device(num_elements);
    };
    inline void alloc(const std::size_t num_elements, const int value) {
      host.resize(num_elements, value);
      alloc_device(num_elements);

      const auto bytes = size * sizeof(T);
      auto bytes_host  = std::vector<byte_type>(bytes, static_cast<byte_type>(value));
      auto host_view   = alpaka::createView(host_device(),
                                          bytes_host.data(),
                                          alpaka::Vec<dim, idx>{static_cast<idx>(bytes)});
      auto device_view = alpaka::createView(q->native_device(),
                                            reinterpret_cast<byte_type*>(ptr),
                                            alpaka::Vec<dim, idx>{static_cast<idx>(bytes)});
      alpaka::memcpy(q->native_queue(), device_view, host_view);
    };

    inline void copy_host2device(const std::size_t copy_size = 0) {
      valid        = true;
      const auto n = copy_size ? copy_size : host.size();
      alloc_device(n);
      auto host_view   = alpaka::createView(host_device(), host.data(), extent(n));
      auto device_view = alpaka::createView(q->native_device(), ptr, extent(n));
      alpaka::memcpy(q->native_queue(), device_view, host_view);
    };
    inline void copy_device2host() {
      valid            = false;
      auto device_view = alpaka::createView(q->native_device(), ptr, extent(size));
      auto host_view   = alpaka::createView(host_device(), host.data(), extent(size));
      alpaka::memcpy(q->native_queue(), host_view, device_view);
    };
    inline void copy_device2device(const buffer_impl<container_type, T, queue_alpaka, args...>& other,
                                   const int n = 0) {
      valid = true;
      alloc_device(other.num_elements());
      host.resize(other.num_elements());

      const auto copy_size = n ? static_cast<std::size_t>(n) : size;
      auto src_view        = alpaka::createView(q->native_device(), other.ptr, extent(copy_size));
      auto dest_view       = alpaka::createView(q->native_device(), ptr, extent(copy_size));
      alpaka::memcpy(q->native_queue(), dest_view, src_view);
    }

    [[nodiscard]] inline auto dev_pointer() { return ptr; }

    [[nodiscard]] T** dev_pointer_ref() { return &ptr; }
  };
#endif

  template<typename T, typename queue_t>
  using buffer_vector = buffer_impl<std::vector, T, queue_t>;

} // namespace mudock
