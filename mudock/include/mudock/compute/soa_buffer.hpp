#pragma once

#include <array>
#include <cstddef>
#include <cstring>
#include <memory>
#include <mudock/compute/object.hpp>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace mudock {

  namespace detail {

    template<typename SoA>
    struct soa_layout_traits {
      using size_type = typename SoA::size_type;
      static constexpr bool is_blocks = requires { SoA::blocksNumber; };

      template<typename... SizeArgs>
      static auto compute_data_size(const SizeArgs&... sizes) {
        if constexpr (is_blocks) {
          static_assert(sizeof...(SizeArgs) == static_cast<std::size_t>(SoA::blocksNumber));
          return SoA::computeDataSize(std::array<size_type, SoA::blocksNumber>{sizes...});
        } else {
          static_assert(sizeof...(SizeArgs) == 1);
          return SoA::computeDataSize(sizes...);
        }
      }

      template<typename... SizeArgs>
      static SoA make_layout(std::byte* data, const SizeArgs&... sizes) {
        if constexpr (is_blocks) {
          static_assert(sizeof...(SizeArgs) == static_cast<std::size_t>(SoA::blocksNumber));
          return SoA(data, std::array<size_type, SoA::blocksNumber>{sizes...});
        } else {
          static_assert(sizeof...(SizeArgs) == 1);
          return SoA(data, sizes...);
        }
      }
    };

  }  // namespace detail

  // Owns the host allocation and its layout/view.
  template<typename SoA>
  class soa_host_buffer {
  public:
    using size_type = typename SoA::size_type;
    using view_type = typename SoA::View;
    using size_result_type = decltype(std::declval<SoA&>().metadata().size());

    soa_host_buffer() = default;
    soa_host_buffer(const soa_host_buffer&) = delete;
    soa_host_buffer& operator=(const soa_host_buffer&) = delete;
    soa_host_buffer(soa_host_buffer&&) = default;
    soa_host_buffer& operator=(soa_host_buffer&&) = default;

    template<typename... SizeArgs>
    void allocate(const SizeArgs&... sizes) {
      storage_.resize(detail::soa_layout_traits<SoA>::compute_data_size(sizes...));
      soa_ = detail::soa_layout_traits<SoA>::make_layout(storage_.data(), sizes...);
      view_ = view_type{soa_};
    }

    [[nodiscard]] std::byte* data() { return storage_.data(); }
    [[nodiscard]] const std::byte* data() const { return storage_.data(); }
    [[nodiscard]] size_result_type size() const { return soa_.metadata().size(); }
    [[nodiscard]] SoA& soa() { return soa_; }
    [[nodiscard]] const SoA& soa() const { return soa_; }
    [[nodiscard]] view_type& view() { return view_; }
    [[nodiscard]] const view_type& view() const { return view_; }

  private:
    std::vector<std::byte> storage_;
    SoA soa_;
    view_type view_;
  };

  // Owns the device allocation and its layout/view. The queue supplies the backend operations.
  template<typename SoA, typename queue_type>
  class soa_device_buffer {
  public:
    using size_type = typename SoA::size_type;
    using view_type = typename SoA::View;
    using size_result_type = decltype(std::declval<SoA&>().metadata().size());

    explicit soa_device_buffer(std::shared_ptr<queue_type> queue)
    {
      if (!queue)
        throw std::invalid_argument("A SOA device buffer requires a queue");
      if (queue->obj_required())
        storage_ = std::make_unique<object<std::byte>>(std::move(queue));
    }

    soa_device_buffer(const soa_device_buffer&) = delete;
    soa_device_buffer& operator=(const soa_device_buffer&) = delete;
    soa_device_buffer(soa_device_buffer&&) = default;
    soa_device_buffer& operator=(soa_device_buffer&&) = default;

    template<typename... SizeArgs>
    void allocate(const SizeArgs&... sizes) {
      if (storage_)
        storage_->alloc(detail::soa_layout_traits<SoA>::compute_data_size(sizes...));
      else
        fallback_.resize(detail::soa_layout_traits<SoA>::compute_data_size(sizes...));
      soa_ = detail::soa_layout_traits<SoA>::make_layout(data(), sizes...);
      view_ = view_type{soa_};
    }

    template<typename SoAHost>
    void copy_from_host(const soa_host_buffer<SoAHost>& host) {
      static_assert(std::is_same_v<SoAHost, SoA>);
      if (size() != host.size())
        throw std::invalid_argument("Host and device SOA buffers have different sizes");
      if (storage_)
        storage_->copy_host2device(host.data(), soa_.metadata().byteSize());
      else
        std::memcpy(fallback_.data(), host.data(), fallback_.size());
    }

    template<typename SoAHost>
    void copy_to_host(soa_host_buffer<SoAHost>& host) const {
      static_assert(std::is_same_v<SoAHost, SoA>);
      if (size() != host.size())
        throw std::invalid_argument("Host and device SOA buffers have different sizes");
      if (storage_)
        storage_->copy_device2host(host.data());
      else
        std::memcpy(host.data(), fallback_.data(), fallback_.size());
    }

    [[nodiscard]] std::byte* data() { return storage_ ? storage_->dev_pointer() : fallback_.data(); }
    [[nodiscard]] const std::byte* data() const { return storage_ ? storage_->dev_pointer() : fallback_.data(); }
    [[nodiscard]] size_result_type size() const { return soa_.metadata().size(); }
    [[nodiscard]] SoA& soa() { return soa_; }
    [[nodiscard]] const SoA& soa() const { return soa_; }
    [[nodiscard]] view_type& view() { return view_; }
    [[nodiscard]] const view_type& view() const { return view_; }

  private:
    std::unique_ptr<object<std::byte>> storage_;
    std::vector<std::byte> fallback_;
    SoA soa_;
    view_type view_;
  };

  // Combined host/device storage used by scratchpad and stage entries.
  template<typename SoA, typename queue_type>
  class soa_buffer {
  public:
    using size_type = typename SoA::size_type;
    using view_type = typename SoA::View;

    explicit soa_buffer(std::shared_ptr<queue_type> queue)
        : device_(std::move(queue)) {}

    soa_buffer(const soa_buffer&) = delete;
    soa_buffer& operator=(const soa_buffer&) = delete;
    soa_buffer(soa_buffer&&) = default;
    soa_buffer& operator=(soa_buffer&&) = default;

    template<typename... SizeArgs>
    void alloc(const SizeArgs&... sizes) {
      host_.allocate(sizes...);
      device_.allocate(sizes...);
      host_view_ = typename soa_host_buffer<SoA>::view_type{host_.soa()};
      device_view_ = typename soa_device_buffer<SoA, queue_type>::view_type{device_.soa()};
    }

    void copy_host2device() {
      device_.copy_from_host(host_);
      valid_ = true;
    }

    void copy_device2host() {
      device_.copy_to_host(host_);
      valid_ = false;
    }

    bool is_valid() const { return valid_; }
    void set_valid() { valid_ = true; }
    void set_not_valid() { valid_ = false; }

    [[nodiscard]] auto host_pointer() { return host_.data(); }
    [[nodiscard]] auto host_pointer() const { return host_.data(); }
    [[nodiscard]] auto operator()() { return host_pointer(); }
    [[nodiscard]] auto num_elements() const { return host_.size(); }
    [[nodiscard]] auto dev_pointer() { return device_.data(); }
    [[nodiscard]] auto dev_pointer() const { return device_.data(); }
    [[nodiscard]] auto& host_view() { return host_view_; }
    [[nodiscard]] const auto& host_view() const { return host_view_; }
    [[nodiscard]] auto& device_view() { return device_view_; }
    [[nodiscard]] const auto& device_view() const { return device_view_; }

  private:
    soa_host_buffer<SoA> host_;
    soa_device_buffer<SoA, queue_type> device_;
    typename soa_host_buffer<SoA>::view_type host_view_;
    typename soa_device_buffer<SoA, queue_type>::view_type device_view_;
    bool valid_ = false;
  };

}  // namespace mudock
