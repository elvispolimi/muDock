#pragma once

#include <concepts>
#include <memory>
#include <mutex>
#include <mudock/compute/queue.hpp>
#include <mudock/data_formats/soa_buffer.hpp>
#include <mudock/data_formats/data_layout_sets.hpp>
#include <mudock/devices.hpp>
#include <mudock/knobs.hpp>

namespace mudock {

  template<typename Layout, typename queue_type>
  struct soa_resource {
    void initialize(const std::shared_ptr<queue_type>& queue) {
      buffer = std::make_shared<soa_buffer<Layout, queue_type>>(queue);
    }

    std::shared_ptr<soa_buffer<Layout, queue_type>> buffer;
    std::mutex mutex;
  };

  template<typename queue_type, typename SoABlock = all_data_layouts>
    requires std::derived_from<queue_type, queue>
  class data_store_soa;

  template<typename queue_type, typename SoABlock>
    requires std::derived_from<queue_type, queue>
  class data_store_soa {
  public:
    data_store_soa(const knobs& conf,
                   const int id,
                   const device_type dev_type,
                   std::shared_ptr<device_memory_tracker> tracker = {})
        : configuration(conf),
          q(std::make_shared<queue_type>(id, dev_type, std::move(tracker))),
          blocks_() {
      blocks_.initialize(q);
    }

    const knobs configuration;

    auto get_queue() { return q; }
    auto get_memory_tracker() const { return q->get_memory_tracker(); }

    template<typename... SizeArgs>
    void alloc(const SizeArgs&... sizes) { blocks_.buffer->alloc(sizes...); }

    bool is_valid() const { return blocks_.buffer->is_valid(); }
    void set_valid() { blocks_.buffer->set_valid(); }
    void set_not_valid() { blocks_.buffer->set_not_valid(); }
    void copy_host2device() { blocks_.buffer->copy_host2device(); }
    void copy_device2host() { blocks_.buffer->copy_device2host(); }
    void copy_device2host_scores() {
      blocks_.buffer->copy_device2host_block([](auto& data) -> decltype(auto) { return data.scores(); });
    }
    void copy_device2host_best() {
      blocks_.buffer->copy_device2host_block([](auto& data) -> decltype(auto) { return data.best(); });
    }

    auto& host_view() { return blocks_.buffer->host_view(); }
    auto& device_view() { return blocks_.buffer->device_view(); }
    const auto& host_view() const { return blocks_.buffer->host_view(); }
    const auto& device_view() const { return blocks_.buffer->device_view(); }

    void invalidate() {
      blocks_.buffer->set_not_valid();
    }

  private:
    std::shared_ptr<queue_type> q;
    soa_resource<SoABlock, queue_type> blocks_;
  };

}  // namespace mudock
