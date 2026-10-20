#pragma once

#include <concepts>
#include <memory>
#include <mutex>
#include <tuple>
#include <mudock/compute/queue.hpp>
#include <mudock/compute/soa_buffer.hpp>
#include <mudock/compute/scratchpad_layout_sets.hpp>
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

  template<typename queue_type, typename LayoutList = all_scratch_layouts>
    requires std::derived_from<queue_type, queue>
  class scratchpad_soa;

  template<typename queue_type, typename... Layouts>
    requires std::derived_from<queue_type, queue>
  class scratchpad_soa<queue_type, type_list<Layouts...>> {
  public:
    scratchpad_soa(const knobs& conf,
                   const int id,
                   const device_type dev_type,
                   std::shared_ptr<device_memory_tracker> tracker = {})
        : configuration(conf),
          q(std::make_shared<queue_type>(id, dev_type, std::move(tracker))),
          resources_() {
      initialize_resources(std::index_sequence_for<Layouts...>{});
    }

    const knobs configuration;

    auto get_queue() { return q; }
    auto get_memory_tracker() const { return q->get_memory_tracker(); }

    template<typename Layout>
      requires layout_in_list<Layout, type_list<Layouts...>>::value
    auto& get() {
      auto& resource = std::get<soa_resource<Layout, queue_type>>(resources_);
      std::lock_guard lock(resource.mutex);
      return *resource.buffer;
    }

    auto& coordinates() { return get<ScratchInputCoordinates<>>(); }
    auto& scratch_coordinates() { return get<ScratchCoordinates<>>(); }
    auto& batch_info() { return get<ScratchBatchInfo<>>(); }
    auto& genetic() { return get<ScratchGenetic<>>(); }
    auto& scores() { return get<ScratchScores<>>(); }
    auto& best() { return get<ScratchBest<>>(); }
    auto& protein_bounds() { return get<ScratchProteinBounds<>>(); }
    auto& protein_dimensions() { return get<ScratchProteinDimensions<>>(); }
    auto& protein_maps() { return get<ScratchProteinMaps<>>(); }

    void invalidate() {
      (invalidate_resource<Layouts>(), ...);
    }

  private:
    template<typename Layout>
    void invalidate_resource() {
      get<Layout>().set_not_valid();
    }

    template<std::size_t... Indices>
    void initialize_resources(std::index_sequence<Indices...>) {
      (std::get<Indices>(resources_).initialize(q), ...);
    }

    std::shared_ptr<queue_type> q;
    std::tuple<soa_resource<Layouts, queue_type>...> resources_;
  };

}  // namespace mudock
