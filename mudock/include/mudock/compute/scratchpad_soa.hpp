#pragma once

#include <concepts>
#include <memory>
#include <mudock/compute/queue.hpp>
#include <mudock/compute/soa_buffer.hpp>
#include <mudock/compute/scratchpad_data.hpp>
#include <mudock/devices.hpp>
#include <mudock/knobs.hpp>

namespace mudock {

  template<typename queue_type>
    requires std::derived_from<queue_type, queue>
  class scratchpad_soa {
  public:
    using coordinates_type = soa_buffer<ScratchCoordinates<>, queue_type>;
    using batch_info_type = soa_buffer<ScratchBatchInfo<>, queue_type>;
    using genetic_type = soa_buffer<ScratchGenetic<>, queue_type>;
    using scores_type = soa_buffer<ScratchScores<>, queue_type>;
    using best_type = soa_buffer<ScratchBest<>, queue_type>;
    using protein_bounds_type = soa_buffer<ScratchProteinBounds<>, queue_type>;
    using protein_dimensions_type = soa_buffer<ScratchProteinDimensions<>, queue_type>;
    using protein_maps_type = soa_buffer<ScratchProteinMaps<>, queue_type>;

    scratchpad_soa(const knobs& conf,
                   const int id,
                   const device_type dev_type,
                   std::shared_ptr<device_memory_tracker> tracker = {})
        : configuration(conf),
          q(make_queue(id, dev_type, std::move(tracker))),
          coordinates_(q),
          scratch_coordinates_(q),
          batch_info_(q),
          genetic_(q),
          scores_(q),
          best_(q),
          protein_bounds_(q),
          protein_dimensions_(q),
          protein_maps_(q) {}

    const knobs configuration;

    auto get_queue() { return q; }
    auto get_memory_tracker() const { return q->get_memory_tracker(); }

    coordinates_type& coordinates() { return coordinates_; }
    coordinates_type& scratch_coordinates() { return scratch_coordinates_; }
    batch_info_type& batch_info() { return batch_info_; }
    genetic_type& genetic() { return genetic_; }
    scores_type& scores() { return scores_; }
    best_type& best() { return best_; }
    protein_bounds_type& protein_bounds() { return protein_bounds_; }
    protein_dimensions_type& protein_dimensions() { return protein_dimensions_; }
    protein_maps_type& protein_maps() { return protein_maps_; }

    void invalidate() {
      coordinates_.set_not_valid();
      scratch_coordinates_.set_not_valid();
      batch_info_.set_not_valid();
      genetic_.set_not_valid();
      scores_.set_not_valid();
      best_.set_not_valid();
      protein_bounds_.set_not_valid();
      protein_dimensions_.set_not_valid();
      protein_maps_.set_not_valid();
    }

  private:
    static std::shared_ptr<queue_type> make_queue(const int id,
                                                   const device_type dev_type,
                                                   std::shared_ptr<device_memory_tracker> tracker) {
      return std::make_shared<queue_type>(id, dev_type, std::move(tracker));
    }

    std::shared_ptr<queue_type> q;
    coordinates_type coordinates_;
    coordinates_type scratch_coordinates_;
    batch_info_type batch_info_;
    genetic_type genetic_;
    scores_type scores_;
    best_type best_;
    protein_bounds_type protein_bounds_;
    protein_dimensions_type protein_dimensions_;
    protein_maps_type protein_maps_;
  };

}  // namespace mudock
