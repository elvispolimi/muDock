#pragma once

#include <concepts>
#include <mudock/compute/buffer.hpp>
#include <mudock/devices.hpp>
#include <mudock/knobs.hpp>
#include <mudock/log.hpp>
#include <mutex>
#include <unordered_map>

namespace mudock {
  template<template<typename, typename> typename buffer_type, typename T, typename queue_t>
  struct buffer_entry {
    std::unordered_map<buffer_data_type, std::shared_ptr<buffer_type<T, queue_t>>> map;
    std::mutex mtx;
  };

  template<template<typename, typename> typename buffer_type,
           typename queue_t,
           // template<typename> typename object_type,
           typename TypeList>
  struct scratchpad_impl;

  // Specialization when TypeList is a std::tuple<...>
  // template<template<typename> typename buffer_type, template<typename> typename object_type, typename... Ts>
  // struct scratchpad_impl<buffer_type, object_type, std::tuple<Ts...>> {
  template<template<typename, typename> typename buffer_type, typename queue_t, typename... Ts>
  struct scratchpad_impl<buffer_type, queue_t, std::tuple<Ts...>> {
    scratchpad_impl(const knobs& conf, const int id, const device_type dev_type)
        : configuration(conf), q(std::make_shared<queue_t>(id, dev_type)) {};

    template<buffer_data_type bdt>
    buffer_type<typename buffer_type_traits<bdt>::type, queue_t>& get(const int dim = 0) {
      using bdt_type = typename buffer_type_traits<bdt>::type;
      auto& entry    = std::get<buffer_entry<buffer_type, bdt_type, queue_t>>(buffers);
      auto dlck      = std::lock_guard<std::mutex>{entry.mtx};
      if (!exists_int<bdt>()) {
        add<bdt>(dim);
      } else if (dim) {
        mudock::info("Requested a scratchpad buffer with size different than zero, but it already exists");
      }
      return *entry.map.at(bdt);
    }

    template<buffer_data_type bdt>
    bool exists(const int dim = 0) {
      using bdt_type = typename buffer_type_traits<bdt>::type;
      auto& entry    = std::get<buffer_entry<buffer_type, bdt_type, queue_t>>(buffers);
      auto dlck      = std::lock_guard<std::mutex>{entry.mtx};
      bool exists    = exists_int<bdt>();
      if (!exists) {
        add<bdt>(dim);
      }
      return exists;
    }

    const knobs configuration;

    auto get_queue() { return q; }

    void invalidate() { invalidate_impl(std::index_sequence_for<Ts...>{}); }

  private:
    std::tuple<buffer_entry<buffer_type, Ts, queue_t>...> buffers;
    std::shared_ptr<queue_t> q;

    template<buffer_data_type bdt>
    void add(const int dim) {
      using bdt_type = typename buffer_type_traits<bdt>::type;
      auto& entry    = std::get<buffer_entry<buffer_type, bdt_type, queue_t>>(buffers);
      entry.map.emplace(bdt, std::make_shared<buffer_type<bdt_type, queue_t>>(q, dim));
    }

    template<buffer_data_type bdt>
    bool exists_int() {
      using bdt_type = typename buffer_type_traits<bdt>::type;
      auto& entry    = std::get<buffer_entry<buffer_type, bdt_type, queue_t>>(buffers);
      auto search    = entry.map.find(bdt);
      return search != entry.map.end();
    }

    template<std::size_t... Is>
    void invalidate_impl(std::index_sequence<Is...>) {
      (invalidate_entry(std::get<Is>(buffers)), ...);
    }

    template<typename Entry>
    static void invalidate_entry(Entry& entry) {
      std::lock_guard<std::mutex> lock(entry.mtx);
      for (auto& [bdt, buf_ptr]: entry.map) {
        if (buf_ptr) {
          buf_ptr->set_not_valid();
        }
      }
    }
  };
  // TODO fix me for all buffers
  template<typename queue_t>
    requires std::derived_from<queue_t, queue>
  using scratchpad = scratchpad_impl<buffer_vector, queue_t, buffer_type_list>;

  // using scratchpad_cpp = scratchpad<queue_cpp>;
} // namespace mudock
