#pragma once

#include <cstddef>
#include <mudock/compute/parse_ids.hpp>
#include <mudock/devices.hpp>
#include <mudock/implementation_types.hpp>
#include <mudock/implementations.hpp>
#include <ranges>
#include <stdexcept>
#include <string>
#include <vector>

namespace mudock::stream_engine {

  /**
   * Parsed form of one backend-spec string (one `--use IMPL:DEV:IDS[:WORKERS][:MEMORY_BYTES]` entry).
   *
   * `ids` are physical device ids on GPU; on CPU they only set the degree of parallelism (no core pinning).
   * `workers_per_device` and `memory_bytes` are consumed only by GPU builders; CPU builders ignore them.
   */
  struct backend_spec {
    static constexpr std::size_t default_workers_per_device = 2;
    static constexpr std::size_t default_memory_bytes       = 1'000'000'000;

    implementation_type impl;
    device_type dev;
    std::vector<int> ids;
    std::size_t workers_per_device;
    std::size_t memory_bytes;
  };

  /**
   * Basically copied from legacy `manager::parse_worker_configuration`,
   * so the graph path can eventually be independent of `mudock/compute/manager.hpp`.
   */
  inline std::size_t parse_positive_size_field(const std::vector<std::string>& parts,
                                               const std::size_t index,
                                               const char* field_name,
                                               const std::size_t default_value) {
    if (parts.size() <= index) {
      return default_value;
    }

    std::size_t value = 0;
    try {
      value = std::stoull(parts[index]);
    } catch (const std::exception&) {
      throw std::runtime_error(std::string{"Invalid "} + field_name + " value '" + parts[index] + "'.");
    }

    if (value == 0) {
      throw std::runtime_error(std::string{field_name} + " must be greater than zero.");
    }

    return value;
  }

  /**
   * Tokenize and validate a backend spec string of the form
   * `IMPLEMENTATION:DEVICE:IDS[:WORKERS][:MEMORY_BYTES]`, returning a `backend_spec`.
   *
   * Basically copied from legacy `manager::parse_worker_configuration`,
   * so the graph path can eventually be independent of `mudock/compute/manager.hpp`.
   */
  inline backend_spec parse_backend_spec(const std::string& configuration) {
    std::vector<std::string> parts;
    for (auto v: configuration | std::views::split(':')) parts.emplace_back(v.begin(), v.end());

    if (parts.size() < 3 || parts.size() > 5) {
      throw std::runtime_error("Invalid backend spec '" + configuration +
                               "'. Expected IMPLEMENTATION:DEVICE:IDS[:WORKERS][:MEMORY_BYTES].");
    }
    if (parts[0].empty() || parts[1].empty() || parts[2].empty()) {
      throw std::runtime_error("Invalid backend spec '" + configuration +
                               "'. Implementation, device, and ids must be non-empty.");
    }

    return backend_spec{
        .impl               = get_impl_type(parts[0]),
        .dev                = get_device_type(parts[1]),
        .ids                = parse_ids(parts[2]),
        .workers_per_device = parse_positive_size_field(parts,
                                                        3,
                                                        "workers_per_device",
                                                        backend_spec::default_workers_per_device),
        .memory_bytes =
            parse_positive_size_field(parts, 4, "memory_bytes", backend_spec::default_memory_bytes),
    };
  }

} // namespace mudock::stream_engine
