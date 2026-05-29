#pragma once

#include <concepts>
#include <memory>
#include <mudock/compute/buffer.hpp>
#include <mudock/compute/parse_ids.hpp>
#include <mudock/compute/pipeline.hpp>
#include <mudock/compute/queue.hpp>
#include <mudock/compute/safe_queue.hpp>
#include <mudock/compute/scratchpad.hpp>
#include <mudock/compute/stage.hpp>
#include <mudock/compute/threadpool.hpp>
#include <mudock/compute/worker.hpp>
#include <mudock/devices.hpp>
#include <mudock/grid.hpp>
#include <mudock/implementations.hpp>
#include <mudock/knobs.hpp>
#include <mudock/log.hpp>
#include <mudock/molecule.hpp>
#include <mudock/utils.hpp>
#include <ranges>
#include <stdexcept>
#include <string>
#include <vector>

namespace mudock {
  inline std::vector<std::string> parse_worker_configuration(const std::string& configuration) {
    std::vector<std::string> parts;
    for (auto v: configuration | std::views::split(':')) parts.emplace_back(v.begin(), v.end());

    if (parts.size() < 3 || parts.size() > 5) {
      throw std::runtime_error("Invalid device configuration '" + configuration +
                               "'. Expected IMPLEMENTATION:DEVICE:IDS[:WORKERS][:MEMORY_BYTES].");
    }
    if (parts[0].empty() || parts[1].empty() || parts[2].empty()) {
      throw std::runtime_error("Invalid device configuration '" + configuration +
                               "'. Implementation, device, and ids must be non-empty.");
    }
    return parts;
  }

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

  template<typename queue_type, typename pipeline_t>
    requires std::derived_from<queue_type, queue>
  inline void launch_worker_cpu(const knobs& knobs,
                                const std::vector<std::string>& parts,
                                threadpool& pool,
                                std::shared_ptr<safe_queue<static_molecule>>& input_molecules,
                                std::shared_ptr<safe_queue<static_molecule>>& output_molecules,
                                pipeline_t& pipe,
                                std::atomic<std::size_t>* in_flight_ligands = nullptr) {
    auto device_scratch = std::make_shared<scratchpad<queue_type>>(knobs, 0, device_type::CPU);
    auto q_b            = device_scratch->get_queue();
    std::function<int(const int)> get_size = [q_b, &knobs](const int x) {
      return pipeline_t::template get_batch_size<queue_type>(x, q_b, knobs);
    };
    auto rob = std::make_shared<reorder_buffer<static_molecule>>(get_size);

    for (const auto id: parse_ids(parts[2])) {
      mudock::info("Starting CPU device ", id, ".");
      pool.add_worker(
          worker(input_molecules,
                 output_molecules,
                 rob,
                 pipe.template get_pipeline<queue_type>(knobs, id, device_type::CPU, device_scratch),
                 in_flight_ligands));
    }
  };

  template<typename queue_type, typename pipeline_t>
    requires std::derived_from<queue_type, queue>
  inline void launch_worker_gpu(const knobs& knobs,
                                const std::vector<std::string>& parts,
                                threadpool& pool,
                                std::shared_ptr<safe_queue<static_molecule>>& input_molecules,
                                std::shared_ptr<safe_queue<static_molecule>>& output_molecules,
                                pipeline_t& pipe,
                                std::atomic<std::size_t>* in_flight_ligands = nullptr) {
    const std::size_t workers_per_device =
        parse_positive_size_field(parts, 3, "workers_per_device", static_cast<std::size_t>(2));
    const std::size_t mem_per_device =
        parse_positive_size_field(parts, 4, "memory_bytes", static_cast<std::size_t>(1000000000));

    for (const auto id: parse_ids(parts[2])) {
      auto q_b                               = std::make_shared<queue_type>(id, device_type::GPU);
      std::function<int(const int)> get_size = [q_b, &knobs, mem_per_device](const int x) {
        return pipeline_t::template get_batch_size<queue_type>(x, q_b, knobs, mem_per_device);
      };
      auto rob = std::make_shared<reorder_buffer<static_molecule>>(get_size);
      mudock::info("Starting GPU device ",
                   id,
                   " with ",
                   workers_per_device,
                   " per device, with ",
                   mem_per_device,
                   " bytes each.");
      auto device_scratch = std::make_shared<scratchpad<queue_type>>(knobs, id, device_type::GPU);
      for (std::size_t i = 0; i < workers_per_device; ++i) {
        pool.add_worker(
            worker(input_molecules,
                   output_molecules,
                   rob,
                   pipe.template get_pipeline<queue_type>(knobs, id, device_type::GPU, device_scratch),
                   in_flight_ligands));
      }
    }
  };

  // this function will configure and create (if needed) cpp workers to the threadpool
  template<typename pipeline_t>
  void manager(const std::vector<std::string>& configurations,
               threadpool& pool,
               const knobs knobs,
               std::shared_ptr<safe_queue<static_molecule>>& input_molecules,
               std::shared_ptr<safe_queue<static_molecule>>& output_molecules,
               pipeline_t& pipe,
               std::atomic<std::size_t>* in_flight_ligands = nullptr) {
    for (auto& configuration: configurations) {
      const auto parts = parse_worker_configuration(configuration);
      auto dev_t       = get_device_type(parts[1]);
      auto impl_t      = get_impl_type(parts[0]);
      switch (dev_t) {
        case device_type::CPU: {
          constexpr_for<0, num_cpu_kernel_type(), 1>([&](const auto kernel) {
            constexpr auto kernel_type = cpu_kernel_type[kernel];
            if (kernel_type == impl_t) {
              using k_t = typename kernel_type_traits<kernel_type>::type;
              launch_worker_cpu<k_t, pipeline_t>(knobs,
                                                 parts,
                                                 pool,
                                                 input_molecules,
                                                 output_molecules,
                                                 pipe,
                                                 in_flight_ligands);
            }
          });
          break;
        }
        case device_type::GPU: {
          constexpr_for<0, num_gpu_kernel_type(), 1>([&](const auto kernel) {
            constexpr auto kernel_type = gpu_kernel_type[kernel];
            if (kernel_type == impl_t) {
              using k_t = typename kernel_type_traits<kernel_type>::type;
              launch_worker_gpu<k_t, pipeline_t>(knobs,
                                                 parts,
                                                 pool,
                                                 input_molecules,
                                                 output_molecules,
                                                 pipe,
                                                 in_flight_ligands);
            }
          });
          break;
        }
        default: throw std::runtime_error("Not supported device type"); break;
      }
    }
  }
} // namespace mudock
