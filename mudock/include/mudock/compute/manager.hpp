#pragma once

#include <array>
#include <concepts>
#include <memory>
#include <mudock/compute/buffer.hpp>
#include <mudock/compute/parse_ids.hpp>
#include <mudock/compute/pipeline.hpp>
#include <mudock/compute/queue.hpp>
#include <mudock/compute/safe_stack.hpp>
#include <mudock/compute/scratchpad.hpp>
#include <mudock/compute/stage.hpp>
#include <mudock/compute/threadpool.hpp>
#include <mudock/compute/worker.hpp>
#include <mudock/devices.hpp>
#include <mudock/grid.hpp>
#include <mudock/implementations.hpp>
#include <mudock/knobs.hpp>
#include <mudock/molecule.hpp>
#include <mudock/utils.hpp>
#include <ranges>
#include <stdexcept>

namespace mudock {

  template<typename queue_type, typename pipeline_t>
    requires std::derived_from<queue_type, queue>
  inline void launch_worker_cpu(const knobs& knobs,
                                const std::vector<std::string>& parts,
                                threadpool& pool,
                                std::shared_ptr<safe_stack<static_molecule>>& input_molecules,
                                std::shared_ptr<safe_stack<static_molecule>>& output_molecules,
                                pipeline_t& pipe) {
    auto device_scratch = std::make_shared<scratchpad<queue_type>>(knobs, 0, device_type::CPU);
    auto q_b            = device_scratch->get_queue();
    std::function<int(const int)> get_size = [q_b](const int x) {
      return pipeline_t::template get_batch_size<queue_type>(x, q_b);
    };
    auto rob = std::make_shared<reorder_buffer<static_molecule>>(get_size);

    for (const auto id: parse_ids(parts[2])) {
      pool.add_worker(
          worker(input_molecules,
                 output_molecules,
                 rob,
                 pipe.template get_pipeline<queue_type>(knobs, id, device_type::CPU, device_scratch)));
    }
  };

  template<typename queue_type, typename pipeline_t>
    requires std::derived_from<queue_type, queue>
  inline void launch_worker_gpu(const knobs& knobs,
                                const std::vector<std::string>& parts,
                                threadpool& pool,
                                std::shared_ptr<safe_stack<static_molecule>>& input_molecules,
                                std::shared_ptr<safe_stack<static_molecule>>& output_molecules,
                                pipeline_t& pipe) {
    auto q_b                               = std::make_shared<queue_type>(0, device_type::GPU);
    std::function<int(const int)> get_size = [q_b](const int x) {
      return pipeline_t::template get_batch_size<queue_type>(x, q_b);
    };
    auto rob = std::make_shared<reorder_buffer<static_molecule>>(get_size);

    for (const auto id: parse_ids(parts[2])) {
      auto device_scratch = std::make_shared<scratchpad<queue_type>>(knobs, id, device_type::GPU);
      pool.add_worker(
          worker(input_molecules,
                 output_molecules,
                 rob,
                 pipe.template get_pipeline<queue_type>(knobs, id, device_type::GPU, device_scratch)));
    }
  };

  // this function will configure and create (if needed) cpp workers to the threadpool
  template<typename pipeline_t>
  void manager(const std::vector<std::string>& configurations,
               threadpool& pool,
               const knobs knobs,
               std::shared_ptr<safe_stack<static_molecule>>& input_molecules,
               std::shared_ptr<safe_stack<static_molecule>>& output_molecules,
               pipeline_t& pipe) {
    for (auto& configuration: configurations) {
      std::vector<std::string> parts;

      for (auto v: configuration | std::views::split(':')) parts.emplace_back(v.begin(), v.end());

      // parts[0], parts[1], parts[2]
      auto dev_t  = get_device_type(parts[1]);
      auto impl_t = get_impl_type(parts[0]);
      switch (dev_t) {
        case device_type::CPU: {
          constexpr_for<0, num_cpu_kernel_type(), 1>([&](const auto kernel) {
            constexpr auto kernel_type = cpu_kernel_type[kernel];
            if (kernel_type == impl_t) {
              using k_t = typename kernel_type_traits<kernel_type>::type;
              launch_worker_cpu<k_t, pipeline_t>(knobs, parts, pool, input_molecules, output_molecules, pipe);
            }
          });
          break;
        }
        case device_type::GPU: {
          constexpr_for<0, num_gpu_kernel_type(), 1>([&](const auto kernel) {
            constexpr auto kernel_type = gpu_kernel_type[kernel];
            if (kernel_type == impl_t) {
              using k_t = typename kernel_type_traits<kernel_type>::type;
              launch_worker_gpu<k_t, pipeline_t>(knobs, parts, pool, input_molecules, output_molecules, pipe);
            }
          });
          break;
        }
        default: throw std::runtime_error("Not supported device type"); break;
      }
    }
  }
} // namespace mudock
