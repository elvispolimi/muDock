#pragma once

#include "backend_spec.hpp"
#include "batcher_body.hpp"
#include "screener_body.hpp"
#include "stream_engine_types.hpp"

#include <concepts>
#include <functional>
#include <iterator>
#include <memory>
#include <mudock/compute/queue.hpp>
#include <mudock/compute/scratchpad.hpp>
#include <mudock/devices.hpp>
#include <mudock/implementations.hpp>
#include <mudock/knobs.hpp>
#include <mudock/log.hpp>
#include <mudock/utils.hpp>
#include <oneapi/tbb/flow_graph.h>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mudock::stream_engine {

  /**
   * The abstraction of a CPU-group (one per backend-spec) or GPU-device (one per id in the backend-spec)
   * backend subgraph: each owns a batcher node feeding a batch buffer that fans out to N screener nodes;
   * the subgraph also owns the batcher state and all screener states backing those nodes.
   *
   * CPU: a number of hardware threads specified in the backend-spec shares host memory,
   * resulting in a shared ROB-holding batcher state,
   * and a screener node per thread sharing a single device scratchpad.
   *
   * GPU: each GPU has isolated memory, resulting in a dedicated ROB-holding batcher state,
   * and a screener node per worker sharing a single device scratchpad.
   *
   * Note: screener nodes have `concurrency = 1` because each screener_state owns its own stage instance
   * that is not thread-safe;
   * for now, parallelism inside a backend-spec is achieved by spawning multiple screener nodes,
   * rather than by making a single node concurrent.
   */
  struct backend_subgraph {
    std::unique_ptr<batcher_state> batching;
    std::vector<std::unique_ptr<screener_state_base>> screening;

    std::unique_ptr<batcher_mfn> batcher;
    std::unique_ptr<batch_buf> batch_buffer;
    std::vector<std::unique_ptr<screener_mfn>> screeners;

    /**
     * Connect this subgraph into the enclosing graph.
     */
    void wire(parsed_ligands_buf& upstream, writer_fn& downstream) {
      tbb_flow::make_edge(upstream, *batcher);
      tbb_flow::make_edge(tbb_flow::output_port<0>(*batcher), *batch_buffer);
      for (auto& screener: screeners) {
        tbb_flow::make_edge(*batch_buffer, *screener);
        tbb_flow::make_edge(tbb_flow::output_port<0>(*screener), downstream);
      }
    }

    /**
     * Flush the ROB's pending batches into the subgraph's batch_buffer.
     */
    void flush_batches() {
      batching->flush_all([this](const lig_batch_ptr& batch) { batch_buffer->try_put(batch); });
    }
  };

  /**
   * Build the subgraph for a CPU backend-spec (see more on `backend_subgraph` above).
   *
   * Mirrors `manager::launch_worker_cpu`, so the graph path can eventually be independent of
   * `mudock/compute/manager.hpp`.
   */
  template<typename queue_type, typename pipeline_t>
    requires std::derived_from<queue_type, queue>
  backend_subgraph build_cpu_subgraph(tbb_flow::graph& g,
                                      const knobs& knobs,
                                      pipeline_t& pipeline,
                                      const backend_spec& spec) {
    auto device_scratch = std::make_shared<scratchpad<queue_type>>(knobs, 0, device_type::CPU);
    auto q_b            = device_scratch->get_queue();
    std::function<int(const int)> get_size = [q_b, &knobs](const int atoms) {
      return pipeline_t::template get_batch_size<queue_type>(atoms, q_b, knobs);
    };
    using stage_t =
        decltype(pipeline.template get_pipeline<queue_type>(knobs, 0, device_type::CPU, device_scratch));

    backend_subgraph cpu_subgraph;
    cpu_subgraph.batching = std::make_unique<batcher_state>(std::move(get_size));
    cpu_subgraph.batcher =
        std::make_unique<batcher_mfn>(g, tbb_flow::serial, batcher_body{*cpu_subgraph.batching});
    cpu_subgraph.batch_buffer = std::make_unique<batch_buf>(g);

    cpu_subgraph.screening.reserve(spec.ids.size());
    cpu_subgraph.screeners.reserve(spec.ids.size());

    for (const auto id: spec.ids) {
      info(std::format("Setting up CPU device {} ...", id));
      auto cpu_screener_stage =
          pipeline.template get_pipeline<queue_type>(knobs, id, device_type::CPU, device_scratch);
      auto& cpu_screener_state = static_cast<screener_state<stage_t>&>(*cpu_subgraph.screening.emplace_back(
          std::make_unique<screener_state<stage_t>>(std::move(cpu_screener_stage))));
      cpu_subgraph.screeners.push_back(
          std::make_unique<screener_mfn>(g, /*concurrency=*/1, screener_body<stage_t>{cpu_screener_state}));
    }

    return cpu_subgraph;
  }

  /**
   * Build the subgraphs for a GPU backend-spec (see more on `backend_subgraph` above).
   *
   * Mirrors `manager::launch_worker_gpu`, so the graph path can eventually be independent of
   * `mudock/compute/manager.hpp`.
   */
  template<typename queue_type, typename pipeline_t>
    requires std::derived_from<queue_type, queue>
  std::vector<backend_subgraph> build_gpu_subgraph(tbb_flow::graph& g,
                                                   const knobs& knobs,
                                                   pipeline_t& pipeline,
                                                   const backend_spec& spec) {
    std::vector<backend_subgraph> gpu_subgraphs;
    gpu_subgraphs.reserve(spec.ids.size());

    for (const auto id: spec.ids) {
      auto device_scratch = std::make_shared<scratchpad<queue_type>>(knobs, id, device_type::GPU);
      auto q_b            = device_scratch->get_queue();
      std::function<int(const int)> get_size =
          [q_b, &knobs, memory_bytes = spec.memory_bytes](const int atoms) {
            return pipeline_t::template get_batch_size<queue_type>(atoms, q_b, knobs, memory_bytes);
          };
      using stage_t =
          decltype(pipeline.template get_pipeline<queue_type>(knobs, id, device_type::GPU, device_scratch));

      info(std::format("Setting up GPU device {} with {} workers, {} bytes each ...",
                       id,
                       spec.workers_per_device,
                       spec.memory_bytes));

      backend_subgraph gpu_subgraph;
      gpu_subgraph.batching = std::make_unique<batcher_state>(std::move(get_size));
      gpu_subgraph.batcher =
          std::make_unique<batcher_mfn>(g, tbb_flow::serial, batcher_body{*gpu_subgraph.batching});
      gpu_subgraph.batch_buffer = std::make_unique<batch_buf>(g);

      gpu_subgraph.screening.reserve(spec.workers_per_device);
      gpu_subgraph.screeners.reserve(spec.workers_per_device);

      for (std::size_t i = 0; i < spec.workers_per_device; ++i) {
        auto gpu_screener_stage =
            pipeline.template get_pipeline<queue_type>(knobs, id, device_type::GPU, device_scratch);
        auto& gpu_screener_state = static_cast<screener_state<stage_t>&>(*gpu_subgraph.screening.emplace_back(
            std::make_unique<screener_state<stage_t>>(std::move(gpu_screener_stage))));
        gpu_subgraph.screeners.push_back(
            std::make_unique<screener_mfn>(g, /*concurrency=*/1, screener_body<stage_t>{gpu_screener_state}));
      }

      gpu_subgraphs.push_back(std::move(gpu_subgraph));
    }

    return gpu_subgraphs;
  }

  /**
   * Dispatch one CPU backend-spec to the queue_type-specialized constructor matching `spec.impl`
   * (CPP / GH / XSIMD).
   */
  template<typename pipeline_t>
  backend_subgraph dispatch_cpu_subgraph(tbb_flow::graph& g,
                                         const knobs& knobs,
                                         pipeline_t& pipeline,
                                         const backend_spec& spec) {
    std::optional<backend_subgraph> cpu_subgraph;

    constexpr_for<0, num_cpu_kernel_type(), 1>([&](const auto kernel) {
      constexpr auto kernel_type = cpu_kernel_type[kernel];
      if (kernel_type == spec.impl) {
        using queue_type_t = typename kernel_type_traits<kernel_type>::type;
        cpu_subgraph.emplace(build_cpu_subgraph<queue_type_t>(g, knobs, pipeline, spec));
      }
    });

    if (!cpu_subgraph) {
      throw std::runtime_error("Requested CPU implementation not compiled in.");
    }

    return std::move(*cpu_subgraph);
  }

  /**
   * Dispatch one GPU backend-spec to the queue_type-specialized constructor matching `spec.impl`
   * (CUDA / HIP / SYCL).
   */
  template<typename pipeline_t>
  std::vector<backend_subgraph> dispatch_gpu_subgraph(tbb_flow::graph& g,
                                                      const knobs& knobs,
                                                      pipeline_t& pipeline,
                                                      const backend_spec& spec) {
    std::vector<backend_subgraph> gpu_subgraphs;

    constexpr_for<0, num_gpu_kernel_type(), 1>([&](const auto kernel) {
      constexpr auto kernel_type = gpu_kernel_type[kernel];
      if (kernel_type == spec.impl) {
        using queue_type_t = typename kernel_type_traits<kernel_type>::type;
        gpu_subgraphs      = build_gpu_subgraph<queue_type_t>(g, knobs, pipeline, spec);
      }
    });

    if (gpu_subgraphs.empty()) {
      throw std::runtime_error("Requested GPU implementation not compiled in.");
    }

    return gpu_subgraphs;
  }

  /**
   * Build the `backend_subgraph`s for every backend-spec string, dispatching each by device type.
   *
   * Mirrors the `manager::manager` pattern, so the graph path can eventually be independent of
   * `mudock/compute/manager.hpp`.
   */
  template<typename pipeline_t>
  std::vector<backend_subgraph> build_backend_subgraphs(tbb_flow::graph& g,
                                                        const knobs& knobs,
                                                        pipeline_t& pipeline,
                                                        const std::vector<std::string>& configurations) {
    std::vector<backend_subgraph> backend_subgraphs;

    for (const auto& configuration: configurations) {
      const auto spec = parse_backend_spec(configuration);

      info(std::format("Setting up backend-spec '{}' ...", configuration));

      switch (spec.dev) {
        case device_type::CPU: {
          backend_subgraphs.push_back(dispatch_cpu_subgraph(g, knobs, pipeline, spec));
          break;
        }
        case device_type::GPU: {
          auto gpu_subgraphs = dispatch_gpu_subgraph(g, knobs, pipeline, spec);
          backend_subgraphs.insert(backend_subgraphs.end(),
                                   std::make_move_iterator(gpu_subgraphs.begin()),
                                   std::make_move_iterator(gpu_subgraphs.end()));
          break;
        }
        default: throw std::runtime_error("Not supported device type.");
      }
    }

    return backend_subgraphs;
  }

} // namespace mudock::stream_engine
