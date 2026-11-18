#pragma once

#include "backend_spec.hpp"
#include "batcher_body.hpp"
#include "screener_body.hpp"
#include "stream_engine_types.hpp"

#include <algorithm>
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
  class backend_subgraph {
    std::unique_ptr<batcher_state> batching;
    std::vector<std::unique_ptr<screener_state_base>> screening;

    std::unique_ptr<batcher_mfn> batcher;
    std::unique_ptr<batch_buf> batch_buffer;
    std::vector<std::unique_ptr<screener_mfn>> screeners;

    observability::backend_metrics* metrics;

  public:
    explicit backend_subgraph(observability::backend_metrics& _metrics): metrics(&_metrics) {}

    /// Install the batcher stage: ROB-holding state, serial batcher node, and downstream batch buffer.
    void set_batcher(tbb_flow::graph& g, std::function<int(const int)> get_size) {
      batching     = std::make_unique<batcher_state>(std::move(get_size));
      batcher      = std::make_unique<batcher_mfn>(g, tbb_flow::serial, batcher_body{*batching, *metrics});
      batch_buffer = std::make_unique<batch_buf>(g);
    }

    /// Reserve capacity for `n` screeners across both the state and node vectors.
    void reserve_screeners(std::size_t n) {
      screening.reserve(n);
      screeners.reserve(n);
    }

    /// Append one screener.
    template<typename stage_t>
    void add_screener(tbb_flow::graph& g, stage_t stage) {
      auto state_ptr = std::make_unique<screener_state<stage_t>>(std::move(stage));
      auto& state    = *state_ptr;
      screening.push_back(std::move(state_ptr));
      screeners.push_back(
          std::make_unique<screener_mfn>(g, tbb_flow::serial, screener_body<stage_t>{state, *metrics}));
    }

    /// Connect this subgraph into the enclosing graph.
    void wire(parsed_ligands_buf& upstream, writer_fn& downstream) {
      tbb_flow::make_edge(upstream, *batcher);
      tbb_flow::make_edge(tbb_flow::output_port<0>(*batcher), *batch_buffer);
      for (auto& screener: screeners) {
        tbb_flow::make_edge(*batch_buffer, *screener);
        tbb_flow::make_edge(tbb_flow::output_port<0>(*screener), downstream);
      }
    }

    /// Flush the ROB's pending batches into the subgraph's batch_buffer.
    void flush_batches() {
      std::size_t count = 0;
      batching->flush_all([this, &count](const lig_batch_ptr& batch) {
        ++count;
        batch_buffer->try_put(batch);
      });
      metrics->batches_emitted.increase(count);
    }
  };

  namespace setup {

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
                                        const backend_spec& spec,
                                        observability::run_metrics& metrics) {
      auto& backend_metrics = *metrics.backends.emplace_back(
          std::make_unique<observability::backend_metrics>(spec.configuration, spec.ids.size()));

      auto device_scratch = std::make_shared<scratchpad<queue_type>>(knobs, 0, device_type::CPU);
      auto q_b            = device_scratch->get_queue();
      std::function<int(const int)> get_size = [q_b, &knobs](const int atoms) {
        return pipeline_t::template get_batch_size<queue_type>(atoms, q_b, knobs);
      };

      backend_subgraph cpu_subgraph(backend_metrics);
      cpu_subgraph.set_batcher(g, std::move(get_size));
      cpu_subgraph.reserve_screeners(spec.ids.size());

      for (const auto id: spec.ids) {
        info(std::format("Setting up CPU device {} ...", id));
        auto stage = pipeline.template get_pipeline<queue_type>(knobs, id, device_type::CPU, device_scratch);
        cpu_subgraph.add_screener(g, std::move(stage));
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
    std::vector<backend_subgraph> build_gpu_subgraphs(tbb_flow::graph& g,
                                                      const knobs& knobs,
                                                      pipeline_t& pipeline,
                                                      const backend_spec& spec,
                                                      observability::run_metrics& metrics) {
      std::vector<backend_subgraph> gpu_subgraphs;
      gpu_subgraphs.reserve(spec.ids.size());

      for (const auto id: spec.ids) {
        auto& backend_metrics = *metrics.backends.emplace_back(
            std::make_unique<observability::backend_metrics>(std::format("{}#{}", spec.configuration, id),
                                                             spec.workers_per_device));

        auto device_scratch = std::make_shared<scratchpad<queue_type>>(knobs, id, device_type::GPU);
        auto q_b            = device_scratch->get_queue();
        std::function<int(const int)> get_size =
            [q_b, &knobs, memory_bytes = spec.memory_bytes](const int atoms) {
              return pipeline_t::template get_batch_size<queue_type>(atoms, q_b, knobs, memory_bytes);
            };

        info(std::format("Setting up GPU device {} with {} workers, {} bytes each ...",
                         id,
                         spec.workers_per_device,
                         spec.memory_bytes));

        backend_subgraph gpu_subgraph(backend_metrics);
        gpu_subgraph.set_batcher(g, std::move(get_size));
        gpu_subgraph.reserve_screeners(spec.workers_per_device);

        for (std::size_t i = 0; i < spec.workers_per_device; ++i) {
          auto stage =
              pipeline.template get_pipeline<queue_type>(knobs, id, device_type::GPU, device_scratch);
          gpu_subgraph.add_screener(g, std::move(stage));
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
                                           const backend_spec& spec,
                                           observability::run_metrics& metrics) {
      std::optional<backend_subgraph> cpu_subgraph;

      constexpr_for<0, num_cpu_kernel_type(), 1>([&](const auto kernel) {
        constexpr auto kernel_type = cpu_kernel_type[kernel];
        if (kernel_type == spec.impl) {
          using queue_type_t = typename kernel_type_traits<kernel_type>::type;
          cpu_subgraph.emplace(build_cpu_subgraph<queue_type_t>(g, knobs, pipeline, spec, metrics));
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
    std::vector<backend_subgraph> dispatch_gpu_subgraphs(tbb_flow::graph& g,
                                                         const knobs& knobs,
                                                         pipeline_t& pipeline,
                                                         const backend_spec& spec,
                                                         observability::run_metrics& metrics) {
      std::vector<backend_subgraph> gpu_subgraphs;

      constexpr_for<0, num_gpu_kernel_type(), 1>([&](const auto kernel) {
        constexpr auto kernel_type = gpu_kernel_type[kernel];
        if (kernel_type == spec.impl) {
          using queue_type_t = typename kernel_type_traits<kernel_type>::type;
          gpu_subgraphs      = build_gpu_subgraphs<queue_type_t>(g, knobs, pipeline, spec, metrics);
        }
      });

      if (gpu_subgraphs.empty()) {
        throw std::runtime_error("Requested GPU implementation not compiled in.");
      }

      return gpu_subgraphs;
    }

  } // namespace setup

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
                                                        const std::vector<std::string>& configurations,
                                                        observability::run_metrics& metrics) {
    std::vector<backend_subgraph> backend_subgraphs;

    for (const auto& configuration: configurations) {
      const auto spec = parse_backend_spec(configuration);

      info(std::format("Setting up backend-spec '{}' ...", configuration));

      switch (spec.dev) {
        case device_type::CPU: {
          backend_subgraphs.push_back(setup::dispatch_cpu_subgraph(g, knobs, pipeline, spec, metrics));
          break;
        }
        case device_type::GPU: {
          auto gpu_subgraphs = setup::dispatch_gpu_subgraphs(g, knobs, pipeline, spec, metrics);
          std::ranges::move(gpu_subgraphs, std::back_inserter(backend_subgraphs));
          break;
        }
        default: throw std::runtime_error("Not supported device type.");
      }
    }

    return backend_subgraphs;
  }

} // namespace mudock::stream_engine
