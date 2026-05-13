#include <cstdint>
#include <iostream>
#include <memory>
#include <chrono>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <vector>
#include <mudock/compute/manager.hpp>
#include <mudock/molecule.hpp>
#include <mudock/mudock.hpp>
#include <mudock/tbb_implementation/parser_filter.hpp>
#include <mudock/tbb_implementation/runtime_services.hpp>
#include <mudock/tbb_implementation/stream_filter.hpp>
#include <mudock/tbb_implementation/tbb_pipeline.hpp>
#include <oneapi/tbb/parallel_pipeline.h>

namespace mudock {

    inline void append_ligand(std::string& buf, const static_molecule& ligand) {
        buf += ligand.properties.get(property_type::NAME);
        buf += ' ';
        buf += ligand.properties.get(property_type::SCORE);
        buf += '\n';
    }

    template<supported_format format, typename pipeline_t>
    void run_tbb_pipeline(std::istream& in,
         const std::vector<std::string>& configurations,
         const knobs& knobs,
         pipeline_t& pipeline,
         std::size_t end,
         std::optional<double> time_limit_sec,
         std::optional<double> observer_sec)
    {
        auto input_queue  = std::make_shared<mudock::safe_queue<mudock::static_molecule>>();
        auto output_queue = std::make_shared<mudock::safe_queue<mudock::static_molecule>>();
        input_queue->initialize(knobs.max_tbb_queue_size);
        output_queue->initialize(knobs.max_tbb_queue_size);
        std::atomic<std::size_t> skipped_ligands{0};
        std::atomic<std::size_t> dropped_by_timeout{0};
        std::atomic<std::size_t> in_flight_ligands{0};
        std::atomic<bool> timeout_triggered{false};
        std::atomic<bool> stop_requested{false};
        const auto start = std::chrono::high_resolution_clock::now();

        // asynchronous thread to print the output ligands as soon as they are ready
        std::thread writer([&] {
            std::string buf;
            buf.reserve(1 << 20);

            std::size_t lines = 0;
            for (auto x = output_queue->dequeue(); x; x = output_queue->dequeue()) {
                append_ligand(buf, *x);

                if (++lines % 4096 == 0) {
                    std::cout << buf;
                    buf.clear();
                }
            }
            if (!buf.empty()) std::cout << buf;
        });

        threadpool pool;

        std::size_t prev_processed = output_queue->get_global_counter();
        auto prev_time             = std::chrono::high_resolution_clock::now();
        detail::periodic_observer observer;
        if (observer_sec && *observer_sec > 0.0) {
            info("Observer enabled with period: ", *observer_sec, " s");
            observer.start(observer_sec, [&]() {
                const auto now                  = std::chrono::high_resolution_clock::now();
                const std::size_t now_processed = output_queue->get_global_counter();
                const std::size_t in_backlog    = input_queue->size();
                const std::size_t in_flight     = in_flight_ligands.load(std::memory_order_relaxed);

                const std::chrono::duration<double> dt = now - prev_time;
                const std::size_t delta_processed      = now_processed - prev_processed;
                const double inst_throughput =
                    dt.count() > 0.0 ? static_cast<double>(delta_processed) / dt.count() : 0.0;
                const std::chrono::duration<double> total = now - start;
                const double avg_throughput =
                    total.count() > 0.0 ? static_cast<double>(now_processed) / total.count() : 0.0;

                info("Observer: processed=",
                     now_processed,
                     ", input_backlog=",
                     in_backlog,
                     ", in_flight=",
                     in_flight,
                     ", inst_throughput=",
                     inst_throughput,
                     " ligands/s, avg_throughput=",
                     avg_throughput,
                     " ligands/s");

                prev_processed = now_processed;
                prev_time      = now;
            });
        }

        detail::deadline_timer timer;
        if (time_limit_sec && *time_limit_sec > 0.0) {
            info("Time limit enabled: ", *time_limit_sec, " s");
            timer.start(time_limit_sec, [&]() {
                stop_requested.store(true, std::memory_order_relaxed);
                input_queue->send_terminate_signal();
                dropped_by_timeout.store(input_queue->clear(), std::memory_order_relaxed);
                timeout_triggered.store(true, std::memory_order_relaxed);
                info("Time limit reached: discarded ",
                     dropped_by_timeout.load(std::memory_order_relaxed),
                     " pending ligand(s) from input queue.");
            });
        }

        bool finalized = false;
        auto finalize_runtime = [&]() {
            if (finalized) {
                return;
            }
            finalized = true;

            stop_requested.store(true, std::memory_order_relaxed);
            input_queue->send_terminate_signal();
            pool.wait();
            output_queue->send_terminate_signal();

            if (writer.joinable()) {
                writer.join();
            }
            observer.stop();
            observer.join();
            timer.cancel();
            timer.join();
        };

        try {
            manager(configurations, pool, knobs, input_queue, output_queue, pipeline, &in_flight_ligands);
            info("Manager done: workers created");

            oneapi::tbb::parallel_pipeline(
                knobs.max_tbb_tokens,
                oneapi::tbb::make_filter<void, std::string>(
                    oneapi::tbb::filter_mode::serial_in_order,
                    stream_filter<format>(in, knobs.max_bytes_per_token, end, &stop_requested))
                    & oneapi::tbb::make_filter<std::string, void>(
                        oneapi::tbb::filter_mode::parallel,
                        parser_filter<format>(input_queue, &skipped_ligands, &stop_requested)));

            info("Pipeline done: closing input queue");
            finalize_runtime();
        } catch (...) {
            finalize_runtime();
            throw;
        }

        if (const auto skipped = skipped_ligands.load(std::memory_order_relaxed); skipped > 0) {
            error("Skipped ", skipped, " ligand(s) due to parse errors.");
        }
        if (timeout_triggered.load(std::memory_order_relaxed)) {
            info("Dropped ligands due to timeout: ", dropped_by_timeout.load(std::memory_order_relaxed));
        }
        info("Output drained");
    }

    template void mudock::run_tbb_pipeline<mudock::supported_format::ADTMOL2, mudock::genetic_adt_pipeline>(
        std::istream&,
        const std::vector<std::string>&,
        const mudock::knobs&,
        mudock::genetic_adt_pipeline&,
        std::size_t,
        std::optional<double>,
        std::optional<double>);

    template void mudock::run_tbb_pipeline<mudock::supported_format::MOL2, mudock::genetic_adt_pipeline>(
        std::istream&,
        const std::vector<std::string>&,
        const mudock::knobs&,
        mudock::genetic_adt_pipeline&,
        std::size_t,
        std::optional<double>,
        std::optional<double>);

} // namespace mudock
