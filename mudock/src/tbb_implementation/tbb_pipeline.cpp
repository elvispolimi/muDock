#include <mudock/tbb_implementation/stream_filter.hpp>
#include <mudock/tbb_implementation/parser_filter.hpp>
#include <mudock/tbb_implementation/tbb_pipeline.hpp>
#include <mudock/compute/manager.hpp>
#include <mudock/compute/safe_stack.hpp>
#include <mudock/molecule.hpp>
#include <mudock/mudock.hpp>

#include <oneapi/tbb/parallel_pipeline.h>
#include <memory>
#include <vector>
#include <iostream>

namespace mudock {

    static inline void printLigand(const static_molecule& ligand) {
        std::cout << ligand.properties.get(property_type::NAME) << " "
                  << ligand.properties.get(property_type::SCORE) << "\n";
    }
    
    void run_tbb_pipeline(std::istream& in,
             const std::vector<std::string>& configurations,
             const knobs& knobs,
             genetic_adt_pipeline& pipeline, 
             std::size_t max_tokens) {
        using VecP = parser_filter::molVec;

        auto input_queue  = std::make_shared<safe_stack<static_molecule>>();
        auto output_queue = std::make_shared<safe_stack<static_molecule>>();

        auto drain_ready = [&]() -> std::size_t {
            std::size_t counter = 0;
            while (auto x = output_queue->dequeue()) {
                printLigand(*x);
                ++counter;
            }
            return counter;
        };

        {
            // Call manager just once before the tbb pipeline
            threadpool pool;
            manager(configurations, pool, knobs, input_queue, output_queue, pipeline);
            info("Manager done: workers created");

            oneapi::tbb::parallel_pipeline(
              max_tokens,
              oneapi::tbb::make_filter<void, std::string>(
                oneapi::tbb::filter_mode::serial_in_order,
                stream_filter(in)
              )
              &
              oneapi::tbb::make_filter<std::string, VecP>(
                oneapi::tbb::filter_mode::parallel,
                parser_filter()
              )
              &
              oneapi::tbb::make_filter<VecP, std::size_t>(
                oneapi::tbb::filter_mode::serial_out_of_order,
                [=](VecP molecules) -> std::size_t {
                    std::size_t enq = 0;
                    for (auto& p : molecules) {
                        if (p) { input_queue->enqueue(std::move(p)); ++enq; }
                    }
                    return enq;
                }
              )
              &
              oneapi::tbb::make_filter<std::size_t, void>(
                oneapi::tbb::filter_mode::serial_out_of_order,
                [&](std::size_t) {
                    drain_ready();
                }
              )
            );

            info("Pipeline done: closing input queue");
            input_queue->close();
        }

        info("Draining output queue");
        drain_ready();
    }
    
} // namespace mudock
