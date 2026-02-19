#include <cstdint>
#include <iostream>
#include <memory>
#include <mudock/tbb_implementation/stream_filter.hpp>
#include <mudock/tbb_implementation/parser_filter.hpp>
#include <mudock/tbb_implementation/tbb_pipeline.hpp>
#include <mudock/compute/manager.hpp>
#include <mudock/compute/safe_stack.hpp>
#include <mudock/molecule.hpp>
#include <mudock/mudock.hpp>
#include <oneapi/tbb/parallel_pipeline.h>
#include <thread>
#include <vector>

namespace mudock {

    static inline void print_ligand(const static_molecule& ligand) {
        std::cout << ligand.properties.get(property_type::NAME) << " "
                  << ligand.properties.get(property_type::SCORE) << "\n";
    }

    void run_tbb_pipeline(std::istream& in,
         const std::vector<std::string>& configurations,
         const knobs& knobs,
         genetic_adt_pipeline& pipeline, 
         std::size_t end)   
    {
        using mol_vec = parser_filter::mol_vec;

        auto input_queue  = std::make_shared<mudock::safe_queue<mudock::static_molecule>>();
        auto output_queue = std::make_shared<mudock::safe_queue<mudock::static_molecule>>();
        input_queue->initialize(1000000);
        output_queue->initialize(1000000);

        // asynchronous thread to print the output ligands as soon as they are ready
        std::thread writer([&]{
            bool ok=false;
            std::string buf;
            buf.reserve(1<<20);

            std::size_t lines = 0;
            for (auto x = output_queue->dequeue(ok); ok; x = output_queue->dequeue(ok)) {
                buf += x->properties.get(property_type::NAME);
                buf += ' ';
                buf += x->properties.get(property_type::SCORE);
                buf += '\n';

                if (++lines % 4096 == 0) {
                    std::cout << buf;
                    buf.clear();
                }
            }
            if (!buf.empty()) std::cout << buf;
        });

        // Helper lambdas to drain results with tbb writer stage (disabled until we integrate the writer stage in the pipeline)
        // auto drain_ready_nb = [&]() {
        //     mudock::safe_queue<mudock::static_molecule>::value_ptr_type x;
        //     while (output_queue->try_dequeue(x)) {     
        //         print_ligand(*x);
        //     }
        // };

        // auto drain_ready = [&]() -> std::size_t {
        //     std::size_t counter = 0;
        //     bool is_retrieved = false;
        //     while (auto x = output_queue->dequeue(is_retrieved)) {
        //         print_ligand(*x);
        //         ++counter;
        //     }
            
        //     return counter;
        // };

        threadpool pool;
        manager(configurations, pool, knobs, input_queue, output_queue, pipeline);
        info("Manager done: workers created");

        const std::size_t effective_max_tokens = knobs.max_tbb_tokens == 0 ? 1 : knobs.max_tbb_tokens;

        oneapi::tbb::parallel_pipeline(
            effective_max_tokens,
            oneapi::tbb::make_filter<void, std::string>(
            oneapi::tbb::filter_mode::serial_in_order,
            stream_filter(in, end)
            )
            &
            oneapi::tbb::make_filter<std::string, mol_vec>(
            oneapi::tbb::filter_mode::parallel,
            parser_filter()
            )
            &
            // A goal is to replace this filter with a tbb stage that directly 
            // computes the score of the ligands and outputs them in parallel, 
            // without using an external thread pool.
            oneapi::tbb::make_filter<mol_vec, void>(
                oneapi::tbb::filter_mode::serial_out_of_order,
            [&](mol_vec molecules) {
                bool is_produced = false;
                for (auto& p : molecules) {
                    if (p) {
                        input_queue->enqueue(p, is_produced);
                    }
                    if (!is_produced) break;
                }
            }
            )
            // &
            // oneapi::tbb::make_filter<std::size_t, void>(
            //         oneapi::tbb::filter_mode::parallel,
            //     [&] (std::size_t) {
            //         ...
            //     }
            // )
        );

        info("Pipeline done: closing input queue");
        input_queue->send_terminate_signal();
        pool.wait();
        output_queue->send_terminate_signal();
        writer.join();

        info("Output drained");
    }

} // namespace mudock
