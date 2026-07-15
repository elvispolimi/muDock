#pragma once

#include <atomic>
#include <exception>
#include <iostream>
#include <memory>
#include <mudock/compute/safe_queue.hpp>
#include <mudock/format/reader.hpp>
#include <mudock/format/supported_format.hpp>
#include <mudock/molecule.hpp>
#include <mudock/mudock.hpp>
#include <string_view>

namespace mudock {

  template<supported_format format>
  class parser_filter {
  private:
    std::shared_ptr<safe_queue<static_molecule>> input_queue;
    std::atomic<std::size_t>* skipped_ligands = nullptr;
    std::atomic<bool>* stop_requested         = nullptr;
    bool parser_debug                         = false;

  public:
    explicit parser_filter(std::shared_ptr<safe_queue<static_molecule>> input,
                           std::atomic<std::size_t>* skipped = nullptr,
                           std::atomic<bool>* stop           = nullptr,
                           bool debug                         = false)
        : input_queue(std::move(input)),
          skipped_ligands(skipped),
          stop_requested(stop),
          parser_debug(debug) {}

    void operator()(std::string_view sv) const {
      type_of_format<format> splitter;

      if (stop_requested != nullptr && stop_requested->load(std::memory_order_relaxed)) {
        return;
      }

      while (!sv.empty()) {
        if (stop_requested != nullptr && stop_requested->load(std::memory_order_relaxed)) {
          return;
        }

        const auto next = splitter.next_molecule_start_index(sv);
        const auto mol  = (next == std::string_view::npos) ? sv : sv.substr(0, next);

        try {
          auto ligand = std::make_unique<static_molecule>(mudock::parser<format, static_molecule>(mol));
          const bool enqueued = input_queue->enqueue(ligand);
          if (!enqueued) {
            return;
          }
        } catch (const std::exception& e) {
          if (parser_debug) {
            std::cerr << "Skipping ligand due to parse error: " << e.what() << '\n';
          }
          if (skipped_ligands != nullptr) {
            skipped_ligands->fetch_add(1, std::memory_order_relaxed);
          }
        }

        if (next == std::string_view::npos)
          break;
        sv.remove_prefix(next);
      }
    }
  };

} // namespace mudock
