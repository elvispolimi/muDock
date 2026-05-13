#include <mudock/tbb_implementation/parser_filter.hpp>
#include <mudock/format/reader.hpp>
#include <exception>
#include <string_view>

namespace mudock {

  template<supported_format format>
  void parser_filter<format>::operator()(std::string_view sv) const {
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
      } catch (const std::exception&) {
        if (skipped_ligands != nullptr) {
          skipped_ligands->fetch_add(1, std::memory_order_relaxed);
        }
      }

      if (next == std::string_view::npos) break;
      sv.remove_prefix(next);
    }
  }

  template class parser_filter<supported_format::ADTMOL2>;
  template class parser_filter<supported_format::MOL2>;

} // namespace mudock
