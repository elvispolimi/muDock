#pragma once

#include <atomic>
#include <memory>
#include <mudock/compute/safe_queue.hpp>
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

  public:
    explicit parser_filter(std::shared_ptr<safe_queue<static_molecule>> input,
                           std::atomic<std::size_t>* skipped = nullptr,
                           std::atomic<bool>* stop = nullptr)
        : input_queue(std::move(input)),
          skipped_ligands(skipped),
          stop_requested(stop) {}

    void operator()(std::string_view sv) const;
  };

} // namespace mudock
