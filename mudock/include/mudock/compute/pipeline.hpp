#pragma once

#include <memory>
#include <mudock/compute/adt_score.hpp>
#include <mudock/compute/genetic.hpp>
#include <mudock/compute/scratchpad.hpp>
#include <mudock/compute/stage.hpp>
#include <mudock/implementations.hpp>
#include <mudock/molecule.hpp>

namespace mudock {
  struct pipeline {
    pipeline(std::shared_ptr<dynamic_molecule> _protein): protein(_protein) {}

    template<typename queue_type>
    stage<queue_type>
        get_pipeline(const knobs&, const int, std::shared_ptr<scratchpad<queue_type>>, dynamic_molecule&) {}

  protected:
    std::shared_ptr<dynamic_molecule> protein;
  };

  struct adt_score_pipeline: pipeline {
    template<typename queue_type>
    adt_score<queue_type> get_pipeline(const knobs& conf,
                                       const int id,
                                       std::shared_ptr<scratchpad<queue_type>> device_scratch) {
      return mudock::adt_score<queue_type>(std::make_shared<mudock::scratchpad<queue_type>>(conf, id),
                                           device_scratch,
                                           *protein);
    }
  };

  struct genetic_adt_pipeline: pipeline {
    template<typename queue_type>
    genetic<queue_type, adt_score> get_pipeline(const knobs& conf,
                                                const int id,
                                                std::shared_ptr<scratchpad<queue_type>> device_scratch) {
      auto q = std::make_shared<mudock::scratchpad<queue_type>>(conf, id);
      return genetic<queue_type, adt_score>(q,
                                            *protein,
                                            mudock::adt_score<queue_type>(q, device_scratch, *protein));
    }
  };
} // namespace mudock
