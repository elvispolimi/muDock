#pragma once

#include <memory>
#include <mudock/compute/adt_score.hpp>
#include <mudock/compute/vina_score.hpp>
#include <mudock/compute/genetic.hpp>
#include <mudock/compute/scratchpad.hpp>
#include <mudock/compute/stage.hpp>
#include <mudock/devices.hpp>
#include <mudock/molecule.hpp>

namespace mudock {
  struct pipeline {
    pipeline(std::shared_ptr<dynamic_molecule> _protein): protein(_protein) {}

    template<typename queue_type>
    stage<queue_type> get_pipeline(const knobs&,
                                   const int,
                                   const device_type,
                                   std::shared_ptr<scratchpad<queue_type>>,
                                   dynamic_molecule&) {}

    template<typename queue_type>
    static int get_batch_size(const int, std::shared_ptr<queue_type>) {
      return 1;
    }

  protected:
    std::shared_ptr<dynamic_molecule> protein;
  };

  struct adt_score_pipeline: pipeline {
    template<typename queue_type>
    adt_score<queue_type> get_pipeline(const knobs& conf,
                                       const int id,
                                       const device_type dev_type,
                                       std::shared_ptr<scratchpad<queue_type>> device_scratch) {
      return mudock::adt_score<queue_type>(
          std::make_shared<mudock::scratchpad<queue_type>>(conf, id, dev_type),
          device_scratch,
          *protein);
    }

    template<typename queue_type>
    static int get_batch_size(const int atoms, std::shared_ptr<queue_type> q) {
      return get_adt_score_batch<queue_type>(atoms, q);
    }
  };

  struct genetic_adt_pipeline: pipeline {
    template<typename queue_type>
    genetic<queue_type, adt_score> get_pipeline(const knobs& conf,
                                                const int id,
                                                const device_type dev_type,
                                                std::shared_ptr<scratchpad<queue_type>> device_scratch) {
      auto q = std::make_shared<mudock::scratchpad<queue_type>>(conf, id, dev_type);
      return genetic<queue_type, adt_score>(q,
                                            *protein,
                                            mudock::adt_score<queue_type>(q, device_scratch, *protein));
    }

    template<typename queue_type>
    static int get_batch_size(const int atoms, std::shared_ptr<queue_type> q) {
      return get_adt_score_batch<queue_type>(atoms, q);
    }
  };

  struct genetic_vina_pipeline : pipeline {
    template<typename queue_type>
    genetic<queue_type, vina_score> get_pipeline(const knobs& conf,
                                                const int id,
                                                const device_type dev_type,
                                                std::shared_ptr<scratchpad<queue_type>> device_scratch) {
      auto q = std::make_shared<mudock::scratchpad<queue_type>>(conf, id, dev_type);
      return genetic<queue_type, vina_score>(q,
                                            *protein,
                                            mudock::vina_score<queue_type>(q, device_scratch, *protein));
    }

    template<typename queue_type>
    static int get_batch_size(const int atoms, std::shared_ptr<queue_type> q) {
      return get_vina_score_batch<queue_type>(atoms, q);
    }
  };

} // namespace mudock
