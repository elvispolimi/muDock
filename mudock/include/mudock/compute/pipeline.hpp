#pragma once

#include <memory>
#include <mudock/compute/adt_score.hpp>
#include <mudock/compute/local_search.hpp>
#include <mudock/compute/adadelta.hpp>
#include <mudock/compute/genetic.hpp>
#include <mudock/compute/lamarckian_genetic.hpp>
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
    static int get_batch_size(const int, std::shared_ptr<queue_type>, const knobs&, const int) {
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
    static int get_batch_size(const int atoms,
                              std::shared_ptr<queue_type> q,
                              const knobs& conf,
                              const size_t max_mem = 1000000000) {
      const int mem = mudock::adt_score<queue_type>::get_ligand_mem(atoms, conf);
      return get_adt_score_batch<queue_type>(atoms, q, max_mem / mem);
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
    static int get_batch_size(const int atoms,
                              std::shared_ptr<queue_type> q,
                              const knobs& conf,
                              const size_t max_mem = 1000000000) {
      const int mem = genetic<queue_type, adt_score>::get_ligand_mem(atoms, conf);
      return get_adt_score_batch<queue_type>(atoms, q, max_mem / mem);
    }
  };

  // TODO i don't know if this is correct
  struct lga_adt_adadelta_pipeline: pipeline {
    template<typename queue_type>
    lamarckian_genetic<queue_type, adt_score, adadelta> get_pipeline(const knobs& conf,
                                                const int id,
                                                const device_type dev_type,
                                                std::shared_ptr<scratchpad<queue_type>> device_scratch) {
      auto q = std::make_shared<mudock::scratchpad<queue_type>>(conf, id, dev_type);
      auto scoring = std::make_shared<mudock::adt_score<queue_type>>(q, device_scratch, *protein);
      return lamarckian_genetic<queue_type, adt_score, adadelta>(q,
                                            *protein,
                                            scoring,
                                            mudock::adadelta<queue_type, adt_score>(q, scoring));
    }

    template<typename queue_type>
    static int get_batch_size(const int atoms,
                              std::shared_ptr<queue_type> q,
                              const knobs& conf,
                              const size_t max_mem = 1000000000) {
      const int mem = lamarckian_genetic<queue_type, adt_score, adadelta>::get_ligand_mem(atoms, conf);
      return get_adt_score_batch<queue_type>(atoms, q, max_mem / mem);
    }
  };

  // TODO L this should be local search generic and be able to accept an implementation. 
  // For the moment adadelta + adt_score is hardcoded.
  struct local_search_pipeline: pipeline {
    template<typename queue_type>
    adadelta<queue_type, adt_score> get_pipeline(const knobs& conf,
                                                      const int id,
                                                      const device_type dev_type,
                                                      std::shared_ptr<scratchpad<queue_type>> device_scratch) {
      auto q       = std::make_shared<mudock::scratchpad<queue_type>>(conf, id, dev_type);
      auto scoring = std::make_shared<mudock::adt_score<queue_type>>(q, device_scratch, *protein);
      return adadelta<queue_type, adt_score>(q, scoring);
    }

    template<typename queue_type>
    static int get_batch_size(const int atoms,
                              std::shared_ptr<queue_type> q,
                              const knobs& conf,
                              const size_t max_mem = 1000000000) {
      const int mem = adadelta<queue_type, adt_score>::get_ligand_mem(atoms, conf);
      return get_adt_score_batch<queue_type>(atoms, q, max_mem / mem);
    }
  };

} // namespace mudock
