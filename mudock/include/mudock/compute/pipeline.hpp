#pragma once

#include <algorithm>
#include <memory>
#include <mudock/compute/adt_score.hpp>
#include <mudock/compute/local_search.hpp>
#include <mudock/compute/adadelta.hpp>
#include <mudock/compute/algorithm.hpp>
#include <mudock/compute/bucket_size.hpp>
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
    static int get_batch_size(const int, std::shared_ptr<queue_type>, const knobs&, const size_t) {
      return 1;
    }

  protected:
    std::shared_ptr<dynamic_molecule> protein;
  };

  template<template<typename> typename scoring_t>
  struct scoring_pipeline: pipeline {
    using pipeline::pipeline;

    template<typename queue_type>
    scoring_t<queue_type> get_pipeline(const knobs& conf,
                                       const int id,
                                       const device_type dev_type,
                                       std::shared_ptr<scratchpad<queue_type>> device_scratch) {
      return scoring_t<queue_type>(std::make_shared<mudock::scratchpad<queue_type>>(conf, id, dev_type),
                                   device_scratch,
                                   *protein);
    }

    template<typename queue_type>
    static int get_batch_size(const int atoms,
                              std::shared_ptr<queue_type> q,
                              const knobs& conf,
                              const size_t max_mem = 1000000000) {
      const size_t mem_per_ligand  = static_cast<size_t>(scoring_t<queue_type>::get_ligand_mem(atoms, conf));
      const size_t max_bucket_size = std::max<size_t>(1, max_mem / mem_per_ligand);
      mudock::stage_bucket_trace("PIPELINE(",
                                 scoring_t<queue_type>::stage_name,
                                 ") pre-resolve for ",
                                 atoms,
                                 " atoms: mem_budget=",
                                 max_mem,
                                 " B, mem_per_ligand=",
                                 mem_per_ligand,
                                 " B, max_bucket_size=",
                                 max_bucket_size);
      return resolve_stage_bucket_size(
          scoring_t<queue_type>::stage_name,
          atoms,
          max_bucket_size,
          mem_per_ligand,
          q->honors_stage_bucket_policy(),
          [&]() { return scoring_t<queue_type>::get_batch_size(atoms, q, conf, max_bucket_size); });
    }
  };

  template<template<typename> typename scoring_t>
  struct genetic_scoring_pipeline: pipeline {
    using pipeline::pipeline;

    template<typename queue_type>
    genetic<queue_type, scoring_t> get_pipeline(const knobs& conf,
                                                const int id,
                                                const device_type dev_type,
                                                std::shared_ptr<scratchpad<queue_type>> device_scratch) {
      auto q = std::make_shared<mudock::scratchpad<queue_type>>(conf, id, dev_type);
      auto scoring = std::make_shared<scoring_t<queue_type>>(q, device_scratch, *protein);
      return genetic<queue_type, scoring_t>(q, *protein, scoring);
    }

    template<typename queue_type>
    static int get_batch_size(const int atoms,
                              std::shared_ptr<queue_type> q,
                              const knobs& conf,
                              const size_t max_mem = 1000000000) {
      const size_t mem_per_ligand =
          static_cast<size_t>(genetic<queue_type, scoring_t>::get_ligand_mem(atoms, conf));
      const size_t max_bucket_size = std::max<size_t>(1, max_mem / mem_per_ligand);
      mudock::stage_bucket_trace("PIPELINE(",
                                 genetic<queue_type, scoring_t>::stage_name,
                                 ") pre-resolve for ",
                                 atoms,
                                 " atoms: mem_budget=",
                                 max_mem,
                                 " B, mem_per_ligand=",
                                 mem_per_ligand,
                                 " B, max_bucket_size=",
                                 max_bucket_size);
      return resolve_stage_bucket_size(
          genetic<queue_type, scoring_t>::stage_name,
          atoms,
          max_bucket_size,
          mem_per_ligand,
          q->honors_stage_bucket_policy(),
          [&]() { return genetic<queue_type, scoring_t>::get_batch_size(atoms, q, conf, max_bucket_size); });
    }
  };

  template<
      template<typename> typename scoring_t,
      template<typename, template<typename> typename> typename local_search_t
  >
  struct lamarckian_genetic_scoring_pipeline: pipeline {
    using pipeline::pipeline;

    template<typename queue_type>
    lamarckian_genetic<
        queue_type, 
        scoring_t, 
        local_search_t
    >
    get_pipeline(const knobs& conf,
                 const int id,
                 const device_type dev_type,
                 std::shared_ptr<scratchpad<queue_type>> device_scratch) {
      auto q = std::make_shared<mudock::scratchpad<queue_type>>(conf, id, dev_type);
      auto scoring = std::make_shared<scoring_t<queue_type>>(q, device_scratch, *protein);
      auto local_search = local_search_t<queue_type, scoring_t>(q, scoring);
      return lamarckian_genetic<queue_type, scoring_t, local_search_t>(q, *protein, scoring, std::move(local_search));
    }

    template<typename queue_type>
    static int get_batch_size(const int atoms,
                              std::shared_ptr<queue_type> q,
                              const knobs& conf,
                              const size_t max_mem = 1000000000) {
      const size_t mem_per_ligand =
          static_cast<size_t>(lamarckian_genetic<queue_type, scoring_t, local_search_t>::get_ligand_mem(atoms, conf));
      const size_t max_bucket_size = std::max<size_t>(1, max_mem / mem_per_ligand);
      mudock::stage_bucket_trace("PIPELINE(",
                                 lamarckian_genetic<queue_type, scoring_t, local_search_t>::stage_name,
                                 ") pre-resolve for ",
                                 atoms,
                                 " atoms: mem_budget=",
                                 max_mem,
                                 " B, mem_per_ligand=",
                                 mem_per_ligand,
                                 " B, max_bucket_size=",
                                 max_bucket_size);
      return resolve_stage_bucket_size(lamarckian_genetic<queue_type, scoring_t, local_search_t>::stage_name,
                                       atoms,
                                       max_bucket_size,
                                       mem_per_ligand,
                                       q->honors_stage_bucket_policy(),
                                       [&]() {
                                         return lamarckian_genetic<queue_type, scoring_t, local_search_t>::get_batch_size(atoms,
                                                                                                q,
                                                                                                conf,
                                                                                                max_bucket_size);
                                       });
    }
  };


  template<
      template<typename> typename scoring_t,
      template<typename, template<typename> typename> typename local_search_t
  >
  struct local_search_pipeline: pipeline {
    using pipeline::pipeline;

    static knobs normalize_knobs(knobs conf) {
      conf.population_number = 1;
      conf.num_generations   = 1;
      conf.lsrate            = 100;
      return conf;
    }

    template<typename queue_type>
    local_search_t<queue_type, scoring_t> get_pipeline(const knobs& conf,
                                                       const int id,
                                                       const device_type dev_type,
                                                       std::shared_ptr<scratchpad<queue_type>> device_scratch) {
      const auto effective_conf = normalize_knobs(conf);
      auto q = std::make_shared<mudock::scratchpad<queue_type>>(effective_conf, id, dev_type);
      auto scoring = std::make_shared<scoring_t<queue_type>>(q, device_scratch, *protein);
      return local_search_t<queue_type, scoring_t>(q, scoring);
    }

    template<typename queue_type>
    static int get_batch_size(const int atoms,
                              std::shared_ptr<queue_type> q,
                              const knobs& conf,
                              const size_t max_mem = 1000000000) {
      const auto effective_conf = normalize_knobs(conf);
      const size_t mem_per_ligand = static_cast<size_t>(local_search_t<queue_type, scoring_t>::get_ligand_mem(atoms, effective_conf));
      const size_t max_bucket_size = std::max<size_t>(1, max_mem / mem_per_ligand);
      mudock::stage_bucket_trace("PIPELINE(",
                                 local_search_t<queue_type, scoring_t>::stage_name,
                                 ") pre-resolve for ",
                                 atoms,
                                 " atoms: mem_budget=",
                                 max_mem,
                                 " B, mem_per_ligand=",
                                 mem_per_ligand,
                                 " B, max_bucket_size=",
                                 max_bucket_size);
      return resolve_stage_bucket_size(local_search_t<queue_type, scoring_t>::stage_name,
                                       atoms,
                                       max_bucket_size,
                                       mem_per_ligand,
                                       q->honors_stage_bucket_policy(),
                                       [&]() {
                                         return local_search_t<queue_type, scoring_t>::get_batch_size(atoms,
                                                                                                q,
                                                                                                effective_conf,
                                                                                                max_bucket_size);
                                       });
    }
  };

  using adt_score_pipeline           = scoring_pipeline<adt_score>;
  using genetic_adt_pipeline         = genetic_scoring_pipeline<adt_score>;
  using lga_adt_adadelta_pipeline    = lamarckian_genetic_scoring_pipeline<adt_score, adadelta>;
  using ls_adt_adadelta_pipeline     = local_search_pipeline<adt_score, adadelta>;
} // namespace mudock
