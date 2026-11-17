#pragma once

#include <algorithm>
#include <memory>
#include <mudock/compute/adt_score.hpp>
#include <mudock/compute/algorithm.hpp>
#include <mudock/compute/bucket_size.hpp>
#include <mudock/compute/genetic.hpp>
#include <mudock/data_formats/data_store.hpp>
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
                                   std::shared_ptr<data_store<queue_type>>,
                                   dynamic_molecule&) {}

    template<typename queue_type>
    static int get_batch_size(const int, std::shared_ptr<queue_type>, const knobs&, const size_t) {
      return 1;
    }

  protected:
    std::shared_ptr<dynamic_molecule> protein;
  };

  template<template<typename, typename> typename scoring_t>
  struct scoring_pipeline: pipeline {
    using pipeline::pipeline;
    using pipeline_data_layouts = scoring_data_layouts;

    template<typename queue_type>
    scoring_t<queue_type, scoring_data_layouts> get_pipeline(const knobs& conf,
                                       const int id,
                                       const device_type dev_type,
                                       std::shared_ptr<data_store<queue_type, protein_data_layouts>> protein_data) {
        return scoring_t<queue_type, scoring_data_layouts>(
          std::make_shared<mudock::data_store<queue_type, scoring_data_layouts>>(
                                       conf, id, dev_type, protein_data->get_memory_tracker()),
                                   protein_data,
                                   *protein);
    }

    template<typename queue_type>
    static int get_batch_size(const int atoms,
                              std::shared_ptr<queue_type> q,
                              const knobs& conf,
                              const size_t max_mem = 1000000000) {
        const size_t mem_per_ligand =
          static_cast<size_t>(scoring_t<queue_type, scoring_data_layouts>::get_ligand_mem(atoms, conf));
      const size_t max_bucket_size = std::max<size_t>(1, max_mem / mem_per_ligand);
      mudock::stage_bucket_trace("PIPELINE(",
                                 scoring_t<queue_type, scoring_data_layouts>::stage_name,
                                 ") pre-resolve for ",
                                 atoms,
                                 " atoms: mem_budget=",
                                 max_mem,
                                 " B, mem_per_ligand=",
                                 mem_per_ligand,
                                 " B, max_bucket_size=",
                                 max_bucket_size,
                                 ", estimated_worker_batch_bytes=",
                                 max_bucket_size * mem_per_ligand);
      return resolve_stage_bucket_size(
          scoring_t<queue_type, scoring_data_layouts>::stage_name,
          atoms,
          max_bucket_size,
          mem_per_ligand,
          q->honors_stage_bucket_policy(),
          [&]() {
            return scoring_t<queue_type, scoring_data_layouts>::get_batch_size(atoms, q, conf, max_bucket_size);
          });
    }
  };

  template<template<typename, typename> typename scoring_t>
  struct genetic_scoring_pipeline: pipeline {
    using pipeline::pipeline;
    using pipeline_data_layouts = genetic_data_layouts;

    template<typename queue_type>
    genetic<queue_type, scoring_t, genetic_data_layouts> get_pipeline(const knobs& conf,
                                                const int id,
                                                const device_type dev_type,
                                                std::shared_ptr<data_store<queue_type, protein_data_layouts>> protein_data) {
        auto q = std::make_shared<mudock::data_store<queue_type, genetic_data_layouts>>(
          conf, id, dev_type, protein_data->get_memory_tracker());
        return genetic<queue_type, scoring_t, genetic_data_layouts>(
          q, *protein, scoring_t<queue_type, genetic_data_layouts>(q, protein_data, *protein));
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
                                 max_bucket_size,
                                 ", estimated_worker_batch_bytes=",
                                 max_bucket_size * mem_per_ligand);
      return resolve_stage_bucket_size(
          genetic<queue_type, scoring_t>::stage_name,
          atoms,
          max_bucket_size,
          mem_per_ligand,
          q->honors_stage_bucket_policy(),
          [&]() { return genetic<queue_type, scoring_t>::get_batch_size(atoms, q, conf, max_bucket_size); });
    }
  };

  using adt_score_pipeline   = scoring_pipeline<adt_score>;
  using genetic_adt_pipeline = genetic_scoring_pipeline<adt_score>;
} // namespace mudock
