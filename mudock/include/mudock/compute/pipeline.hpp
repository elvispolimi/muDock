#pragma once

#include <algorithm>
#include <memory>
#include <mudock/compute/adt_score.hpp>
#include <mudock/compute/algorithm.hpp>
#include <mudock/compute/bucket_size.hpp>
#include <mudock/compute/genetic.hpp>
#include <mudock/compute/ligand_placement.hpp>
#include <mudock/compute/scratchpad.hpp>
#include <mudock/compute/stage.hpp>
#include <mudock/chem/autodock_grid.hpp>
#include <mudock/devices.hpp>
#include <mudock/molecule.hpp>
#include <stdexcept>

namespace mudock {

  struct pipeline {
    pipeline(std::shared_ptr<autodock_grid> _grid, ligand_placement _placement = {})
        : grid(std::move(_grid)), center(grid_center(grid)), placement(_placement) {}
    pipeline(std::shared_ptr<dynamic_molecule> _protein, ligand_placement _placement = {})
        : protein(std::move(_protein)), center(protein_center(protein)), placement(_placement) {}

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
    static point3D grid_center(const std::shared_ptr<autodock_grid>& value) {
      if (!value) throw std::invalid_argument("AutoDock scoring pipeline requires a grid");
      return value->get_center();
    }
    static point3D protein_center(const std::shared_ptr<dynamic_molecule>& value) {
      if (!value) throw std::invalid_argument("Protein scoring pipeline requires a protein");
      return value->get_center();
    }

    std::shared_ptr<dynamic_molecule> protein;
    std::shared_ptr<autodock_grid> grid;
    point3D center;
    ligand_placement placement;
  };

  template<template<typename> typename scoring_t>
  struct scoring_pipeline: pipeline {
    scoring_pipeline(std::shared_ptr<autodock_grid> grid_ptr, ligand_placement placement_config = {})
        : pipeline(std::move(grid_ptr), placement_config) {}
    scoring_pipeline(std::shared_ptr<dynamic_molecule> protein_ptr, ligand_placement placement_config = {})
        : pipeline(std::move(protein_ptr), placement_config) {}

    template<typename queue_type>
    scoring_t<queue_type> get_pipeline(const knobs& conf,
                                       const int id,
                                       const device_type dev_type,
                                       std::shared_ptr<scratchpad<queue_type>> device_scratch) {
      auto scratch = std::make_shared<mudock::scratchpad<queue_type>>(conf, id, dev_type);
      if constexpr (requires { scoring_t<queue_type>{scratch, device_scratch, *grid, placement}; }) {
        if (!grid) throw std::invalid_argument("This scoring function requires AutoGrid maps");
        return scoring_t<queue_type>{scratch, device_scratch, *grid, placement};
      } else {
        if (!protein) throw std::invalid_argument("This scoring function requires a protein structure");
        return scoring_t<queue_type>{scratch, device_scratch, *protein, placement};
      }
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
    genetic_scoring_pipeline(std::shared_ptr<autodock_grid> grid_ptr,
                             ligand_placement placement_config = {})
        : pipeline(std::move(grid_ptr), placement_config) {}
    genetic_scoring_pipeline(std::shared_ptr<dynamic_molecule> protein_ptr,
                             ligand_placement placement_config = {})
        : pipeline(std::move(protein_ptr), placement_config) {}

    template<typename queue_type>
    genetic<queue_type, scoring_t> get_pipeline(const knobs& conf,
                                                const int id,
                                                const device_type dev_type,
                                                std::shared_ptr<scratchpad<queue_type>> device_scratch) {
      auto q = std::make_shared<mudock::scratchpad<queue_type>>(conf, id, dev_type);
      const ligand_placement preserve{ligand_placement_mode::preserve, std::nullopt};
      auto score = [&]() {
        if constexpr (requires { scoring_t<queue_type>{q, device_scratch, *grid, preserve}; }) {
          if (!grid) throw std::invalid_argument("This scoring function requires AutoGrid maps");
          return scoring_t<queue_type>{q, device_scratch, *grid, preserve};
        } else {
          if (!protein) throw std::invalid_argument("This scoring function requires a protein structure");
          return scoring_t<queue_type>{q, device_scratch, *protein, preserve};
        }
      }();
      return genetic<queue_type, scoring_t>(q, center, std::move(score), placement);
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

  using adt_score_pipeline   = scoring_pipeline<adt_score>;
  using genetic_adt_pipeline = genetic_scoring_pipeline<adt_score>;
} // namespace mudock
