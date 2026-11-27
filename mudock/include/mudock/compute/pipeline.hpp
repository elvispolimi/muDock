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
#include <mudock/devices.hpp>
#include <mudock/molecule.hpp>

namespace mudock {

  struct pipeline {
    pipeline(std::shared_ptr<dynamic_molecule> _protein,
             ligand_placement _placement = {}):
        protein(_protein), placement(_placement) {}

    template<typename queue_type>
    stage<queue_type> get_pipeline(const knobs&,
                                   const int,
                                   const device_type,
                                   std::shared_ptr<scratchpad<queue_type>>,
                                   dynamic_molecule&) {}

    template<typename queue_type>
    static int get_batch_size(const int, const knobs&, const size_t = 1000000000) {
      return 1;
    }

  protected:
    std::shared_ptr<dynamic_molecule> protein;
    ligand_placement placement;
  };

  template<scoring_function scoring>
  struct scoring_pipeline: pipeline {
    template<typename queue_type>
    using scoring_stage = typename scoring_backend<queue_type, scoring>::stage_type;

    scoring_pipeline(std::shared_ptr<dynamic_molecule> protein_ptr, ligand_placement placement_config = {})
        : pipeline(std::move(protein_ptr), placement_config) {}

    template<typename queue_type>
    scoring_stage<queue_type> get_pipeline(const knobs& conf,
                                           const int id,
                                           const device_type dev_type,
                                           std::shared_ptr<scratchpad<queue_type>> device_scratch) {
      return scoring_stage<queue_type>(std::make_shared<mudock::scratchpad<queue_type>>(conf, id, dev_type),
                                       device_scratch,
                                       *protein,
                                       placement);
    }

    template<typename queue_type>
    static int get_batch_size(const int atoms,
                              const knobs& conf,
                              const size_t max_mem = 1000000000) {
      using selected_scoring = scoring_stage<queue_type>;
      const size_t mem_per_ligand  = static_cast<size_t>(selected_scoring::get_ligand_mem(atoms, conf));
      const size_t max_bucket_size = std::max<size_t>(1, max_mem / mem_per_ligand);
      mudock::stage_bucket_trace("PIPELINE(",
                                 selected_scoring::stage_name,
                                 ") pre-resolve for ",
                                 atoms,
                                 " atoms: mem_budget=",
                                 max_mem,
                                 " B, mem_per_ligand=",
                                 mem_per_ligand,
                                 " B, max_bucket_size=",
                                 max_bucket_size);
      return resolve_stage_bucket_size(
          selected_scoring::stage_name, atoms, max_bucket_size, mem_per_ligand);
    }
  };

  template<scoring_function scoring>
  struct genetic_scoring_pipeline: pipeline {
    template<typename queue_type>
    using scoring_stage = typename scoring_backend<queue_type, scoring>::stage_type;

    genetic_scoring_pipeline(std::shared_ptr<dynamic_molecule> protein_ptr,
                             ligand_placement placement_config = {})
        : pipeline(std::move(protein_ptr), placement_config) {}

    template<typename queue_type>
    genetic<queue_type, scoring> get_pipeline(const knobs& conf,
                                              const int id,
                                              const device_type dev_type,
                                              std::shared_ptr<scratchpad<queue_type>> device_scratch) {
      auto q = std::make_shared<mudock::scratchpad<queue_type>>(conf, id, dev_type);
      return genetic<queue_type, scoring>(
          q,
          *protein,
          scoring_stage<queue_type>(
              q, device_scratch, *protein, ligand_placement{ligand_placement_mode::preserve, std::nullopt}),
          placement);
    }

    template<typename queue_type>
    static int get_batch_size(const int atoms,
                              const knobs& conf,
                              const size_t max_mem = 1000000000) {
      const size_t mem_per_ligand =
          static_cast<size_t>(genetic<queue_type, scoring>::get_ligand_mem(atoms, conf));
      const size_t max_bucket_size = std::max<size_t>(1, max_mem / mem_per_ligand);
      mudock::stage_bucket_trace("PIPELINE(",
                                 genetic<queue_type, scoring>::stage_name,
                                 ") pre-resolve for ",
                                 atoms,
                                 " atoms: mem_budget=",
                                 max_mem,
                                 " B, mem_per_ligand=",
                                 mem_per_ligand,
                                 " B, max_bucket_size=",
                                 max_bucket_size);
      return resolve_stage_bucket_size(
          genetic<queue_type, scoring>::stage_name, atoms, max_bucket_size, mem_per_ligand);
    }
  };

} // namespace mudock
