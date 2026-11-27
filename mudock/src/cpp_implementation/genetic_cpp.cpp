#include "mudock/type_alias.hpp"

#include <algorithm>
#include <cmath>
#include <mudock/compute/devices_memory.hpp>
#include <mudock/compute/genetic.hpp>
#include <mudock/compute/pipeline_registry.hpp>
#include <mudock/cpp_implementation/mutate.hpp>
#include <mudock/cpp_implementation/queue_cpp.hpp>
#include <mudock/cpp_implementation/scoring_cpp.hpp>
#include <mudock/cpp_implementation/vectorization.hpp>
#include <mudock/utils.hpp>
#include <random>
#include <stdexcept>
#include <tuple>

namespace mudock {
  inline constexpr char genetic_compute_region_name[] = "genetic_compute";
  namespace {
    template<typename queue_t>
    inline constexpr cpu_vectorization cpu_vectorization_for = cpu_vectorization::AUTO;

#ifdef MUDOCK_USE_GH
    template<>
    inline constexpr cpu_vectorization cpu_vectorization_for<queue_gh> = cpu_vectorization::GH;
#endif

#ifdef MUDOCK_USE_XSIMD
    template<>
    inline constexpr cpu_vectorization cpu_vectorization_for<queue_xsimd> = cpu_vectorization::XSIMD;
#endif

    constexpr fp_type coordinate_step = static_cast<fp_type>(0.2);
    constexpr fp_type angle_step      = static_cast<fp_type>(4);
    thread_local device_memory<std::mt19937> random_generator;

    fp_type random_unit(std::mt19937& generator, std::uniform_real_distribution<fp_type>& distribution) {
      if constexpr (is_debug())
        return fp_type{0.4f};
      return distribution(generator);
    }

    fp_type random_range(std::mt19937& generator,
                         std::uniform_real_distribution<fp_type>& distribution,
                         const fp_type minimum,
                         const fp_type maximum) {
      return random_unit(generator, distribution) * (maximum - minimum) + minimum;
    }

    int random_population_index(std::mt19937& generator,
                                std::uniform_real_distribution<fp_type>& distribution,
                                const int population_size) {
      return static_cast<int>(
          random_range(generator, distribution, fp_type{0}, static_cast<fp_type>(population_size - 1)));
    }

    void constrain_translation(chromosome& genes, const fp_type max_half_width) {
      if (max_half_width <= fp_type{0})
        return;
      for (int i = 0; i < 3; ++i) genes[i] = std::clamp(genes[i], -max_half_width, max_half_width);
    }

    template<cpu_vectorization vectorization, typename queue_t, scoring_function scoring>
    void genetic_cpu(const genetic_args& genetic, adt_score_args scoring_args) {
      random_generator.init(static_cast<int>(genetic.seed));
      auto& generator = random_generator();
      std::uniform_real_distribution<fp_type> distribution{fp_type{0}, fp_type{1}};
      const auto& geometry             = genetic.geometry;
      chromosome* population           = genetic.population;
      chromosome* next_population      = genetic.next_population;
      const int genes_per_ligand_limit = static_cast<int>(std::tuple_size_v<chromosome>);

      for (int ligand = 0; ligand < genetic.batch_ligands; ++ligand) {
        const int rotamers   = genetic.num_rotamers[ligand];
        const int gene_count = 6 + rotamers;
        if (gene_count > genes_per_ligand_limit)
          throw std::runtime_error("Ligand chromosome exceeds the fixed chromosome capacity");
        chromosome* const ligand_population = population + ligand * genetic.population_size;
        fp_type* const ligand_scores        = genetic.scores + ligand * genetic.population_size;
        for (int individual_index = 0; individual_index < genetic.population_size; ++individual_index) {
          chromosome& genes               = ligand_population[individual_index];
          ligand_scores[individual_index] = big_bound<fp_type>();
          for (int i = 0; i < 3; ++i)
            genes[i] = random_range(generator, distribution, fp_type{-45}, fp_type{45}) * coordinate_step;
          constrain_translation(genes, genetic.max_translation_half_width);
          for (int i = 3; i < gene_count; ++i)
            genes[i] = random_range(generator, distribution, fp_type{-45}, fp_type{45}) * angle_step;
        }
      }

      for (int generation = 0; generation < genetic.generations; ++generation) {
        for (int ligand = 0; ligand < genetic.batch_ligands; ++ligand) {
          const int num_atoms             = geometry.num_atoms[ligand];
          const int num_rotamers          = genetic.num_rotamers[ligand];
          const int atom_offset           = ligand * geometry.batch_atom_stride;
          const fp_type* const original_x = geometry.original_x + atom_offset;
          const fp_type* const original_y = geometry.original_y + atom_offset;
          const fp_type* const original_z = geometry.original_z + atom_offset;
          const int* const fragment_masks = geometry.fragments + geometry.ligand_fragments_start[ligand];
          const int fragment_index        = geometry.frag_indices_start[ligand];

          for (int pose = 0; pose < genetic.population_size; ++pose) {
            const auto pose_offset =
                (static_cast<std::size_t>(ligand) * static_cast<std::size_t>(genetic.population_size) +
                 static_cast<std::size_t>(pose)) *
                static_cast<std::size_t>(geometry.batch_atom_stride);
            auto* const pose_x = geometry.x_scratch + pose_offset;
            auto* const pose_y = geometry.y_scratch + pose_offset;
            auto* const pose_z = geometry.z_scratch + pose_offset;
            std::copy_n(original_x, num_atoms, pose_x);
            std::copy_n(original_y, num_atoms, pose_y);
            std::copy_n(original_z, num_atoms, pose_z);
            const int* const fragment_starts = geometry.fragments_start_index + fragment_index;
            const int* const fragment_stops  = geometry.fragments_stop_index + fragment_index;
            apply<vectorization>(pose_x,
                                 pose_y,
                                 pose_z,
                                 population[ligand * genetic.population_size + pose],
                                 num_atoms,
                                 num_rotamers,
                                 fragment_masks,
                                 fragment_starts,
                                 fragment_stops);
          }
        }
        scoring_args.scores_per_ligand = genetic.population_size;
        scoring_args.batch_atom_stride = geometry.batch_atom_stride;
        scoring_args.scores            = genetic.scores;
        genetic_scoring<queue_t, scoring>::score(scoring_args);

        if (generation + 1 == genetic.generations)
          break;

        for (int ligand = 0; ligand < genetic.batch_ligands; ++ligand) {
          const int rotamers                 = genetic.num_rotamers[ligand];
          const int gene_count               = 6 + rotamers;
          const chromosome* const parents    = population + ligand * genetic.population_size;
          chromosome* const children         = next_population + ligand * genetic.population_size;
          const fp_type* const ligand_scores = genetic.scores + ligand * genetic.population_size;
          for (int child_index = 0; child_index < genetic.population_size; ++child_index) {
            chromosome& child = children[child_index];
            int parent1       = random_population_index(generator, distribution, genetic.population_size);
            int parent2       = random_population_index(generator, distribution, genetic.population_size);
            for (int tournament = 0; tournament < genetic.tournament_length; ++tournament) {
              const int candidate1 =
                  random_population_index(generator, distribution, genetic.population_size);
              const int candidate2 =
                  random_population_index(generator, distribution, genetic.population_size);
              if (ligand_scores[candidate1] < ligand_scores[parent1])
                parent1 = candidate1;
              if (ligand_scores[candidate2] < ligand_scores[parent2])
                parent2 = candidate2;
            }
            const int split = std::min(
                static_cast<int>(
                    random_range(generator, distribution, fp_type{0}, static_cast<fp_type>(gene_count))),
                gene_count);
            std::copy_n(parents[parent1].begin(), split, child.begin());
            std::copy(parents[parent2].begin() + split,
                      parents[parent2].begin() + gene_count,
                      child.begin() + split);

            for (int gene = 0; gene < 3; ++gene) {
              if (random_unit(generator, distribution) < genetic.mutation_probability)
                child[gene] +=
                    random_range(generator, distribution, fp_type{-10}, fp_type{10}) * coordinate_step;
            }
            constrain_translation(child, genetic.max_translation_half_width);
            for (int gene = 3; gene < gene_count; ++gene) {
              if (random_unit(generator, distribution) < genetic.mutation_probability)
                child[gene] += random_range(generator, distribution, fp_type{-10}, fp_type{10}) * angle_step;
            }
          }
        }
        std::swap(population, next_population);
      }

      for (int ligand = 0; ligand < genetic.batch_ligands; ++ligand) {
        const int rotamers                        = genetic.num_rotamers[ligand];
        const chromosome* const ligand_population = population + ligand * genetic.population_size;
        const fp_type* const ligand_scores        = genetic.scores + ligand * genetic.population_size;
        int best                                  = 0;
        for (int index = 1; index < genetic.population_size; ++index) {
          if (ligand_scores[index] < ligand_scores[best])
            best = index;
        }
        genetic.best_scores[ligand] = ligand_scores[best];
        std::copy_n(ligand_population[best].begin(), 6 + rotamers, genetic.best_chromosomes[ligand].begin());
      }
    }
  } // namespace

  template<typename queue_t>
  template<scoring_function scoring>
  void genetic_backend<queue_t>::run(
      std::shared_ptr<queue_t> q,
      const genetic_args& genetic,
      geometric<queue_t>&,
      typename scoring_backend<queue_t, scoring>::stage_type& scoring_stage_instance) {
    const auto scoring_args = scoring_stage_instance.kernel_args();
    q->template invoke_kernel<genetic_compute_region_name>(
        [&]() { genetic_cpu<cpu_vectorization_for<queue_t>, queue_t, scoring>(genetic, scoring_args); });
  }

  [[maybe_unused]] extern const auto registered_genetic_cpp_backends =
      make_genetic_backend_references<genetic_backend<queue_cpp>>(genetic_pipeline_registry{});
#ifdef MUDOCK_USE_GH
  [[maybe_unused]] extern const auto registered_genetic_gh_backends =
      make_genetic_backend_references<genetic_backend<queue_gh>>(genetic_pipeline_registry{});
#endif
#ifdef MUDOCK_USE_XSIMD
  [[maybe_unused]] extern const auto registered_genetic_xsimd_backends =
      make_genetic_backend_references<genetic_backend<queue_xsimd>>(genetic_pipeline_registry{});
#endif
} // namespace mudock
