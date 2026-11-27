#include <mudock/compute/genetic.hpp>
#include <mudock/compute/geometric_transform.hpp>
#include <mudock/compute/pipeline_registry.hpp>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/sycl_implementation/invoke_kernel_sycl.hpp>
#include <mudock/sycl_implementation/mutate.hpp>
#include <mudock/sycl_implementation/scoring_sycl.hpp>
#include <mudock/sycl_implementation/sycl_random.hpp>
#include <mudock/utils.hpp>
#include <sycl/sycl.hpp>

#ifndef MUDOCK_SYCL_WG_SIZE
  #define MUDOCK_SYCL_WG_SIZE 32
#endif

namespace mudock {

  static constexpr fp_type coordinate_step = static_cast<fp_type>(0.2);
  static constexpr fp_type angle_step      = static_cast<fp_type>(4);

  // TODO check the real randomness
  template<typename T>
  const T random_gen_sycl(XORWOWState& state, const T min, const T max) {
    fp_type value;
    if constexpr (is_debug())
      // TODO value here for debug
      value = static_cast<fp_type>(0.4);
    else {
      value = syclrand_uniform(&state);
    }
    return static_cast<T>((value * static_cast<fp_type>(max - min)) + static_cast<fp_type>(min));
  }

  inline int get_selection_distribution(XORWOWState& state, const int* population_number) {
    return random_gen_sycl<int>(state, 0, *population_number - 1);
  };

  inline fp_type get_init_change_distribution(XORWOWState& state) {
    return random_gen_sycl<fp_type>(state, -45, 45);
  }
  inline fp_type get_mutation_change_distribution(XORWOWState& state) {
    return random_gen_sycl<fp_type>(state, -10, 10);
  };
  inline fp_type get_mutation_coin_distribution(XORWOWState& state) {
    return random_gen_sycl<fp_type>(state, 0, 1);
  };

  template<bool bounded>
  inline void limit_translation(chromosome& chrom, const fp_type half_width) {
    if constexpr (!bounded)
      return;
    chrom[0] = sycl::clamp(chrom[0], -half_width, half_width);
    chrom[1] = sycl::clamp(chrom[1], -half_width, half_width);
    chrom[2] = sycl::clamp(chrom[2], -half_width, half_width);
  }
  inline int get_crossover_distribution(XORWOWState& state, const int* num_rotamers) {
    return random_gen_sycl<int>(state, 0, 6 + *num_rotamers);
  };

  inline int tournament_selection_sycl(XORWOWState& state,
                                       const int tournament_length,
                                       const int chromosome_number,
                                       const fp_type* __restrict__ scores) {
    const int num_iterations = tournament_length;
    int best_individual      = get_selection_distribution(state, &chromosome_number);
    for (int i = 0; i < num_iterations; ++i) {
      auto contended = get_selection_distribution(state, &chromosome_number);
      if (scores[contended] < scores[best_individual]) {
        best_individual = contended;
      }
    }
    return best_individual;
  }

  template<int MAX_ATOMS, bool bounded, class scoring_backend_t>
  struct genetic_gpu {
    using scoring_kernel_t = typename scoring_backend_t::device_kernel;

    void operator()(sycl::nd_item<3> it,
                    genetic_args genetic,
                    typename scoring_backend_t::arguments scoring) const {
      const int ligand_id        = static_cast<int>(it.get_group(2));
      const int lane             = static_cast<int>(it.get_local_id(2));
      const int thread_per_block = static_cast<int>(it.get_local_range(2));
      const int state_index      = ligand_id * thread_per_block + lane;
      XORWOWState local_state;
      syclrand_init(static_cast<unsigned long long>(genetic.seed),
                    static_cast<unsigned long long>(state_index),
                    &local_state);
      const fp_type max_translation_half_width = genetic.max_translation_half_width;
      const int num_rotamers                   = genetic.num_rotamers[ligand_id];
      const int num_genes                      = 6 + num_rotamers;
      chromosome* population                   = genetic.population + ligand_id * genetic.population_size;
      chromosome* next_population = genetic.next_population + ligand_id * genetic.population_size;
      fp_type* scores             = genetic.scores + ligand_id * genetic.population_size;

      for (int index = lane; index < genetic.population_size; index += thread_per_block) {
        scores[index]     = std::numeric_limits<fp_type>::infinity();
        chromosome& genes = population[index];
        MUDOCK_PRAGMA_UNROLL(MUDOCK_UNROLL_FACTOR)
        for (int gene = 0; gene < 3; ++gene)
          genes[gene] = get_init_change_distribution(local_state) * coordinate_step;
        limit_translation<bounded>(genes, max_translation_half_width);
        for (int gene = 3; gene < num_genes; ++gene)
          genes[gene] = get_init_change_distribution(local_state) * angle_step;
      }
      sycl::group_barrier(it.get_group());

      for (int generation = 0; generation < genetic.generations; ++generation) {
        const auto& geometry = genetic.geometry;
        for (int pose = 0; pose < genetic.population_size; ++pose) {
          apply_sycl<MAX_ATOMS>{}(it,
                                  genetic.population_size,
                                  geometry.batch_atom_stride,
                                  geometry.original_x,
                                  geometry.original_y,
                                  geometry.original_z,
                                  geometry.x_scratch,
                                  geometry.y_scratch,
                                  geometry.z_scratch,
                                  population - ligand_id * genetic.population_size,
                                  geometry.fragments,
                                  geometry.ligand_fragments_start,
                                  geometry.fragments_start_index,
                                  geometry.fragments_stop_index,
                                  geometry.frag_indices_start,
                                  geometry.num_rotamers,
                                  geometry.num_atoms,
                                  pose);
          sycl::group_barrier(it.get_group());
        }

        for (int pose = 0; pose < genetic.population_size; ++pose) {
          const int coordinate_offset =
              (ligand_id * genetic.population_size + pose) * geometry.batch_atom_stride;
          auto ligand_scoring              = scoring;
          ligand_scoring.ligand_index      = ligand_id;
          ligand_scoring.x                 = geometry.x_scratch + coordinate_offset;
          ligand_scoring.y                 = geometry.y_scratch + coordinate_offset;
          ligand_scoring.z                 = geometry.z_scratch + coordinate_offset;
          ligand_scoring.batch_atom_stride = geometry.batch_atom_stride;
          const fp_type score              = scoring_kernel_t::template score<MAX_ATOMS>(it, ligand_scoring);
          if (lane == 0)
            scores[pose] = score;
        }
        sycl::group_barrier(it.get_group());

        if (generation + 1 == genetic.generations)
          break;

        for (int index = lane; index < genetic.population_size; index += thread_per_block) {
          chromosome& child = next_population[index];
          const int parent1 = tournament_selection_sycl(local_state,
                                                        genetic.tournament_length,
                                                        genetic.population_size,
                                                        scores);
          const int parent2 = tournament_selection_sycl(local_state,
                                                        genetic.tournament_length,
                                                        genetic.population_size,
                                                        scores);
          const int split   = get_crossover_distribution(local_state, &num_rotamers);
          const fp_type* p1 = population[parent1].data();
          const fp_type* p2 = population[parent2].data();
          for (int gene = 0; gene < num_genes; ++gene) child[gene] = gene < split ? p1[gene] : p2[gene];
          MUDOCK_PRAGMA_UNROLL(MUDOCK_UNROLL_FACTOR)
          for (int gene = 0; gene < 3; ++gene) {
            if (get_mutation_coin_distribution(local_state) < genetic.mutation_probability)
              child[gene] += get_mutation_change_distribution(local_state) * coordinate_step;
          }
          limit_translation<bounded>(child, max_translation_half_width);
          for (int gene = 3; gene < num_genes; ++gene) {
            if (get_mutation_coin_distribution(local_state) < genetic.mutation_probability)
              child[gene] += get_mutation_change_distribution(local_state) * angle_step;
          }
        }
        sycl::group_barrier(it.get_group());
        chromosome* swap = population;
        population       = next_population;
        next_population  = swap;
      }

      int min_index = lane;
      fp_type min_score =
          min_index < genetic.population_size ? scores[min_index] : std::numeric_limits<fp_type>::infinity();
      for (int index = lane + thread_per_block; index < genetic.population_size; index += thread_per_block) {
        if (min_score > scores[index]) {
          min_index = index;
          min_score = scores[index];
        }
      }
      const fp_type best_score = sycl::reduce_over_group(it.get_group(), min_score, sycl::minimum());
      if (min_score == best_score) {
        genetic.best_scores[ligand_id] = min_score;
        const fp_type* source          = population[min_index].data();
        fp_type* target                = genetic.best_chromosomes[ligand_id].data();
        for (int gene = 0; gene < num_genes; ++gene) target[gene] = source[gene];
      }
    }
  };

  template<>
  template<scoring_function scoring>
  void genetic_backend<queue_sycl>::run(
      std::shared_ptr<queue_sycl> q,
      const genetic_args& genetic,
      geometric<queue_sycl>&,
      typename scoring_backend<queue_sycl, scoring>::stage_type& scoring_stage) {
    const auto scoring_args                  = scoring_stage.kernel_args();
    const fp_type max_translation_half_width = genetic.max_translation_half_width;
    constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>(
        [&](const auto atom_index) {
          constexpr int max_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
          using selected_scoring  = scoring_backend<queue_sycl, scoring>;
          if (max_translation_half_width > fp_type{0}) {
            q->invoke_kernel<genetic_gpu<max_atoms, true, selected_scoring>>(genetic.batch_ligands,
                                                                             MUDOCK_SYCL_WG_SIZE,
                                                                             genetic,
                                                                             scoring_args);
          } else {
            q->invoke_kernel<genetic_gpu<max_atoms, false, selected_scoring>>(genetic.batch_ligands,
                                                                              MUDOCK_SYCL_WG_SIZE,
                                                                              genetic,
                                                                              scoring_args);
          }
        },
        genetic.geometry.batch_atom_stride,
        reorder_buffer<static_molecule>::atoms_clusters.data());
  }

  [[maybe_unused]] extern const auto registered_genetic_sycl_backends =
      make_genetic_backend_references<genetic_backend<queue_sycl>>(genetic_pipeline_registry{});

} // namespace mudock
