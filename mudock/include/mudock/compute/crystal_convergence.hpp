#pragma once

#include "mudock/type_alias.hpp"

#include <concepts>
#include <memory>
#include <mudock/chem/geom_ligand.hpp>
#include <mudock/compute/batch_multiple.hpp>
#include <mudock/compute/queue.hpp>
#include <mudock/log.hpp>
#if !defined(__CUDACC__) && !defined(__HIPCC__)
  #include <mudock/compute/buffer_utils.hpp>
  #include <mudock/compute/scratchpad.hpp>
  #include <mudock/compute/convergence.hpp>
#endif
#include <mudock/cpp_implementation/center_of_mass.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/cpp_implementation/mutate_cpp.hpp>
#include <mudock/grid/point3D.hpp>
#include <mudock/molecule.hpp>

namespace mudock {
  template<typename queue_type>
  batch_multiple get_crystal_conv_batch_multiple(const int, std::shared_ptr<queue_type>) {
    return {};
  }

  template<typename queue_type>
    requires std::derived_from<queue_type, queue>
  struct crystal_convergence_kernel {
    static constexpr char crystal_convergence_region_name[] = "crystal_convergence";
    crystal_convergence_kernel(const int batch_ligands_,
                const int batch_atoms_,
                const int population_number_,
                // fp_type* __restrict__ x_scratch_b_,
                // fp_type* __restrict__ y_scratch_b_,
                // fp_type* __restrict__ z_scratch_b_,
                fp_type *__restrict__ scores_b_,
                std::shared_ptr<queue_type> q_)
        : batch_ligands(batch_ligands_),
          batch_atoms(batch_atoms_),
          population_number(population_number_),
          // x_scratch_b(x_scratch_b_),
          // y_scratch_b(y_scratch_b_),
          // z_scratch_b(z_scratch_b_),
          scores_b(scores_b_),
          q(q_) {}

    void operator()();

    crystal_convergence_kernel(const crystal_convergence_kernel&)            = default;
    crystal_convergence_kernel(crystal_convergence_kernel&&)                 = default;
    crystal_convergence_kernel& operator=(const crystal_convergence_kernel&) = delete;
    crystal_convergence_kernel& operator=(crystal_convergence_kernel&&)      = delete;

    ~crystal_convergence_kernel() = default;

  private:
  // scores_b for best score to check crystal
  // x/y/z_scratch_b and ligand_template for rmsd (num_atoms can be retrieved from ligand_template)
    const int batch_ligands;
    const int batch_atoms;
    const int population_number;
    // fp_type* __restrict__ x_scratch_b;
    // fp_type* __restrict__ y_scratch_b;
    // fp_type* __restrict__ z_scratch_b;
    fp_type *__restrict__ scores_b;
    std::shared_ptr<queue_type> q;
  };

#if !defined(__CUDACC__) && !defined(__HIPCC__)
  template<typename queue_t>
    requires std::derived_from<queue_t, queue>
  struct crystal_convergence: public convergence<queue_t> {
    crystal_convergence(std::shared_ptr<scratchpad<queue_t>> _scratch)
        : convergence<queue_t>(_scratch) {};

    void prepare(batch<static_molecule>& batch) {
      batch_ligands                           = batch.num_ligands;
      batch_atoms                             = batch.batch_max_atoms;
      const int tot_atoms_in_batch            = batch_ligands * batch_atoms;
      auto q                                  = (*this->scratch).get_queue();
      const auto population_number            = std::max(1, static_cast<int>((*this->scratch).configuration.population_number));

      load_num_rotamers<queue_t>(batch, this->scratch);
      load_num_atoms<queue_t>(batch, this->scratch);
      // auto& x_scratch_b = (*this->scratch).template get<buffer_data_type::X_SCRATCH>();
      // auto& y_scratch_b = (*this->scratch).template get<buffer_data_type::Y_SCRATCH>();
      // auto& z_scratch_b = (*this->scratch).template get<buffer_data_type::Z_SCRATCH>();
      auto& scores_b    = (*this->scratch).template get<buffer_data_type::SCORES>();
      
      scores_b.alloc(population_number * batch_ligands);

      // Binding
      int* num_atoms_p    = (*this->scratch).template get<buffer_data_type::NUM_ATOMS>().dev_pointer();
      int* num_rotamers_p = (*this->scratch).template get<buffer_data_type::NUM_ROTAMERS>().dev_pointer();

      // fp_type* x_scratch_p           = x_scratch_b.dev_pointer();
      // fp_type* y_scratch_p           = y_scratch_b.dev_pointer();
      // fp_type* z_scratch_p           = z_scratch_b.dev_pointer();
      fp_type* __restrict__ scores_p = scores_b.dev_pointer();


      kernel = std::make_unique<crystal_convergence_kernel<queue_t>>(batch_ligands,
                                                                     batch_atoms,
                                                                     population_number,
                                                                    //  x_scratch_p,
                                                                    //  y_scratch_p,
                                                                    //  z_scratch_p,
                                                                     scores_p,
                                                                     q);
    }

    void operator()() {
      // assert(
      //     ((*this->scratch).template get<buffer_data_type::CHROMOSOMES>().num_elements() % batch_ligands ==
      //      0) &&
      //     "Number of chromosomes per ligand is not a multiple of the number of ligands expected to be docked");

      assert(kernel && "Kernel method not yet prepared");
      (*kernel)();
    };

    // TODO fix function
    static std::size_t get_shared_ligand_mem(const int max_atoms, const knobs conf) {
      (void) max_atoms;
      (void) conf;
      return 0;
    }

    // TODO fix function
    static std::size_t get_private_ligand_mem(const int max_atoms, const knobs) {
      const int batch_rotamers              = max_atoms - 3;
      const int tot_rotamers_atoms_in_batch = max_atoms * batch_rotamers;
      std::size_t mem{0};
      return mem;
    }

    static int get_ligand_mem(const int max_atoms, const knobs conf) {
      return static_cast<int>(get_shared_ligand_mem(max_atoms, conf) +
                              get_private_ligand_mem(max_atoms, conf));
    }

    static batch_multiple get_batch_size(const int atoms,
                                         std::shared_ptr<queue_t> q,
                                         const knobs& conf,
                                         const size_t max_bucket_size) {
      (void) conf;
      const auto plain_multiple_info =
          normalize_batch_multiple(get_crystal_conv_batch_multiple<queue_t>(atoms, q));
      mudock::stage_bucket_trace("CONVERGENCE stage plain multiple for ",
                                 atoms,
                                 " atoms -> total=",
                                 plain_multiple_info.total_multiple(),
                                 " (active_blocks_per_sm=",
                                 plain_multiple_info.active_blocks_per_sm,
                                 ", num_sms=",
                                 plain_multiple_info.num_sms,
                                 ")",
                                 " (max_bucket_size hint=",
                                 max_bucket_size,
                                 ")");
      return plain_multiple_info;
    }

  private:
    int batch_ligands;
    int batch_atoms;
    std::unique_ptr<crystal_convergence_kernel<queue_t>> kernel;

    // TODO fix function
    void teardown_impl(batch<static_molecule>& batch) {  
    };
  };
#endif
} // namespace mudock
