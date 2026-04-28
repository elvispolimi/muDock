#pragma once

#include <cstddef>
#include <cstring>
#include <mudock/batch.hpp>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/compute/adadelta_kernel.hpp>
#ifndef __CUDACC__
  #include <mudock/compute/buffer_utils.hpp>
  #include <mudock/compute/scoring.hpp>
  #include <mudock/compute/scratchpad.hpp>
#endif
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  #define MAX_ADADELTA_ITERATIONS 300
  #define ADADELTA_RHO 0.95f
  #define ADADELTA_EPSILON 1e-6f

#ifndef __CUDACC__
  // TODO check that the object type and the kernel impl are the same
  template<typename queue_type>
  struct adadelta: public local_search<queue_type> {
    adadelta(std::shared_ptr<scratchpad<queue_type>> _scratch)
        : local_search<queue_type>(_scratch)/*,
        score_stage(_scratch)*/ {}

    void prepare(batch<static_molecule> &batch) {
      // Allocate gradient buffer for AdaDelta (one gradient per individual per ligand)
      batch_ligands = batch.num_ligands;
      const int individuals_per_ligand = std::max(1, static_cast<int>((*this->scratch).configuration.population_number));
      auto &gradient_b = (*this->scratch).template get<buffer_data_type::GRADIENTS>();
      const size_t gradient_count = static_cast<size_t>(batch_ligands) * static_cast<size_t>(individuals_per_ligand);
      if (!gradient_b.is_valid() || gradient_b.num_elements() != gradient_count) {
        gradient_b.alloc(gradient_count);
        gradient_b.set_valid();
      }
      gradient *gradients_b = gradient_b.dev_pointer();

      auto &chromosomes_b = (*this->scratch).template get<buffer_data_type::CHROMOSOMES>();
      chromosome *population_b = chromosomes_b.dev_pointer();

      ls_ad_kernel = std::make_unique<adadelta_kernel<queue_type>>(individuals_per_ligand,
                                                                  batch_ligands,
                                                                  gradients_b,
                                                                  population_b,
                                                                  q);
      
    }

    void operator()() {
      assert(
          (((*this->scratch).template get<buffer_data_type::GRADIENTS>().num_elements() % batch_ligands) == 0) &&
          "Number of gradients is not a multiple of ligands in the batch");
      assert(ls_ad_kernel && "Adadelta local search kernel method not yet prepared");
      for (int i = 0; i < MAX_ADADELTA_ITERATIONS; ++i){
        ls_ad_kernel->compute_gradients();
        ls_ad_kernel->apply_adadelta();
      }
    }

    // static int get_ligand_mem(const int max_atoms, const knobs conf) {
    //   // TODO L what should I do here?
    // }

  private:
    int batch_ligands;

    std::unique_ptr<adadelta_kernel<queue_type>> ls_ad_kernel;

    void teardown_impl(batch<static_molecule> &batch) override {
      // TODO L what to do here???
    };
  };
#endif
} // namespace mudock
