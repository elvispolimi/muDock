#pragma once

#include <concepts>
#include <memory>
#include <mudock/compute/queue.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>


namespace mudock {
  
  template<typename queue_type>
    requires std::derived_from<queue_type, queue>
  struct adadelta_kernel {
    static constexpr char adadelta_region_name[] = "adadelta_kernel";
    adadelta_kernel(const int individuals_per_ligand_,
                    const int batch_ligands_,
                    const int batch_atoms_,
                    gradient *__restrict__ gradients_b_,
                    chromosome *__restrict__ chromosomes_b_,
                    int *__restrict__ num_rotamers_b_,
                    chromosome *__restrict__ adadelta_e_g2_b_,
                    chromosome *__restrict__ adadelta_e_dw2_b_,
                    int *__restrict__ active_individuals_b_,
                    int *__restrict__ converged_ligands_b_,
                    std::shared_ptr<queue_type> q_,
                    const fp_type rho_,
                    const fp_type epsilon_)
        : individuals_per_ligand(individuals_per_ligand_),
          batch_ligands(batch_ligands_),
          batch_atoms(batch_atoms_),
          gradients_b(gradients_b_),
          chromosomes_b(chromosomes_b_),
          num_rotamers_b(num_rotamers_b_),
          adadelta_e_g2_b(adadelta_e_g2_b_),
          adadelta_e_dw2_b(adadelta_e_dw2_b_),
          active_individuals_b(active_individuals_b_),
          converged_ligands_b(converged_ligands_b_),
          q(q_),
          rho(rho_),
          epsilon(epsilon_) {}

    void operator()();
    inline void set_chromosomes_buffer(chromosome* chromosomes_b_) { chromosomes_b = chromosomes_b_; }

    adadelta_kernel(const adadelta_kernel &)            = default;
    adadelta_kernel(adadelta_kernel &&)                 = default;
    adadelta_kernel &operator=(const adadelta_kernel &) = delete;
    adadelta_kernel &operator=(adadelta_kernel &&)      = delete;

    ~adadelta_kernel() = default;

  private:
    const int individuals_per_ligand;
    const int batch_ligands;
    const int batch_atoms;
    gradient *__restrict__ gradients_b;
    chromosome* chromosomes_b;
    int *__restrict__ num_rotamers_b;
    chromosome *__restrict__ adadelta_e_g2_b;
    chromosome *__restrict__ adadelta_e_dw2_b;
    int *__restrict__ active_individuals_b;
    int *__restrict__ converged_ligands_b;
    std::shared_ptr<queue_type> q;
    const fp_type rho;
    const fp_type epsilon;
  };

} // namespace mudock
