#pragma once

#include <concepts>
#include <memory>
#include <mudock/compute/queue.hpp>
#include <mudock/compute/scoring.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {
  
  template<typename queue_type>
    requires std::derived_from<queue_type, queue>
  struct adadelta_kernel {
    static constexpr char adadelta_region_name[] = "adadelta_kernel";
    static constexpr char gradient_region_name[] = "adadelta_gradient_kernel";
    adadelta_kernel(const int individuals_per_ligand_,
                    const int batch_ligands_,
                    std::shared_ptr<differentiable_scoring<queue_type>> score_stage_,
                    gradient *__restrict__ gradients_b_,
                    chromosome *__restrict__ population_b_,
                    int *__restrict__ num_rotamers_b_,
                    chromosome *__restrict__ adadelta_e_g2_b_,
                    chromosome *__restrict__ adadelta_e_dw2_b_,
                    int *__restrict__ stall_counter_b_,
                    int *__restrict__ inactive_b_,
                    std::shared_ptr<queue_type> q_)
        : individuals_per_ligand(individuals_per_ligand_),
          batch_ligands(batch_ligands_),
          score_stage(score_stage_),
          gradients_b(gradients_b_),
          population_b(population_b_),
          num_rotamers_b(num_rotamers_b_),
          adadelta_e_g2_b(adadelta_e_g2_b_),
          adadelta_e_dw2_b(adadelta_e_dw2_b_),
          stall_counter_b(stall_counter_b_),
          inactive_b(inactive_b_),
          q(q_) {}

    void operator()();
    void compute_gradients();
    void apply_adadelta(int i);

    adadelta_kernel(const adadelta_kernel &)            = default;
    adadelta_kernel(adadelta_kernel &&)                 = default;
    adadelta_kernel &operator=(const adadelta_kernel &) = delete;
    adadelta_kernel &operator=(adadelta_kernel &&)      = delete;

    ~adadelta_kernel() = default;

  private:
    const int individuals_per_ligand;
    const int batch_ligands;
    std::shared_ptr<differentiable_scoring<queue_type>> score_stage;
    gradient *__restrict__ gradients_b;
    chromosome* population_b;
    int *__restrict__ num_rotamers_b;
    chromosome *__restrict__ adadelta_e_g2_b;
    chromosome *__restrict__ adadelta_e_dw2_b;
    int *__restrict__ stall_counter_b;
    int *__restrict__ inactive_b;
    std::shared_ptr<queue_type> q;
  };

} // namespace mudock
