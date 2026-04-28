#pragma once

#include <concepts>
#include <memory>
#include <mudock/compute/queue.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {
  
  template<typename queue_type>
    requires std::derived_from<queue_type, queue>
  struct adadelta_kernel {
    static constexpr char adadelta_region_name[] = "adadelta_kernel";
    static constexpr char gradient_region_name[] = "adadelta_gradient_kernel";
    adadelta_kernel(const int individuals_per_ligand,
                    const int batch_ligands_,
                    gradient *__restrict__ gradients_b_,
                    chromosome *__restrict__ population_b_,
                    std::shared_ptr<queue_type> q_)
        : scores_per_ligand(individuals_per_ligand),
          batch_ligands(batch_ligands_),
          gradients_b(gradients_b_),
          population_b(population_b_),
          q(q_) {}

    void operator()();
    void compute_gradients();
    void apply_adadelta();

    adadelta_kernel(const adadelta_kernel &)            = default;
    adadelta_kernel(adadelta_kernel &&)                 = default;
    adadelta_kernel &operator=(const adadelta_kernel &) = delete;
    adadelta_kernel &operator=(adadelta_kernel &&)      = delete;

    ~adadelta_kernel() = default;

  private:
    const int batch_ligands;
    const int individuals_per_ligand;
    gradient *__restrict__ gradients_b;
    chromosome* population_b;
    std::shared_ptr<queue_type> q;
  };

} // namespace mudock
