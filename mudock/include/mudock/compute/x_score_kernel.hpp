#pragma once

#include <concepts>
#include <memory>
#include <mudock/compute/queue.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {
  template<typename queue_type>
    requires std::derived_from<queue_type, queue>
  struct x_score_kernel {
    static constexpr char x_region_name[] = "x_score_kernel";

    x_score_kernel(const int batch_ligands_,
                   const int batch_atoms_,
                   const int *__restrict__ num_atoms_b_,
                   const fp_type *__restrict__ lig_x_b_,
                   const fp_type *__restrict__ lig_y_b_,
                   const fp_type *__restrict__ lig_z_b_,
                   const fp_type *__restrict__ lig_vdw_b_,
                   const int *__restrict__ lig_scorable_b_,
                   const int *__restrict__ lig_hb_b_,
                   const fp_type *__restrict__ lig_rt_b_,
                   const fp_type *__restrict__ lig_hbt_b_,
                   const int num_prot_atoms_,
                   const fp_type *__restrict__ prot_x_b_,
                   const fp_type *__restrict__ prot_y_b_,
                   const fp_type *__restrict__ prot_z_b_,
                   const fp_type *__restrict__ prot_vdw_b_,
                   const int *__restrict__ prot_scorable_b_,
                   const int *__restrict__ prot_hb_b_,
                   fp_type *__restrict__ terms_b_,
                   std::shared_ptr<queue_type> q_)
        : batch_ligands(batch_ligands_),
          batch_atoms(batch_atoms_),
          num_atoms_b(num_atoms_b_),
          lig_x_b(lig_x_b_),
          lig_y_b(lig_y_b_),
          lig_z_b(lig_z_b_),
          lig_vdw_b(lig_vdw_b_),
          lig_scorable_b(lig_scorable_b_),
          lig_hb_b(lig_hb_b_),
          lig_rt_b(lig_rt_b_),
          lig_hbt_b(lig_hbt_b_),
          num_prot_atoms(num_prot_atoms_),
          prot_x_b(prot_x_b_),
          prot_y_b(prot_y_b_),
          prot_z_b(prot_z_b_),
          prot_vdw_b(prot_vdw_b_),
          prot_scorable_b(prot_scorable_b_),
          prot_hb_b(prot_hb_b_),
          terms_b(terms_b_),
          q(q_) {}

    void operator()();

    x_score_kernel(const x_score_kernel &)            = default;
    x_score_kernel(x_score_kernel &&)                 = default;
    x_score_kernel &operator=(const x_score_kernel &) = delete;
    x_score_kernel &operator=(x_score_kernel &&)      = delete;

    ~x_score_kernel() = default;

  private:
    const int batch_ligands;
    const int batch_atoms;
    const int *__restrict__ num_atoms_b;

    // ligand atom data (strided per ligand by batch_atoms)
    const fp_type *__restrict__ lig_x_b;
    const fp_type *__restrict__ lig_y_b;
    const fp_type *__restrict__ lig_z_b;
    const fp_type *__restrict__ lig_vdw_b;
    const int *__restrict__ lig_scorable_b;
    const int *__restrict__ lig_hb_b;

    // per-ligand RT and HB Terms (indexed by ligand and not strided by batch_atoms)
    const fp_type *__restrict__ lig_rt_b;
    const fp_type *__restrict__ lig_hbt_b;

    // protein atom data, used for the entire ligand batch
    const int num_prot_atoms;
    const fp_type *__restrict__ prot_x_b;
    const fp_type *__restrict__ prot_y_b;
    const fp_type *__restrict__ prot_z_b;
    const fp_type *__restrict__ prot_vdw_b;
    const int *__restrict__ prot_scorable_b;
    const int *__restrict__ prot_hb_b;

    // per-ligand XScore terms
    fp_type *__restrict__ terms_b;
    std::shared_ptr<queue_type> q;
  };

} // namespace mudock