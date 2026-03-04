#pragma once

#include <concepts>
#include <memory>
#include <mudock/compute/queue.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {
  // TODO check maybe the kernel can be fused togheter with main adt score
  // May become an issue to keep separate the TU and the CUDA/etc dependencies
  template<typename queue_type>
    requires std::derived_from<queue_type, queue>
    struct vina_score_kernel {
      static constexpr char vina_region_name[] = "vina_score_kernel";
      vina_score_kernel(const int scores_per_ligand_,
          const int batch_ligands_,
          const int batch_atoms_,
          //Protein data
          const size_t num_atoms_protein_,
          const fp_type* __restrict__ protein_x_,
          const fp_type* __restrict__ protein_y_,
          const fp_type* __restrict__ protein_z_,
          const int* __restrict__ p_is_hbond_acceptor_,
          const int* __restrict__ p_is_hbond_donor_,
          const int* __restrict__ p_is_hydrophobic_,
          const fp_type* __restrict__ p_vdw_radius_,
          // Ligand data
          const int *__restrict__ num_ligand_atoms_b_, 
          const fp_type *__restrict__ x_scratch_b_,
          const fp_type *__restrict__ y_scratch_b_,
          const fp_type *__restrict__ z_scratch_b_,
          const int *__restrict__ l_is_hbond_acceptor_b_,
          const int *__restrict__ l_is_hbond_donor_b_,
          const int *__restrict__ l_is_hydrophobic_b_,
          const fp_type *__restrict__ l_vdw_radius_b_,
          const int *__restrict__ active_torsions_b_,
          const int *__restrict__ interacting_pairs_first_b_, 
          const int *__restrict__ interacting_pairs_second_b_,
          const int *__restrict__ num_interacting_pairs_b_, 
          fp_type *__restrict__ scores_b_,
          std::shared_ptr<queue_type> q_)
            : scores_per_ligand(scores_per_ligand_),
            batch_ligands(batch_ligands_),
            batch_atoms(batch_atoms_),
            num_atoms_protein(num_atoms_protein_),
            protein_x(protein_x_),
            protein_y(protein_y_),
            protein_z(protein_z_),
            p_is_hbond_acceptor(p_is_hbond_acceptor_),
            p_is_hbond_donor(p_is_hbond_donor_),
            p_is_hydrophobic(p_is_hydrophobic_),
            p_vdw_radius(p_vdw_radius_),
            num_ligand_atoms_b(num_ligand_atoms_b_),
            x_scratch_b(x_scratch_b_),
            y_scratch_b(y_scratch_b_),
            z_scratch_b(z_scratch_b_),
            l_is_hbond_acceptor_b(l_is_hbond_acceptor_b_),
            l_is_hbond_donor_b(l_is_hbond_donor_b_),
            l_is_hydrophobic_b(l_is_hydrophobic_b_),
            l_vdw_radius_b(l_vdw_radius_b_),
            active_torsions_b(active_torsions_b_),
            interacting_pairs_first_b(interacting_pairs_first_b_),
            interacting_pairs_second_b(interacting_pairs_second_b_),
            num_interacting_pairs_b(num_interacting_pairs_b_),
            scores_b(scores_b_),
            q(q_) {}

      void operator()();

      vina_score_kernel(const vina_score_kernel &)            = default;
      vina_score_kernel(vina_score_kernel &&)                 = default;
      vina_score_kernel &operator=(const vina_score_kernel &) = delete;
      vina_score_kernel &operator=(vina_score_kernel &&)      = delete;

      ~vina_score_kernel() = default;

      private:
      const int scores_per_ligand;
      const int batch_ligands;
      const int batch_atoms;   
      //Protein data
      const size_t num_atoms_protein;
      const fp_type* __restrict__ protein_x;
      const fp_type* __restrict__ protein_y;
      const fp_type* __restrict__ protein_z;
      const int* __restrict__ p_is_hbond_acceptor;
      const int* __restrict__ p_is_hbond_donor;
      const int* __restrict__ p_is_hydrophobic;
      const fp_type* __restrict__ p_vdw_radius;
      // Ligand data
      const int *__restrict__ num_ligand_atoms_b; 
      const fp_type *__restrict__ x_scratch_b;
      const fp_type *__restrict__ y_scratch_b;
      const fp_type *__restrict__ z_scratch_b;
      const int* __restrict__ l_is_hbond_acceptor_b;
      const int* __restrict__ l_is_hbond_donor_b;
      const int* __restrict__ l_is_hydrophobic_b;
      const fp_type* __restrict__ l_vdw_radius_b;
      const int *__restrict__ active_torsions_b;
      const int* __restrict__ interacting_pairs_first_b; 
      const int* __restrict__ interacting_pairs_second_b;
      const int* __restrict__ num_interacting_pairs_b; 

      fp_type *__restrict__ scores_b;
      std::shared_ptr<queue_type> q;
    };

}
