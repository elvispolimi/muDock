#include <cmath>

#include <mudock/format.hpp>
#include <mudock/molecule.hpp>
#include <mudock/cpp_implementation/vina_score_cpp.hpp>

#define GAUSS1_COEFF        (- 0.035579f)
#define GAUSS2_COEFF        (- 0.005156f)
#define REPULSION_COEFF     (0.840245f)
#define HYDROPHOBIC_COEFF   (- 0.035069f)
#define H_BOND_COEFF        (- 0.587439f)
#define NROT_COEFF          (0.05846f)

#define IS_DIFF_FROM_ZERO(dst) (std::fabs(dst) > 1e-6f)

namespace mudock { 

  inline fp_type distance(fp_type x, fp_type y, fp_type z) {
    return sqrt( x*x + y*y + z*z );
  }

  inline fp_type gauss1(const fp_type dst) {
    fp_type gauss1 = 0;
    if(IS_DIFF_FROM_ZERO(dst)) gauss1 = exp(- powf(dst / 0.5f, 2));
    return gauss1;
  }

  inline fp_type gauss2(const fp_type dst) {
    fp_type gauss2 = 0;
    if(IS_DIFF_FROM_ZERO(dst)) gauss2 = expf(- powf((dst - 3) / 2, 2));
    return gauss2;
  }

  inline fp_type repulsion(const fp_type dst) {
    return powf((dst < 0) * dst, 2);
  }

  inline fp_type hydrophobic(const fp_type dst, const int rec_lig_is_hydrophobic) {
    // if(rec_lig_is_hydrophobic < 0) return 0; // Sentinel value, no interaction
    bool hydro_1 = rec_lig_is_hydrophobic && (dst <= 0.5f);
    bool hydro_2_cond = rec_lig_is_hydrophobic && (dst > 0.5f) && (dst < 1.5f);
    fp_type hydro_2 = 1.5f * hydro_2_cond - hydro_2_cond * dst;
    return hydro_1 + hydro_2;
  }

  inline fp_type hbonding(const fp_type dst, const int rec_lig_is_hb) {
    // if(rec_lig_is_hb < 0) return 0; // Sentinel value, no interaction
    bool h_bond_1 = rec_lig_is_hb && (dst <= -0.7f);
    bool h_bond_2_cond = rec_lig_is_hb && (dst < 0) && (dst > -0.7f);
    fp_type h_bond_2 = h_bond_2_cond * (- dst) / 0.7f;
    return h_bond_1 + h_bond_2;
  }

  inline fp_type compute_pair_energy(
      fp_type dx, fp_type dy, fp_type dz,
      fp_type vdw1, fp_type vdw2,
      bool is_hba1, bool is_hbd1, bool is_hba2, bool is_hbd2,
      bool is_hydro1, bool is_hydro2
      ) {
    fp_type dst = distance(dx, dy, dz);
    if (dst > 8) return 0;
    
    dst -= (vdw1 + vdw2);
    int is_h = (is_hba1 && is_hbd2) || (is_hba2 && is_hbd1);
    int is_hydro = is_hydro1 && is_hydro2;

    return GAUSS1_COEFF * gauss1(dst) +
      GAUSS2_COEFF * gauss2(dst) +
      REPULSION_COEFF * repulsion(dst) +
      HYDROPHOBIC_COEFF * hydrophobic(dst, is_hydro) +
      H_BOND_COEFF * hbonding(dst, is_h);
  }

  inline fp_type score_inter(
      /// Protein data
      const int num_atoms_protein,
      const fp_type* __restrict__ protein_x,
      const fp_type* __restrict__ protein_y,
      const fp_type* __restrict__ protein_z,
      const int* __restrict__ p_is_hbond_acceptor,
      const int* __restrict__ p_is_hbond_donor,
      const int* __restrict__ p_is_hydrophobic,
      const fp_type* __restrict__ p_vdw_radius,

      /// Ligand data
      const int num_atoms_ligand,
      const fp_type* __restrict__ ligand_x,
      const fp_type* __restrict__ ligand_y,
      const fp_type* __restrict__ ligand_z,
      const int* __restrict__ l_is_hbond_acceptor,
      const int* __restrict__ l_is_hbond_donor,
      const int* __restrict__ l_is_hydrophobic,
      const fp_type* __restrict__ l_vdw_radius
      ) {
    fp_type total = 0;

#pragma omp simd
    for (int pIdx = 0; pIdx < num_atoms_protein; pIdx++) {

      const fp_type px = protein_x[pIdx];
      const fp_type py = protein_y[pIdx];
      const fp_type pz = protein_z[pIdx];
      const fp_type pvdw = p_vdw_radius[pIdx];
      const bool phba = p_is_hbond_acceptor[pIdx];
      const bool phbd = p_is_hbond_donor[pIdx];
      const bool phydro = p_is_hydrophobic[pIdx];

      for (int lIdx = 0; lIdx < num_atoms_ligand; lIdx++) {
        total += compute_pair_energy(
            px - ligand_x[lIdx],
            py - ligand_y[lIdx],
            pz - ligand_z[lIdx],
            pvdw, l_vdw_radius[lIdx],
            phba, phbd,
            l_is_hbond_acceptor[lIdx], l_is_hbond_donor[lIdx],
            phydro, l_is_hydrophobic[lIdx]
            );
      }
    }    
    return total;        
  }

  inline fp_type score_intra(
      const fp_type* __restrict__ ligand_x,
      const fp_type* __restrict__ ligand_y,
      const fp_type* __restrict__ ligand_z,
      const int* __restrict__ l_is_hbond_acceptor,
      const int* __restrict__ l_is_hbond_donor,
      const int* __restrict__ l_is_hydrophobic,
      const fp_type* __restrict__ l_vdw_radius,
      const int* __restrict__ interacting_pairs_first,
      const int* __restrict__ interacting_pairs_second,
      const int num_interacting_pairs
      ) {
    fp_type total = 0;

#pragma omp simd
  for (int i = 0; i < num_interacting_pairs; i++) {
      const int a1 = interacting_pairs_first[i];
      const int a2 = interacting_pairs_second[i];
      total += compute_pair_energy(
          ligand_x[a1] - ligand_x[a2],
          ligand_y[a1] - ligand_y[a2],
          ligand_z[a1] - ligand_z[a2],
          l_vdw_radius[a1], l_vdw_radius[a2],
          l_is_hbond_acceptor[a1], l_is_hbond_donor[a1],
          l_is_hbond_acceptor[a2], l_is_hbond_donor[a2],
          l_is_hydrophobic[a1], l_is_hydrophobic[a2]
          );
    }
    return total;
  }

#define print_matrix(namematrix, size, msg, mat) do{      \
  printf(namematrix);                             \
  printf("[");                                    \
  for(int i = 0; i < (size); i++){             \
    if(i != 0) printf(", ");                      \
    printf(msg, (mat)[i]);                        \
  }                                               \
  printf("]\n");                                  \
}while(false)                                             \

inline fp_type vina_scoring(  
    /// Protein data
    const int num_atoms_protein,
    const fp_type* __restrict__ protein_x,
    const fp_type* __restrict__ protein_y,
    const fp_type* __restrict__ protein_z,
    const int* __restrict__ p_is_hbond_acceptor,
    const int* __restrict__ p_is_hbond_donor,
    const int* __restrict__ p_is_hydrophobic,
    const fp_type* __restrict__ p_vdw_radius,

    ///Ligand data
    const int num_atoms_ligand,
    const fp_type* __restrict__ ligand_x,
    const fp_type* __restrict__ ligand_y,
    const fp_type* __restrict__ ligand_z,
    const int* __restrict__ l_is_hbond_acceptor,
    const int* __restrict__ l_is_hbond_donor,
    const int* __restrict__ l_is_hydrophobic,
    const fp_type* __restrict__ l_vdw_radius,
    const int active_torsions,
    const int* __restrict__ interacting_pairs_first,
    const int* __restrict__ interacting_pairs_second,
    const int num_interacting_pairs
    ){


#if 0
      printf("num_atoms_protein: %i\n", num_atoms_protein);
      print_matrix("protein_x", num_atoms_protein, "%f", protein_x);
      print_matrix("protein_y", num_atoms_protein, "%f", protein_y);
      print_matrix("protein_z", num_atoms_protein, "%f", protein_z);
      print_matrix("p_is_hbond_acceptor", num_atoms_protein, "%i", p_is_hbond_acceptor);
      print_matrix("p_is_hbond_donor", num_atoms_protein, "%i", p_is_hbond_donor);
      print_matrix("p_is_hydrophobic", num_atoms_protein, "%i", p_is_hydrophobic);
      print_matrix("p_vdw_radius", num_atoms_protein, "%f", p_vdw_radius);
#endif

#if 0
      printf("num_atoms_ligand: %i\n", num_atoms_ligand);
      print_matrix("ligand_x", num_atoms_ligand, "%f", ligand_x);
      print_matrix("ligand_y", num_atoms_ligand, "%f", ligand_y);
      print_matrix("ligand_z", num_atoms_ligand, "%f", ligand_z);
      print_matrix("l_is_hbond_acceptor", num_atoms_ligand, "%i", l_is_hbond_acceptor);
      print_matrix("l_is_hbond_donor", num_atoms_ligand, "%i",l_is_hbond_donor);
      print_matrix("l_is_hydrophobic", num_atoms_ligand, "%i", l_is_hydrophobic);
      print_matrix("l_vdw_radius", num_atoms_ligand, "%f", l_vdw_radius);

      printf("active torsions: %i\n", active_torsions);
      print_matrix("interacting_pairs_first", num_interacting_pairs, "%i", interacting_pairs_first);
      print_matrix("interacting_pairs_second", num_interacting_pairs, "%i", interacting_pairs_second);
#endif
      fp_type inter_score = score_inter(
          num_atoms_protein,
          protein_x,
          protein_y,
          protein_z,
          p_is_hbond_acceptor,
          p_is_hbond_donor,
          p_is_hydrophobic,
          p_vdw_radius,
          num_atoms_ligand,
          ligand_x,
          ligand_y,
          ligand_z,
          l_is_hbond_acceptor,
          l_is_hbond_donor,
          l_is_hydrophobic,
          l_vdw_radius
          );

      fp_type intra_score = score_intra(
          ligand_x,
          ligand_y,
          ligand_z,
          l_is_hbond_acceptor,
          l_is_hbond_donor,
          l_is_hydrophobic,
          l_vdw_radius,
          interacting_pairs_first,
          interacting_pairs_second,
          num_interacting_pairs
          );


      fp_type score = (inter_score + intra_score) / ( 1 + NROT_COEFF * active_torsions);

      // printf("Score inter %f, Score intra %f, Score %f\n", inter_score, intra_score, score); 

      return score;
    }

inline void calc_energy(
    const int batch_atoms,
    const int batch_ligands,
    const int scores_per_ligand,
    const int *__restrict__ num_atoms_b, 
    const fp_type *__restrict__ x_scratch_b,
    const fp_type *__restrict__ y_scratch_b,
    const fp_type *__restrict__ z_scratch_b,
    const int* __restrict__ l_is_hbond_acceptor_b,
    const int* __restrict__ l_is_hbond_donor_b,
    const int* __restrict__ l_is_hydrophobic_b,
    const fp_type* __restrict__ l_vdw_radius_b,
    const int* __restrict__ active_torsions_b,
    const int* __restrict__ interacting_pairs_first_b, 
    const int* __restrict__ interacting_pairs_second_b,
    const int* __restrict__ num_interacting_pairs_b, 
    const int* __restrict__ interacting_pairs_offset_b, 

    /// Protein data
    const int num_atoms_protein,
    const fp_type* __restrict__ protein_x,
    const fp_type* __restrict__ protein_y,
    const fp_type* __restrict__ protein_z,
    const int* __restrict__ p_is_hbond_acceptor,
    const int* __restrict__ p_is_hbond_donor,
    const int* __restrict__ p_is_hydrophobic,
    const fp_type* __restrict__ p_vdw_radius,
    fp_type *__restrict__ scores_b
    )
{

  for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
    int offset_interacting_pairs = interacting_pairs_offset_b[ligand_index];
    const int atom_stride  = ligand_index * batch_atoms;
    const int num_atoms_ligand    = num_atoms_b[ligand_index];

    const fp_type *__restrict__ scratch_x = x_scratch_b + atom_stride * scores_per_ligand;
    const fp_type *__restrict__ scratch_y = y_scratch_b + atom_stride * scores_per_ligand;
    const fp_type *__restrict__ scratch_z = z_scratch_b + atom_stride * scores_per_ligand;


    ///Ligand data
    const int* __restrict__ l_is_hbond_acceptor = l_is_hbond_acceptor_b + atom_stride;
    const int* __restrict__ l_is_hbond_donor = l_is_hbond_donor_b + atom_stride;
    const int* __restrict__ l_is_hydrophobic = l_is_hydrophobic_b + atom_stride;
    const fp_type* __restrict__ l_vdw_radius = l_vdw_radius_b + atom_stride;
    const int active_torsions = active_torsions_b[ligand_index];
    const int* __restrict__ interacting_pairs_first =  interacting_pairs_first_b + offset_interacting_pairs;
    const int* __restrict__ interacting_pairs_second =  interacting_pairs_second_b + offset_interacting_pairs;
    const int num_interacting_pairs = num_interacting_pairs_b[ligand_index];


    fp_type *__restrict__ scores_l = scores_b + ligand_index * scores_per_ligand;
    for (int scores_index = 0; scores_index < scores_per_ligand; ++scores_index) {

        const fp_type *__restrict__ scratch_x_l = scratch_x + scores_index * batch_atoms;
        const fp_type *__restrict__ scratch_y_l = scratch_y + scores_index * batch_atoms;
        const fp_type *__restrict__ scratch_z_l = scratch_z + scores_index * batch_atoms;

      fp_type score = vina_scoring(num_atoms_protein,
          protein_x,
          protein_y,
          protein_z,
          p_is_hbond_acceptor,
          p_is_hbond_donor, 
          p_is_hydrophobic,
          p_vdw_radius,
          num_atoms_ligand,
          scratch_x_l,
          scratch_y_l,
          scratch_z_l,
          l_is_hbond_acceptor,
          l_is_hbond_donor,
          l_is_hydrophobic,
          l_vdw_radius,
          active_torsions,
          interacting_pairs_first,
          interacting_pairs_second,
          num_interacting_pairs
          );

      scores_l[scores_index] = score;
    }
  }

}

template<>
void vina_score_kernel<queue_cpp>::operator()() {
  q->invoke_kernel<this->vina_region_name>(
      calc_energy,
      batch_atoms,
      batch_ligands,
      scores_per_ligand,
      num_ligand_atoms_b, 
      x_scratch_b,
      y_scratch_b,
      z_scratch_b,
      l_is_hbond_acceptor_b,
      l_is_hbond_donor_b,
      l_is_hydrophobic_b,
      l_vdw_radius_b,
      active_torsions_b,
      interacting_pairs_first_b, 
      interacting_pairs_second_b,
      num_interacting_pairs_b,
      interacting_pairs_offset_b,
      num_atoms_protein,
      protein_x,
      protein_y,
      protein_z,
      p_is_hbond_acceptor,
      p_is_hbond_donor,
      p_is_hydrophobic,
      p_vdw_radius,
      scores_b
        );
}
}
