#include <cstdio>
#include <vector>
#include <cmath>

#include <mudock/format.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/cpp_implementation/vina_score_cpp.hpp>

#define GAUSS1_COEFF        (- 0.035579f)
#define GAUSS2_COEFF        (- 0.005156f)
#define REPULSION_COEFF     (0.840245f)
#define HYDROPHOBIC_COEFF   (- 0.035069f)
#define H_BOND_COEFF        (- 0.587439f)
#define NROT_COEFF          (0.05846f)

namespace mudock {

  std::unordered_map<int, std::vector<int>> get_atoms_in_frag(
      const std::span<const bond>& bonds, 
      const std::size_t num_atom
      ){
    auto graph = make_graph(bonds, num_atom);
    const auto ligand_fragments =
      std::make_unique<fragments<static_containers>>(graph,
          bonds,
          num_atom);

    auto rigid_pieces = ligand_fragments.get()->get_rigid_pieces();

    std::unordered_map<int, std::vector<int>> atoms_in_fragment;

    for (size_t i = 0; i < num_atom; ++i) {
      atoms_in_fragment[rigid_pieces[i]].emplace_back(i);
    }

    return atoms_in_fragment;
  }


  std::pair<std::vector<int>, std::vector<int>> get_interactive_pairs(const static_molecule& ligand){

    std::pair<std::vector<int>, std::vector<int>> out;

    const std::span<const bond>& bonds = ligand.get_bonds(); 
    const size_t num_atom = ligand.num_atoms();

    std::unordered_map<int, std::vector<int>> atoms_in_fragment = get_atoms_in_frag(bonds, num_atom);

    const auto num_rotamers = atoms_in_fragment.size();

    for (size_t rot1 = 0; rot1 < num_rotamers; ++rot1) {

      /// const auto* bitmask_rot1 = frag_masks + rot1 * num_atoms;
      std::vector<int> atoms_rot1 = atoms_in_fragment[rot1];
      for(size_t i = 0; i < atoms_rot1.size(); ++i){

        int atom1 = atoms_rot1[i];

        for(size_t rot2 = rot1 + 1; rot2 < num_rotamers; ++rot2){

          /// const auto* bitmask_rot2 = frag_masks + rot2 * num_atoms;
          std::vector<int> atoms_rot2 = atoms_in_fragment[rot2];

          for(size_t j = 0; j < atoms_rot2.size(); ++j){

            int atom2 = atoms_rot2[j];

            /// Search for the atom2 in the neighbors of atom1
            bool found = false;
            for(size_t nb = 0; nb < max_static_neighbors() && !found; nb++) {
              int atom1_nb = ligand.neighbors(atom1, nb);
              if (atom1_nb == atom2) found = true;
              else if (atom1_nb == -1) break;
            }
            if (found) continue;

            /// Opendock: if [i, j] in self.torsion_bond_index or [j, i] in self.torsion_bond_index: continue

            //int mask_1 = bitmask_rot1[atom1];
            //int mask_2 = bitmask_rot2[atom2];
            ///if(mask_1 == 3 || mask_2 == 3 || mask_1 == 2 || mask_2 == 2) continue;

            /// Order pair before adding it to the output
            int a = std::min(atom1, atom2);
            int b = std::max(atom1, atom2);
            out.first.emplace_back(a);
            out.second.emplace_back(b);
          }
        }
      }
    }

    info("Total interacting pairs: ", out.first.size());
    return out;
  }


  inline fp_type distance(fp_type x, fp_type y, fp_type z) {
    return sqrt( x*x + y*y + z*z );
  }

  inline fp_type gauss1(const size_t idx, const fp_type* __restrict__ dst_mtx) {
    fp_type gauss1 = 0;
    if(dst_mtx[idx] != 0) gauss1 = exp(- pow(dst_mtx[idx] / 0.5f, 2));
    return gauss1;
  }

  inline fp_type gauss2(const size_t idx, const fp_type* __restrict__ dst_mtx) {
    fp_type gauss2 = 0;
    if(dst_mtx[idx] != 0) gauss2 = exp(- pow((dst_mtx[idx] - 3) / 2, 2));
    return gauss2;
  }

  inline fp_type repulsion(const size_t idx, const fp_type* __restrict__ dst_mtx) {
    return pow((dst_mtx[idx] < 0) * dst_mtx[idx], 2);
  }

  inline fp_type hydrophobic(const size_t idx, const fp_type* __restrict__ dst_mtx, const int* __restrict__ rec_lig_is_hydrophobic) {
    bool hydro_1 = rec_lig_is_hydrophobic[idx] && (dst_mtx[idx] <= 0.5f);
    bool hydro_2_cond = rec_lig_is_hydrophobic[idx] && (dst_mtx[idx] > 0.5f) && (dst_mtx[idx] < 1.5f);
    fp_type hydro_2 = 1.5f * hydro_2_cond - hydro_2_cond * dst_mtx[idx];
    return hydro_1 + hydro_2;
  }

  inline fp_type hbonding(const size_t idx, const fp_type* __restrict__ dst_mtx, const int* __restrict__ rec_lig_is_hb) {
    bool h_bond_1 = rec_lig_is_hb[idx] && (dst_mtx[idx] <= -0.7f);
    bool h_bond_2_cond = rec_lig_is_hb[idx] && (dst_mtx[idx] < 0.f) && (dst_mtx[idx] > -0.7f);
    fp_type h_bond_2 = h_bond_2_cond * (- dst_mtx[idx]) / 0.7f;
    return h_bond_1 + h_bond_2;
  }

  fp_type score_function(
      const fp_type* __restrict__ dst_mtx, 
      const int* __restrict__ ij_is_hydrophobic,
      const int* __restrict__ ij_is_hbond,
      const size_t size
      ) {

    fp_type g1 = 0; 
    fp_type g2 = 0; 
    fp_type rep = 0;
    fp_type hydro = 0;
    fp_type hbond = 0;

    for(size_t i = 0; i < size; i++){
      g1  += gauss1(i, dst_mtx);
      g2  += gauss2(i, dst_mtx);
      rep += repulsion(i, dst_mtx); 
      hydro += hydrophobic(i, dst_mtx, ij_is_hydrophobic);
      hbond += hbonding(i, dst_mtx, ij_is_hbond);
    }

    return GAUSS1_COEFF * g1 + GAUSS2_COEFF * g2 + REPULSION_COEFF * rep + HYDROPHOBIC_COEFF * hydro + H_BOND_COEFF * hbond;
  }

  void parse_data(
      /// Protein data
      const size_t num_atoms_protein,
      const fp_type* __restrict__ protein_x,
      const fp_type* __restrict__ protein_y,
      const fp_type* __restrict__ protein_z,
      const int* __restrict__ p_is_hbond_acceptor,
      const int* __restrict__ p_is_hbond_donor,
      const int* __restrict__ p_is_hydrophobic,
      const fp_type* __restrict__ p_vdw_radius,

      ///Ligand data
      const size_t num_atoms_ligand,
      const fp_type* __restrict__ ligand_x,
      const fp_type* __restrict__ ligand_y,
      const fp_type* __restrict__ ligand_z,
      const int* __restrict__ l_is_hbond_acceptor,
      const int* __restrict__ l_is_hbond_donor,
      const int* __restrict__ l_is_hydrophobic,
      const fp_type* __restrict__ l_vdw_radius,

      std::vector<fp_type>& dst_mtx,
      std::vector<int>& rec_lig_is_hbond, 
      std::vector<int>& rec_lig_is_hydrophobic   
        ){

          for(size_t proteinIdx = 0; proteinIdx < num_atoms_protein; proteinIdx++){
            for(size_t ligandIdx = 0; ligandIdx < num_atoms_ligand; ligandIdx++){
              fp_type dst = distance(
                  protein_x[proteinIdx] - ligand_x[ligandIdx],
                  protein_y[proteinIdx] - ligand_y[ligandIdx],
                  protein_z[proteinIdx] - ligand_z[ligandIdx]
                  );

              if(dst > 8) continue;

              dst -= p_vdw_radius[proteinIdx] + l_vdw_radius[ligandIdx];

              dst_mtx.push_back(dst);
              rec_lig_is_hbond.push_back(
                  (p_is_hbond_acceptor[proteinIdx] && l_is_hbond_donor[ligandIdx]) || (l_is_hbond_acceptor[ligandIdx] && p_is_hbond_donor[proteinIdx])
                  );
              rec_lig_is_hydrophobic.push_back(
                  p_is_hydrophobic[proteinIdx] && l_is_hydrophobic[ligandIdx]
                  );
            }
          }
        }

  void parse_intra_data(
      const fp_type* __restrict__ ligand_x,
      const fp_type* __restrict__ ligand_y,
      const fp_type* __restrict__ ligand_z,
      const int* __restrict__ l_is_hbond_acceptor,
      const int* __restrict__ l_is_hbond_donor,
      const int* __restrict__ l_is_hydrophobic,
      const fp_type* __restrict__ l_vdw_radius,
      const int* __restrict__ interacting_pairs_first,
      const int* __restrict__ interacting_pairs_second,
      const size_t num_interacting_pairs,

      std::vector<fp_type>& intra_dst_mtx,
      std::vector<int>& intra_rec_lig_is_hbond, 
      std::vector<int>& intra_rec_lig_is_hydrophobic   
      ){

    for(size_t i = 0; i < num_interacting_pairs; i++){
      int atom_1 = interacting_pairs_first[i];
      int atom_2 = interacting_pairs_second[i];

      fp_type dst = distance(
          ligand_x[atom_1] - ligand_x[atom_2],
          ligand_y[atom_1] - ligand_y[atom_2],
          ligand_z[atom_1] - ligand_z[atom_2]
          );

      if(dst > 8) continue;

      dst -= l_vdw_radius[atom_1] + l_vdw_radius[atom_2];

      intra_dst_mtx.push_back(dst);
      intra_rec_lig_is_hbond.push_back(
          (l_is_hbond_acceptor[atom_1] && l_is_hbond_donor[atom_2]) || (l_is_hbond_acceptor[atom_2] && l_is_hbond_donor[atom_1])
          );
      intra_rec_lig_is_hydrophobic.push_back(
          l_is_hydrophobic[atom_1] && l_is_hydrophobic[atom_2]
          );
    }
  }

#define print_matrix(namematrix, size, msg, mat) do{      \
          printf(namematrix);                             \
          printf("[");                                    \
          for(size_t i = 0; i < (size); i++){             \
            if(i != 0) printf(", ");                      \
            printf(msg, (mat)[i]);                        \
          }                                               \
          printf("]\n");                                  \
}while(false)                                             \

  inline fp_type vina_scoring(  
      /// Protein data
      const size_t num_atoms_protein,
      const fp_type* __restrict__ protein_x,
      const fp_type* __restrict__ protein_y,
      const fp_type* __restrict__ protein_z,
      const int* __restrict__ p_is_hbond_acceptor,
      const int* __restrict__ p_is_hbond_donor,
      const int* __restrict__ p_is_hydrophobic,
      const fp_type* __restrict__ p_vdw_radius,

      ///Ligand data
      const size_t num_atoms_ligand,
      const fp_type* __restrict__ ligand_x,
      const fp_type* __restrict__ ligand_y,
      const fp_type* __restrict__ ligand_z,
      const int* __restrict__ l_is_hbond_acceptor,
      const int* __restrict__ l_is_hbond_donor,
      const int* __restrict__ l_is_hydrophobic,
      const fp_type* __restrict__ l_vdw_radius,
      const size_t active_torsions,
      const int* __restrict__ interacting_pairs_first,
      const int* __restrict__ interacting_pairs_second,
      const size_t num_interacting_pairs
        ){

#if 0
          print_matrix("protein_x", num_atoms_protein, "%f", protein_x);
          print_matrix("protein_y", num_atoms_protein, "%f", protein_y);
          print_matrix("protein_z", num_atoms_protein, "%f", protein_z);
          print_matrix("p_is_hbond_acceptor", num_atoms_protein, "%i", p_is_hbond_acceptor);
          print_matrix("p_is_hbond_donor", num_atoms_protein, "%i", p_is_hbond_donor);
          print_matrix("p_is_hydrophobic", num_atoms_protein, "%i", p_is_hydrophobic);
          print_matrix("p_vdw_radius", num_atoms_protein, "%f", p_vdw_radius);

          print_matrix("ligand_x", num_atoms_ligand, "%f", ligand_x);
          print_matrix("ligand_y", num_atoms_ligand, "%f", ligand_y);
          print_matrix("ligand_z", num_atoms_ligand, "%f", ligand_z);
          print_matrix("l_is_hbond_acceptor", num_atoms_ligand, "%i", l_is_hbond_acceptor);
          print_matrix("l_is_hbond_donor", num_atoms_ligand, "%i", l_is_hbond_donor);
          print_matrix("l_is_hydrophobic", num_atoms_ligand, "%i", l_is_hydrophobic);
          print_matrix("l_vdw_radius", num_atoms_ligand, "%f", l_vdw_radius);

          printf("active torsions: %zu\n", active_torsions);
          print_matrix("interacting_pairs_first", num_interacting_pairs, "%i", interacting_pairs_first);
          print_matrix("interacting_pairs_second", num_interacting_pairs, "%i", interacting_pairs_second);
#endif

          std::vector<fp_type> dst_mtx                    = std::vector<fp_type>();
          std::vector<int> rec_lig_is_hbond               = std::vector<int>();
          std::vector<int> rec_lig_is_hydrophobic         = std::vector<int>();

          parse_data(
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
              l_vdw_radius,

              dst_mtx, 
              rec_lig_is_hbond, 
              rec_lig_is_hydrophobic
                );

          fp_type inter_score = score_function(dst_mtx.data(), rec_lig_is_hydrophobic.data(), rec_lig_is_hbond.data(), dst_mtx.size());

          std::vector<fp_type> intra_dst_mtx                  = std::vector<fp_type>();
          std::vector<int> intra_rec_lig_is_hbond             = std::vector<int>();
          std::vector<int> intra_rec_lig_is_hydrophobic       = std::vector<int>();

          parse_intra_data(
              ligand_x,
              ligand_y,
              ligand_z,
              l_is_hbond_acceptor,
              l_is_hbond_donor,
              l_is_hydrophobic,
              l_vdw_radius,
              interacting_pairs_first,
              interacting_pairs_second,
              num_interacting_pairs,

              intra_dst_mtx,
              intra_rec_lig_is_hbond, 
              intra_rec_lig_is_hydrophobic   
              );

          fp_type intra_score = score_function(intra_dst_mtx.data(), intra_rec_lig_is_hydrophobic.data(), intra_rec_lig_is_hbond.data(), intra_dst_mtx.size());
          fp_type score = (inter_score + intra_score) / ( 1 + NROT_COEFF * active_torsions);

          printf("Score inter %f, Score intra %f, Score %f\n", inter_score, intra_score, score); 

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

      /// Protein data
      const size_t num_atoms_protein,
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

          size_t offset_interacting_pairs = 0;
          for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
            const int atom_stride  = ligand_index * batch_atoms;
            const int num_atoms_ligand    = num_atoms_b[ligand_index];

            const fp_type *__restrict__ scratch_x = x_scratch_b + atom_stride;
            const fp_type *__restrict__ scratch_y = y_scratch_b + atom_stride;
            const fp_type *__restrict__ scratch_z = z_scratch_b + atom_stride;

            ///Ligand data
            const int* __restrict__ l_is_hbond_acceptor = l_is_hbond_acceptor_b + atom_stride;
            const int* __restrict__ l_is_hbond_donor = l_is_hbond_donor_b + atom_stride;
            const int* __restrict__ l_is_hydrophobic = l_is_hydrophobic_b + atom_stride;
            const fp_type* __restrict__ l_vdw_radius = l_vdw_radius_b + atom_stride;
            const int active_torsions = active_torsions_b[ligand_index];
            const int* __restrict__ interacting_pairs_first =  interacting_pairs_first_b + offset_interacting_pairs;
            const int* __restrict__ interacting_pairs_second =  interacting_pairs_second_b + offset_interacting_pairs;
            const size_t num_interacting_pairs = num_interacting_pairs_b[ligand_index];
            offset_interacting_pairs += num_interacting_pairs;

            fp_type score = vina_scoring(num_atoms_protein,
                protein_x,
                protein_y,
                protein_z,
                p_is_hbond_acceptor,
                p_is_hbond_donor, 
                p_is_hydrophobic,
                p_vdw_radius,
                num_atoms_ligand,
                scratch_x,
                scratch_y,
                scratch_z,
                l_is_hbond_acceptor,
                l_is_hbond_donor,
                l_is_hydrophobic,
                l_vdw_radius,
                active_torsions,
                interacting_pairs_first,
                interacting_pairs_second,
                num_interacting_pairs
                );

            scores_b[ligand_index] = score;
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
