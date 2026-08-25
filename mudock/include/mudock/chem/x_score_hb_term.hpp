#pragma once

#include <mudock/chem/x_score_hb.hpp>
#include <mudock/chem/x_score_ligand.hpp>
#include <mudock/chem/x_score_protein.hpp>
#include <mudock/type_alias.hpp>
#include <vector>

namespace mudock {

  // Hydrogen-bond (HB) term

  struct x_score_hb_atoms {
    std::vector<fp_type> x, y, z;                // atom coordinates
    std::vector<fp_type> root_x, root_y, root_z; // HB root
    std::vector<fp_type> radius;                 // vdw radius
    std::vector<x_score_hb> hb;                  // HB
    std::vector<int> donor_type;                 // XScore Get_Donor_Type (0 none / 1 straight / 2 angled)
    std::vector<int> acceptor_type;              // XScore Get_Acceptor_Type (0 none / 1 straight / 2 angled)
    std::vector<int> donor_limit;                // Sum_HBonds step 3 saturation limit
    std::vector<int> acceptor_limit;             // Sum_HBonds step 4 lone-pair limit`
    std::vector<char> has_root;                  // whether a heavy-neighbour root could be computed

    [[nodiscard]] int size() const { return static_cast<int>(x.size()); }
  };

  // Step 1

  // build the H-bond atom data for the protein
  [[nodiscard]] x_score_hb_atoms build_protein_hb_atoms(const x_score_protein& prot);

  // build the H-bond atom data for the ligand
  [[nodiscard]] x_score_hb_atoms build_ligand_hb_atoms(const x_score_ligand& lig);

  // One candidate H-bond between a ligand atom and a protein atom
  struct x_score_hb_candidate {
    int li;       // ligand index
    int pj;       // protein index
    int type;     // 1: lig donor; 2: lig acceptor / prot donor; 3: lig acceptor / metal
    fp_type score;
  };

  // Step 2 
  
  // Every ligand/protein pair that can form an H-bond
  [[nodiscard]] std::vector<x_score_hb_candidate> get_hbond_pair_pl(const x_score_hb_atoms& lig_atoms,
                                                                    const x_score_hb_atoms& prot_atoms);

  // Step 3 

  // filtering candidates and calculation
  [[nodiscard]] fp_type sum_hbonds(std::vector<x_score_hb_candidate>& candidates,
                                   const x_score_hb_atoms& lig_atoms,
                                   const x_score_hb_atoms& prot_atoms);

  // Compute HB term for this ligand-protein pair
  [[nodiscard]] fp_type compute_x_score_hb(const x_score_hb_atoms& lig_atoms,
                                           const x_score_hb_atoms& prot_atoms);

} // namespace mudock