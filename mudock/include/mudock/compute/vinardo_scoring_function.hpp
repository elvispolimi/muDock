#pragma once

#include <mudock/chem/vinardo_preprocessing.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

struct vinardo_pair_terms {
  fp_type gauss;
  fp_type repulsion;
  fp_type hydrophobic;
  fp_type hbond;
  fp_type weighted_gauss;
  fp_type weighted_repulsion;
  fp_type weighted_hydrophobic;
  fp_type weighted_hbond;
  fp_type weighted_energy;
};

struct vinardo_score_breakdown {
  fp_type protein_ligand;
  fp_type ligand_ligand;
  fp_type total;
};

//Debug API's
vinardo_pair_terms compute_vinardo_pair_terms(fp_type surface_distance,
                                              bool hydrophobic_possible,
                                              bool hbond_possible);

vinardo_score_breakdown compute_vinardo_score_breakdown(
    vinardo_layer<dynamic_containers>& protein,
    vinardo_layer<static_containers>& ligand,
    std::vector<vinardo_protein_ligand_pair>& protein_ligand_pairs,
    std::vector<vinardo_ligand_ligand_pair>& ligand_ligand_pairs);

//Fast API's
fp_type compute_vinardo_pair_energy(fp_type surface_distance,
                                    bool hydrophobic_possible,
                                    bool hbond_possible);

fp_type compute_vinardo_score(vinardo_layer<dynamic_containers>& protein,
                              vinardo_layer<static_containers>& ligand,
                              std::vector<vinardo_protein_ligand_pair>& protein_ligand_pairs,
                              std::vector<vinardo_ligand_ligand_pair>& ligand_ligand_pairs);

} // namespace mudock
