#include <mudock/chem/vinardo_layer.hpp>
#include <mudock/chem/vinardo_preprocessing.hpp>
#include <mudock/compute/vinardo_scoring_function.hpp>
#include <mudock/type_alias.hpp>

#include <algorithm>
#include <cmath>

namespace mudock {

inline fp_type gauss(const fp_type d) {
  constexpr fp_type s1 = fp_type{0.8};

  const fp_type x = d / s1;
  return std::exp(-(x * x));
}

inline fp_type repulsion(const fp_type d) {
  //If it's greater than 0 there is no repulsion, otherwise we have a quadratic repulsion
  const fp_type clamped = std::min(d, fp_type{0});
  return clamped * clamped;
}

inline fp_type hydrophobic(const fp_type d) {
  constexpr fp_type p1 = fp_type{0};
  constexpr fp_type p2 = fp_type{2.5};

  if (d <= p1) {
    return fp_type{1};
  } else if (d < p2) {
    // Paper formula: p2 - d. Smina implements the same interval as a normalized slope_step.
    // return fp_type{p2 - d};
    return (p2 - d) / (p2 - p1);
  } else {
    return fp_type{0};
  }
}

inline fp_type hbond(const fp_type d) {
  constexpr fp_type h1 = -0.6;

  if (d <= h1) {
    return fp_type{1};
  } else if (d < fp_type{0}) {
    return d / (h1);
  } else {
    return fp_type{0};
  }
}

fp_type compute_vinardo_pair_energy(const fp_type surface_distance,
                                    const bool hydrophobic_possible,
                                    const bool hbond_possible) {
  //Here we define the constant weights of the scoring function
  constexpr fp_type w1 = -0.045;
  constexpr fp_type w2 = 0.800;
  constexpr fp_type w3 = -0.035;
  constexpr fp_type w4 = -0.600;

  //We use this to avoid the if in the inner loop
  const fp_type hydro_mask = static_cast<fp_type>(hydrophobic_possible);
  const fp_type hbond_mask = static_cast<fp_type>(hbond_possible);

  return w1 * gauss(surface_distance) + w2 * repulsion(surface_distance) +
         hydro_mask * w3 * hydrophobic(surface_distance) + hbond_mask * w4 * hbond(surface_distance);
}

vinardo_pair_terms compute_vinardo_pair_terms(const fp_type surface_distance,
                                              const bool hydrophobic_possible,
                                              const bool hbond_possible) {
  //Here we define the constant weights of the scoring function
  constexpr fp_type w1 = -0.045;
  constexpr fp_type w2 = 0.800;
  constexpr fp_type w3 = -0.035;
  constexpr fp_type w4 = -0.600;

  //We use this to avoid the if in the inner loop
  const fp_type hydro_mask = static_cast<fp_type>(hydrophobic_possible);
  const fp_type hbond_mask = static_cast<fp_type>(hbond_possible);

  vinardo_pair_terms terms{};
  terms.gauss                = gauss(surface_distance);
  terms.repulsion            = repulsion(surface_distance);
  terms.hydrophobic          = hydro_mask * hydrophobic(surface_distance);
  terms.hbond                = hbond_mask * hbond(surface_distance);
  terms.weighted_gauss       = w1 * terms.gauss;
  terms.weighted_repulsion   = w2 * terms.repulsion;
  terms.weighted_hydrophobic = w3 * terms.hydrophobic;
  terms.weighted_hbond       = w4 * terms.hbond;
  terms.weighted_energy =
      terms.weighted_gauss + terms.weighted_repulsion + terms.weighted_hydrophobic + terms.weighted_hbond;
  return terms;
}

vinardo_score_breakdown compute_vinardo_score_breakdown(vinardo_layer<dynamic_containers>& protein,
                                                        vinardo_layer<static_containers>& ligand,
                                                        std::vector<vinardo_protein_ligand_pair>& pl_pairs,
                                                        std::vector<vinardo_ligand_ligand_pair>& ll_pairs) {

  fp_type pl_score = fp_type{0};
  fp_type ll_score = fp_type{0};

  const auto protein_x = protein().x();
  const auto protein_y = protein().y();
  const auto protein_z = protein().z();
  const auto ligand_x  = ligand().x();
  const auto ligand_y  = ligand().y();
  const auto ligand_z  = ligand().z();

  //Now we have to process all the couples, the vinardo scoring function is characterized by being
  //an averaged sum of four terms, more precisely given a couple i,j we have:
  // f_{i,j} = w1 * gauss(d) + w2 * repulsion(d) + w3 * hydrophobic + w4 * hbond
  // where d is defined as surface distance and can be obtained by:
  // d = r_i_j - (R_i + R_j)
  //Where the sum R_i and R_j is already precomputed for every couple and the distance between the two atoms has to be calculated

  //Now we can process the protein-ligand pairs

  for (const auto& pair: pl_pairs) {
    const fp_type dx    = protein_x[pair.protein_atom_idx] - ligand_x[pair.ligand_atom_idx];
    const fp_type dy    = protein_y[pair.protein_atom_idx] - ligand_y[pair.ligand_atom_idx];
    const fp_type dz    = protein_z[pair.protein_atom_idx] - ligand_z[pair.ligand_atom_idx];
    const fp_type r_i_j = std::sqrt(dx * dx + dy * dy + dz * dz);

    const fp_type d = r_i_j - pair.radius_sum;

    //This skip is done is smina, but it's not grounded on the paper
    if (r_i_j >= fp_type{8}) {
      continue;
    }

    pl_score += compute_vinardo_pair_energy(d, pair.hydrophobic_possible, pair.hbond_possible);
  }

  for (const auto& pair: ll_pairs) {
    const fp_type dx    = ligand_x[pair.ligand_atom_i_idx] - ligand_x[pair.ligand_atom_j_idx];
    const fp_type dy    = ligand_y[pair.ligand_atom_i_idx] - ligand_y[pair.ligand_atom_j_idx];
    const fp_type dz    = ligand_z[pair.ligand_atom_i_idx] - ligand_z[pair.ligand_atom_j_idx];
    const fp_type r_i_j = std::sqrt(dx * dx + dy * dy + dz * dz);

    const fp_type d = r_i_j - pair.radius_sum;

    //This skip is done is smina, but it's not grounded on the paper
    if (r_i_j >= fp_type{8}) {
      continue;
    }

    ll_score += compute_vinardo_pair_energy(d, pair.hydrophobic_possible, pair.hbond_possible);
  }

  return vinardo_score_breakdown{pl_score, ll_score, pl_score + ll_score};
}

fp_type compute_vinardo_score(vinardo_layer<dynamic_containers>& protein,
                              vinardo_layer<static_containers>& ligand,
                              std::vector<vinardo_protein_ligand_pair>& pl_pairs,
                              std::vector<vinardo_ligand_ligand_pair>& ll_pairs) {
  fp_type score = fp_type{0};

  const auto protein_x = protein().x();
  const auto protein_y = protein().y();
  const auto protein_z = protein().z();
  const auto ligand_x  = ligand().x();
  const auto ligand_y  = ligand().y();
  const auto ligand_z  = ligand().z();

  for (const auto& pair: pl_pairs) {
    const fp_type dx    = protein_x[pair.protein_atom_idx] - ligand_x[pair.ligand_atom_idx];
    const fp_type dy    = protein_y[pair.protein_atom_idx] - ligand_y[pair.ligand_atom_idx];
    const fp_type dz    = protein_z[pair.protein_atom_idx] - ligand_z[pair.ligand_atom_idx];
    const fp_type r_i_j = std::sqrt(dx * dx + dy * dy + dz * dz);

    if (r_i_j >= fp_type{8}) {
      continue;
    }

    const fp_type d = r_i_j - pair.radius_sum;
    score += compute_vinardo_pair_energy(d, pair.hydrophobic_possible, pair.hbond_possible);
  }

  for (const auto& pair: ll_pairs) {
    const fp_type dx    = ligand_x[pair.ligand_atom_i_idx] - ligand_x[pair.ligand_atom_j_idx];
    const fp_type dy    = ligand_y[pair.ligand_atom_i_idx] - ligand_y[pair.ligand_atom_j_idx];
    const fp_type dz    = ligand_z[pair.ligand_atom_i_idx] - ligand_z[pair.ligand_atom_j_idx];
    const fp_type r_i_j = std::sqrt(dx * dx + dy * dy + dz * dz);

    if (r_i_j >= fp_type{8}) {
      continue;
    }

    const fp_type d = r_i_j - pair.radius_sum;
    score += compute_vinardo_pair_energy(d, pair.hydrophobic_possible, pair.hbond_possible);
  }

  return score;
}

} // namespace mudock
