#pragma once

#include <mudock/chem/vinardo_layer.hpp>
#include <mudock/chem/vinardo_type.hpp>
#include <mudock/molecule/graph.hpp>
#include <mudock/type_alias.hpp>

#include <cstdint>
#include <queue>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace mudock {

struct vinardo_protein_ligand_pair {
  int protein_atom_idx;
  int ligand_atom_idx;
  std::uint8_t hydrophobic_possible;
  std::uint8_t hbond_possible;
  mudock::fp_type radius_sum;
};

struct vinardo_ligand_ligand_pair {
  int ligand_atom_i_idx;
  int ligand_atom_j_idx;
  std::uint8_t hydrophobic_possible;
  std::uint8_t hbond_possible;
  mudock::fp_type radius_sum;
};

inline void append_protein_ligand_pairs(std::vector<vinardo_protein_ligand_pair>& pl_pairs,
                                        const vinardo_layer<mudock::dynamic_containers>& protein_layer,
                                        const vinardo_layer<mudock::static_containers>& ligand_layer,
                                        const auto& protein_radii,
                                        const auto& ligand_radii,
                                        const auto& protein_types,
                                        const auto& ligand_types,
                                        const auto& protein_donor,
                                        const auto& protein_acceptor,
                                        const auto& protein_hydro,
                                        const auto& ligand_donor,
                                        const auto& ligand_acceptor,
                                        const auto& ligand_hydro) {
  const auto protein_atoms = static_cast<std::size_t>(protein_layer.num_atoms());
  const auto ligand_atoms  = static_cast<std::size_t>(ligand_layer.num_atoms());

  pl_pairs.reserve(protein_atoms * ligand_atoms);
  for (std::size_t i = 0; i < protein_atoms; ++i) {
    for (std::size_t j = 0; j < ligand_atoms; ++j) {
      if (is_hydrogen(protein_types[i]) || is_hydrogen(ligand_types[j])) {
        continue;
      }
      vinardo_protein_ligand_pair pair;
      pair.protein_atom_idx     = static_cast<int>(i);
      pair.ligand_atom_idx      = static_cast<int>(j);
      pair.hydrophobic_possible = protein_hydro[i] && ligand_hydro[j];
      pair.hbond_possible       = (protein_donor[i] && ligand_acceptor[j]) ||
                            (protein_acceptor[i] && ligand_donor[j]);
      pair.radius_sum           = protein_radii[i] + ligand_radii[j];
      pl_pairs.push_back(pair);
    }
  }
}

inline void append_ligand_ligand_pairs(std::vector<vinardo_ligand_ligand_pair>& ll_pairs,
                                       const vinardo_layer<mudock::static_containers>& ligand_layer,
                                       const auto& ligand_radii,
                                       const auto& ligand_types,
                                       const auto& ligand_donor,
                                       const auto& ligand_acceptor,
                                       const auto& ligand_hydro,
                                       std::span<const std::uint8_t> relatively_movable,
                                       const std::vector<std::uint8_t>& within_three_bonds) {
  const auto num_atoms = static_cast<std::size_t>(ligand_layer.num_atoms());

  ll_pairs.reserve((num_atoms * (num_atoms - 1)) / 2);
  for (std::size_t i = 0; i < num_atoms; ++i) {
    for (std::size_t j = i + 1; j < num_atoms; ++j) {
      if (is_hydrogen(ligand_types[i]) || is_hydrogen(ligand_types[j])) {
        continue;
      }
      if (relatively_movable[i * num_atoms + j] && !(within_three_bonds[i * num_atoms + j])) {
        vinardo_ligand_ligand_pair pair;
        pair.ligand_atom_i_idx    = static_cast<int>(i);
        pair.ligand_atom_j_idx    = static_cast<int>(j);
        pair.hydrophobic_possible = ligand_hydro[i] && ligand_hydro[j];
        pair.hbond_possible       = (ligand_donor[i] && ligand_acceptor[j]) ||
                              (ligand_acceptor[i] && ligand_donor[j]);
        pair.radius_sum           = ligand_radii[i] + ligand_radii[j];
        ll_pairs.push_back(pair);
      }
    }
  }
}

inline std::vector<std::uint8_t> precompute_within_n_bonds(const auto& g,
                                                           const std::size_t num_atoms,
                                                           const std::size_t n) {
  using vertex_type = typename std::decay_t<decltype(g)>::vertex_descriptor;
  std::vector<std::uint8_t> within_n_bonds(num_atoms * num_atoms, 0);

  for (size_t i = 0; i < num_atoms; ++i) {
    //For each atom i, we perform a breadth first search up to depth n to find all atoms that are within n bonds from i.
    std::queue<vertex_type> current_queue;
    std::queue<vertex_type> next_queue;

    current_queue.push(static_cast<vertex_type>(i));
    std::vector<std::uint8_t> visited(num_atoms, 0);
    visited[i] = 1;
    // Note: this BFS assumes atom indices match graph vertex descriptors.
    // This is currently true because make_graph adds one vertex per atom in index order, but imo it's fragile.
    //TODO: consider finding another way to build the graph that guarantees this property.
    for (size_t depth = 0; depth < n; ++depth) {
      while (!current_queue.empty()) {
        auto visiting_node = current_queue.front();  //This return the value
        current_queue.pop();                         //This remove the value
        const auto [begin, end] = boost::adjacent_vertices(visiting_node, g);

        for (auto iterator = begin; iterator != end; ++iterator) {
          //The iterator points towards a vertex descriptr
          auto neighbor     = *iterator;
          auto neighbor_idx = g[neighbor].atom_index;
          if (visited[neighbor_idx] != 1) {
            next_queue.push(neighbor);
            visited[neighbor_idx] = 1;
          }
        }
      }
      //Swap the queue
      std::swap(current_queue, next_queue);
    }

    // visited[j] == 1 means that atom j is reachable from source atom i
    // within at most n BFS levels, i.e. within n covalent bonds.

    for (size_t j = i + 1; j < num_atoms; ++j) {
      if (visited[j] == 1) {
        within_n_bonds[i * num_atoms + j] = 1;
        within_n_bonds[j * num_atoms + i] = 1;
      }
    }
  }
  return within_n_bonds;
}

inline std::vector<vinardo_ligand_ligand_pair>
preprocess_ligand_vinardo(vinardo_layer<mudock::static_containers>& ligand_layer,
                          std::span<const std::uint8_t> relatively_movable_matrix) {
  const auto ligand_radii    = ligand_layer.get_radius();
  const auto ligand_donor    = ligand_layer.get_is_hbond_donor();
  const auto ligand_acceptor = ligand_layer.get_is_hbond_acceptor();
  const auto ligand_hydro    = ligand_layer.get_is_hydrophobic();
  const auto ligand_types    = ligand_layer.get_vinardo_type();

  std::vector<vinardo_ligand_ligand_pair> ll_pairs;
  auto& ligand         = ligand_layer.get_base_molecule();
  auto graph           = mudock::make_graph(ligand.get_bonds(), ligand.num_atoms());
  const auto num_atoms = static_cast<std::size_t>(ligand.num_atoms());

  if (relatively_movable_matrix.size() != num_atoms * num_atoms) {
    throw std::runtime_error("Vinardo preprocessing requires a ligand mobility matrix");
  }
  auto within_three_bonds = precompute_within_n_bonds(graph, num_atoms, 3);

  append_ligand_ligand_pairs(ll_pairs,
                             ligand_layer,
                             ligand_radii,
                             ligand_types,
                             ligand_donor,
                             ligand_acceptor,
                             ligand_hydro,
                             relatively_movable_matrix,
                             within_three_bonds);
  return ll_pairs;
}

inline std::vector<vinardo_protein_ligand_pair>
preprocess_protein_ligand_vinardo(vinardo_layer<mudock::static_containers>& ligand_layer,
                                  vinardo_layer<mudock::dynamic_containers>& protein_layer) {
  const auto protein_radii    = protein_layer.get_radius();
  const auto ligand_radii     = ligand_layer.get_radius();
  const auto protein_donor    = protein_layer.get_is_hbond_donor();
  const auto protein_acceptor = protein_layer.get_is_hbond_acceptor();
  const auto protein_hydro    = protein_layer.get_is_hydrophobic();
  const auto protein_types    = protein_layer.get_vinardo_type();
  const auto ligand_donor     = ligand_layer.get_is_hbond_donor();
  const auto ligand_acceptor  = ligand_layer.get_is_hbond_acceptor();
  const auto ligand_hydro     = ligand_layer.get_is_hydrophobic();
  const auto ligand_types     = ligand_layer.get_vinardo_type();

  std::vector<vinardo_protein_ligand_pair> pl_pairs;
  append_protein_ligand_pairs(pl_pairs,
                              protein_layer,
                              ligand_layer,
                              protein_radii,
                              ligand_radii,
                              protein_types,
                              ligand_types,
                              protein_donor,
                              protein_acceptor,
                              protein_hydro,
                              ligand_donor,
                              ligand_acceptor,
                              ligand_hydro);

  return pl_pairs;
}

} // namespace mudock
