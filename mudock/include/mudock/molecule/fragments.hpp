#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <mudock/grid/mdindex.hpp>
#include <mudock/molecule/bond.hpp>
#include <mudock/molecule/containers.hpp>
#include <mudock/molecule/graph.hpp>
#include <mudock/type_alias.hpp>
#include <span>
#include <utility>

namespace mudock {
  //===------------------------------------------------------------------------------------------------------
  // Utility function to find the rotatable edge on a graph (using a consistent container type)
  //===------------------------------------------------------------------------------------------------------
  // NOTE: we define an edge with their two edges, to avoid problems in graph update
  struct edge_description {
    using vertex_type = typename molecule_graph_type::vertex_descriptor;
    vertex_type source;
    vertex_type dest;
  };

  using vertex_type = typename molecule_graph_type::vertex_descriptor;
  template<class container_aliases>
  std::pair<typename container_aliases::template bonds_size<edge_description>, int>
      get_rotatable_edges(const std::span<const bond>& bonds, const molecule_graph_type& g);

  template<class container_aliases>
  requires is_container_specification<container_aliases>
  class fragments {
    using fragment_mask_type    = int;
    using mask_container_type   = container_aliases::template fragments_size<fragment_mask_type>;
    using pieces_container_type = container_aliases::template fragments_size<fragment_mask_type>;
    using index_container_type  = container_aliases::template fragments_size<int>;

    // this is a 2D grid that keep track of which atom belong to the related fragment. Each atom mask has the
    // following meaning:
    //  0 - the atom does not belong to the rotatable fragment
    //  1 - the atom belongs to the rotatable fragment
    //  2 - the atom is part of the rotatable bond and it is in the first fragment
    //  3 - the atom is part of the rotatable bond and it is in the second fragment
    mask_container_type mask;
    pieces_container_type rigid_pieces;
    index2D index;

    // store which is the atom index that define the rotatable bonds
    index_container_type start_atom_indices;
    index_container_type stop_atom_indices;

  public:
    using value_type = fragment_mask_type;
    fragments(molecule_graph_type& graph, const std::span<const bond>& bonds, const int num_atoms);

    // utility functions to get the whole container
    [[nodiscard]] inline auto get_mask(const int bond_index) {
      assert(bond_index < index.size_y());
      return std::span(std::begin(mask) + index.to1D(0, bond_index), index.size_x());
    }
    [[nodiscard]] inline auto get_mask(const int bond_index) const {
      assert(bond_index < index.size_y());
      return std::span(std::cbegin(mask) + index.to1D(0, bond_index), index.size_x());
    }
    [[nodiscard]] inline auto get_rigid_pieces() const { return std::span(rigid_pieces); }

    // utility functions to access the data
    [[nodiscard]] inline int& get_mask(const int bond_index, const int atom_index) {
      assert(bond_index < index.size_y());
      assert(atom_index < index.size_x());
      return mask[index.to1D(atom_index, bond_index)];
    }
    [[nodiscard]] inline const int& get_mask(const int bond_index, const int atom_index) const {
      assert(bond_index < index.size_y());
      assert(atom_index < index.size_x());
      return mask[index.to1D(atom_index, bond_index)];
    }

    // utility function to get the indices of the atoms related to the rotatable bonds
    [[nodiscard]] inline std::pair<int, int> get_rotatable_atoms(const int bond_index) const {
      assert(bond_index < index.size_y());
      return std::make_pair(start_atom_indices[bond_index], stop_atom_indices[bond_index]);
    }

    // utility function to get the number of rotatable bonds
    [[nodiscard]] inline auto get_num_rotatable_bonds() const { return index.size_y(); }
  };

  // //===------------------------------------------------------------------------------------------------------
  // // Out-of-class method definitions
  // //===------------------------------------------------------------------------------------------------------

  template<>
  fragments<static_containers>::fragments(molecule_graph_type& graph,
                                          const std::span<const bond>& bonds,
                                          const int num_atoms);

  template<>
  fragments<dynamic_containers>::fragments(molecule_graph_type& graph,
                                           const std::span<const bond>& bonds,
                                           const int num_atoms);

  //===------------------------------------------------------------------------------------------------------
  // Helper method definitions
  //===------------------------------------------------------------------------------------------------------

  inline auto same_fragment(const std::span<const int>& mask, const int& atom_id1, const int& atom_id2) {
    return mask[atom_id1] == mask[atom_id2] || (mask[atom_id1] == 0 && mask[atom_id2] == 2) ||
           (mask[atom_id2] == 0 && mask[atom_id1] == 2) || (mask[atom_id1] == 1 && mask[atom_id2] == 3) ||
           (mask[atom_id2] == 1 && mask[atom_id1] == 3);
  }

  inline auto rotable_atoms(const std::span<const int>& mask, const int& atom_id1, const int& atom_id2) {
    return (mask[atom_id1] == 2 && mask[atom_id2] == 3) || (mask[atom_id2] == 2 && mask[atom_id1] == 3);
  }

} // namespace mudock
