#pragma once

#include <cstdint>
#include <mudock/grid/mdindex.hpp>
#include <mudock/molecule/bond.hpp>
#include <mudock/molecule/containers.hpp>
#include <mudock/type_alias.hpp>
#include <span>
#include <utility>

namespace mudock {

  template<class container_aliases>
    requires is_container_specification<container_aliases>
  class fragments {
    using fragment_mask_type   = int;
    using mask_container_type  = container_aliases::template fragments_size<fragment_mask_type>;
    using index_container_type = container_aliases::template fragments_size<std::size_t>;

    // this is a 2D grid that keep track of which atom belong to the related fragment
    mask_container_type mask;
    index2D index;

    // store which is the atom index that define the rotatable bonds
    index_container_type start_atom_indices;
    index_container_type stop_atom_indices;

  public:
    using value_type = fragment_mask_type;
    fragments(const std::span<const bond>& bonds, const std::size_t num_atoms);

    // utility functions to get the whole container
    [[nodiscard]] inline auto get_mask(const std::size_t bond_index) {
      return std::span(std::begin(mask) + index.to1D(0, bond_index), index.size_x());
    }
    [[nodiscard]] inline auto get_mask(const std::size_t bond_index) const {
      return std::span(std::cbegin(mask) + index.to1D(0, bond_index), index.size_x());
    }

    // utility functions to access the data
    [[nodiscard]] inline int& get_mask(const std::size_t bond_index, const std::size_t atom_index) {
      return mask[index.to1D(atom_index, bond_index)];
    }
    [[nodiscard]] inline const int& get_mask(const std::size_t bond_index,
                                             const std::size_t atom_index) const {
      return mask[index.to1D(atom_index, bond_index)];
    }

    // utility function to get the indices of the atoms related to the rotatable bonds
    [[nodiscard]] inline std::pair<std::size_t, std::size_t>
        get_rotatable_atoms(const std::size_t bond_index) const {
      return std::make_pair(start_atom_indices[bond_index], stop_atom_indices[bond_index]);
    }

    // utility function to get the number of rotatable bonds
    [[nodiscard]] inline auto get_num_rotatable_bonds() const { return index.size_y(); }
  };

  //===------------------------------------------------------------------------------------------------------
  // Out-of-class method definitions
  //===------------------------------------------------------------------------------------------------------

  template<>
  fragments<static_containers>::fragments(const std::span<const bond>& bonds, const std::size_t num_atoms);

  template<>
  fragments<dynamic_containers>::fragments(const std::span<const bond>& bonds, const std::size_t num_atoms);

} // namespace mudock
