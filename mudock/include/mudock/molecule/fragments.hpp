#pragma once

#include <cstdint>
#include <mudock/grid/mdindex.hpp>
#include <mudock/molecule/bond.hpp>
#include <mudock/molecule/containers.hpp>
#include <mudock/type_alias.hpp>
#include <span>

namespace mudock {

  template<class container_aliases>
  class fragments {
    template<typename T>
    using array_type = container_aliases::template fragments_size<T>;

    array_type<coordinate_type> storage;
    index2D index;

  public:
    fragments(const std::size_t num_atoms, const std::size_t num_bonds);
    [[nodiscard]] inline auto get_num_rotatable_bonds() const { return index.size_y(); }

    // utility functions to get the whole container
    [[nodiscard]] inline std::span<coordinate_type> get_mask(const std::size_t bond_index) {
      const auto begin = storage.begin() + index.to1D(0, bond_index);
      const auto end   = begin + index.size_x();
      assert(begin != std::end(storage) && end != std::end(storage));
      return std::span(begin, end);
    }
    [[nodiscard]] inline std::span<const coordinate_type> get_mask(const std::size_t bond_index) const {
      const auto begin = storage.cbegin() + index.to1D(0, bond_index);
      const auto end   = begin + index.size_x();
      assert(begin != std::end(storage) && end != std::end(storage));
      return std::span(begin, end);
    }

    // utility functions to access the data
    [[nodiscard]] inline coordinate_type& get_mask(const auto bond_index, const auto atom_index) {
      return storage[index.to1D(atom_index, bond_index)];
    }
    [[nodiscard]] inline const coordinate_type& get_mask(const auto bond_index, const auto atom_index) const {
      return storage[index.to1D(atom_index, bond_index)];
    }
  };

  template<class container_aliases>
  [[nodiscard]] fragments<container_aliases> make_fragments(const std::span<const bond>& bonds,
                                                            const std::size_t num_atoms,
                                                            const std::size_t num_bonds);

  //===------------------------------------------------------------------------------------------------------
  // Out-of-class method definitions
  //===------------------------------------------------------------------------------------------------------

  template<>
  fragments<static_containers>::fragments(const std::size_t num_atoms, const std::size_t num_bonds);

  template<>
  fragments<dynamic_containers>::fragments(const std::size_t num_atoms, const std::size_t num_bonds);

  template<>
  [[nodiscard]] fragments<static_containers>
      make_fragments<static_containers>(const std::span<const bond>& bonds,
                                        const std::size_t num_atoms,
                                        const std::size_t num_bonds);
  template<>
  [[nodiscard]] fragments<dynamic_containers>
      make_fragments<dynamic_containers>(const std::span<const bond>& bonds,
                                         const std::size_t num_atoms,
                                         const std::size_t num_bonds);

} // namespace mudock
