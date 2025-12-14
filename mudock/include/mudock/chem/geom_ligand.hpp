#pragma once

#include <functional>
#include <mudock/chem/molecule_layer.hpp>
#include <mudock/molecule/containers.hpp>
#include <mudock/molecule/fragments.hpp>

namespace mudock {
  struct geom_ligand: public molecule_layer<static_containers> {
    geom_ligand(static_molecule& _molecule): molecule_layer<static_containers>(_molecule) { prepare(); };
    [[nodiscard]] inline auto* fragments_masks() const { return frag_masks.data(); }
    [[nodiscard]] inline auto* fragmets_starts() const { return frag_start_indexes.data(); }
    [[nodiscard]] inline auto* fragments_stops() const { return frag_stop_indexes.data(); }

  private:
    // TODO probably you can optimize this, and write directly where they are needed on the buffers
    std::vector<int> frag_masks;
    std::vector<int> frag_start_indexes;
    std::vector<int> frag_stop_indexes;

    void prepare() {
      get_linearized_fragments_mask((*this).get_base_molecule().num_atoms(),
                                    (*this).get_base_molecule().num_rotamers(),
                                    frag_masks,
                                    frag_start_indexes,
                                    frag_stop_indexes,
                                    (*this).get_base_molecule());
    };
  };
} // namespace mudock
