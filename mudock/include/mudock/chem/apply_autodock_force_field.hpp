#pragma once

#include "mudock/chem/autodock_parameters.hpp"
#include "mudock/molecule/graph.hpp"

#include <mudock/chem/assign_autodock_babel_types.hpp>
#include <mudock/chem/assign_autodock_types.hpp>
#include <mudock/chem/autodock_babel_types.hpp>
#include <mudock/chem/autodock_types.hpp>
#include <mudock/molecule.hpp>

namespace mudock {
  template<class container_aliases>
    requires is_container_specification<container_aliases>
  void apply_autodock_forcefield(molecule<container_aliases>& molecule) {
    // allocate memory for the support vectors required to allocate the atoms type
    const std::size_t num_atoms = molecule.num_atoms();
    typename container_aliases::template atoms_size<autodock_babel_ff> mol_autodock_babel_types;
    typename container_aliases::template atoms_size<autodock_ff> mol_autodock_types;
    mudock::resize(mol_autodock_babel_types, num_atoms);
    mudock::resize(mol_autodock_types, num_atoms);

    // create the graph of the molecule
    const auto graph = make_graph(molecule.get_bonds());

    // assign the autodock babel type
    auto babel_type_span = make_span(mol_autodock_babel_types, num_atoms);
    assign_autodock_babel_types(babel_type_span,
                                molecule.get_x(),
                                molecule.get_y(),
                                molecule.get_z(),
                                molecule.get_elements(),
                                graph);

    // assign the autodock type
    assign_autodock_types(make_span(mol_autodock_types, num_atoms),
                          molecule.get_elements(),
                          molecule.get_is_aromatic(),
                          babel_type_span,
                          graph);

    // fill the atom properties using the autodock force field
    for (std::size_t index{0}; index < num_atoms; ++index) {
      const auto& ff_entry          = get_description(mol_autodock_types[index]);
      molecule.autodock_type(index) = ff_entry.value;
      molecule.Rii(index)           = ff_entry.Rii;
      molecule.epsii(index)         = ff_entry.epsii * autodock_parameters::coeff_vdW;
      molecule.vol(index)           = ff_entry.vol;
      molecule.solpar(index)        = ff_entry.solpar;
      molecule.Rij_hb(index)        = ff_entry.Rij_hb;
      molecule.epsij_hb(index)      = ff_entry.epsij_hb * autodock_parameters::coeff_hbond;
      molecule.num_hbond(index)     = ff_entry.hbond;
    }
  }
} // namespace mudock
