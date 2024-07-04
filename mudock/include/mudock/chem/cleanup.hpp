#pragma once

#include <mudock/molecule.hpp>

namespace mudock {
  // TODO clean up functions
  // - mergeLPS (lonepairs)
  // - cleanUpResidues (waters and nonstdres)

  // From AutoDockTools mergeNPHS 
  template<class container_aliases>
  void cleanup_nonpolar_hydrogen(molecule<container_aliases>& molecule, const molecule_graph_type& graph) {
    const auto& elements = molecule.elements();

    size_t hydrogen_number = 0;
    size_t no_bonded_hs    = 0;
    size_t bonded_hs       = 0;
    std::vector<size_t> nphs;
    std::vector<size_t> nphs_f;
    for (size_t index = 0; index < molecule.num_atoms(); index++)
      if (elements[index] == element::H) {
        // Count hydrogens
        hydrogen_number++;

        // TODO out_edges is bidirectional?
        const auto [begin, end] = boost::out_edges(index, graph);
        if (begin == end)
          // Count the hydrogens with no bonds
          no_bonded_hs++;
        else {
          // Count the hydrogens with at least one bond
          bonded_hs++;
          // Check whether there are any nphs hydrogens
          if (std::any_of(begin, end, [&](const auto bond) {
                return elements[graph[boost::target(bond, graph)].atom_index] == element::C;
              }))
            nphs.push_back(index);
          if (std::any_of(begin, end, [&](const auto bond) {
                return elements[graph[boost::target(bond, graph)].atom_index] == element::C;
              }))
            nphs.push_back(index);
        }
      }
    // MAKE sure there are some hydrogens
    if (!hydrogen_number)
      return;
    // Check whether there are any hydrogens with no bonds
    if (!no_bonded_hs)
      return;
    //  next check is superfluous
    if (!bonded_hs)
      return;

    // Check whether there are any nphs hydrogens
    // ASSUMPTION everything is on the same molecule
    // Remove NHPS atoms
    for (auto i: nphs) { molecule.remove_atom(i); }

    // Remove bonds connected to them
  }
} // namespace mudock
