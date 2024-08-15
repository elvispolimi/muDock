#include "mudock/grid/grid_map.hpp"
#include "mudock/grid/mdindex.hpp"
#include "mudock/molecule.hpp"
#include "mudock/molecule/containers.hpp"

#include <cmath>
#include <cstdint>
#include <mudock/grid.hpp>
#include <mudock/kernels/calc_energy_cpp.hpp>

namespace mudock {

  fp_type trilinear_interpolation(const grid_map& map, const point3D& coord) {
    const point3D p1 = map.get_index_from_coordinates(coord);
    const point3D p2 = map.get_index_from_coordinates(add(coord, point3D{grid_spacing}));

    // TODO autodock use a different method, not much readable
    const point3D xyz_d   = divide(difference(coord, p1), difference(p2, p1));
    const point3D xyz_d_o = difference(static_cast<const point3D>(point3D{1}), xyz_d);

    const fp_type c00 = map.at(p1) * xyz_d_o.x + map.at({p2.x, p1.y, p1.z}) * xyz_d.x;
    const fp_type c10 = map.at({p1.x, p2.y, p1.z}) * xyz_d_o.x + map.at({p2.x, p2.y, p1.z}) * xyz_d.x;
    const fp_type c01 = map.at({p1.x, p1.y, p2.z}) * xyz_d_o.x + map.at({p2.x, p1.y, p2.z}) * xyz_d.x;
    const fp_type c11 = map.at({p1.x, p2.y, p2.z}) * xyz_d_o.x + map.at({p1.x, p1.y, p1.z}) * xyz_d.x;

    const fp_type c0 = c00 * xyz_d_o.y + c10 * xyz_d.y;
    const fp_type c1 = c01 * xyz_d_o.y + c11 * xyz_d.y;

    return c0 * xyz_d_o.z + c1 * xyz_d.z;
  }

  // nonbonds.cc for nbmatrix required by weed_bonds
  void nonbonds(grid<uint_fast8_t, index2D>& nbmatrix, const dynamic_molecule& ligand) {
    //
    // in "nbmatrix", the values 1 (and 4) mean this pair of atoms will be included in the internal, non-bonded list
    //                           0                                         ignored
    const size_t num_atoms = ligand.num_atoms();

    // set all nonbonds in nbmatrix to 1, except "1-1 interactions" (self interaction)
    for (size_t i = 0; i < num_atoms; i++) {
      for (size_t j = 0; j < num_atoms; j++) { nbmatrix.at(i, j) = 1; } // j
      nbmatrix.at(i, i) = 0;                                            /* 2005-01-10 RH & GMM */
    }
    for (auto& bond: ligand.get_bonds()) {
      // Ignore 1-2 Interactions
      nbmatrix.at(bond.source, bond.dest) = 0;
      nbmatrix.at(bond.dest, bond.source) = 0;
    }

    for (auto& bond_1: ligand.get_bonds())
      for (auto& bond_2: ligand.get_bonds()) { // loop over each atom "k" bonded to the current atom "j"
        if (bond_1.dest != bond_2.source)
          continue;

        // Ignore "1-3 Interactions"
        nbmatrix.at(bond_2.dest, bond_1.source) = 0;
        nbmatrix.at(bond_1.source, bond_2.dest) = 0;

        for (auto& bond_3: ligand.get_bonds()) {
          if (bond_2.dest != bond_3.source)
            continue;
          nbmatrix.at(bond_1.source, bond_3.dest) = 0;
          nbmatrix.at(bond_3.dest, bond_1.source) = 0;
        }
      }
  }

  // weedbonds.cc for nonbondlist only for the first group
  /*___________________________________________________________________________
  |    ENDBRANCH---TORS---BRANCH---R O O T---BRANCH---ENDBRANCH                |
  |                                  /              \                          |
  |                                BRANCH            BRANCH--TORS---ENDBRANCH  |
  |                                /                  \                        |
  |                               ENDBRANCH            ENDBRANCH               |
  |____________________________________________________________________________|
  |  Eliminate all rigidly bonded atoms:                                       |
  |                                     * those atoms which are at the ROOT;   |
  |                                     * atoms between TORS and BRANCH;       |
  |                                     * atoms between BRANCH and ENDBRANCH.  |
  |  This is necessary for internal energy calculations.                       |
  |____________________________________________________________________________|
  | Weed out bonds in rigid pieces,                                            |
  |____________________________________________________________________________|
  */
  void weed_bonds(grid<uint_fast8_t, index2D>& nbmatrix,
                  const static_molecule& ligand,
                  const fragments<dynamic_containers>& ligand_fragments) {
    for (size_t i = 0; i < ligand.num_bonds(); i++) {
      const auto& bond = ligand.get_bonds()[i];
      // TODO @Davide check this
      if (!bond.can_rotate)
        continue;
      const auto mask = ligand_fragments.get_mask(i);
      for (size_t m = 0; m < mask.size(); m++)
        for (size_t n = 0; n < mask.size(); n++) {
          // TODO @Davide means that they are on the same rigid piece?
          // Check readPDBQT.cc rigid_pieces
          // TODO what happens when they have type 2?
          if (mask[n] == mask[m]) {
            // Set the entry for atoms "i" and "j" in the
            //    nonbond matrix to 0 ;
            //
            // Later on, we will not calculate the interaction energy
            //    between atoms "i" and "j"
            nbmatrix.at(m, n) = 0;
            nbmatrix.at(n, m) = 0;
          }
        }
    }
  }

  fp_type calc_energy(const dynamic_molecule& receptor,
                      const static_molecule& ligand,
                      const fragments<static_containers>& ligand_fragments,
                      const grid_atom_mapper& grid_maps,
                      const grid_map& electro_map,
                      const grid_map& desolv_map) {
    fp_type elect_total = 0;
    fp_type emap_total  = 0;
    fp_type dmap_total  = 0;
    for (size_t index = 0; index < ligand.num_atoms(); ++index) {
      const point3D coord{ligand.x(index), ligand.y(index), ligand.z(index)};
      const auto& atom_charge = ligand.charge(index);
      const auto& atom_map    = grid_maps.get_atom_map(receptor.autodock_type(index));

      if (atom_map.outside_grid(coord)) {
        // TODO outside grid
      }

      // Trilinear Interpolation
      elect_total += trilinear_interpolation(electro_map, coord) * atom_charge;
      emap_total += trilinear_interpolation(atom_map, coord);
      dmap_total += trilinear_interpolation(desolv_map, coord) * std::fabs(atom_charge);
    }

    const int n_torsions = ligand_fragments.get_num_rotatable_bonds();
    if (n_torsions > 0) {
      const size_t num_atoms = ligand.num_atoms();
      // TODO @Davide suppose that the receptor does not have Flexible residues eintcal.cc:147
      // TODO
      grid<uint_fast8_t, index2D> nbmatrix{{num_atoms, num_atoms}};
    }
    return elect_total + emap_total + dmap_total;
  }

} // namespace mudock
