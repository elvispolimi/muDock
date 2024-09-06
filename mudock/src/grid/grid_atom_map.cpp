#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <mudock/chem.hpp>
#include <mudock/grid.hpp>
#include <mudock/grid/grid_map.hpp>
#include <mudock/log.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <ostream>
#include <unordered_map>
#include <vector>

namespace mudock {

  struct vdw_hb_energy_table {
    vdw_hb_energy_table(const fp_type cA, const fp_type cB, const fp_type dxA, const fp_type dxB) {
      e_vdW_Hb[0] = EINTCLAMP;
      for (size_t indx_r = 1; indx_r < NEINT - 1; ++indx_r) {
        const fp_type r  = indx_r / A_DIV;
        const fp_type rA = std::pow(r, dxA);
        const fp_type rB = std::pow(r, dxB);

        e_vdW_Hb[indx_r] = std::min(EINTCLAMP, (cA / rA - cB / rB));
      }
      e_vdW_Hb[NEINT - 1] = 0;

      /* smooth with min function */ /* GPF_MAP */

      /* Angstrom is divided by A_DIV in look-up table. */
      /* Typical value of r_smooth is 0.5 Angstroms  */
      /* so i_smooth = 0.5 * 100. / 2 = 25 */
      const size_t i_smooth = std::floor(r_smooth * A_DIV / fp_type{2});
      std::vector<fp_type> energy_smooth;
      energy_smooth.resize(NEINT, EINTCLAMP);
      if (i_smooth > 0) {
        for (size_t indx_r = 0; indx_r < NEINT; ++indx_r) {
          const auto temp = indx_r - i_smooth;
          size_t j        = temp <= indx_r ? temp : size_t{0};
          for (; j < std::min(size_t{NEINT}, indx_r + i_smooth + 1); j++)
            energy_smooth[indx_r] = std::min(energy_smooth[indx_r], e_vdW_Hb[j]);
        }
        for (size_t indx_r = 0; indx_r < NEINT; indx_r++) { e_vdW_Hb[indx_r] = energy_smooth[indx_r]; }
      } /* endif smoothing */
    };

    fp_type operator[](const size_t index) const { return e_vdW_Hb[index]; }

  private:
    static constexpr fp_type r_smooth{0.5}; //Angstrom

    std::vector<fp_type> e_vdW_Hb = std::vector<fp_type>(NEINT); // vdW & Hb energies
  };

  struct interaction {
    const fp_type cA;
    const fp_type cB;      // coefficients if specified in gpf
    const fp_type nbp_r;   // radius of energy-well minimum
    const fp_type nbp_eps; // depth of energy-well minimum
    const size_t xA;       // generally 12
    const size_t xB;       // 6 for non-hbonders 10 for h-bonders
    const size_t hbonder;
    const autodock_ff receptor_type;
    const vdw_hb_energy_table vdw_hb_table;

    interaction(const fp_type _cA,
                const fp_type _cB,
                const fp_type _nbp_r,
                const fp_type _nbp_eps,
                const size_t _xA,
                const size_t _xB,
                const size_t _hbonder,
                const autodock_ff _receptor_type)
        : cA(_cA),
          cB(_cB),
          nbp_r(_nbp_r),
          nbp_eps(_nbp_eps),
          xA(_xA),
          xB(_xB),
          hbonder(_hbonder),
          receptor_type(_receptor_type),
          vdw_hb_table(cA, cB, xA, xB) {};
  };

  struct scratchpad {
    grid_atom_map atom_map;
    fp_type hbondmin{999999};
    fp_type hbondmax{-999999};
    fp_type hbondflag{false};
    fp_type energy;
    const autodock_ff_description& grid_type_desc;

    scratchpad(const std::set<autodock_ff> receptor_types,
               const autodock_ff& ligand_type,
               const index3D& npts,
               const point3D& min,
               const point3D& max)
        : atom_map(ligand_type, npts, min, max), grid_type_desc(get_description(atom_map.get_atom_type())) {
      // Initialize data structure
      for (auto& receptor_t: receptor_types) {
        const auto& receptor_type_desc = get_description(receptor_t);

        fp_type nbp_r = (grid_type_desc.Rii + receptor_type_desc.Rii) / fp_type{2.};
        // TODO check if it is ok to use the floating point version
        fp_type nbp_eps = sqrtf(grid_type_desc.epsii * autodock_parameters::coeff_vdW *
                                receptor_type_desc.epsii * autodock_parameters::coeff_vdW);
        // TODO probably they are constant
        size_t xA      = 12;
        size_t xB      = 6;
        size_t hbonder = 0;
        if (grid_type_desc.hbond > 2 && (receptor_type_desc.hbond == 1 ||
                                         receptor_type_desc.hbond == 2)) { /*AS,A1,A2 map vs DS,D1 probe*/
          xB                  = 10;
          hbonder             = 1;
          atom_map.is_hbonder = true;
          nbp_r               = grid_type_desc.Rij_hb;
          nbp_eps             = grid_type_desc.epsij_hb * autodock_parameters::coeff_hbond;
        } else if ((grid_type_desc.hbond == 1 || grid_type_desc.hbond == 2) &&
                   (receptor_type_desc.hbond > 2)) { /*DS,D1 map vs AS,A1,A2 probe*/
          xB                  = 10;
          hbonder             = 1;
          atom_map.is_hbonder = true;
          nbp_r               = receptor_type_desc.Rij_hb;
          nbp_eps             = receptor_type_desc.epsij_hb * autodock_parameters::coeff_hbond;
        }; /*initialize energy parms for each possible receptor type*/
        // TODO to be checked later in line 1707 from main.cpp
        const fp_type cA = (nbp_eps / (xA - xB)) * std::pow(nbp_r, static_cast<fp_type>(xA)) * xB;
        const fp_type cB = nbp_eps / (xA - xB) * std::pow(nbp_r, static_cast<fp_type>(xB)) * xA;

        interactions.emplace(receptor_t, interaction{cA, cB, nbp_r, nbp_eps, xA, xB, hbonder, receptor_t});
      }
    }

    inline void reset() {
      hbondmin  = 999999;
      hbondmax  = -999999;
      hbondflag = false;
      energy    = 0;
      return;
    }

    inline void write_to_map(const size_t coord_x, const size_t coord_y, const size_t coord_z) {
      energy += hbondmin + hbondmax;
      /*
      * O U T P U T . . .
      *
      * Now output this grid point's energies to the maps:
      *
      */
      if (std::fabs(energy) < precision)
        energy = 0;
      atom_map.at(coord_x, coord_y, coord_z) = energy;
      // TODO
      // gridmap[k].energy_max = max(gridmap[k].energy_max, gridmap[k].energy);
      // gridmap[k].energy_min = min(gridmap[k].energy_min, gridmap[k].energy);
    }

    [[nodiscard]] inline auto& get_interaction(const autodock_ff r_t) { return interactions.at(r_t); }

  private:
    // TODO @Davide
    std::unordered_map<autodock_ff, interaction> interactions;
  };

  grid_atom_mapper generate_atom_grid_maps(dynamic_molecule& receptor) {
    // Allocate the correspondent scratchpads
    std::vector<scratchpad> scratchpads;

    // Get autodock types found in the receptor
    const auto& receptor_autodock_types = receptor.get_autodock_type();
    std::set<autodock_ff> receptor_types{receptor_autodock_types.begin(), receptor_autodock_types.end()};
    // TODO check if H should be removed
    // receptor_types.erase(std::remove(receptor_types.begin(), receptor_types.end(), mudock::autodock_ff::H));

    // Define the autodock ligand types
    constexpr std::array<autodock_ff, 12> ligand_types{autodock_ff::A,
                                                       autodock_ff::C,
                                                       autodock_ff::H,
                                                       autodock_ff::HD,
                                                       autodock_ff::N,
                                                       autodock_ff::NA,
                                                       autodock_ff::OA,
                                                       autodock_ff::SA,
                                                       autodock_ff::Cl,
                                                       autodock_ff::F,
                                                       autodock_ff::S,
                                                       autodock_ff::Br};
    //  Get maximum and minimum of the bounding box around the receptor
    const fp_type receptor_max_x = std::ranges::max(receptor.get_x());
    const fp_type receptor_max_y = std::ranges::max(receptor.get_y());
    const fp_type receptor_max_z = std::ranges::max(receptor.get_z());
    const fp_type receptor_min_x = std::ranges::min(receptor.get_x());
    const fp_type receptor_min_y = std::ranges::min(receptor.get_y());
    const fp_type receptor_min_z = std::ranges::min(receptor.get_z());

    const point3D grid_minimum{std::floor((receptor_min_x - cutoff_distance) * fp_type{2}) / fp_type{2},
                               std::floor((receptor_min_y - cutoff_distance) * fp_type{2}) / fp_type{2},
                               std::floor((receptor_min_z - cutoff_distance) * fp_type{2}) / fp_type{2}};
    const point3D grid_maximum{std::ceil((receptor_max_x + cutoff_distance) * fp_type{2}) / fp_type{2},
                               std::ceil((receptor_max_y + cutoff_distance) * fp_type{2}) / fp_type{2},
                               std::ceil((receptor_max_z + cutoff_distance) * fp_type{2}) / fp_type{2}};
    const index3D npts{static_cast<size_t>((grid_maximum.x - grid_minimum.x) / grid_spacing),
                       static_cast<size_t>((grid_maximum.y - grid_minimum.y) / grid_spacing),
                       static_cast<size_t>((grid_maximum.z - grid_minimum.z) / grid_spacing)};

    for (auto ligand_type: ligand_types) {
      // grid_atom_maps.push_back({ligand_type, npts});
      scratchpads.emplace_back(receptor_types, ligand_type, npts, grid_minimum, grid_maximum);
    }

    /* exponential function for receptor and ligand desolvation */
    /* note: the solvation term ranges beyond the non-bond cutoff
    * and will not be smoothed
    */
    std::array<fp_type, NDIEL> sol_fn;
    for (size_t indx_r = 1; indx_r < NDIEL; ++indx_r) {
      const fp_type r = indx_r / A_DIV;
      sol_fn[indx_r] =
          autodock_parameters::coeff_desolv * std::exp(-(r * r) / (fp_type{2} * (sigma * sigma)));
    }

    // TODO what are these???
    // TODO check it rvector and rvector2 are really needed
    std::vector<point3D> rvector(receptor.num_atoms());
    // TODO why rvector2? Seems to be like the square of the previous one
    std::vector<point3D> rvector2(receptor.num_atoms());
    // TODO seems to be often 0....
    std::vector<size_t> rexp(receptor.num_atoms());
    // TODO all scale with inv_rd should become a method of point3D for normalizing? Vector stuff?
    for (size_t index = 0; index < receptor.num_atoms(); ++index) {
      const auto receptor_type = receptor_autodock_types[index];
      /*
      * If 'ia' is a hydrogen atom, it could be a
      * RECEPTOR hydrogen-BOND DONOR,
      */
      // TODO we should create an enum for hbond types
      const auto temp          = index - range_near_atom_receptor;
      const auto from_neighbor = temp <= index ? temp : size_t{0};
      const auto to_neighbor   = std::min(index + range_near_atom_receptor, receptor.num_atoms());
      if (receptor.num_hbond(index) == 2) {
        for (size_t other_index = from_neighbor; other_index < to_neighbor; ++other_index)
          if (index != other_index) {
            /*
            * =>  NH-> or OH->
            */
            const point3D diff = difference(
                point3D{receptor.x(index), receptor.y(index), receptor.z(index)},
                point3D{receptor.x(other_index), receptor.y(other_index), receptor.z(other_index)});
            fp_type square_distance = sum_components(square(diff));
            /*
            * If ia & ib are less than 1.3 A apart -- they are covalently bonded,
            */
            if (square_distance < covalence_distance) { /*INCREASED for H-S bonds*/
              // TODO check this, ask @Davide
              if (square_distance == fp_type{0}) {
                error(
                    "While calculating an H-O or H-N bond vector, attempt to divide by zero was just prevented.");
                square_distance = std::numeric_limits<fp_type>::epsilon();
              }

              const fp_type inv_rd = fp_type{1} / sqrtf(square_distance);
              /*
              * N-H: Set exponent rexp to 2 for m/m H-atom,
              */
              if (receptor_type != autodock_ff::OA && receptor_autodock_types[other_index] != autodock_ff::SA)
                rexp[index] = 2;
              /*
              * O-H: Set exponent rexp to 4 for m/m H-atom,
              * and flag disordered hydroxyls
              */
              if (receptor_type == autodock_ff::OA ||
                  receptor_autodock_types[other_index] == autodock_ff::SA) {
                rexp[index] = 4;
                // TODO what is disorder_h
                // if (disorder_h == TRUE)
                //   disorder[ia] = TRUE;
              }
              /*
              * Normalize the vector from ib to ia, N->H or O->H...
              */
              rvector[index] = scale(diff, inv_rd);
              /*
              * First O-H/N-H H-bond-donor found; Go on to next atom,
              */
              break;
            }
          }
      } else if (receptor.num_hbond(index) == 5) { /*A2*/
        /*
         * Scan from at most, (ia-20)th m/m atom, or ia-th (if ia<20)
         *        to (ia + 5)th m/m-atom
         * determine number of atoms bonded to the oxygen
         */
        // TODO check these index_1 and _2, seems odd to me
        size_t nbond = 0, index_1 = 0, index_2 = 0;
        for (size_t other_index = from_neighbor; other_index < to_neighbor; ++other_index)
          if (index != other_index) {
            const fp_type square_distance =
                distance2(point3D{receptor.x(index), receptor.y(index), receptor.z(index)},
                          point3D{receptor.x(other_index), receptor.y(other_index), receptor.z(other_index)});

            if ((square_distance < unknown_distance_1 &&
                 (receptor_autodock_types[other_index] != autodock_ff::HD &&
                  receptor_autodock_types[other_index] != autodock_ff::H)) ||
                (square_distance < unknown_distance_2 &&
                 (receptor_autodock_types[other_index] == autodock_ff::HD ||
                  receptor_autodock_types[other_index] == autodock_ff::H))) {
              if (nbond == 2) {
                error("Found an H-bonding atom with three bonded atoms.");
              }
              if (nbond == 1) {
                nbond   = 2;
                index_2 = other_index;
              }
              if (nbond == 0) {
                nbond   = 1;
                index_1 = other_index;
              }
            }
          } /* ( ib != ia ) */
        /* if no bonds, something is wrong */
        if (nbond == 0) {
          error("Oxygen atom found with no bonded atoms.");
        }

        /* one bond: Carbonyl Oxygen O=C-X */

        if (nbond == 1) {
          /* calculate normalized carbonyl bond vector rvector[ia][] */

          rvector[index] = normalize(point3D{receptor.x(index), receptor.y(index), receptor.z(index)},
                                     point3D{receptor.x(index_1), receptor.y(index_1), receptor.z(index_1)});

          /* find a second atom (i2) bonded to carbonyl carbon (i1) */
          for (index_2 = from_neighbor; index_2 < to_neighbor; ++index_2) {
            if ((index_2 != index_1) && (index_2 != index)) {
              const fp_type square_distance =
                  distance2(point3D{receptor.x(index_1), receptor.y(index_1), receptor.z(index_1)},
                            point3D{receptor.x(index_2), receptor.y(index_2), receptor.z(index_2)});
              if ((square_distance < unknown_distance_3 &&
                   receptor_autodock_types[index_2] != autodock_ff::HD) ||
                  (square_distance < unknown_distance_4 &&
                   receptor_autodock_types[index_2] == autodock_ff::HD)) {
                /* found one */
                /* d[i] vector from carbon to second atom */
                // TODO @Davide it appears the same to me as the previous one
                // rd2 = 0.;
                // for (i = 0; i < XYZ; i++) {
                //   d[i] = index[i2][iindexoord[i1][i];
                //   rd2 += sq(d[i]);
                // }
                const point3D d =
                    normalize(point3D{receptor.x(index_2), receptor.y(index_2), receptor.z(index_2)},
                              point3D{receptor.x(index_1), receptor.y(index_1), receptor.z(index_1)});

                /* C=O cross C-X gives the lone pair plane normal */
                rvector2[index] = normalize(point3D{
                    rvector[index].y * d.z - rvector[index].z * d.y,
                    rvector[index].z * d.x - rvector[index].x * d.z,
                    rvector[index].x * d.y - rvector[index].y * d.x,
                });
              }
            }
          } /*i2-loop*/
        } /* endif nbond==1 */

        /* two bonds: Hydroxyl or Ether Oxygen X1-O-X2 */
        if (nbond == 2) {
          /* not a disordered hydroxyl */
          /* normalized X1 to X2 vector, defines lone pair plane */

          rvector2[index] = normalize(point3D{receptor.x(index_2), receptor.y(index_2), receptor.z(index_2)},
                                      point3D{receptor.x(index_1), receptor.y(index_1), receptor.z(index_1)});

          /* vector pointing between the lone pairs:
                ** front of the vector is the oxygen atom,
                ** X1->O vector dotted with normalized X1->X2 vector plus
                ** coords of X1 gives the point on the X1-X2 line for the
                ** back of the vector.
                */
          const point3D diff =
              difference(point3D{receptor.x(index), receptor.y(index), receptor.z(index)},
                         point3D{receptor.x(index_1), receptor.y(index_1), receptor.z(index_1)});
          //  TODO check this static const to const
          const fp_type rdot = inner_product(diff, static_cast<const point3D&>(rvector2[index]));

          rvector[index] = normalize(
              difference(point3D{receptor.x(index), receptor.y(index), receptor.z(index)},
                         add(scale(rvector2[index], rdot),
                             point3D{receptor.x(index_1), receptor.y(index_1), receptor.z(index_1)})));
        } /* end two bonds to Oxygen */
        /* NEW Directional N Acceptor */
      } else if (receptor.num_hbond(index) == 4) { /*A1*/
        /*
        ** Scan from at most, (ia-20)th m/m atom, or ia-th (if ia<20)
        **        to (ia+5)th m/m-atom
        ** determine number of atoms bonded to the oxygen
        */
        size_t nbond = 0, index_1 = 0, index_2 = 0, index_3 = 0;
        for (size_t other_index = from_neighbor; other_index < to_neighbor; ++other_index)
          if (index != other_index) {
            const fp_type square_distance =
                distance2(point3D{receptor.x(index), receptor.y(index), receptor.z(index)},
                          point3D{receptor.x(other_index), receptor.y(other_index), receptor.z(other_index)});

            if ((square_distance < unknown_distance_1 &&
                 (receptor_autodock_types[other_index] != autodock_ff::HD &&
                  receptor_autodock_types[other_index] != autodock_ff::H)) ||
                (square_distance < unknown_distance_2 &&
                 (receptor_autodock_types[other_index] == autodock_ff::HD ||
                  receptor_autodock_types[other_index] == autodock_ff::H))) {
              if (nbond == 2) {
                nbond   = 3;
                index_3 = other_index;
              }
              if (nbond == 1) {
                nbond   = 2;
                index_2 = other_index;
              }
              if (nbond == 0) {
                nbond   = 1;
                index_1 = other_index;
              }
            }
          } /* ( ib != ia ) */
        /* if no bonds, something is wrong */
        if (nbond == 0) {
          error("Nitrogen atom found with no bonded atoms.");
        }

        /* one bond: Azide Nitrogen :N=C-X */

        if (nbond == 1) {
          /* calculate normalized N=C bond vector rvector[ia][] */
          rvector[index] = normalize(point3D{receptor.x(index), receptor.y(index), receptor.z(index)},

                                     point3D{receptor.x(index_1), receptor.y(index_1), receptor.z(index_1)});
        } /* endif nbond==1 */

        /* two bonds: X1-N=X2 */
        if (nbond == 2) {
          /* normalized vector from Nitrogen to midpoint between X1 and X2 */
          rvector[index] = normalize(difference(
              point3D{receptor.x(index), receptor.y(index), receptor.z(index)},

              scalar_division(add(point3D{receptor.x(index_2), receptor.y(index_2), receptor.z(index_2)},
                                  point3D{receptor.x(index_1), receptor.y(index_1), receptor.z(index_1)}),
                              2)));
        } /* end two bonds for nitrogen*/
        /* three bonds: X1,X2,X3 */
        if (nbond == 3) {
          /* normalized vector from Nitrogen to midpoint between X1, X2, and X3 */
          rvector[index] = normalize(difference(
              point3D{receptor.x(index), receptor.y(index), receptor.z(index)},

              scalar_division(add(add(point3D{receptor.x(index_1), receptor.y(index_1), receptor.z(index_1)},
                                      point3D{receptor.x(index_2), receptor.y(index_2), receptor.z(index_2)}),
                                  point3D{receptor.x(index_3), receptor.y(index_3), receptor.z(index_3)}),
                              3)));
        } /* end three bonds for Nitrogen */
        /* endNEW directional N Acceptor */
      }
    }

    /*
    * Iterate over all grid points, Z( Y ( X ) ) (X is fastest)...
    */
    for (size_t index_z = 0; index_z < npts.size_z(); ++index_z) {
      /*
      *  c[0:2] contains the current grid point.
      */
      const fp_type coord_z = grid_minimum.z + index_z * grid_spacing;
      for (size_t index_y = 0; index_y < npts.size_y(); ++index_y) {
        const fp_type coord_y = grid_minimum.y + index_y * grid_spacing;
        for (size_t index_x = 0; index_x < npts.size_x(); ++index_x) {
          const fp_type coord_x = grid_minimum.x + index_x * grid_spacing;

          /* Initialize Min Hbond variables  for each new point*/
          for (auto& scracth: scratchpads) scracth.reset();

          /* NEW2: Find Closest Hbond */
          fp_type rmin    = std::numeric_limits<fp_type>::max();
          size_t closestH = 0;
          for (size_t index = 0; index < receptor.num_atoms(); ++index) {
            if (receptor.num_hbond(index) == 1 || receptor.num_hbond(index) == 2) { /*DS or D1*/
              const fp_type d = distance(point3D{receptor.x(index), receptor.y(index), receptor.z(index)},
                                         point3D{coord_x, coord_y, coord_z});
              if (d < rmin) {
                rmin     = d;
                closestH = index;
              }
            } /* Hydrogen test */
          } /* ia loop */
          /* END NEW2: Find Min Hbond */

          for (size_t index = 0; index < receptor.num_atoms(); ++index) {
            const auto receptor_type       = receptor_autodock_types[index];
            const auto& receptor_type_desc = get_description(receptor_type);
            const auto receptor_hbond      = receptor.num_hbond(index);
            point3D dist = difference(point3D{receptor.x(index), receptor.y(index), receptor.z(index)},
                                      point3D{coord_x, coord_y, coord_z});
            fp_type d    = sqrtf(sum_components(square(dist)));
            dist         = normalize(dist);

            const size_t indx_n = std::min<size_t>(std::floor(d * A_DIV), NEINT - 1);
            const size_t indx_r = std::min<size_t>(std::floor(d * A_DIV), NDIEL - 1);

            /*
            * If distance from grid point to atom ia is too large,
            * or if atom is a disordered hydrogen,
            *   add nothing to the grid-point's non-bond energy;
            *   just continue to next atom...
            */
            if (d > NBC)
              continue; /* onto the next atom... */

            fp_type racc{1};
            fp_type rdon{1};
            /* NEW2 Hramp ramps in Hbond acceptor probes */
            fp_type Hramp{1};
            /* END NEW2 Hramp ramps in Hbond acceptor probes */

            if (receptor_hbond == 2) { /*D1*/
              /*
              *  ia-th receptor atom = Hydrogen ( 4 = H )
              *  => receptor H-bond donor, OH or NH.
              *  calculate racc for H-bond ACCEPTOR PROBES at this grid pt.
              *            ====     ======================
              */

              /*
              *  d[] = Unit vector from current grid pt to ia_th m/m atom.
              *  cos_theta = d dot rvector == cos(angle) subtended.
              */
              fp_type cos_theta = -sum_components(product(dist, rvector[index]));
              if (cos_theta <= 0) {
                /*
                         *  H->current-grid-pt vector >= 90 degrees from
                         *  N->H or O->H vector,
                         */
                racc = 0.;
              } else {
                /*
                         *  racc = [cos(theta)]^2.0 for N-H
                         *  racc = [cos(theta)]^4.0 for O-H,
                         */
                switch (rexp[index]) {
                  case 1:
                  default: racc = cos_theta; break;
                  case 2: racc = cos_theta * cos_theta; break;
                  case 4:
                    fp_type tmp = cos_theta * cos_theta;
                    racc        = tmp * tmp;
                    break;
                }
                /* racc = pow( cos_theta, (double)rexp[ia]); */
                /* NEW2 calculate dot product of bond vector with bond vector of best hbond */
                if (index == closestH) {
                  Hramp = fp_type{1};
                } else {
                  cos_theta     = sum_components(product(rvector[closestH], rvector[index]));
                  cos_theta     = std::min(cos_theta, fp_type{1});
                  cos_theta     = std::max(cos_theta, fp_type{-1});
                  fp_type theta = std::acos(cos_theta);
                  Hramp         = fp_type{0.5} - fp_type{0.5} * std::cos(theta * fp_type{120} / fp_type{90});
                } /* ia test for closestH */
                /* END NEW2 calculate dot product of bond vector with bond vector of best hbond */
              }
              /* endif (atom_type[ia] == hydrogen) */
            } else if (receptor_hbond == 4) { /*A1*/
              /* NEW Directional N acceptor */
              /*
              **  ia-th macromolecule atom = Nitrogen ( 4 = H )
              **  calculate rdon for H-bond Donor PROBES at this grid pt.
              **            ====     ======================
              */

              /*
              **  d[] = Unit vector from current grid pt to ia_th m/m atom.
              **  cos_theta = d dot rvector == cos(angle) subtended.
              */
              fp_type cos_theta = -sum_components(product(dist, rvector[index]));
              if (cos_theta <= 0) {
                /*
                **  H->current-grid-pt vector >= 90 degrees from
                **  X->N vector,
                */
                rdon = 0.;
              } else {
                /*
                **  racc = [cos(theta)]^2.0 for H->N
                */
                rdon = cos_theta * cos_theta;
              }
              /* endif (atom_type[ia] == nitrogen) */
              /* end NEW Directional N acceptor */
            } else if (receptor_hbond == 5) { /*A2*/
              /*
              **  ia-th receptor atom = Oxygen
              **  => receptor H-bond acceptor, oxygen.
              */
              rdon = fp_type{0};

              /* check to see that probe is in front of oxygen, not behind */
              fp_type cos_theta = -sum_components(product(dist, rvector[index]));
              /*
                    ** t0 is the angle out of the lone pair plane, calculated
                    ** as 90 deg - acos (vector to grid point DOT lone pair
                    ** plane normal)
                    */
              fp_type t0 = sum_components(product(dist, rvector2[index]));
              t0         = std::clamp(t0, fp_type{-1}, fp_type{1});

              t0 = math::pi_halved - std::acos(t0);
              /*
              ** ti is the angle in the lone pair plane, away from the
              ** vector between the lone pairs,
              ** calculated as (grid vector CROSS lone pair plane normal)
              ** DOT C=O vector - 90 deg
              */
              const point3D cross =
                  normalize(point3D{dist.y * rvector2[index].z - dist.z * rvector2[index].y,
                                    dist.z * rvector2[index].x - dist.x * rvector2[index].z,
                                    dist.x * rvector2[index].y - dist.y * rvector2[index].x});
              fp_type ti = sum_components(product(cross, static_cast<const point3D&>(rvector[index])));

              /* rdon expressions from Goodford */
              rdon = 0.;
              if (cos_theta >= 0) {
                ti = std::clamp(ti, fp_type{-1}, fp_type{1});
                ti = std::acos(ti) - math::pi_halved;
                if (ti < 0) {
                  ti = -ti;
                }
                /* the 2.0*ti can be replaced by (ti + ti) in: rdon = (0.9 + 0.1*sin(2.0*ti))*cos(t0);*/
                rdon = (fp_type{0.9} + fp_type{0.1} * std::sin(ti + ti)) * std::cos(t0);
              } else if (cos_theta >= fp_type{-0.34202}) {
                /* 0.34202 = cos (100 deg) */
                // TODO @Davide ok here the fp_type?
                rdon = fp_type{562.25} * std::pow(fp_type{0.116978} - cos_theta * cos_theta, fp_type{3}) *
                       std::cos(t0);
              }
              /* endif atom_type == oxygen, not disordered */
            }

            /*
              * For each probe atom-type,
              * Sum pairwise interactions between each probe
              * at this grid point (c[0:2])
              * and the current receptor atom, ia...
            */
            for (scratchpad& scratch: scratchpads) {
              grid_atom_map& atom_map    = scratch.atom_map;
              const auto& grid_type_desc = get_description(atom_map.get_atom_type());
              const auto& inter          = scratch.get_interaction(receptor_type);
              const auto& vdw_hb_value   = inter.vdw_hb_table[indx_n];

              if (scratch.atom_map.is_hbonder) {
                fp_type rsph = vdw_hb_value / fp_type{100};
                rsph         = std::clamp(rsph, fp_type{0}, fp_type{1});
                if ((grid_type_desc.hbond == 3 || grid_type_desc.hbond == 5) /*AS or A2*/
                    && (receptor_hbond == 1 || receptor_hbond == 2)) {       /*DS or D1*/
                  scratch.energy += vdw_hb_value * Hramp * (racc + (fp_type{1} - racc) * rsph);
                } else if ((grid_type_desc.hbond == 4)                        /*A1*/
                           && (receptor_hbond == 1 || receptor_hbond == 2)) { /*DS,D1*/
                  scratch.hbondmin =
                      std::min(scratch.hbondmin, vdw_hb_value * (racc + (fp_type{1} - racc) * rsph));
                  scratch.hbondmax =
                      std::max(scratch.hbondmax, vdw_hb_value * (racc + (fp_type{1} - racc) * rsph));
                  scratch.hbondflag = true;
                } else if ((grid_type_desc.hbond == 1 || grid_type_desc.hbond == 2) &&
                           (receptor_hbond > 2)) { /*DS,D1 vs AS,A1,A2*/
                  /*  PROBE is H-BOND DONOR, */
                  fp_type temp_hbond_enrg = vdw_hb_value * (rdon + (fp_type{1} - rdon) * rsph);
                  scratch.hbondmin        = std::min(scratch.hbondmin, temp_hbond_enrg);
                  scratch.hbondmax        = std::max(scratch.hbondmax, temp_hbond_enrg);
                  scratch.hbondflag       = true;
                } else { /*end of is_hbonder*/
                  /*  hbonder PROBE-ia cannot form a H-bond..., */
                  scratch.energy += vdw_hb_value;
                }
              } else {
                /*  PROBE does not form H-bonds..., */
                scratch.energy += vdw_hb_value;
              } /* end hbonder tests */

              /* add desolvation energy  */
              /* forcefield desolv coefficient/weight in sol_fn*/
              // TODO Check charges they seems to be different
              scratch.energy +=
                  scratch.grid_type_desc.solpar * get_description(receptor_type).vol * sol_fn[indx_r] +
                  (receptor_type_desc.solpar + solpar_q * std::fabs(receptor.charge(index))) *
                      scratch.grid_type_desc.vol * sol_fn[indx_r];
            }
          } /* ia loop, over all receptor atoms... */

          /* adjust maps of hydrogen-bonding atoms by adding largest and
          * smallest interaction of all 'pair-wise' interactions with receptor atoms
          */
          for (scratchpad& scratch: scratchpads) { scratch.write_to_map(index_x, index_y, index_z); }
        }
      }
    }
    // Allocate the results grid maps
    // TODO check if it works with the move
    std::vector<grid_atom_map> grid_atom_maps;
    for (auto& scratch: scratchpads) grid_atom_maps.push_back(std::move(scratch.atom_map));
    return grid_atom_maps;
  }
} // namespace mudock
