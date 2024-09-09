#pragma once

#include "mudock/chem/autodock_types.hpp"

#include <cmath>
#include <mudock/chem.hpp>
#include <mudock/type_alias.hpp>
#include <set>
#include <utility>
#include <vector>

namespace mudock {
  inline fp_type calc_ddd_Mehler_Solmajer(fp_type distance) {
    /*____________________________________________________________________________
     * Distance-dependent dielectric ewds: Mehler and Solmajer, Prot Eng 4, 903-910.
     *____________________________________________________________________________*/
    const fp_type lambda{0.003627};
    const fp_type epsilon0{78.4};
    const fp_type A{-8.5525};
    const fp_type B = epsilon0 - A;
    const fp_type rk{7.7839};
    const fp_type lambda_B = -lambda * B;

    fp_type epsilon = A + B / (fp_type{1} + rk * std::exp(lambda_B * distance));

    if (epsilon < std::numeric_limits<fp_type>::epsilon()) {
      epsilon = 1.0L;
    }
    return epsilon;
  }

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
          size_t j        = temp <= size_t{0} ? size_t{0} : temp;
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
          vdw_hb_table(cA, cB, xA, xB){};
  };

  // Custom hash function for std::pair
  struct pair_hash {
    template<class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& p) const {
      auto hash1 = std::hash<T1>{}(p.first);
      auto hash2 = std::hash<T2>{}(p.second);
      return hash1 ^ (hash2 << 1); // Combine the two hash values
    }
  };

  struct precompute {
    std::array<fp_type, NDIEL> sol_fn;
    std::array<fp_type, NDIEL> epsilon_fn;
    std::array<fp_type, NDIEL> r_epsilon_fn;

    precompute(const std::set<autodock_ff> types) {
      // Initialize data structure
      for (auto& ligand_t: types) {
        const autodock_ff_description& grid_type_desc = get_description(ligand_t);
        for (auto& receptor_t: types) {
          const auto& receptor_type_desc = get_description(receptor_t);

          fp_type nbp_r = (grid_type_desc.Rii + receptor_type_desc.Rii) / fp_type{2.};
          // TODO check if it is ok to use the floating point version
          fp_type nbp_eps = std::sqrt(grid_type_desc.epsii * autodock_parameters::coeff_vdW *
                                      receptor_type_desc.epsii * autodock_parameters::coeff_vdW);
          // TODO probably they are constant
          size_t xA      = 12;
          size_t xB      = 6;
          size_t hbonder = 0;
          if (grid_type_desc.hbond > 2 && (receptor_type_desc.hbond == 1 ||
                                           receptor_type_desc.hbond == 2)) { /*AS,A1,A2 map vs DS,D1 probe*/
            xB      = 10;
            hbonder = 1;
            nbp_r   = grid_type_desc.Rij_hb;
            nbp_eps = grid_type_desc.epsij_hb * autodock_parameters::coeff_hbond;
          } else if ((grid_type_desc.hbond == 1 || grid_type_desc.hbond == 2) &&
                     (receptor_type_desc.hbond > 2)) { /*DS,D1 map vs AS,A1,A2 probe*/
            xB      = 10;
            hbonder = 1;
            nbp_r   = receptor_type_desc.Rij_hb;
            nbp_eps = receptor_type_desc.epsij_hb * autodock_parameters::coeff_hbond;
          }; /*initialize energy parms for each possible receptor type*/
          // TODO to be checked later in line 1707 from main.cpp
          const fp_type cA = (nbp_eps / (xA - xB)) * std::pow(nbp_r, static_cast<fp_type>(xA)) * xB;
          const fp_type cB = nbp_eps / (xA - xB) * std::pow(nbp_r, static_cast<fp_type>(xB)) * xA;
          interactions.emplace(receptor_t, interaction{cA, cB, nbp_r, nbp_eps, xA, xB, hbonder, receptor_t});
        }
      }

      /* exponential function for receptor and ligand desolvation */
      /* note: the solvation term ranges beyond the non-bond cutoff 
      * and will not be smoothed 
      */
      for (size_t indx_r = 1; indx_r < NDIEL; ++indx_r) {
        const fp_type r = indx_r / A_DIV;
        sol_fn[indx_r] =
            autodock_parameters::coeff_desolv * std::exp(-(r * r) / (fp_type{2} * (sigma * sigma)));
      }

      epsilon_fn[0] = 1.0;
      for (size_t indx_r = 1; indx_r < NDIEL; ++indx_r) {
        epsilon_fn[indx_r] = calc_ddd_Mehler_Solmajer(indx_r / A_DIV);
      }
      /* convert epsilon to factor / epsilon */
      for (size_t i = 0; i < NDIEL; ++i) { r_epsilon_fn[i] = factor / epsilon_fn[i]; }
    }

    [[nodiscard]] inline auto& get_interaction(const autodock_ff l_t, const autodock_ff r_t) {
      return interactions.at(std::pair{l_t, r_t});
    }

  private:
    // TODO @Davide
    std::unordered_map<std::pair<autodock_ff, autodock_ff>, interaction, pair_hash> interactions;
  };
} // namespace mudock
