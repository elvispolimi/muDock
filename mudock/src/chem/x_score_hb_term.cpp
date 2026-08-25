#include <algorithm>
#include <cmath>
#include <concepts>
#include <mudock/chem/x_score_hb_term.hpp>
#include <mudock/chem/x_score_validity.hpp>
#include <mudock/chem/x_score_xtool_types.hpp>
#include <string_view>
#include <type_traits>
#include <vector>

namespace mudock {

  namespace {

    // geometry helpers

    // angle between two vectors
    fp_type angle_of_two_vectors(const fp_type v1[3], const fp_type v2[3]) {
      const double l1 = std::sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
      const double l2 = std::sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
      const double dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
      double c         = dot / (l1 * l2);
      if (c > 1.0)
        c = 1.0; 
      else if (c < -1.0)
        c = -1.0;
      return static_cast<fp_type>(std::acos(c) / M_PI * 180.0);
    }

    // angle a-b-c
    fp_type angle_abc(const fp_type a[3], const fp_type b[3], const fp_type c[3]) {
      const fp_type v1[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
      const fp_type v2[3] = {b[0] - c[0], b[1] - c[1], b[2] - c[2]};
      return angle_of_two_vectors(v1, v2);
    }

    // ramp
    fp_type ramp(const fp_type v, const fp_type a1, const fp_type a2, const fp_type a3, const fp_type a4) {
      if (v < a1)
        return fp_type{0};
      if (v < a2)
        return (v - a1) / (a2 - a1);
      if (v < a3)
        return fp_type{1};
      if (v < a4)
        return (a4 - v) / (a4 - a3);
      return fp_type{0};
    }

    // chem helpers
    
    bool is_hydrogen(const xtool_ff t) {
      return t == xtool_ff::H || t == xtool_ff::Hhb || t == xtool_ff::Hg;
    }
    
    bool is_N4(const xtool_ff t) { return t == xtool_ff::N4; }
    bool is_N3(const xtool_ff t) { return t == xtool_ff::N3 || t == xtool_ff::N3h || t == xtool_ff::N3un; }
    bool is_Npl3(const xtool_ff t) {
      return t == xtool_ff::Npl3 || t == xtool_ff::Npl3h || t == xtool_ff::Npl3un;
    }
    bool is_N2(const xtool_ff t) { return t == xtool_ff::N2 || t == xtool_ff::N2h || t == xtool_ff::N2un; }
    bool is_Nar(const xtool_ff t) { return t == xtool_ff::Nar || t == xtool_ff::Narh || t == xtool_ff::Narun; }

    char element_of(const xtool_ff basic) {
      switch (basic) {
        case xtool_ff::O3:
        case xtool_ff::O3h:
        case xtool_ff::O3un:
        case xtool_ff::O2:
        case xtool_ff::O2un:
        case xtool_ff::Oco2:
        case xtool_ff::Ow: return 'O';
        case xtool_ff::N4:
        case xtool_ff::N3:
        case xtool_ff::N3h:
        case xtool_ff::N3un:
        case xtool_ff::Npl3:
        case xtool_ff::Npl3h:
        case xtool_ff::Npl3un:
        case xtool_ff::N2:
        case xtool_ff::N2h:
        case xtool_ff::N2un:
        case xtool_ff::Nar:
        case xtool_ff::Narh:
        case xtool_ff::Narun:
        case xtool_ff::N1:
        case xtool_ff::N1un:
        case xtool_ff::Nam: return 'N';
        case xtool_ff::S3:
        case xtool_ff::S3h:
        case xtool_ff::S3un:
        case xtool_ff::S2:
        case xtool_ff::S2un:
        case xtool_ff::So:
        case xtool_ff::So2: return 'S';
        default: return 'X';
      }
    }

    
    int get_donor_type(const int origin,
                       const xtool_ff basic,
                       const xtool_ff xtype,
                       const int num_nonh,
                       const x_score_hb hb,
                       const std::string_view name,
                       const residue res) {
      if (hb != x_score_hb::D && hb != x_score_hb::DA && hb != x_score_hb::M)
        return 0;

      if (origin == 2) { // protein
        if (hb == x_score_hb::M)
          return 1;
        if (basic == xtool_ff::O3 || basic == xtool_ff::Oco2 || basic == xtool_ff::Ow)
          return 2;
        if (is_N4(xtype) || is_N3(xtype))
          return 2;
        if (is_Npl3(xtype) || is_N2(xtype) || is_Nar(xtype)) {
          if (name.find("NH") != std::string_view::npos && res == residue::ARG)
            return 2;
          if (name.find("ND") != std::string_view::npos && res == residue::ASN)
            return 2;
          if (name.find("NE") != std::string_view::npos && res == residue::GLN)
            return 2;
          return 1;
        }
        if (basic == xtool_ff::S3)
          return 2;
        return 2;
      }

      // ligand
      if (basic == xtool_ff::O3 || basic == xtool_ff::Oco2)
        return 2;
      if (is_N4(xtype) || is_N3(xtype))
        return (num_nonh <= 2) ? 2 : 1;
      if (is_Npl3(xtype) || is_N2(xtype) || is_Nar(xtype))
        return (num_nonh <= 1) ? 2 : 1;
      if (basic == xtool_ff::S3)
        return 2;
      return 2;
    }

    
    int get_acceptor_type(const int origin,
                          const xtool_ff basic,
                          const xtool_ff xtype,
                          const int num_nonh,
                          const x_score_hb hb) {
      if (hb != x_score_hb::A && hb != x_score_hb::DA)
        return 0;

      if (origin == 2) { // protein
        if (basic == xtool_ff::O3 || basic == xtool_ff::O2 || basic == xtool_ff::Oco2 ||
            basic == xtool_ff::Ow)
          return 2;
        if (is_Npl3(xtype) || is_N2(xtype) || is_Nar(xtype))
          return 1;
        return 2;
      }

      // ligand
      if (basic == xtool_ff::O3 || basic == xtool_ff::O2 || basic == xtool_ff::Oco2)
        return 2;
      if (is_Npl3(xtype) || is_N2(xtype) || is_Nar(xtype))
        return (num_nonh <= 1) ? 2 : 1;
      return 2;
    }

    
    template<typename layer_t>
    std::vector<std::vector<int>> build_adjacency(const layer_t& layer, const int num_atoms) {
      const auto& mol  = layer.get_base_molecule();
      const auto bonds = mol.get_bonds();
      std::vector<std::vector<int>> adj(num_atoms);
      for (const auto& b: bonds) {
        adj[b.source].push_back(b.dest);
        adj[b.dest].push_back(b.source);
      }
      return adj;
    }

    
    bool hb_eligible(const x_score_hb hb, const bool is_protein) {
      if (hb == x_score_hb::D || hb == x_score_hb::A || hb == x_score_hb::DA)
        return true;
      return is_protein && hb == x_score_hb::M;
    }

    int donor_limit_of(const x_score_hb hb, const int num_h) {
      return (hb == x_score_hb::DA) ? 1 : num_h;
    }

    int acceptor_limit_of(const xtool_ff basic) { return (element_of(basic) == 'N') ? 1 : 2; }

    // Distinguishes ligand vs protein for build_hb_atoms behavior
    template<typename layer_t>
    constexpr bool is_protein_layer_v = std::same_as<std::remove_cvref_t<layer_t>, x_score_protein>;

    template<typename layer_t>
    x_score_hb_atoms build_hb_atoms(const layer_t& layer) {
      constexpr bool is_protein = is_protein_layer_v<layer_t>;
      constexpr int origin      = is_protein ? 2 : 1;

      const auto& mol     = layer.get_base_molecule();
      const int num_atoms = static_cast<int>(mol.num_atoms());
      const auto adj      = build_adjacency(layer, num_atoms);

      x_score_hb_atoms atoms;
      atoms.x.reserve(num_atoms);
      atoms.y.reserve(num_atoms);
      atoms.z.reserve(num_atoms);
      atoms.root_x.reserve(num_atoms);
      atoms.root_y.reserve(num_atoms);
      atoms.root_z.reserve(num_atoms);
      atoms.radius.reserve(num_atoms);
      atoms.hb.reserve(num_atoms);
      atoms.donor_type.reserve(num_atoms);
      atoms.acceptor_type.reserve(num_atoms);
      atoms.donor_limit.reserve(num_atoms);
      atoms.acceptor_limit.reserve(num_atoms);
      atoms.has_root.reserve(num_atoms);

      // filtering out atoms not useful for Hydrogen Bonding

      for (int i = 0; i < num_atoms; ++i) {
        if (layer.valid(i) == x_score_validity::invalid)
          continue;
        const xtool_ff basic = mol.atom_type(i);
        if (is_hydrogen(basic))
          continue;
        if (basic == xtool_ff::Ow)
          continue;
        const x_score_hb hb = layer.hb(i);
        if (!hb_eligible(hb, is_protein))
          continue;

        fp_type rx = 0, ry = 0, rz = 0;
        int num_nonh = 0, num_h = 0;
        for (const int nb: adj[i]) {
          if (is_hydrogen(mol.atom_type(nb))) {
            ++num_h;
          } else if (layer.hb(nb) == x_score_hb::M) {
            continue;
          } else {
            rx += mol.x(nb);
            ry += mol.y(nb);
            rz += mol.z(nb);
            ++num_nonh;
          }
        }

        // Ligand-only
        if constexpr (!is_protein) {
          if (num_nonh == 0)
            continue;
        }

        atoms.x.push_back(mol.x(i));
        atoms.y.push_back(mol.y(i));
        atoms.z.push_back(mol.z(i));
        atoms.radius.push_back(layer.vdw_radius(i));
        atoms.hb.push_back(hb);

        if (num_nonh > 0) {
          const fp_type n = static_cast<fp_type>(num_nonh);
          atoms.root_x.push_back(rx / n);
          atoms.root_y.push_back(ry / n);
          atoms.root_z.push_back(rz / n);
          atoms.has_root.push_back(1);
        } else {
          atoms.root_x.push_back(mol.x(i));
          atoms.root_y.push_back(mol.y(i));
          atoms.root_z.push_back(mol.z(i));
          atoms.has_root.push_back(0);
        }

        const xtool_ff xtype = layer.x_score_xtool_type(i);
        std::string_view name;
        residue res = residue::UNKNOWN;
        if constexpr (is_protein) {
          name = mol.atom_name(i);
          res  = mol.residue_types(i);
        }
        atoms.donor_type.push_back(get_donor_type(origin, basic, xtype, num_nonh, hb, name, res));
        atoms.acceptor_type.push_back(get_acceptor_type(origin, basic, xtype, num_nonh, hb));

        atoms.donor_limit.push_back(donor_limit_of(hb, num_h));
        atoms.acceptor_limit.push_back(acceptor_limit_of(basic));
      }

      return atoms;
    }


    fp_type value_hbond_2(const x_score_hb_atoms& donor,
                          const int di,
                          const x_score_hb_atoms& acceptor,
                          const int ai) {
      const fp_type dco[3] = {donor.x[di], donor.y[di], donor.z[di]};
      const fp_type aco[3] = {acceptor.x[ai], acceptor.y[ai], acceptor.z[ai]};
      const fp_type dro[3] = {donor.root_x[di], donor.root_y[di], donor.root_z[di]};
      const fp_type aro[3] = {acceptor.root_x[ai], acceptor.root_y[ai], acceptor.root_z[ai]};

      const fp_type dx = dco[0] - aco[0], dy = dco[1] - aco[1], dz = dco[2] - aco[2];
      const fp_type d  = std::sqrt(dx * dx + dy * dy + dz * dz);

      const bool donor_is_metal    = donor.hb[di] == x_score_hb::M;
      const bool acceptor_is_metal = acceptor.hb[ai] == x_score_hb::M;

      // DR-D-A angle
      bool mark1 = false;
      fp_type a1 = 0;
      if (!donor_is_metal && donor.has_root[di]) {
        a1    = fp_type{180} - angle_abc(dro, dco, aco);
        mark1 = true;
      }

      // D-A-AR angle
      bool mark2 = false;
      fp_type a2 = 0;
      if (!acceptor_is_metal && acceptor.has_root[ai]) {
        a2    = fp_type{180} - angle_abc(dco, aco, aro);
        mark2 = true;
      }

      const int d_type = donor.donor_type[di];
      const int a_type = acceptor.acceptor_type[ai];

      const fp_type d0 = donor.radius[di] + acceptor.radius[ai];
      const fp_type d1 = 0, d2 = 1, d3 = d0 - fp_type{0.4f}, d4 = d0 + fp_type{0.2f};

      fp_type tmp1;
      if (d < d1)
        tmp1 = 0;
      else if (d < d2)
        tmp1 = (d - d1) / (d2 - d1);
      else if (d < d3)
        tmp1 = 1;
      else if (d < d4)
        tmp1 = (d4 - d) / (d4 - d3);
      else
        tmp1 = 0;

      fp_type tmp3 = 1;
      if (mark1) {
        if (d_type == 1)
          tmp3 = ramp(a1, 0, fp_type{0.001f}, 25, 50);
        else
          tmp3 = ramp(a1, 25, 50, 75, 100);
      }

      fp_type tmp4 = 1;
      if (mark2) {
        if (a_type == 1)
          tmp4 = ramp(a2, 0, fp_type{0.001f}, 30, 55);
        else
          tmp4 = ramp(a2, 0, 5, 70, 95);
      }

      return (tmp3 >= tmp4) ? tmp1 * tmp4 : tmp1 * tmp3;
    }

  } // namespace


  x_score_hb_atoms build_protein_hb_atoms(const x_score_protein& prot) { return build_hb_atoms(prot); }

  x_score_hb_atoms build_ligand_hb_atoms(const x_score_ligand& lig) { return build_hb_atoms(lig); }


  std::vector<x_score_hb_candidate> get_hbond_pair_pl(const x_score_hb_atoms& lig_atoms,
                                                      const x_score_hb_atoms& prot_atoms) {
    constexpr fp_type cutoff = fp_type{5.0};

    std::vector<x_score_hb_candidate> candidates;
    const int num_lig  = lig_atoms.size();
    const int num_prot = prot_atoms.size();

    for (int li = 0; li < num_lig; ++li) {
      const x_score_hb lhb = lig_atoms.hb[li];
      const fp_type lx = lig_atoms.x[li], ly = lig_atoms.y[li], lz = lig_atoms.z[li];

      for (int pj = 0; pj < num_prot; ++pj) {
        const x_score_hb phb = prot_atoms.hb[pj];

        // determine the H-bond type from the donor/acceptor characters
        int type = 0;
        if (lhb == x_score_hb::D) {
          if (phb == x_score_hb::A || phb == x_score_hb::DA)
            type = 1;
        } else if (lhb == x_score_hb::A) {
          if (phb == x_score_hb::D || phb == x_score_hb::DA)
            type = 2;
          else if (phb == x_score_hb::M)
            type = 3;
        } else if (lhb == x_score_hb::DA) {
          if (phb == x_score_hb::A || phb == x_score_hb::DA)
            type = 1;
          else if (phb == x_score_hb::D)
            type = 2;
          else if (phb == x_score_hb::M)
            type = 3;
        }
        if (type == 0)
          continue;

        const fp_type dx = lx - prot_atoms.x[pj];
        const fp_type dy = ly - prot_atoms.y[pj];
        const fp_type dz = lz - prot_atoms.z[pj];
        const fp_type d  = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (d > cutoff)
          continue;

        // donor/acceptor assignment then geometric strength
        fp_type score;
        if (type == 1)
          score = value_hbond_2(lig_atoms, li, prot_atoms, pj); // ligand donor, protein acceptor
        else
          score = value_hbond_2(prot_atoms, pj, lig_atoms, li); // protein donor (or metal), ligand acceptor

        if (std::fabs(score) > fp_type{0})
          candidates.push_back({li, pj, type, score});
      }
    }

    return candidates;
  }


  fp_type sum_hbonds(std::vector<x_score_hb_candidate>& candidates,
                     const x_score_hb_atoms& lig_atoms,
                     const x_score_hb_atoms& prot_atoms) {
    const int n = static_cast<int>(candidates.size());

    // Step 1: decreasing |score|
    std::stable_sort(candidates.begin(),
                     candidates.end(),
                     [](const x_score_hb_candidate& a, const x_score_hb_candidate& b) {
                       return std::fabs(a.score) > std::fabs(b.score);
                     });

    // Step 2: two H-bonds on the same ligand atom must be at least 45 degrees apart
    for (int i = 0; i < n - 1; ++i)
      for (int j = i + 1; j < n; ++j) {
        if (candidates[i].li != candidates[j].li)
          continue;
        const int li = candidates[i].li, pi = candidates[i].pj;
        const int lj = candidates[j].li, pj = candidates[j].pj;
        const fp_type v1[3] = {prot_atoms.x[pi] - lig_atoms.x[li],
                               prot_atoms.y[pi] - lig_atoms.y[li],
                               prot_atoms.z[pi] - lig_atoms.z[li]};
        const fp_type v2[3] = {prot_atoms.x[pj] - lig_atoms.x[lj],
                               prot_atoms.y[pj] - lig_atoms.y[lj],
                               prot_atoms.z[pj] - lig_atoms.z[lj]};
        const fp_type angle = std::fabs(angle_of_two_vectors(v1, v2));
        if (angle < fp_type{45})
          candidates[j].score = 0;
      }

    // Step 3: a donor ligand atom forms no more H-bonds than the hydrogens it carries
    for (int i = 0; i < n - 1; ++i) {
      if (candidates[i].type != 1)
        continue;
      int count       = 1;
      const int limit = lig_atoms.donor_limit[candidates[i].li];
      for (int j = i + 1; j < n; ++j) {
        if (candidates[j].type != 1 || candidates[i].li != candidates[j].li)
          continue;
        ++count;
        if (count > limit)
          candidates[j].score = 0;
      }
    }

    // Step 4: an acceptor ligand atom forms no more H-bonds than its lone pairs
    for (int i = 0; i < n - 1; ++i) {
      if (candidates[i].type != 2 && candidates[i].type != 3)
        continue;
      int count       = 1;
      const int limit = lig_atoms.acceptor_limit[candidates[i].li];
      for (int j = i + 1; j < n; ++j) {
        if ((candidates[j].type != 2 && candidates[j].type != 3) ||
            candidates[i].li != candidates[j].li)
          continue;
        ++count;
        if (count > limit)
          candidates[j].score = 0;
      }
    }

    // sum the survived contributions
    fp_type sum = 0;
    for (const auto& c: candidates)
      if (std::fabs(c.score) >= fp_type{0.01f})
        sum += c.score;

    return sum;
  }

  fp_type compute_x_score_hb(const x_score_hb_atoms& lig_atoms, const x_score_hb_atoms& prot_atoms) {
    auto candidates = get_hbond_pair_pl(lig_atoms, prot_atoms);
    return sum_hbonds(candidates, lig_atoms, prot_atoms);
  }

} // namespace mudock