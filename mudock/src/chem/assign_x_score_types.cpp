#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include <cctype>
#include <mudock/chem/assign_x_score_types.hpp>
#include <mudock/chem/sybyl_atom_types.hpp>
#include <mudock/chem/x_score_xtool_types.hpp>
#include <mudock/molecule/graph.hpp>
#include <string>
#include <vector>

namespace mudock {

  // Assigns specific X-Score atom types
  // 
  // There are two distinct specializations:
  // 
  // 1. Ligands (x_score_static_layer): 
  //    Performed on the molecule's graph. Detects 5- and 6-membered 
  //    aromatic rings and assigns the correct X-Score type (x_score_xtool_types)
  // 
  // 2. Proteins (x_score_dynamic_layer): 
  //    Uses a dictionary-based approach. It looks up the residue and atom name in the dedicated dictionary (x_score_residue_xtool_types)
  //    to assign predefined parameters.
  //    It also includes a workaround to normalize 
  //    non-standard hydrogen names often found in older PDB files.


  static bool is_hydrogen_type(xtool_ff t) { return t == xtool_ff::H || t == xtool_ff::Hhb; }

  static bool is_heteroatom_type(xtool_ff t) {
    switch (t) {
      case xtool_ff::F:
      case xtool_ff::Cl:
      case xtool_ff::Br:
      case xtool_ff::I:
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
      case xtool_ff::Nam:
      case xtool_ff::O3:
      case xtool_ff::O3h:
      case xtool_ff::O3un:
      case xtool_ff::O2:
      case xtool_ff::O2un:
      case xtool_ff::Oco2:
      case xtool_ff::Ow:
      case xtool_ff::P3:
      case xtool_ff::S3:
      case xtool_ff::S3h:
      case xtool_ff::S3un:
      case xtool_ff::S2:
      case xtool_ff::S2un:
      case xtool_ff::So:
      case xtool_ff::So2: return true;
      default: return false;
    }
  }

  static bool is_oxygen_or_nitrogen_type(xtool_ff t) {
    switch (t) {
      case xtool_ff::O3:
      case xtool_ff::O3h:
      case xtool_ff::O3un:
      case xtool_ff::O2:
      case xtool_ff::O2un:
      case xtool_ff::Oco2:
      case xtool_ff::Ow:
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
      case xtool_ff::Nam: return true;
      default: return false;
    }
  }


  struct atom_environment {
    int num_h      = 0;
    int num_nonh   = 0;
    int num_hetero = 0;
    bool is_bonded_to_ON = false;
  };


  static void assign_xtool_types_from_sybyl(x_score_static_layer& layer) {
    auto& mol           = layer.get_base_molecule();
    const int num_atoms = static_cast<int>(mol.num_atoms());

    for (int index = 0; index < num_atoms; ++index)
      mol.atom_type(index) = xtool_type_from_sybyl(mol.sybyl_type(index));
  }


  static bool is_6ring_pi_atom(xtool_ff t) {
    return t == xtool_ff::Car || t == xtool_ff::Nar || t == xtool_ff::C2 || t == xtool_ff::N2;
  }

  static int get_5ring_pi_count(xtool_ff t) {
    switch (t) {
      case xtool_ff::Car:
      case xtool_ff::Nar:
      case xtool_ff::C2:
      case xtool_ff::N2: return 1;
      case xtool_ff::Npl3:
      case xtool_ff::Nam:
      case xtool_ff::O3:
      case xtool_ff::O3h:
      case xtool_ff::S3:
      case xtool_ff::S3h: return 2;
      default: return 0;
    }
  }

  static bool is_single_or_amide_bond(bond_type bt) {
    return bt == bond_type::SINGLE || bt == bond_type::AMIDE;
  }

  static std::vector<bool> detect_xscore_aromaticity(const molecule_graph_type& graph,
                                                     const std::span<const xtool_ff> atom_types,
                                                     const std::span<const bond> bonds,
                                                     const std::size_t num_atoms) {
    std::vector<bool> xscore_aromatic(num_atoms, false);

    std::vector<int> atom_ring(num_atoms, -1);

    for (std::size_t i = 0; i < num_atoms; ++i) {
      if (boost::degree(i, graph) <= 1) {
        atom_ring[i] = 0;
      }
    }

    bool changed = true;
    while (changed) {
      changed = false;
      for (std::size_t i = 0; i < num_atoms; ++i) {
        if (atom_ring[i] != -1)
          continue;

        int ring_neighbor_count = 0;
        auto [vi, vi_end]       = boost::adjacent_vertices(i, graph);
        for (; vi != vi_end; ++vi) {
          if (atom_ring[*vi] != 0) {
            ring_neighbor_count++;
          }
        }
        if (ring_neighbor_count <= 1) {
          atom_ring[i] = 0;
          changed      = true;
        }
      }
    }

    for (std::size_t i = 0; i < num_atoms; ++i) {
      if (atom_ring[i] == -1)
        atom_ring[i] = 1;
    }

    auto find_ring_of_size = [&](int src, int dst, int target_size) -> std::vector<int> {
      struct dfs_state {
        int atom;
        std::vector<int> path;
      };

      std::vector<dfs_state> stack;
      stack.push_back({dst, {src, dst}});

      while (!stack.empty()) {
        auto [current, path] = std::move(stack.back());
        stack.pop_back();

        if (static_cast<int>(path.size()) > target_size)
          continue;

        auto [ni, ni_end] = boost::adjacent_vertices(current, graph);
        for (; ni != ni_end; ++ni) {
          int next = static_cast<int>(*ni);

          if (next == src && static_cast<int>(path.size()) == target_size) {
            return path;
          }

          if (atom_ring[next] != 1)
            continue;

          if (std::find(path.begin(), path.end(), next) != path.end())
            continue;

          if (next == src)
            continue;

          if (static_cast<int>(path.size()) < target_size) {
            auto new_path = path;
            new_path.push_back(next);
            stack.push_back({next, std::move(new_path)});
          }
        }
      }
      return {};
    };

    auto get_ring_bond_types = [&](const std::vector<int>& ring_atoms) -> std::vector<bond_type> {
      std::vector<bond_type> ring_bonds;
      int n = static_cast<int>(ring_atoms.size());
      for (int j = 0; j < n; ++j) {
        int a1 = ring_atoms[j];
        int a2 = ring_atoms[(j + 1) % n];
        auto [ei, ei_end] = boost::out_edges(a1, graph);
        for (; ei != ei_end; ++ei) {
          int target = static_cast<int>(boost::target(*ei, graph));
          if (target == a2) {
            ring_bonds.push_back(bonds[graph[*ei].bond_index].type);
            break;
          }
        }
      }
      return ring_bonds;
    };

    auto make_bond_key = [](int a, int b) -> long long {
      if (a > b)
        std::swap(a, b);
      return static_cast<long long>(a) * 1000000 + b;
    };

    {
      std::vector<long long> checked_bonds;
      auto ei_pair = boost::edges(graph);
      for (auto ei = ei_pair.first; ei != ei_pair.second; ++ei) {
        int a1 = static_cast<int>(boost::source(*ei, graph));
        int a2 = static_cast<int>(boost::target(*ei, graph));

        if (atom_ring[a1] != 1 || atom_ring[a2] != 1)
          continue;

        long long bk = make_bond_key(a1, a2);
        if (std::find(checked_bonds.begin(), checked_bonds.end(), bk) != checked_bonds.end())
          continue;
        checked_bonds.push_back(bk);

        auto ring_atoms = find_ring_of_size(a1, a2, 6);
        if (ring_atoms.empty())
          continue;

        bool all_pi  = true;
        int pi_count = 0;
        for (int atom_idx: ring_atoms) {
          if (is_6ring_pi_atom(atom_types[atom_idx])) {
            pi_count++;
          } else {
            all_pi = false;
            break;
          }
        }
        if (!all_pi || pi_count != 6)
          continue;

        auto ring_bond_types = get_ring_bond_types(ring_atoms);
        if (static_cast<int>(ring_bond_types.size()) != 6)
          continue;

        bool has_consecutive_single = false;
        for (int j = 0; j < 5; ++j) {
          if (is_single_or_amide_bond(ring_bond_types[j]) &&
              is_single_or_amide_bond(ring_bond_types[j + 1])) {
            has_consecutive_single = true;
            break;
          }
        }
        if (!has_consecutive_single) {
          if (is_single_or_amide_bond(ring_bond_types[5]) && is_single_or_amide_bond(ring_bond_types[0])) {
            has_consecutive_single = true;
          }
        }
        if (has_consecutive_single)
          continue;

        for (int atom_idx: ring_atoms) { xscore_aromatic[atom_idx] = true; }
      }
    }

    {
      std::vector<long long> checked_bonds;
      auto ei_pair = boost::edges(graph);
      for (auto ei = ei_pair.first; ei != ei_pair.second; ++ei) {
        int a1 = static_cast<int>(boost::source(*ei, graph));
        int a2 = static_cast<int>(boost::target(*ei, graph));

        if (atom_ring[a1] != 1 || atom_ring[a2] != 1)
          continue;
        if (xscore_aromatic[a1] && xscore_aromatic[a2])
          continue;

        long long bk = make_bond_key(a1, a2);
        if (std::find(checked_bonds.begin(), checked_bonds.end(), bk) != checked_bonds.end())
          continue;
        checked_bonds.push_back(bk);

        auto ring_atoms = find_ring_of_size(a1, a2, 5);
        if (ring_atoms.empty())
          continue;

        bool valid   = true;
        int total_pi = 0;
        std::vector<int> pi_per_atom(5);
        for (int j = 0; j < 5; ++j) {
          int pi = get_5ring_pi_count(atom_types[ring_atoms[j]]);
          if (pi == 0) {
            valid = false;
            break;
          }
          pi_per_atom[j] = pi;
          total_pi += pi;
        }
        if (!valid || total_pi != 6)
          continue;

        auto ring_bond_types = get_ring_bond_types(ring_atoms);
        if (static_cast<int>(ring_bond_types.size()) != 5)
          continue;

        bool has_bad_consecutive = false;
        for (int j = 0; j < 4; ++j) {
          if (is_single_or_amide_bond(ring_bond_types[j]) &&
              is_single_or_amide_bond(ring_bond_types[j + 1])) {
            if (pi_per_atom[j + 1] == 2)
              continue;
            has_bad_consecutive = true;
            break;
          }
        }
        if (!has_bad_consecutive) {
          if (is_single_or_amide_bond(ring_bond_types[4]) && is_single_or_amide_bond(ring_bond_types[0])) {
            if (pi_per_atom[0] != 2) {
              has_bad_consecutive = true;
            }
          }
        }
        if (has_bad_consecutive)
          continue;

        for (int atom_idx: ring_atoms) { xscore_aromatic[atom_idx] = true; }
      }
    }

    return xscore_aromatic;
  }


  // LIGAND SPECIALIZATION

  template<>
  void assign_x_score_types(x_score_static_layer& layer) {
    assign_xtool_types_from_sybyl(layer);

    const auto& mol             = layer.get_base_molecule();
    const int num_atoms = mol.num_atoms();
    const auto sybyl_types = mol.get_atom_type();

    const auto graph = make_graph(mol.get_bonds(), num_atoms);

    const auto xscore_aromatic = detect_xscore_aromaticity(graph, sybyl_types, mol.get_bonds(), num_atoms);

    for (int i = 0; i < num_atoms; ++i) {
      atom_environment env;

      auto [vi, vi_end] = boost::adjacent_vertices(i, graph);
      for (; vi != vi_end; ++vi) {
        const auto neighbor_idx   = *vi;
        const xtool_ff neib_type = sybyl_types[neighbor_idx];

        if (is_hydrogen_type(neib_type)) {
          env.num_h++;
        } else {
          env.num_nonh++;
          if (is_heteroatom_type(neib_type)) {
            env.num_hetero++;
          }
        }
      }

      if (is_hydrogen_type(sybyl_types[i])) {
        auto [vi2, vi2_end] = boost::adjacent_vertices(i, graph);
        for (; vi2 != vi2_end; ++vi2) {
          if (is_oxygen_or_nitrogen_type(sybyl_types[*vi2])) {
            env.is_bonded_to_ON = true;
            break;
          }
        }
      }

      xtool_ff assigned_type     = xtool_ff::Un;
      const xtool_ff center_type = sybyl_types[i];
      const bool is_aromatic = xscore_aromatic[i];

      if (is_hydrogen_type(center_type)) {
        if (env.is_bonded_to_ON)
          assigned_type = xtool_ff::Hhb;
        else
          assigned_type = xtool_ff::H;
      }

      if (center_type == xtool_ff::C3) {
        if (env.num_hetero == 0)
          assigned_type = xtool_ff::C3;
        else if (env.num_hetero > 0)
          assigned_type = xtool_ff::C3x;
        else
          assigned_type = xtool_ff::C3un;
      }

      if (center_type == xtool_ff::C2 && !is_aromatic) {
        if (env.num_hetero == 0)
          assigned_type = xtool_ff::C2;
        else if (env.num_hetero > 0)
          assigned_type = xtool_ff::C2x;
        else
          assigned_type = xtool_ff::C2un;
      }

      if (center_type == xtool_ff::Car || (center_type == xtool_ff::C2 && is_aromatic)) {
        if (env.num_hetero == 0)
          assigned_type = xtool_ff::Car;
        else if (env.num_hetero > 0)
          assigned_type = xtool_ff::Carx;
        else
          assigned_type = xtool_ff::Carun;
      }

      if (center_type == xtool_ff::C1) {
        if (env.num_hetero == 0)
          assigned_type = xtool_ff::C1;
        else if (env.num_hetero > 0)
          assigned_type = xtool_ff::C1x;
        else
          assigned_type = xtool_ff::C1un;
      }

      if (center_type == xtool_ff::Ccat) {
        assigned_type = xtool_ff::Ccat;
      }

      if (center_type == xtool_ff::N4 || center_type == xtool_ff::N3) {
        if (env.num_nonh <= 2)
          assigned_type = xtool_ff::N4;
        else if (env.num_nonh == 3)
          assigned_type = xtool_ff::N3;
        else
          assigned_type = xtool_ff::N3un;
      }

      if ((center_type == xtool_ff::Nam && !is_aromatic) || (center_type == xtool_ff::Npl3 && !is_aromatic)) {
        if (env.num_nonh == 1)
          assigned_type = xtool_ff::Npl3h;
        else if (env.num_nonh == 2)
          assigned_type = xtool_ff::Npl3h;
        else if (env.num_nonh == 3)
          assigned_type = xtool_ff::Npl3;
        else
          assigned_type = xtool_ff::Npl3un;
      }

      if (center_type == xtool_ff::N2 && !is_aromatic) {
        if (env.num_nonh == 1)
          assigned_type = xtool_ff::N2h;
        else if (env.num_nonh == 2)
          assigned_type = xtool_ff::N2;
        else
          assigned_type = xtool_ff::N2un;
      }

      if (center_type == xtool_ff::Nar || (center_type == xtool_ff::N2 && is_aromatic) ||
          (center_type == xtool_ff::Npl3 && is_aromatic) || (center_type == xtool_ff::Nam && is_aromatic)) {
        if (env.num_h == 1)
          assigned_type = xtool_ff::Narh;
        else if (env.num_h == 0)
          assigned_type = xtool_ff::Nar;
        else
          assigned_type = xtool_ff::Narun;
      }

      if (center_type == xtool_ff::N1) {
        if (env.num_nonh == 1)
          assigned_type = xtool_ff::N1;
        else
          assigned_type = xtool_ff::N1un;
      }

      if (center_type == xtool_ff::O3) {
        if (env.num_nonh == 1)
          assigned_type = xtool_ff::O3h;
        else if (env.num_nonh == 2)
          assigned_type = xtool_ff::O3;
        else
          assigned_type = xtool_ff::O3un;
      }

      if (center_type == xtool_ff::O2) {
        assigned_type = xtool_ff::O2;
      }

      if (center_type == xtool_ff::Oco2) {
        assigned_type = xtool_ff::Oco2;
      }

      if (center_type == xtool_ff::S3) {
        if (env.num_nonh == 1)
          assigned_type = xtool_ff::S3h;
        else if (env.num_nonh == 2)
          assigned_type = xtool_ff::S3;
        else
          assigned_type = xtool_ff::S3un;
      }

      if (center_type == xtool_ff::S2) {
        assigned_type = xtool_ff::S2;
      }

      if (center_type == xtool_ff::So) {
        assigned_type = xtool_ff::So;
      }

      if (center_type == xtool_ff::So2) {
        assigned_type = xtool_ff::So;
      }

      if (center_type == xtool_ff::P3) {
        assigned_type = xtool_ff::P3;
      }

      if (center_type == xtool_ff::F)
        assigned_type = xtool_ff::F;

      if (center_type == xtool_ff::Cl)
        assigned_type = xtool_ff::Cl;

      if (center_type == xtool_ff::Br)
        assigned_type = xtool_ff::Br;

      if (center_type == xtool_ff::I)
        assigned_type = xtool_ff::I;

      if (center_type == xtool_ff::Si)
        assigned_type = xtool_ff::Si;

      // DATA SAVING
      layer.x_score_xtool_type(i) = assigned_type;
      layer.vdw_radius(i)         = get_description(assigned_type).vdw_radius;
      layer.hb(i) = get_description(assigned_type).hbond;
      layer.valid(i) =
          (assigned_type == xtool_ff::Un) ? x_score_validity::invalid : x_score_validity::valid;
    }
  }


  static std::string normalize_old_pdb_hydrogen_name(std::string_view name) {
    if (name.size() >= 2 && std::isdigit(static_cast<unsigned char>(name[0])) &&
        !std::isdigit(static_cast<unsigned char>(name[1]))) {
      std::string normalized(name.substr(1));
      normalized += name[0];
      return normalized;
    }
    return {};
  }

  // PROTEIN SPECIALIZATION
  template<>
  void assign_x_score_types(x_score_dynamic_layer& layer) {
    auto& mol = layer.get_base_molecule();

    for (int i = 0; i < mol.num_atoms(); ++i) {
      residue res_type           = mol.residue_types(i);
      std::string_view atom_name = mol.atom_name(i);

      int res_enum_val = static_cast<int>(res_type);

      bool found = false;

      // if your residue is not UNKNOWN (27)
      if (res_enum_val >= 0 && res_enum_val < mudock::num_residues()) {
        const auto& res_desc = get_description(res_type);

        for (const auto& atom_tmpl: res_desc.atoms) {
          if (atom_tmpl.name == atom_name) {
            mol.atom_type(i)            = atom_tmpl.basic_atom_type;
            layer.x_score_xtool_type(i) = atom_tmpl.x_tool_atom_type;
            layer.vdw_radius(i)         = atom_tmpl.vdw_radius;
            layer.hb(i) = atom_tmpl.hbond;
            //logp type can also be assigned here if needed
            found = true;
            break;
          }
        }

        std::string normalized_name;
        if (!found) {
          normalized_name = normalize_old_pdb_hydrogen_name(atom_name);
          if (!normalized_name.empty()) {
            for (const auto& atom_tmpl: res_desc.atoms) {
              if (atom_tmpl.name == normalized_name) {
                mol.atom_type(i)            = atom_tmpl.basic_atom_type;
                layer.x_score_xtool_type(i) = atom_tmpl.x_tool_atom_type;
                layer.vdw_radius(i)         = atom_tmpl.vdw_radius;
                layer.hb(i)                 = atom_tmpl.hbond;
                found                       = true;
                break;
              }
            }
          }
        }
      }

      if (found) {
        layer.valid(i) = x_score_validity::valid;
        continue;
      }

      if (!found) {
        const auto& ter_desc = get_description(residue::TER);
        for (const auto& atom_tmpl: ter_desc.atoms) {
          if (atom_tmpl.name == atom_name) {
            mol.atom_type(i)            = atom_tmpl.basic_atom_type;
            layer.x_score_xtool_type(i) = atom_tmpl.x_tool_atom_type;
            layer.vdw_radius(i)         = atom_tmpl.vdw_radius;
            layer.hb(i) = atom_tmpl.hbond;
            found                       = true;
            break;
          }
        }
        std::string normalized_name;
        if (!found) {
          normalized_name = normalize_old_pdb_hydrogen_name(atom_name);
          if (!normalized_name.empty()) {
            for (const auto& atom_tmpl: ter_desc.atoms) {
              if (atom_tmpl.name == normalized_name) {
                mol.atom_type(i)            = atom_tmpl.basic_atom_type;
                layer.x_score_xtool_type(i) = atom_tmpl.x_tool_atom_type;
                layer.vdw_radius(i)         = atom_tmpl.vdw_radius;
                layer.hb(i)                 = atom_tmpl.hbond;
                found                       = true;
                break;
              }
            }
          }
        }
      }

      if (!found) {
        const auto& het_desc = get_description(residue::HET);
        for (const auto& atom_tmpl: het_desc.atoms) {
          if (atom_tmpl.name == atom_name) {
            mol.atom_type(i)            = atom_tmpl.basic_atom_type;
            layer.x_score_xtool_type(i) = atom_tmpl.x_tool_atom_type;
            layer.vdw_radius(i)         = atom_tmpl.vdw_radius;
            layer.hb(i) = atom_tmpl.hbond;
            found                       = true;
            break;
          }
        }
        std::string normalized_name;
        if (!found) {
          normalized_name = normalize_old_pdb_hydrogen_name(atom_name);
          if (!normalized_name.empty()) {
            for (const auto& atom_tmpl: het_desc.atoms) {
              if (atom_tmpl.name == normalized_name) {
                mol.atom_type(i)            = atom_tmpl.basic_atom_type;
                layer.x_score_xtool_type(i) = atom_tmpl.x_tool_atom_type;
                layer.vdw_radius(i)         = atom_tmpl.vdw_radius;
                layer.hb(i)                 = atom_tmpl.hbond;
                found                       = true;
                break;
              }
            }
          }
        }
      }

      if (!found) {
        mol.atom_type(i)            = xtool_ff::Un;
        layer.x_score_xtool_type(i) = xtool_ff::Un;
        layer.vdw_radius(i)         = get_description(xtool_ff::Un).vdw_radius;
        layer.hb(i)         = get_description(xtool_ff::Un).hbond;
      }

      layer.valid(i) = found ? x_score_validity::valid : x_score_validity::invalid;
    }
  }

} // namespace mudock