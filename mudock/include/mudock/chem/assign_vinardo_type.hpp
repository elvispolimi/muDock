#pragma once

#include <cstdint>
#include <mudock/chem/autodock_types.hpp>
#include <mudock/chem/vinardo_type.hpp>
#include <mudock/molecule.hpp>
#include <mudock/molecule/graph.hpp>
#include <mudock/type_alias.hpp>
#include <stdexcept>

namespace mudock {
   // Convert muDock AutoDock atom types into Vinardo atom types in three phases:
   // 1. normalize muDock autodock_ff values to the Vinardo/Smina lookup domain;
   // 2. assign an initial Vinardo type from the lookup table;
   // 3. resolve ambiguous cases using neighboring temporary Vinardo types.

  template<class container_aliases, typename T>
  using atoms_array_type = typename container_aliases::template atoms_size<T>;

   [[nodiscard]] inline ad_lookup_table_input_type normalize_autodock_to_vinardo_lookup(const autodock_ff type) {
      switch (type) {
         case autodock_ff::H:  return ad_lookup_table_input_type::H;
         case autodock_ff::HD: return ad_lookup_table_input_type::HD;
         case autodock_ff::HS: return ad_lookup_table_input_type::HD;
         case autodock_ff::C:  return ad_lookup_table_input_type::C;
         case autodock_ff::A:  return ad_lookup_table_input_type::A;
         case autodock_ff::N:  return ad_lookup_table_input_type::N;
         case autodock_ff::NS: return ad_lookup_table_input_type::NA;
         case autodock_ff::NA: return ad_lookup_table_input_type::NA;
         case autodock_ff::OA: return ad_lookup_table_input_type::OA;
         case autodock_ff::OS: return ad_lookup_table_input_type::OA;
         case autodock_ff::S:  return ad_lookup_table_input_type::S;
         case autodock_ff::SA: return ad_lookup_table_input_type::SA;
         case autodock_ff::P:  return ad_lookup_table_input_type::P;
         case autodock_ff::F:  return ad_lookup_table_input_type::F;
         case autodock_ff::Cl: return ad_lookup_table_input_type::Cl;
         case autodock_ff::Br: return ad_lookup_table_input_type::Br;
         case autodock_ff::I:  return ad_lookup_table_input_type::I;
         case autodock_ff::Mg: return ad_lookup_table_input_type::Mg;
         case autodock_ff::Mn: return ad_lookup_table_input_type::Mn;
         case autodock_ff::Zn: return ad_lookup_table_input_type::Zn;
         case autodock_ff::Ca: return ad_lookup_table_input_type::Ca;
         case autodock_ff::Fe: return ad_lookup_table_input_type::Fe;

         case autodock_ff::Se:
            return ad_lookup_table_input_type::S;

         // fallback to GenericMetal for metal-like types
         case autodock_ff::Cu:
         case autodock_ff::Na:
         case autodock_ff::K:
         case autodock_ff::Hg:
         case autodock_ff::Co:
         case autodock_ff::U:
         case autodock_ff::Cd:
         case autodock_ff::Ni:
            return ad_lookup_table_input_type::GenericMetal;

         // Unsupported Types
         // case autodock_ff::He:
         // case autodock_ff::Ne:
         // case autodock_ff::Al:
         // case autodock_ff::Si:
         // case autodock_ff::Z:
         // case autodock_ff::G:
         // case autodock_ff::GA:
         // case autodock_ff::J:
         // case autodock_ff::Q:
         //   return ad_lookup_table_input_type::Unsupported;

         default:
            return ad_lookup_table_input_type::Unsupported;
      }
   }

   [[nodiscard]] inline bool is_polar_hydrogen(const vinardo_atom_type type) {
      return type == vinardo_atom_type::PolarHydrogen;
   }

   [[nodiscard]] inline bool is_heteroatom(const vinardo_atom_type type) {
      switch (type) {
         case vinardo_atom_type::Nitrogen:
         case vinardo_atom_type::NitrogenXSDonor:
         case vinardo_atom_type::NitrogenXSDonorAcceptor:
         case vinardo_atom_type::NitrogenXSAcceptor:
         case vinardo_atom_type::Oxygen:
         case vinardo_atom_type::OxygenXSDonor:
         case vinardo_atom_type::OxygenXSDonorAcceptor:
         case vinardo_atom_type::OxygenXSAcceptor:
         case vinardo_atom_type::Sulfur:
         case vinardo_atom_type::SulfurAcceptor:
         case vinardo_atom_type::Phosphorus:
         case vinardo_atom_type::Fluorine:
         case vinardo_atom_type::Chlorine:
         case vinardo_atom_type::Bromine:
         case vinardo_atom_type::Iodine:
         case vinardo_atom_type::Magnesium:
         case vinardo_atom_type::Manganese:
         case vinardo_atom_type::Zinc:
         case vinardo_atom_type::Calcium:
         case vinardo_atom_type::Iron:
         case vinardo_atom_type::GenericMetal:
            return true;
         default:
            return false;
      }
   }

   // Resolve ambiguous Vinardo atom types using local topology and the
   // temporary Vinardo types already assigned to neighboring atoms.
   [[nodiscard]] inline vinardo_atom_type adjust_vinardo_type(const auto& atom_vinardo_type,
                                                              const auto& graph,
                                                              std::size_t atom_idx) {
      bool hbonded = false;
      bool hetero_bonded = false;
      const auto type = atom_vinardo_type[atom_idx];
      const auto [begin, end] = boost::adjacent_vertices(atom_idx, graph);
      for (auto vertex_iterator = begin; vertex_iterator != end; ++vertex_iterator) {
         const auto neighb_vertex = *vertex_iterator;
         const auto neighb_idx = graph[neighb_vertex].atom_index;
         const auto neighb_type = atom_vinardo_type[neighb_idx];

         if (is_polar_hydrogen(neighb_type)) {
            hbonded = true;
         }

         if (is_heteroatom(neighb_type)) {
            hetero_bonded = true;
         }
         //Early Stop if possible
         if (hbonded && hetero_bonded) {
            break;
         }
      }

      //Now I need to resolve the ambiguity using the graph

      switch (type) {
         case vinardo_atom_type::AliphaticCarbonXSHydrophobe:
         case vinardo_atom_type::AliphaticCarbonXSNonHydrophobe:
            return hetero_bonded
               ? vinardo_atom_type::AliphaticCarbonXSNonHydrophobe
               : vinardo_atom_type::AliphaticCarbonXSHydrophobe;

         case vinardo_atom_type::AromaticCarbonXSHydrophobe:
         case vinardo_atom_type::AromaticCarbonXSNonHydrophobe:
            return hetero_bonded
               ? vinardo_atom_type::AromaticCarbonXSNonHydrophobe
               : vinardo_atom_type::AromaticCarbonXSHydrophobe;

         case vinardo_atom_type::NitrogenXSDonor:
         case vinardo_atom_type::Nitrogen:
            return hbonded
               ? vinardo_atom_type::NitrogenXSDonor
               : vinardo_atom_type::Nitrogen;

         case vinardo_atom_type::NitrogenXSDonorAcceptor:
         case vinardo_atom_type::NitrogenXSAcceptor:
            return hbonded
               ? vinardo_atom_type::NitrogenXSDonorAcceptor
               : vinardo_atom_type::NitrogenXSAcceptor;

         case vinardo_atom_type::OxygenXSDonor:
         case vinardo_atom_type::Oxygen:
            return hbonded
               ? vinardo_atom_type::OxygenXSDonor
               : vinardo_atom_type::Oxygen;

         case vinardo_atom_type::OxygenXSDonorAcceptor:
         case vinardo_atom_type::OxygenXSAcceptor:
            return hbonded
               ? vinardo_atom_type::OxygenXSDonorAcceptor
               : vinardo_atom_type::OxygenXSAcceptor;

         default:
            return type;
      }
   }
   [[nodiscard]] inline const vinardo_type_info& get_vinardo_info(const vinardo_atom_type type) {
      for (std::size_t i = 0; i < static_cast<std::size_t>(vinardo_atom_type::NumTypes); ++i) {
         if (vinardo_data[i].vd_type == type) {
            return vinardo_data[i];
         }
      }
      throw std::runtime_error("get_vinardo_info: no entry found for vinardo_atom_type " + std::to_string(static_cast<int>(type)));
   }

   template<class container_aliases>
    requires is_container_specification<container_aliases>
   void assign_vinardo_type(
      molecule<container_aliases>& molecule,
      atoms_array_type<container_aliases, vinardo_atom_type>& atom_vinardo_type,
      atoms_array_type<container_aliases, fp_type>& atom_radius,
      atoms_array_type<container_aliases, std::uint8_t>& atom_is_hydrophobic,
      atoms_array_type<container_aliases, std::uint8_t>& atom_is_hbond_donor,
      atoms_array_type<container_aliases, std::uint8_t>& atom_is_hbond_acceptor) {

      //Compute the graph for resolving ambiguity
      auto graph = mudock::make_graph(molecule.get_bonds(), molecule.num_atoms());

      //Here i need the code logic to convert an autodock atom type to a vinardo atom type

      const auto types = molecule.get_autodock_type();
      for (std::size_t i = 0; i < static_cast<std::size_t>(molecule.num_atoms()); ++i) {
         const autodock_ff ad_type = types[i];
         const ad_lookup_table_input_type lookup_type = normalize_autodock_to_vinardo_lookup(ad_type);
         if (lookup_type == ad_lookup_table_input_type::Unsupported) {
            throw std::runtime_error("Unsupported autodock_ff for Vinardo conversion: " + std::to_string(static_cast<int>(ad_type)));
         }
         for (std::size_t j = 0; j < static_cast<std::size_t>(vinardo_atom_type::NumTypes); j++) {
            if (vinardo_data[j].ad_type == lookup_type) {
               atom_vinardo_type[i] = vinardo_data[j].vd_type;
               break; //I take only the first match
            }
         }
      }
      //Now i need to resolve the ambiguity using the graph

      for (std::size_t i = 0; i < static_cast<std::size_t>(molecule.num_atoms()); ++i) {
         atom_vinardo_type[i] = adjust_vinardo_type(atom_vinardo_type, graph, i);
         //Now I can compile the field in the layer
         const auto& info = get_vinardo_info(atom_vinardo_type[i]);
         atom_radius[i] = info.xs_radius;
         atom_is_hydrophobic[i] = info.xs_hydrophobe;
         atom_is_hbond_donor[i] = info.xs_donor;
         atom_is_hbond_acceptor[i] = info.xs_acceptor;
      }

   }

} // namespace mudock
