#pragma once

#include <cassert>
#include <filesystem>
#include <fstream>
#include <mudock/chem/autodock_layer.hpp>
#include <mudock/chem/autodock_parameters.hpp>
#include <mudock/chem/autodock_types.hpp>
#include <sstream>
#include <mudock/utils.hpp>

namespace mudock {
  struct pdbqt {
    static constexpr auto PDBQT_ATOM_TOKEN    = "ATOM";
    static constexpr auto PDBQT_HETATOM_TOKEN = "HETATOM";
    static constexpr auto PDBQT_TORSDOF_TOKEN = "TORSDOF";
    static constexpr auto PDBQT_START_TOKEN   = "REMARK  Name =";

    std::string_view::size_type next_molecule_start_index(std::string_view text) const;
  };
  // The following function relies on the same order of atom loading and file
  // FIXME use the same approach as for the check rotor
  template<class molecule_type>
    requires std::derived_from<molecule_type, autodock_static_layer> ||
             std::derived_from<molecule_type, autodock_dynamic_layer>
  void apply_autodock_forcefield_pdbqt(molecule_type& molecule, const std::filesystem::path input_path) {
    const auto desc = read_from_stream(std::ifstream(input_path));
    std::stringstream desc_s{desc};

    [[maybe_unused]] const std::size_t num_atoms = molecule.num_atoms();

    int index = 0;
    std::string line;
    while (std::getline(desc_s, line)) {
      if (line.find(pdbqt::PDBQT_ATOM_TOKEN) != std::string::npos ||
          line.find(pdbqt::PDBQT_HETATOM_TOKEN) != std::string::npos) {
        assert(index < static_cast<int>(num_atoms));
        // FIXMED Really bad, at the moment we rely on OpenBabel PDBQT structure
        // What if PDBQT is standardized..
        if (line.size() < 79)
          line += " ";
        std::string adt_value = line.substr(77, 2);
        assert(adt_value.size() == 2);
        if (adt_value[1] == ' ')
          adt_value.pop_back();
        // FIX ME add check that the order of atoms is the same

        const auto adt                  = parse_autodock_type(adt_value);
        molecule().autodock_type(index) = adt;
        const auto& ff_entry            = get_description(adt);
        // molecule.autodock_type(index)   = ff_entry.value;
        molecule.Rii(index)         = ff_entry.Rii;
        molecule.epsii(index)       = ff_entry.epsii * autodock_parameters::coeff_vdW;
        molecule.vol(index)         = ff_entry.vol;
        molecule.solpar(index)      = ff_entry.solpar;
        molecule.Rij_hb(index)      = ff_entry.Rij_hb;
        molecule.epsij_hb(index)    = ff_entry.epsij_hb * autodock_parameters::coeff_hbond;
        molecule().num_hbond(index) = ff_entry.hbond;
        ++index;
      }
    }
  }
} // namespace mudock
