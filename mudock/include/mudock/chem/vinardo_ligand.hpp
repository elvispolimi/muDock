#pragma once

#include <mudock/chem/vinardo_layer.hpp>
#include <mudock/chem/vinardo_preprocessing.hpp>
#include <mudock/chem/vinardo_smina_helpers.hpp>
#include <mudock/molecule.hpp>

#include <stdexcept>

namespace mudock {

  struct vinardo_ligand: public vinardo_layer<static_containers> {
    vinardo_ligand(static_molecule& ligand): vinardo_layer<static_containers>(ligand) {
      prepare();
    }

    [[nodiscard]] const auto& get_ligand_ligand_pairs() const {
      return ll_pairs;
    }

    [[nodiscard]] unsigned get_num_tors() const {
      return num_tors;
    }

  private:
    std::vector<vinardo_ligand_ligand_pair> ll_pairs;
    unsigned int num_tors = 0;

    void prepare() {
      auto& ligand = this->get_base_molecule();
      if (!ligand.pdbqt_ligand_data.valid) {
        throw std::runtime_error("Missing PDBQT ligand data");
      }

      ll_pairs = preprocess_ligand_vinardo(*this, ligand.pdbqt_ligand_data.mobility_matrix);
      num_tors = smina_num_tors(ligand, ligand.pdbqt_ligand_data.rotors, this->get_vinardo_type());
    }
  };

} // namespace mudock
