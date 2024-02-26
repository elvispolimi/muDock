#pragma once

#include "mudock/molecule/graph.hpp"
#include "mudock/type_alias.hpp"

#include <mudock/chem.hpp>
#include <mudock/molecule.hpp>
#include <ostream>

namespace mudock {

  class dot {
  public:
    template<class molecule_type>
      requires is_molecule<molecule_type>
    void print(molecule_type&& molecule, std::ostream& out) const {
      out << "graph G {" << std::endl;
      out << "  node [ style=rounded ]" << std::endl;
      for (std::size_t i{0}; i < molecule.num_atoms(); ++i) {
        out << "  " << i << " [ label=\"" << i << "." << get_description(molecule.elements(i)).symbol
            << "\" ]" << std::endl;
      }
      for (std::size_t i{0}; i < molecule.num_bonds(); ++i) {
        const auto& bond = molecule.bonds(i);
        out << "  " << bond.source << " -- " << bond.dest << " [";
        if (bond.can_rotate) {
          out << " style=\"bold\" color=\"red\" ";
        } else {
          out << " style=\"dotted\" ";
        }
        out << "label=\"" << get_description(bond.type).name << "\"]" << std::endl;
      }
      out << "}" << std::endl;
    }

    void print(const molecule_graph_type& g, std::ostream& out) const;
  };

} // namespace mudock
