#include <mudock/format/dot.hpp>

namespace mudock {
  void dot::print(const molecule_graph_type& g, std::ostream& out) const {
    out << "graph G {" << std::endl;
    out << "  node [ style=rounded ]" << std::endl;
    const auto [atom_begin, atom_end] = boost::vertices(g);
    for (auto it = atom_begin; it != atom_end; ++it) {
      out << "  " << g[*it].atom_index << " [ label=\"" << g[*it].atom_index << "\" ]" << std::endl;
    }

    const auto [edge_begin, edge_end] = boost::edges(g);
    for (auto it = edge_begin; it != edge_end; ++it) {
      const auto source = boost::source(*it, g);
      const auto dest   = boost::target(*it, g);
      out << "  " << g[source].atom_index << " -- " << g[dest].atom_index << std::endl;
    }
    out << "}" << std::endl;
  }
} // namespace mudock
