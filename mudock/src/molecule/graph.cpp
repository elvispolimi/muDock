#include <mudock/molecule/graph.hpp>
#include <mudock/type_alias.hpp>
#include <unordered_map>

namespace mudock {
  molecule_graph_type make_graph(const std::span<const bond>& bonds) {
    using vertex_type = typename molecule_graph_type::vertex_descriptor;

    // describe support data structures that we need to compute the fragment mask
    molecule_graph_type g;
    std::unordered_map<int, vertex_type> index2vertex;
    const auto get_vertex = [&index2vertex, &g](const auto index_value) {
      auto vertex_it = index2vertex.find(index_value);
      if (vertex_it == std::end(index2vertex)) { // it's a new vertex
        const auto [new_vertex_it, is_inserted] = index2vertex.emplace(index_value, boost::add_vertex(g));
        assert(is_inserted);
        g[new_vertex_it->second].atom_index = index_value;
        vertex_it                           = new_vertex_it;
      }
      assert(vertex_it != std::end(index2vertex));
      return vertex_it->second;
    };

    // populate the graph with the molecule topology
    for (std::size_t i{0}; i < bonds.size(); ++i) {
      const auto& bond_description = bonds[i];
      const auto [edge, is_inserted] =
          boost::add_edge(get_vertex(bond_description.source), get_vertex(bond_description.dest), g);
      assert(is_inserted);
      g[edge].bond_index = i;
    }
    return g;
  }
} // namespace mudock
