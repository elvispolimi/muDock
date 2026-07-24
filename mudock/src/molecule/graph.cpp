#include <cstddef>
#include <iterator>
#include <mudock/molecule/graph.hpp>
#include <mudock/type_alias.hpp>
#include <unordered_map>

namespace mudock {
  molecule_graph_type make_graph(const std::span<const bond>& bonds, const std::size_t num_atoms) {
    using vertex_type = typename molecule_graph_type::vertex_descriptor;

    // describe support data structures that we need to compute the fragment mask
    molecule_graph_type g;
    std::unordered_map<int, vertex_type> index2vertex;
    for (std::size_t i = 0; i < num_atoms; ++i) {
      const auto [new_vertex_it, is_inserted] = index2vertex.emplace(i, boost::add_vertex(g));
      assert(is_inserted);
      g[new_vertex_it->second].atom_index = static_cast<int>(i);
    }

    // populate the graph with the molecule topology
    for (std::size_t i{0}; i < bonds.size(); ++i) {
      const auto& bond_description = bonds[i];
      const auto source            = index2vertex.find(bond_description.source);
      const auto dest              = index2vertex.find(bond_description.dest);
      assert(source != std::end(index2vertex) && dest != std::end(index2vertex));
      const auto [edge, is_inserted] = boost::add_edge(source->second, dest->second, g);
      assert(is_inserted);
      g[edge].bond_index = static_cast<int>(i);
    }

    return g;
  }
} // namespace mudock
