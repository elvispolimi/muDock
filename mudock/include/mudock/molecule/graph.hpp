#pragma once

#include <boost/graph/adjacency_list.hpp>
#include <mudock/molecule/bond.hpp>
#include <mudock/type_alias.hpp>
#include <span>

namespace mudock {

  // graph description
  struct molecule_edge_type {
    index_type bond_index;
  };
  struct molecule_vertex_type {
    index_type atom_index;
  };
  using molecule_graph_type = boost::
      adjacency_list<boost::setS, boost::vecS, boost::undirectedS, molecule_vertex_type, molecule_edge_type>;

  // utility function that build a graph from the molecules bond
  [[nodiscard]] molecule_graph_type make_graph(const std::span<const bond>& bonds,
                                               const index_type num_bonds);

} // namespace mudock
