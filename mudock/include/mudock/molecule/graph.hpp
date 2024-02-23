#pragma once

#include <boost/graph/adjacency_list.hpp>
#include <mudock/molecule/bond.hpp>
#include <mudock/molecule/fragments.hpp>
#include <mudock/type_alias.hpp>

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
  template<template<typename> class container_type>
  molecule_graph_type make_graph(const container_type<bond>& bonds, const index_type num_bonds) {
    using vertex_type = typename molecule_graph_type::vertex_descriptor;

    // describe support data structures that we need to compute the fragment mask
    molecule_graph_type g;
    std::unordered_map<index_type, vertex_type> index2vertex;
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
    for (index_type i{0}; i < num_bonds; ++i) {
      const auto& bond_description = bonds[i];
      const auto [edge, is_inserted] =
          boost::add_edge(get_vertex(bond_description.source), get_vertex(bond_description.dest), g);
      assert(is_inserted);
      g[edge].bond_index = static_cast<index_type>(i);
    }
    return g;
  }

} // namespace mudock
