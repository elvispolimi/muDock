#include <algorithm>
#include <boost/graph/breadth_first_search.hpp>
#include <cassert>
#include <cstdint>
#include <iterator>
#include <mudock/grid.hpp>
#include <mudock/molecule/constraints.hpp>
#include <mudock/molecule/containers.hpp>
#include <mudock/molecule/fragments.hpp>
#include <mudock/molecule/graph.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  template<>
  fragments<static_container_type>::fragments(const std::size_t num_atoms, const std::size_t num_bonds)
      : index(num_atoms, num_bonds) {
    assert(num_atoms < max_static_atoms());
    assert(num_bonds < max_static_bonds());
    storage.fill(coordinate_type{0});
    index = index2D{num_atoms, num_bonds};
  }

  template<>
  fragments<dynamic_container_type>::fragments(const std::size_t num_atoms, const std::size_t num_bonds)
      : index(num_atoms, num_bonds) {
    storage.clear();
    storage.resize(num_atoms * num_bonds, coordinate_type{0});
    index = index2D{num_atoms, num_bonds};
  }

  //===------------------------------------------------------------------------------------------------------
  // Implementation fo the methods that build the fragment mask
  //===------------------------------------------------------------------------------------------------------

  // this is a simple functor that set to one all the atoms that it visits
  class bitmask_setter: public boost::default_bfs_visitor {
    using graph_type  = molecule_graph_type;
    using vertex_type = typename boost::graph_traits<graph_type>::vertex_descriptor;

    std::span<coordinate_type> bitmask;

  public:
    bitmask_setter(std::span<coordinate_type> mask): bitmask(mask) {}
    void discover_vertex(vertex_type u, const graph_type &g) {
      bitmask[g[u].atom_index] = static_cast<coordinate_type>(1);
    }
  };

  // this is a simple functor that count how many atoms it visits
  class atom_counter: public boost::default_bfs_visitor {
    using graph_type  = molecule_graph_type;
    using vertex_type = typename boost::graph_traits<graph_type>::vertex_descriptor;

    std::size_t counter = 0;

  public:
    void discover_vertex([[maybe_unused]] vertex_type u, [[maybe_unused]] const graph_type &g) { ++counter; }
    inline auto get_counter() const { return counter; }
  };

  // find all the rotatable bonds in the bonds
  using edge_type = typename molecule_graph_type::edge_descriptor;
  template<template<typename> class container_type>
  std::pair<container_type<edge_type>, index_type> get_rotatable_edges(const container_type<bond> &bonds,
                                                                       const molecule_graph_type &g);

  template<>
  std::pair<static_container_type<edge_type>, index_type>
      get_rotatable_edges<static_container_type>(const static_container_type<bond> &bonds,
                                                 const molecule_graph_type &g) {
    assert(boost::num_edges(g) < max_static_bonds());
    static_container_type<edge_type> result;
    auto rotatable_bond_counter   = index_type{0};
    const auto [begin_it, end_it] = boost::edges(g);
    for (auto it = begin_it; it != end_it; ++it) {
      if (bonds[g[*it].bond_index].can_rotate) {
        result[rotatable_bond_counter] = *it;
        ++rotatable_bond_counter;
      }
    }
    return std::make_pair(result, rotatable_bond_counter);
  }

  template<>
  std::pair<dynamic_container_type<edge_type>, index_type>
      get_rotatable_edges<dynamic_container_type>(const dynamic_container_type<bond> &bonds,
                                                  const molecule_graph_type &g) {
    dynamic_container_type<edge_type> result;
    result.reserve(static_cast<std::size_t>(boost::num_edges(g)));
    const auto [begin_it, end_it] = boost::edges(g);
    for (auto it = begin_it; it != end_it; ++it) {
      if (bonds[g[*it].bond_index].can_rotate) {
        result.emplace_back(*it);
      }
    }
    return std::make_pair(result, static_cast<index_type>(result.size()));
  }

  // this is an internal templatized version to fill the fragment mask that is agnostic with respect
  // to the actual container type
  template<template<typename> class container_type>
  fragments<container_type> make_fragments_internal(const container_type<bond> &bonds,
                                                    const index_type num_atoms,
                                                    const index_type num_bonds) {
    auto g                                            = make_graph<container_type>(bonds, num_bonds);
    const auto [rotatable_edges, num_rotatable_edges] = get_rotatable_edges<container_type>(bonds, g);
    auto result = fragments<container_type>(num_atoms, num_rotatable_edges);
    for (index_type i{0}; i < num_rotatable_edges; ++i) {
      const auto edge          = rotatable_edges[i];
      auto mask                = result.get_mask(i);
      const auto source_vertex = boost::source(edge, g);
      const auto dest_vertex   = boost::target(edge, g);

      boost::remove_edge(source_vertex, dest_vertex, g);
      atom_counter counter_source, counter_dest;
      boost::breadth_first_search(g, source_vertex, boost::visitor(counter_source));
      boost::breadth_first_search(g, dest_vertex, boost::visitor(counter_dest));
      if (counter_source.get_counter() > counter_dest.get_counter()) {
        boost::breadth_first_search(g, dest_vertex, boost::visitor(bitmask_setter{mask}));
      } else {
        boost::breadth_first_search(g, source_vertex, boost::visitor(bitmask_setter{mask}));
      }
      boost::add_edge(source_vertex, dest_vertex, g);
    }
    return result;
  }

  // specialization for static containers
  template<>
  fragments<static_container_type>
      make_fragments<static_container_type>(const static_container_type<bond> &bonds,
                                            const index_type num_atoms,
                                            const index_type num_bonds) {
    return make_fragments_internal<static_container_type>(bonds, num_atoms, num_bonds);
  }

  // specialization for dynamic containers
  template<>
  fragments<dynamic_container_type>
      make_fragments<dynamic_container_type>(const dynamic_container_type<bond> &bonds,
                                             const index_type num_atoms,
                                             const index_type num_bonds) {
    return make_fragments_internal<dynamic_container_type>(bonds, num_atoms, num_bonds);
  }

} // namespace mudock
