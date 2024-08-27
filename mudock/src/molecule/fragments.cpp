#include <algorithm>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/connected_components.hpp>
#include <cassert>
#include <cstdint>
#include <gsl/pointers>
#include <iterator>
#include <mudock/grid.hpp>
#include <mudock/molecule/constraints.hpp>
#include <mudock/molecule/containers.hpp>
#include <mudock/molecule/fragments.hpp>
#include <mudock/type_alias.hpp>
#include <utility>

namespace mudock {

  //===------------------------------------------------------------------------------------------------------
  // Utility function to find the rotatable edge on a graph (using a consistent container type)
  //===------------------------------------------------------------------------------------------------------
  // NOTE: we define an edge with their two edges, to avoid problems in graph update
  template<>
  std::pair<typename static_containers::template bonds_size<edge_description>, std::size_t>
      get_rotatable_edges<static_containers>(const std::span<const bond> &bonds,
                                             const molecule_graph_type &g) {
    assert(boost::num_edges(g) < max_static_bonds());
    static_containers::template bonds_size<edge_description> result;
    auto rotatable_bond_counter   = std::size_t{0};
    const auto [begin_it, end_it] = boost::edges(g);
    for (auto it = begin_it; it != end_it; ++it) {
      if (bonds[g[*it].bond_index].can_rotate) {
        result[rotatable_bond_counter] = edge_description{boost::source(*it, g), boost::target(*it, g)};
        ++rotatable_bond_counter;
      }
    }
    return std::make_pair(result, rotatable_bond_counter);
  }

  template<>
  std::pair<typename dynamic_containers::template bonds_size<edge_description>, std::size_t>
      get_rotatable_edges<dynamic_containers>(const std::span<const bond> &bonds,
                                              const molecule_graph_type &g) {
    dynamic_containers::template bonds_size<edge_description> result;
    result.reserve(static_cast<std::size_t>(boost::num_edges(g)));
    const auto [begin_it, end_it] = boost::edges(g);
    for (auto it = begin_it; it != end_it; ++it) {
      if (bonds[g[*it].bond_index].can_rotate) {
        result.emplace_back(boost::source(*it, g), boost::target(*it, g));
      }
    }
    return std::make_pair(result, result.size());
  }

  //===------------------------------------------------------------------------------------------------------
  // Functor that work with graph (to count and find the bitmask)
  //===------------------------------------------------------------------------------------------------------

  // this is a simple functor that set to one all the atoms that it visits
  class bitmask_setter: public boost::default_bfs_visitor {
    using graph_type      = molecule_graph_type;
    using vertex_type     = typename boost::graph_traits<graph_type>::vertex_descriptor;
    using mask_value_type = typename fragments<static_containers>::value_type;

    std::span<mask_value_type> bitmask;

  public:
    bitmask_setter(std::span<mask_value_type> mask): bitmask(mask) {}
    void discover_vertex(vertex_type u, const graph_type &g) {
      bitmask[g[u].atom_index] = static_cast<mask_value_type>(1);
    }
  };

  // this is a simple functor that count how many atoms it visits
  class atom_counter: public boost::default_bfs_visitor {
    using graph_type  = molecule_graph_type;
    using vertex_type = typename boost::graph_traits<graph_type>::vertex_descriptor;

    std::size_t &counter;

  public:
    atom_counter(std::size_t &c): counter(c) {}
    void discover_vertex([[maybe_unused]] vertex_type u, [[maybe_unused]] const graph_type &g) { ++counter; }
  };

  //===------------------------------------------------------------------------------------------------------
  // Utility function that fill the information of the fragments
  //===------------------------------------------------------------------------------------------------------
  template<class container_aliases>
    requires is_container_specification<container_aliases>
  void fragments<container_aliases>::fill_fragment_mask(const std::size_t index_mask,
                                                        gsl::not_null<std::size_t *> start_index,
                                                        gsl::not_null<std::size_t *> stop_index,
                                                        const edge_description &rotatable_bond,
                                                        molecule_graph_type &g) {
    std::span<fragments<static_containers>::value_type> tmp_mask = get_mask(index_mask);
    const auto source_vertex                                     = rotatable_bond.source;
    const auto dest_vertex                                       = rotatable_bond.dest;
    boost::remove_edge(source_vertex, dest_vertex, g);
    std::size_t counter_source{0}, counter_dest{0};
    boost::breadth_first_search(g, source_vertex, boost::visitor(atom_counter{counter_source}));
    boost::breadth_first_search(g, dest_vertex, boost::visitor(atom_counter{counter_dest}));
    if (counter_source > counter_dest) {
      boost::breadth_first_search(g, dest_vertex, boost::visitor(bitmask_setter{tmp_mask}));
      *start_index.get() = g[source_vertex].atom_index;
      *stop_index.get()  = g[dest_vertex].atom_index;
    } else {
      boost::breadth_first_search(g, source_vertex, boost::visitor(bitmask_setter{tmp_mask}));
      *start_index.get() = g[dest_vertex].atom_index;
      *stop_index.get()  = g[source_vertex].atom_index;
    }
    tmp_mask[g[source_vertex].atom_index] = fragments<static_containers>::value_type{2};
    tmp_mask[g[dest_vertex].atom_index]   = fragments<static_containers>::value_type{3};
    boost::add_edge(source_vertex, dest_vertex, g);
  }

  template<class container_aliases>
    requires is_container_specification<container_aliases>
  void fragments<container_aliases>::fill_rigid_pieces(molecule_graph_type &g) {
    boost::connected_components(g, rigid_pieces.data());
  }
} // namespace mudock
