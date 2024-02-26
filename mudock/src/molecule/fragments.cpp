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
#include <utility>

namespace mudock {

  template<>
  fragments<static_containers>::fragments(const std::size_t num_atoms, const std::size_t num_bonds)
      : index(num_atoms, num_bonds) {
    assert(num_atoms < max_static_atoms());
    assert(num_bonds < max_static_bonds());
    storage.fill(coordinate_type{0});
    index = index2D{num_atoms, num_bonds};
  }

  template<>
  fragments<dynamic_containers>::fragments(const std::size_t num_atoms, const std::size_t num_bonds)
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

    std::size_t counter = 0;

  public:
    void discover_vertex([[maybe_unused]] vertex_type u, [[maybe_unused]] const graph_type &g) { ++counter; }
    inline auto get_counter() const { return counter; }
  };

  // find all the rotatable bonds in the bonds
  // NOTE: we define an edge with their two edges, to avoid problems in graph update
  struct edge_description {
    using vertex_type = typename molecule_graph_type::vertex_descriptor;
    vertex_type source;
    vertex_type dest;
  };

  using vertex_type = typename molecule_graph_type::vertex_descriptor;
  template<class container_aliases>
  std::pair<typename container_aliases::template bonds_size<edge_description>, std::size_t>
      get_rotatable_edges(const std::span<const bond> &bonds, const molecule_graph_type &g);

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

  // this is an internal templatized version to fill the fragment mask that is agnostic with respect
  // to the actual container type
  template<class container_aliases>
  fragments<container_aliases> make_fragments_internal(const std::span<const bond> &bonds,
                                                       const std::size_t num_atoms) {
    auto g                                            = make_graph(bonds);
    const auto [rotatable_edges, num_rotatable_edges] = get_rotatable_edges<container_aliases>(bonds, g);
    auto result = fragments<container_aliases>(num_atoms, num_rotatable_edges);
    for (std::size_t i{0}; i < num_rotatable_edges; ++i) {
      const auto edge          = rotatable_edges[i];
      auto mask                = result.get_mask(i);
      const auto source_vertex = edge.source;
      const auto dest_vertex   = edge.dest;

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
  fragments<static_containers> make_fragments<static_containers>(const std::span<const bond> &bonds,
                                                                 const std::size_t num_atoms) {
    return make_fragments_internal<static_containers>(bonds, num_atoms);
  }

  // specialization for dynamic containers
  template<>
  fragments<dynamic_containers> make_fragments<dynamic_containers>(const std::span<const bond> &bonds,
                                                                   const std::size_t num_atoms) {
    return make_fragments_internal<dynamic_containers>(bonds, num_atoms);
  }

} // namespace mudock
