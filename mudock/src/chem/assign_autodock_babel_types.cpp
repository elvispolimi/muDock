#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <mudock/chem/assign_autodock_babel_types.hpp>
#include <mudock/grid.hpp>
#include <mudock/type_alias.hpp>
#include <stdexcept>
#include <utility>

namespace mudock {

  //===------------------------------------------------------------------------------------------------------
  // Utility functions to compute properties useful for determine the atom type
  //===------------------------------------------------------------------------------------------------------

  template<class edge_operator>
  static auto count_neighbors_if(const molecule_graph_type::vertex_descriptor v,
                                 const molecule_graph_type& graph,
                                 edge_operator&& op) {
    const auto [begin, end] = boost::out_edges(v, graph);
    return std::count_if(begin, end, op);
  }

  class is_heavy_edge {
    const std::span<element>& elements;
    const molecule_graph_type& graph;

  public:
    inline is_heavy_edge(const std::span<element>& e, const molecule_graph_type& g): elements(e), graph(g) {}

    inline auto operator()(const molecule_graph_type::edge_descriptor edge) {
      return elements[graph[boost::target(edge, graph)].atom_index] == element::H;
    }
  };

  class is_free_ox {
    const std::span<element>& elements;
    const molecule_graph_type& graph;

  public:
    inline is_free_ox(const std::span<element>& e, const molecule_graph_type& g): elements(e), graph(g) {}

    inline auto operator()(const molecule_graph_type::edge_descriptor edge) {
      return count_neighbors_if(boost::target(edge, graph), graph, is_heavy_edge{elements, graph}) == 1;
    }
  };

  static auto find_closest_neighbor(const molecule_graph_type::vertex_descriptor v,
                                    const molecule_graph_type& graph,
                                    const std::span<coordinate_type>& x,
                                    const std::span<coordinate_type>& y,
                                    const std::span<coordinate_type>& z) {
    const auto index_source = graph[v].atom_index;
    const auto source_point = point3D{x[index_source], y[index_source], z[index_source]};
    const auto [begin, end] = boost::out_edges(v, graph);
    auto min_distance       = std::numeric_limits<coordinate_type>::max();
    auto min_index          = std::numeric_limits<std::size_t>::max();
    for (auto it = begin; it != end; ++it) {
      const auto index_target = graph[boost::target(*it, graph)].atom_index;
      const auto target_point = point3D{x[index_target], y[index_target], z[index_target]};
      const auto distance     = std::sqrt(distance2(source_point, target_point));
      if (distance < min_distance) {
        min_distance = distance;
        min_index    = index_target;
      }
    }
    return std::make_pair(min_distance, min_index);
  }

  static auto find_mean_angle_3_neighbors(const molecule_graph_type::vertex_descriptor v,
                                          const molecule_graph_type& graph,
                                          const std::span<coordinate_type>& x,
                                          const std::span<coordinate_type>& y,
                                          const std::span<coordinate_type>& z) {
    const auto origin_index = graph[v].atom_index;
    std::array<std::size_t, 3> neighbors_index;
    const auto [begin, end] = boost::out_edges(v, graph);
    std::transform(begin, end, std::begin(neighbors_index), [&](const auto bond) {
      return graph[boost::target(bond, graph)].atom_index;
    });
    const auto origin    = point3D{x[origin_index], y[origin_index], z[origin_index]};
    const auto k         = point3D{x[neighbors_index[0]], y[neighbors_index[0]], z[neighbors_index[0]]};
    const auto l         = point3D{x[neighbors_index[1]], y[neighbors_index[1]], z[neighbors_index[1]]};
    const auto m         = point3D{x[neighbors_index[2]], y[neighbors_index[2]], z[neighbors_index[2]]};
    const auto angle_kal = angle(origin, k, l);
    const auto angle_kam = angle(origin, k, m);
    const auto angle_mal = angle(origin, l, m);
    return rad_to_deg((angle_kal + angle_kam + angle_mal) / coordinate_type{3});
  }

  static auto find_angle_2_neighbors(const molecule_graph_type::vertex_descriptor v,
                                     const molecule_graph_type& graph,
                                     const std::span<coordinate_type>& x,
                                     const std::span<coordinate_type>& y,
                                     const std::span<coordinate_type>& z) {
    const auto origin_index = graph[v].atom_index;
    std::array<std::size_t, 2> neighbors_index;
    const auto [begin, end] = boost::out_edges(v, graph);
    std::transform(begin, end, std::begin(neighbors_index), [&](const auto bond) {
      return graph[boost::target(bond, graph)].atom_index;
    });
    const auto origin    = point3D{x[origin_index], y[origin_index], z[origin_index]};
    const auto k         = point3D{x[neighbors_index[0]], y[neighbors_index[0]], z[neighbors_index[0]]};
    const auto m         = point3D{x[neighbors_index[1]], y[neighbors_index[1]], z[neighbors_index[1]]};
    const auto angle_kam = angle(origin, k, m);
    return angle_kam < coordinate_type{180} ? angle_kam : coordinate_type{360} - angle_kam;
  }

  static constexpr auto square(const coordinate_type x) { return x * x; }

  //===------------------------------------------------------------------------------------------------------
  // Utility function to set the initial atom type
  //===------------------------------------------------------------------------------------------------------

  static auto handleH(const std::span<element>& elements,
                      const molecule_graph_type& graph,
                      const molecule_graph_type::vertex_descriptor v) -> autodock_babel_ff {
    assert(elements[graph[v].atom_index] == element::H);
    const auto is_not_carbon = [&](const molecule_graph_type::edge_descriptor edge) {
      return elements[graph[boost::target(edge, graph)].atom_index] != element::C;
    };
    const auto [begin, end] = boost::out_edges(v, graph);
    return std::all_of(begin, end, is_not_carbon) ? autodock_babel_ff::H : autodock_babel_ff::HC;
  }

  static auto handleC(const std::span<element>& elements,
                      const std::span<coordinate_type>& x,
                      const std::span<coordinate_type>& y,
                      const std::span<coordinate_type>& z,
                      const molecule_graph_type& graph,
                      const molecule_graph_type::vertex_descriptor v) -> autodock_babel_ff {
    const auto num_neighbors = boost::out_degree(v, graph);
    if (num_neighbors == 4) {
      return autodock_babel_ff::C3;
    } else if (num_neighbors == 3) {
      // compute the average angles with the neighbors k,m,l
      const auto mean_angle = find_mean_angle_3_neighbors(v, graph, x, y, z);

      // assign the correct type based on the average angle and the neighbors
      if (mean_angle < coordinate_type{114.8}) {
        return autodock_babel_ff::C3;
      } else {
        const auto num_free_ox = count_neighbors_if(v, graph, is_free_ox{elements, graph});
        if (num_free_ox >= 2) {
          return autodock_babel_ff::Cac;
        } else {
          return autodock_babel_ff::C2;
        }
      }
    } else if (num_neighbors == 2) {
      // compute the angle with the neighbors k,m
      const auto angle = find_angle_2_neighbors(v, graph, x, y, z);

      // find out the closest atom
      const auto [min_distance, min_index] = find_closest_neighbor(v, graph, x, y, z);

      // assign the correct type based on angle and distance
      if (angle < coordinate_type{114.8}) {
        if ((min_distance < coordinate_type{1.42} && elements[min_index] == element::C) ||
            (min_distance < coordinate_type{1.41} && elements[min_index] == element::N)) {
          return autodock_babel_ff::C2;
        } else {
          return autodock_babel_ff::C3;
        }
      } else if (angle < coordinate_type{122}) {
        if ((min_distance > coordinate_type{1.41} && elements[min_index] == element::C) ||
            (min_distance > coordinate_type{1.46} && elements[min_index] == element::N) ||
            (min_distance > coordinate_type{1.44} && elements[min_index] == element::O)) {
          return autodock_babel_ff::C3;
        } else {
          return autodock_babel_ff::C2;
        }
      } else if (angle < coordinate_type{160}) {
        return autodock_babel_ff::C2;
      } else {
        return autodock_babel_ff::C1;
      }
    }
    return autodock_babel_ff::C;
  }

  static auto handleN(const std::span<element>& elements,
                      const std::span<coordinate_type>& x,
                      const std::span<coordinate_type>& y,
                      const std::span<coordinate_type>& z,
                      const molecule_graph_type& graph,
                      const molecule_graph_type::vertex_descriptor v) -> autodock_babel_ff {
    const auto num_neighbors = boost::out_degree(v, graph);
    if (num_neighbors == 4) {
      const auto num_free_ox = count_neighbors_if(v, graph, is_free_ox{elements, graph});
      if (num_free_ox >= 1) {
        return autodock_babel_ff::Nox;
      } else {
        return autodock_babel_ff::N3_plus;
      }
    } else if (num_neighbors == 3) {
      // compute the average angles with the neighbors k,m,l
      const auto mean_angle = find_mean_angle_3_neighbors(v, graph, x, y, z);

      // assign the correct type based on the average angle and the neighbors
      if (mean_angle < coordinate_type{114.8}) {
        return autodock_babel_ff::N3;
      } else {
        const auto num_free_ox = count_neighbors_if(v, graph, is_free_ox{elements, graph});
        if (num_free_ox >= 2) {
          return autodock_babel_ff::Ntr;
        } else {
          return autodock_babel_ff::Npl;
        }
      }
    } else if (num_neighbors == 2) {
      // compute the angle with the neighbors k,m
      const auto angle = find_angle_2_neighbors(v, graph, x, y, z);

      // find out the closest atom
      const auto [min_distance, min_index] = find_closest_neighbor(v, graph, x, y, z);

      // assign the correct type based on angle and distance
      if (angle < coordinate_type{114.8}) {
        if ((min_distance < coordinate_type{1.38} && elements[min_index] == element::C) ||
            (min_distance < coordinate_type{1.32} && elements[min_index] == element::N)) {
          return autodock_babel_ff::Npl;
        } else {
          return autodock_babel_ff::N3;
        }
      } else if (angle < coordinate_type{160}) {
        return autodock_babel_ff::Npl;
      } else {
        return autodock_babel_ff::N1;
      }
    }
    return autodock_babel_ff::N;
  }

  static auto handleO(const molecule_graph_type& graph, const molecule_graph_type::vertex_descriptor v)
      -> autodock_babel_ff {
    const auto num_neighbors = boost::out_degree(v, graph);
    if (num_neighbors == 2) {
      return autodock_babel_ff::O3;
    }
    return autodock_babel_ff::O;
  }

  static auto handleP(const std::span<element>& elements,
                      const molecule_graph_type& graph,
                      const molecule_graph_type::vertex_descriptor v) -> autodock_babel_ff {
    const auto num_neighbors = boost::out_degree(v, graph);
    if (num_neighbors == 4) {
      const auto num_free_ox = count_neighbors_if(v, graph, is_free_ox{elements, graph});
      if (num_free_ox >= 2) {
        return autodock_babel_ff::Pac;
      } else if (num_free_ox == 1) {
        return autodock_babel_ff::Pox;
      } else {
        return autodock_babel_ff::P3_plus;
      }
    }
    return autodock_babel_ff::P;
  }

  static auto handleS(const std::span<element>& elements,
                      const molecule_graph_type& graph,
                      const molecule_graph_type::vertex_descriptor v) -> autodock_babel_ff {
    const auto num_neighbors = boost::out_degree(v, graph);
    const auto num_free_ox   = count_neighbors_if(v, graph, is_free_ox{elements, graph});
    if (num_neighbors == 4) {
      if (num_free_ox >= 3) {
        return autodock_babel_ff::Sac;
      } else if (num_free_ox >= 1) {
        return autodock_babel_ff::Sox;
      } else {
        return autodock_babel_ff::S;
      }
    } else if (num_neighbors == 3) {
      if (num_free_ox >= 1) {
        return autodock_babel_ff::Sox;
      } else {
        return autodock_babel_ff::S3_plus;
      }
    } else if (num_neighbors == 2) {
      return autodock_babel_ff::S3;
    }
    return autodock_babel_ff::S;
  }

  //===------------------------------------------------------------------------------------------------------
  // The implementation of the actual function that assign the atom typing
  //===------------------------------------------------------------------------------------------------------

  void assign_autodock_babel_types(std::span<autodock_babel_ff> types,
                                   const std::span<coordinate_type> x,
                                   const std::span<coordinate_type> y,
                                   const std::span<coordinate_type> z,
                                   const std::span<element> elements,
                                   const molecule_graph_type& graph) {
    // the sensing scheme is not straightforward. It starts by assign to each atom its element type, then it
    // has some rules to further specify atom types. It basically is a porting from:
    // AutoDockTools/MGLToolsPckgs/PyBabel/atomTypes.py
    assert(types.size() == elements.size());

    // perform the initial assignment of the neighbors when they have more than one neighbor
    const auto elements2babel = [&](const molecule_graph_type::vertex_descriptor v) {
      switch (elements[graph[v].atom_index]) {
        case (element::H): return handleH(elements, graph, v);
        case (element::He): return autodock_babel_ff::He;
        case (element::Li): return autodock_babel_ff::Li;
        case (element::Be): return autodock_babel_ff::Be;
        case (element::B): return autodock_babel_ff::B; // the original code is dead code
        case (element::C): return handleC(elements, x, y, z, graph, v);
        case (element::N): return handleN(elements, x, y, z, graph, v);
        case (element::O): return handleO(graph, v);
        case (element::F): return autodock_babel_ff::F;
        case (element::Ne): return autodock_babel_ff::Ne;
        case (element::Na): return autodock_babel_ff::Na;
        case (element::Mg): return autodock_babel_ff::Mg;
        case (element::Al): return autodock_babel_ff::Al;
        case (element::Si): return autodock_babel_ff::Si;
        case (element::P): return handleP(elements, graph, v);
        case (element::S): return handleS(elements, graph, v);
        case (element::Cl): return autodock_babel_ff::Cl;
        case (element::Ar): return autodock_babel_ff::Ar;
        case (element::K): return autodock_babel_ff::K;
        case (element::Ca): return autodock_babel_ff::Ca;
        case (element::Sc): return autodock_babel_ff::Sc;
        case (element::Ti): return autodock_babel_ff::Ti;
        case (element::V): return autodock_babel_ff::V;
        case (element::Cr): return autodock_babel_ff::Cr;
        case (element::Mn): return autodock_babel_ff::Mn;
        case (element::Fe): return autodock_babel_ff::Fe;
        case (element::Co): return autodock_babel_ff::Co;
        case (element::Ni): return autodock_babel_ff::Ni;
        case (element::Cu): return autodock_babel_ff::Cu;
        case (element::Zn): return autodock_babel_ff::Zn;
        case (element::Ga): return autodock_babel_ff::Ga;
        case (element::Ge): return autodock_babel_ff::Ge;
        case (element::As): return autodock_babel_ff::As;
        case (element::Se): return autodock_babel_ff::Se;
        case (element::Br): return autodock_babel_ff::Br;
        case (element::Kr): return autodock_babel_ff::Kr;
        case (element::Rb): return autodock_babel_ff::Rb;
        case (element::Sr): return autodock_babel_ff::Sr;
        case (element::Y): return autodock_babel_ff::Y;
        case (element::Zr): return autodock_babel_ff::Zr;
        case (element::Nb): return autodock_babel_ff::Nb;
        case (element::Mo): return autodock_babel_ff::Mo;
        case (element::Tc): return autodock_babel_ff::Tc;
        case (element::Ru): return autodock_babel_ff::Ru;
        case (element::Rh): return autodock_babel_ff::Rh;
        case (element::Pd): return autodock_babel_ff::Pd;
        case (element::Ag): return autodock_babel_ff::Ag;
        case (element::Cd): return autodock_babel_ff::Cd;
        case (element::In): return autodock_babel_ff::In;
        case (element::Sn): return autodock_babel_ff::Sn;
        case (element::Sb): return autodock_babel_ff::Sb;
        case (element::Te): return autodock_babel_ff::Te;
        case (element::I): return autodock_babel_ff::I;
        case (element::Xe): return autodock_babel_ff::Xe;
        case (element::Cs): return autodock_babel_ff::Cs;
        case (element::Ba): return autodock_babel_ff::Ba;
        case (element::La): return autodock_babel_ff::La;
        case (element::Ce): return autodock_babel_ff::Ce;
        case (element::Pr): return autodock_babel_ff::Pr;
        case (element::Nd): return autodock_babel_ff::Nd;
        case (element::Pm): return autodock_babel_ff::Pm;
        case (element::Sm): return autodock_babel_ff::Sm;
        case (element::Eu): return autodock_babel_ff::Eu;
        case (element::Gd): return autodock_babel_ff::Gd;
        case (element::Tb): return autodock_babel_ff::Tb;
        case (element::Dy): return autodock_babel_ff::Dy;
        case (element::Ho): return autodock_babel_ff::Ho;
        case (element::Er): return autodock_babel_ff::Er;
        case (element::Tm): return autodock_babel_ff::Tm;
        case (element::Yb): return autodock_babel_ff::Yb;
        case (element::Lu): return autodock_babel_ff::Lu;
        case (element::Hf): return autodock_babel_ff::Hf;
        case (element::Ta): return autodock_babel_ff::Ta;
        case (element::W): return autodock_babel_ff::W;
        case (element::Re): return autodock_babel_ff::Re;
        case (element::Os): return autodock_babel_ff::Os;
        case (element::Ir): return autodock_babel_ff::Ir;
        case (element::Pt): return autodock_babel_ff::Pt;
        case (element::Au): return autodock_babel_ff::Au;
        case (element::Hg): return autodock_babel_ff::Hg;
        case (element::Tl): return autodock_babel_ff::Tl;
        case (element::Pb): return autodock_babel_ff::Pb;
        case (element::Bi): return autodock_babel_ff::Bi;
        case (element::Po): return autodock_babel_ff::Po;
        case (element::At): return autodock_babel_ff::At;
        case (element::Rn): return autodock_babel_ff::Rn;
        case (element::Fr): return autodock_babel_ff::Fr;
        case (element::Ra): return autodock_babel_ff::Ra;
        case (element::Ac): return autodock_babel_ff::Ac;
        case (element::Th): return autodock_babel_ff::Th;
        case (element::Pa): return autodock_babel_ff::Pa;
        case (element::U): return autodock_babel_ff::U;
        case (element::Np): return autodock_babel_ff::Np;
        case (element::Pu): return autodock_babel_ff::Pu;
        case (element::Am): return autodock_babel_ff::Am;
        case (element::Cm): return autodock_babel_ff::Cm;
        case (element::Bk): return autodock_babel_ff::Bk;
        case (element::Cf): return autodock_babel_ff::Cf;
        case (element::Es): return autodock_babel_ff::Es;
        case (element::Fm): return autodock_babel_ff::Fm;
        case (element::Md): return autodock_babel_ff::Md;
        case (element::No): return autodock_babel_ff::No;
        case (element::Lr): return autodock_babel_ff::Lr;
        case (element::Rf): return autodock_babel_ff::Rf;
        case (element::Db): return autodock_babel_ff::Db;
        case (element::Sg): return autodock_babel_ff::Sg;
        case (element::Bh): return autodock_babel_ff::Bh;
        case (element::Hs): return autodock_babel_ff::Hs;
        case (element::Mt): return autodock_babel_ff::Mt;
        case (element::Ds): return autodock_babel_ff::Ds;
        case (element::Rg): return autodock_babel_ff::Rg;
        case (element::Cn): return autodock_babel_ff::Cn;
        case (element::Nh): return autodock_babel_ff::Nh;
        case (element::Fl): return autodock_babel_ff::Fl;
        case (element::Mc): return autodock_babel_ff::Mc;
        case (element::Lv): return autodock_babel_ff::Lv;
        case (element::Ts): return autodock_babel_ff::Ts;
        case (element::Og): return autodock_babel_ff::Og;
        case (element::Uue): return autodock_babel_ff::Uue;
        default: throw std::runtime_error("Unknown element, internal error");
      };
    };
    const auto [vertex_begin, vertex_end] = boost::vertices(graph);
    std::transform(vertex_begin, vertex_end, std::begin(types), elements2babel);

    // handle the atoms with "valence_one"
    // NOTE: we need to do it after the main switch since it depends on types computed previously
    for (auto it = vertex_begin; it != vertex_end; ++it) {
      const auto vertex           = *it;
      const auto vertex_index     = graph[vertex].atom_index;
      const auto vertex_atom_type = elements[vertex_index];
      const auto vertex_point     = point3D{x[vertex_index], y[vertex_index], z[vertex_index]};

      const auto num_edges = boost::out_degree(vertex, graph);
      if (num_edges == 1) {
        const auto [neigh_begin, neigh_end] = boost::out_edges(vertex, graph);
        const auto neighbor_vertex          = boost::target(*neigh_begin, graph);
        const auto neighbor_index           = graph[neighbor_vertex].atom_index;
        const auto neighbor_point = point3D{x[neighbor_index], y[neighbor_index], z[neighbor_index]};

        switch (vertex_atom_type) {
          case (element::C):
            if (types[neighbor_index] == autodock_babel_ff::C1 &&
                distance2(vertex_point, neighbor_point) < square(coordinate_type{1.22}))
              types[vertex_index] = autodock_babel_ff::C1;
            else if (elements[neighbor_index] == element::C &&
                     distance2(vertex_point, neighbor_point) < square(coordinate_type{1.41}))
              types[vertex_index] = autodock_babel_ff::C1;
            else if (elements[neighbor_index] == element::N &&
                     distance2(vertex_point, neighbor_point) < square(coordinate_type{1.37}))
              types[vertex_index] = autodock_babel_ff::C2;
            else
              types[vertex_index] = autodock_babel_ff::C3;
            break;

          case (element::N):
            if (types[neighbor_index] == autodock_babel_ff::C1 &&
                distance2(vertex_point, neighbor_point) < square(coordinate_type{1.20}))
              types[vertex_index] = autodock_babel_ff::N1;
            else if ((types[neighbor_index] == autodock_babel_ff::C2 ||
                      types[neighbor_index] == autodock_babel_ff::C3) &&
                     distance2(vertex_point, neighbor_point) > square(coordinate_type{1.38}))
              types[vertex_index] = autodock_babel_ff::N3;
            else
              types[vertex_index] = autodock_babel_ff::Npl;
            break;

          case (element::O):
            if (types[neighbor_index] == autodock_babel_ff::Cac ||
                types[neighbor_index] == autodock_babel_ff::Pac ||
                types[neighbor_index] == autodock_babel_ff::Sac ||
                types[neighbor_index] == autodock_babel_ff::Ntr)
              types[vertex_index] = autodock_babel_ff::O_minus;
            else if (types[neighbor_index] == autodock_babel_ff::Nox ||
                     types[neighbor_index] == autodock_babel_ff::Pox ||
                     types[neighbor_index] == autodock_babel_ff::Sox)
              types[vertex_index] = autodock_babel_ff::O2;
            else if (elements[neighbor_index] == element::C &&
                     distance2(vertex_point, neighbor_point) > square(coordinate_type{1.30})) {
              types[vertex_index]   = autodock_babel_ff::O2;
              types[neighbor_index] = autodock_babel_ff::C2;
            } else if (elements[neighbor_index] == element::As &&
                       distance2(vertex_point, neighbor_point) > square(coordinate_type{1.685}))
              types[vertex_index] = autodock_babel_ff::O2;
            else
              types[vertex_index] = autodock_babel_ff::O3;
            break;

          case (element::S):
            if (elements[neighbor_index] == element::P)
              types[vertex_index] = autodock_babel_ff::S2;
            else if (elements[neighbor_index] == element::C &&
                     distance2(vertex_point, neighbor_point) > square(coordinate_type{1.76})) {
              types[vertex_index]   = autodock_babel_ff::S2;
              types[neighbor_index] = autodock_babel_ff::C2;
            } else if (elements[neighbor_index] == element::As &&
                       distance2(vertex_point, neighbor_point) > square(coordinate_type{2.11}))
              types[vertex_index] = autodock_babel_ff::S2;
            else
              types[vertex_index] = autodock_babel_ff::S3;
            break;

          default: break; // do nothing
        };
      }
    }

    // post process the C2 to make sure that they are not C3
    static constexpr std::array<autodock_babel_ff, 12> phase5_types = {{autodock_babel_ff::C3,
                                                                        autodock_babel_ff::HC,
                                                                        autodock_babel_ff::N3,
                                                                        autodock_babel_ff::N3_plus,
                                                                        autodock_babel_ff::O3,
                                                                        autodock_babel_ff::Pac,
                                                                        autodock_babel_ff::Sac,
                                                                        autodock_babel_ff::Sox,
                                                                        autodock_babel_ff::C1,
                                                                        autodock_babel_ff::S3,
                                                                        autodock_babel_ff::Cac}};
    const auto is_phase5_interesting                                = [&](const std::size_t i) {
      return std::any_of(std::begin(phase5_types), std::end(phase5_types), [&](const auto ff_type) {
        return types[i] == ff_type;
      });
    };
    for (auto it = vertex_begin; it != vertex_end; ++it) {
      const auto vertex       = *it;
      const auto vertex_index = graph[vertex].atom_index;
      if (types[vertex_index] == autodock_babel_ff::C2) {
        const auto [edge_begin, edge_end] = boost::out_edges(vertex, graph);
        if (std::all_of(edge_begin, edge_end, [&](const auto edge) {
              return is_phase5_interesting(graph[boost::target(edge, graph)].atom_index);
            })) {
          types[vertex_index] = autodock_babel_ff::C3;
        }
      }
    }

    // post-process the nitrogen and carbon to see if they are some special case
    for (auto it = vertex_begin; it != vertex_end; ++it) {
      const auto vertex       = *it;
      const auto vertex_index = graph[vertex].atom_index;
      if (types[vertex_index] == autodock_babel_ff::N3) {
        auto protonated                   = true;
        const auto [edge_begin, edge_end] = boost::out_edges(vertex, graph);
        for (auto neigh_it = edge_begin; neigh_it != edge_end; ++neigh_it) {
          const auto neigh       = *neigh_it;
          const auto neigh_index = graph[boost::target(neigh, graph)].atom_index;
          const auto neigh_type  = types[neigh_index];
          if (boost::out_degree(vertex, graph) == 2 &&
              (neigh_type == autodock_babel_ff::C2 || neigh_type == autodock_babel_ff::Sox ||
               neigh_type == autodock_babel_ff::Sac || neigh_type == autodock_babel_ff::Pac)) {
            protonated          = false;
            types[vertex_index] = autodock_babel_ff::Npl;
            break;
          } else if (neigh_type != autodock_babel_ff::C3 && elements[neigh_index] != element::H) {
            protonated = false; // we cannot break otherwise we can alter the outcome
          }
        }
        if (protonated) {
          types[vertex_index] = autodock_babel_ff::N3_plus;
        }
      } else if (types[vertex_index] == autodock_babel_ff::C2) {
        const auto [edge_begin, edge_end] = boost::out_edges(vertex, graph);
        const auto m                      = std::count_if(edge_begin, edge_end, [&](const auto edge) {
          const auto neigh_type = types[graph[boost::target(edge, graph)].atom_index];
          return neigh_type == autodock_babel_ff::Npl || neigh_type == autodock_babel_ff::N2 ||
                 neigh_type == autodock_babel_ff::Ng_plus;
        });
        if (m == 3) {
          types[vertex_index] = autodock_babel_ff::C_plus;
          std::for_each(edge_begin, edge_end, [&](const auto edge) {
            types[boost::target(edge, graph)] = autodock_babel_ff::Ng_plus;
          });
        }
      } else if (types[vertex_index] == autodock_babel_ff::Cac) {
        const auto [edge_begin, edge_end] = boost::out_edges(vertex, graph);
        for (auto neigh_it = edge_begin; neigh_it != edge_end; ++neigh_it) {
          const auto neigh       = boost::target(*neigh_it, graph);
          const auto neigh_index = graph[neigh].atom_index;
          if (elements[neigh_index] == element::O &&
              count_neighbors_if(neigh, graph, is_heavy_edge{elements, graph}) == 1) {
            types[neigh_index] = autodock_babel_ff::O_minus;
          }
        }
      }
    }

    // look for amid nitrogens
    for (auto vertex_it = vertex_begin; vertex_it != vertex_end; ++vertex_it) {
      const auto vertex       = *vertex_it;
      const auto vertex_index = graph[vertex].atom_index;
      if (types[vertex_index] == autodock_babel_ff::Npl) {
        const auto [edge_begin, edge_end] = boost::out_edges(vertex, graph);
        for (auto neigh_it = edge_begin; neigh_it != edge_end; ++neigh_it) {
          const auto neigh       = *neigh_it;
          const auto neigh_index = graph[boost::target(neigh, graph)].atom_index;
          const auto neigh_type  = types[neigh_index];
          if (neigh_type == autodock_babel_ff::Cac || neigh_type == autodock_babel_ff::Sox) {
            types[vertex_index] = autodock_babel_ff::Nam;
            break;
          } else if (neigh_type == autodock_babel_ff::C2) {
            const auto [neigh_edge_begin, neigh_edge_end] = boost::out_edges(vertex, graph);
            auto is_interesting                           = false;
            for (auto it = neigh_edge_begin; it != neigh_edge_end; ++it) {
              const auto type = types[graph[boost::target(*it, graph)].atom_index];
              if (type == autodock_babel_ff::O2 || type == autodock_babel_ff::S2) {
                is_interesting = true;
                break;
              }
            }
            if (is_interesting) {
              types[vertex_index] = autodock_babel_ff::Nam;
              break;
            }
          }
        }
      }
    }

    // if we reach this point it means that we have correctly categorized the atoms in the autodock babel type
  }

} // namespace mudock
