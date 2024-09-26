#include <mudock/chem/assign_autodock_types.hpp>

namespace mudock {

  static auto handleH(const std::span<element>& elements,
                      const molecule_graph_type& graph,
                      const molecule_graph_type::vertex_descriptor v) -> autodock_ff {
    assert(elements[graph[v].atom_index] == element::H);
    const auto num_neighbors = boost::out_degree(v, graph);
    if (num_neighbors) {
      molecule_graph_type::adjacency_iterator ai, ai_end;
      boost::tie(ai, ai_end) = adjacent_vertices(v, graph);
      // num_neighbors greater than 0, no needs to check if ai!=ai_end
      // IS AN HYDROGEN ATOM -> ONLY ON BOND
      // Element type A -> Should be an aromatic carbon -> babel_type does not exists -> ASSUMPTION no check required
      // AutoDockTools/atomTypeTools.py:378
      if (elements[graph[*ai].atom_index] != element::C)
        return autodock_ff::HD;
    } else {
      return autodock_ff::HD;
    }
    return autodock_ff::H;
    // TODO missing HS case
  }

  static auto handleN([[maybe_unused]] const std::span<element>& elements,
                      const std::span<autodock_babel_ff>& babel_type,
                      const molecule_graph_type& graph,
                      const molecule_graph_type::vertex_descriptor v) -> autodock_ff {
    assert(elements[graph[v].atom_index] == element::N);
    const auto num_neighbors = boost::out_degree(v, graph);
    if (babel_type[graph[v].atom_index] == autodock_babel_ff::N3 ||
        babel_type[graph[v].atom_index] == autodock_babel_ff::N2 ||
        babel_type[graph[v].atom_index] == autodock_babel_ff::N1 ||
        ((babel_type[graph[v].atom_index] == autodock_babel_ff::Nam ||
          babel_type[graph[v].atom_index] == autodock_babel_ff::Npl ||
          babel_type[graph[v].atom_index] == autodock_babel_ff::Ng_plus) &&
         num_neighbors == 2)) {
      return autodock_ff::NA;
    }
    return autodock_ff::N;
    // TODO missing NS case
  }

  static auto handleS([[maybe_unused]] const std::span<element>& elements,
                      const std::span<autodock_babel_ff>& babel_type,
                      const molecule_graph_type& graph,
                      const molecule_graph_type::vertex_descriptor v) -> autodock_ff {
    assert(elements[graph[v].atom_index] == element::S);
    if (babel_type[graph[v].atom_index] != autodock_babel_ff::Sox ||
        babel_type[graph[v].atom_index] != autodock_babel_ff::Sac) {
      return autodock_ff::SA;
    }
    return autodock_ff::S;
  }

  static auto handleC([[maybe_unused]] const std::span<element>& elements,
                      const molecule_graph_type& graph,
                      const std::span<int> is_aromatic,
                      const molecule_graph_type::vertex_descriptor v) -> autodock_ff {
    assert(elements[graph[v].atom_index] == element::C);
    if (is_aromatic[graph[v].atom_index])
      return autodock_ff::A;
    return autodock_ff::C;
  }

  void assign_autodock_types(std::span<autodock_ff> types,
                             const std::span<element> elements,
                             const std::span<int> is_aromatic,
                             const std::span<autodock_babel_ff> babel_type,
                             const molecule_graph_type& graph) {
    // It basically is a porting from:
    // AutoDockTools/atomTypeTools.py
    assert(types.size() == elements.size());

    const auto elements2autodock = [&](const molecule_graph_type::vertex_descriptor v) {
      switch (elements[graph[v].atom_index]) {
        case (element::H): return handleH(elements, graph, v);
        case (element::He): return autodock_ff::He;
        case (element::Li): return autodock_ff::Li;
        case (element::Be): return autodock_ff::Be;
        case (element::B): return autodock_ff::B; // the original code is dead code
        case (element::C): return handleC(elements, graph, is_aromatic, v);
        case (element::N): return handleN(elements, babel_type, graph, v);
        case (element::O): return autodock_ff::OA;
        // TODO missing OS case
        case (element::F): return autodock_ff::F;
        case (element::Ne): return autodock_ff::Ne;
        case (element::Na): return autodock_ff::Na;
        case (element::Mg): return autodock_ff::Mg;
        case (element::Al): return autodock_ff::Al;
        case (element::Si): return autodock_ff::Si;
        case (element::P): return autodock_ff::P;
        case (element::S): return handleS(elements, babel_type, graph, v);
        case (element::Cl): return autodock_ff::Cl;
        case (element::Ar): return autodock_ff::Ar;
        case (element::K): return autodock_ff::K;
        case (element::Ca): return autodock_ff::Ca;
        case (element::Sc): return autodock_ff::Sc;
        case (element::Ti): return autodock_ff::Ti;
        case (element::V): return autodock_ff::V;
        case (element::Cr): return autodock_ff::Cr;
        case (element::Mn): return autodock_ff::Mn;
        case (element::Fe): return autodock_ff::Fe;
        case (element::Co): return autodock_ff::Co;
        case (element::Ni): return autodock_ff::Ni;
        case (element::Cu): return autodock_ff::Cu;
        case (element::Zn): return autodock_ff::Zn;
        case (element::Ga): return autodock_ff::Ga;
        case (element::Ge): return autodock_ff::Ge;
        case (element::As): return autodock_ff::As;
        case (element::Se): return autodock_ff::Se;
        case (element::Br): return autodock_ff::Br;
        case (element::Kr): return autodock_ff::Kr;
        case (element::Rb): return autodock_ff::Rb;
        case (element::Sr): return autodock_ff::Sr;
        case (element::Y): return autodock_ff::Y;
        case (element::Zr): return autodock_ff::Zr;
        case (element::Nb): return autodock_ff::Nb;
        case (element::Mo): return autodock_ff::Mo;
        case (element::Tc): return autodock_ff::Tc;
        case (element::Ru): return autodock_ff::Ru;
        case (element::Rh): return autodock_ff::Rh;
        case (element::Pd): return autodock_ff::Pd;
        case (element::Ag): return autodock_ff::Ag;
        case (element::Cd): return autodock_ff::Cd;
        case (element::In): return autodock_ff::In;
        case (element::Sn): return autodock_ff::Sn;
        case (element::Sb): return autodock_ff::Sb;
        case (element::Te): return autodock_ff::Te;
        case (element::I): return autodock_ff::I;
        case (element::Xe): return autodock_ff::Xe;
        case (element::Cs): return autodock_ff::Cs;
        case (element::Ba): return autodock_ff::Ba;
        case (element::La): return autodock_ff::La;
        case (element::Ce): return autodock_ff::Ce;
        case (element::Pr): return autodock_ff::Pr;
        case (element::Nd): return autodock_ff::Nd;
        case (element::Pm): return autodock_ff::Pm;
        case (element::Sm): return autodock_ff::Sm;
        case (element::Eu): return autodock_ff::Eu;
        case (element::Gd): return autodock_ff::Gd;
        case (element::Tb): return autodock_ff::Tb;
        case (element::Dy): return autodock_ff::Dy;
        case (element::Ho): return autodock_ff::Ho;
        case (element::Er): return autodock_ff::Er;
        case (element::Tm): return autodock_ff::Tm;
        case (element::Yb): return autodock_ff::Yb;
        case (element::Lu): return autodock_ff::Lu;
        case (element::Hf): return autodock_ff::Hf;
        case (element::Ta): return autodock_ff::Ta;
        case (element::W): return autodock_ff::W;
        case (element::Re): return autodock_ff::Re;
        case (element::Os): return autodock_ff::Os;
        case (element::Ir): return autodock_ff::Ir;
        case (element::Pt): return autodock_ff::Pt;
        case (element::Au): return autodock_ff::Au;
        case (element::Hg): return autodock_ff::Hg;
        case (element::Tl): return autodock_ff::Tl;
        case (element::Pb): return autodock_ff::Pb;
        case (element::Bi): return autodock_ff::Bi;
        case (element::Po): return autodock_ff::Po;
        case (element::At): return autodock_ff::At;
        case (element::Rn): return autodock_ff::Rn;
        case (element::Fr): return autodock_ff::Fr;
        case (element::Ra): return autodock_ff::Ra;
        case (element::Ac): return autodock_ff::Ac;
        case (element::Th): return autodock_ff::Th;
        case (element::Pa): return autodock_ff::Pa;
        case (element::U): return autodock_ff::U;
        case (element::Np): return autodock_ff::Np;
        case (element::Pu): return autodock_ff::Pu;
        case (element::Am): return autodock_ff::Am;
        case (element::Cm): return autodock_ff::Cm;
        case (element::Bk): return autodock_ff::Bk;
        case (element::Cf): return autodock_ff::Cf;
        case (element::Es): return autodock_ff::Es;
        case (element::Fm): return autodock_ff::Fm;
        case (element::Md): return autodock_ff::Md;
        case (element::No): return autodock_ff::No;
        case (element::Lr): return autodock_ff::Lr;
        case (element::Rf): return autodock_ff::Rf;
        case (element::Db): return autodock_ff::Db;
        case (element::Sg): return autodock_ff::Sg;
        case (element::Bh): return autodock_ff::Bh;
        case (element::Hs): return autodock_ff::Hs;
        case (element::Mt): return autodock_ff::Mt;
        case (element::Ds): return autodock_ff::Ds;
        case (element::Rg): return autodock_ff::Rg;
        case (element::Cn): return autodock_ff::Cn;
        case (element::Nh): return autodock_ff::Nh;
        case (element::Fl): return autodock_ff::Fl;
        case (element::Mc): return autodock_ff::Mc;
        case (element::Lv): return autodock_ff::Lv;
        case (element::Ts): return autodock_ff::Ts;
        case (element::Og): return autodock_ff::Og;
        case (element::Uue): return autodock_ff::Uue;
        // TODO missing case for Z, G, GA, J, Q
        default: throw std::runtime_error("Unknown element, internal error");
      };
    };
    const auto [vertex_begin, vertex_end] = boost::vertices(graph);
    for (auto it = vertex_begin; it != vertex_end; ++it) {
      types[graph[*it].atom_index] = elements2autodock(*it);
    }
  }

} // namespace mudock
