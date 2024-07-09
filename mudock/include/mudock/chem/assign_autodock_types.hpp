#pragma once

#include <mudock/chem/autodock_babel_types.hpp>
#include <mudock/chem/autodock_types.hpp>
#include <mudock/chem/elements.hpp>
#include <mudock/molecule.hpp>
#include <mudock/molecule/graph.hpp>
#include <mudock/type_alias.hpp>
#include <span>

namespace mudock {

  void assign_autodock_types(std::span<autodock_ff> types,
                             const std::span<element> elements,
                             const std::span<int> is_aromatic,
                             const std::span<autodock_babel_ff> babel_type,
                             const molecule_graph_type& graph);

}
