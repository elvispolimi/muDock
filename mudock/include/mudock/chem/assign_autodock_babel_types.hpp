#pragma once

#include <mudock/chem/autodock_babel_types.hpp>
#include <mudock/chem/elements.hpp>
#include <mudock/molecule.hpp>
#include <mudock/molecule/graph.hpp>
#include <mudock/type_alias.hpp>
#include <span>

namespace mudock {

  void assign_autodock_babel_types(std::span<autodock_babel_ff> types,
                                   const std::span<coordinate_type> x,
                                   const std::span<coordinate_type> y,
                                   const std::span<coordinate_type> z,
                                   const std::span<element> elements,
                                   const molecule_graph_type& graph);

}
