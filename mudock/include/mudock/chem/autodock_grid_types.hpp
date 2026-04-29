#pragma once

#include <mudock/chem/autodock_types.hpp>
#include <stdexcept>
#include <string>

namespace mudock {
  enum class autodock_grid_type : int {
    A = 0,
    C,
    H,
    HD,
    N,
    NA,
    OA,
    SA,
    Cl,
    F,
    S,
    Si,
    Br,
    P,
    I,
    ELEC,
    DESOLV
  };

  static constexpr int num_autodock_ff_grids() { return 15; }
  static constexpr int num_autodock_grids() { return num_autodock_ff_grids() + 2; }

  inline autodock_grid_type autodock_grid_from_ff(const autodock_ff& autodock_type) {
    switch (autodock_type) {
      case autodock_ff::A: return autodock_grid_type::A;
      case autodock_ff::C: return autodock_grid_type::C;
      case autodock_ff::H: return autodock_grid_type::H;
      case autodock_ff::HD: return autodock_grid_type::HD;
      case autodock_ff::N: return autodock_grid_type::N;
      case autodock_ff::NA: return autodock_grid_type::NA;
      case autodock_ff::OA: return autodock_grid_type::OA;
      case autodock_ff::SA: return autodock_grid_type::SA;
      case autodock_ff::Cl: return autodock_grid_type::Cl;
      case autodock_ff::F: return autodock_grid_type::F;
      case autodock_ff::S: return autodock_grid_type::S;
      case autodock_ff::Si: return autodock_grid_type::Si;
      case autodock_ff::Br: return autodock_grid_type::Br;
      case autodock_ff::P: return autodock_grid_type::P;
      case autodock_ff::I: return autodock_grid_type::I;
      default:
        const auto error_msg =
            std::string("Missing Autodock Grid Type ") + std::string(get_description(autodock_type).name);
        throw std::runtime_error(error_msg);
    }
  }

  inline autodock_ff autodock_ff_from_grid(const autodock_grid_type& ligand_map) {
    switch (ligand_map) {
      case autodock_grid_type::A: return autodock_ff::A;
      case autodock_grid_type::C: return autodock_ff::C;
      case autodock_grid_type::H: return autodock_ff::H;
      case autodock_grid_type::HD: return autodock_ff::HD;
      case autodock_grid_type::N: return autodock_ff::N;
      case autodock_grid_type::NA: return autodock_ff::NA;
      case autodock_grid_type::OA: return autodock_ff::OA;
      case autodock_grid_type::SA: return autodock_ff::SA;
      case autodock_grid_type::Cl: return autodock_ff::Cl;
      case autodock_grid_type::F: return autodock_ff::F;
      case autodock_grid_type::S: return autodock_ff::S;
      case autodock_grid_type::Si: return autodock_ff::Si;
      case autodock_grid_type::Br: return autodock_ff::Br;
      case autodock_grid_type::P: return autodock_ff::P;
      case autodock_grid_type::I: return autodock_ff::I;
      default: throw std::runtime_error("Missing Autodock Force Field Type");
    }
  }

  struct map_description {
    autodock_grid_type value;
    std::string_view name;
  };
  extern const std::array<map_description, num_autodock_grids()> MAP_DICTIONARY;

  inline const map_description& get_description(const autodock_grid_type e) {
    assert(MAP_DICTIONARY[static_cast<size_t>(e)].value == e);
    return MAP_DICTIONARY[static_cast<size_t>(e)];
  }
  autodock_grid_type parse_map_symbol(const std::string_view symbol);
} // namespace mudock
