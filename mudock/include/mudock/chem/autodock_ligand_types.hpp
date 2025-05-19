#pragma once

#include <mudock/chem/autodock_types.hpp>
#include <stdexcept>

namespace mudock {
  enum class autodock_ligand_ff : int { A = 0, C, H, HD, N, NA, OA, SA, Cl, F, S, Br, P, I };

  static constexpr int num_autodock_ligand_types() { return 14; }

  inline autodock_ligand_ff autodock_ligand_from_type(const autodock_ff& autodock_type) {
    switch (autodock_type) {
      case autodock_ff::A: return autodock_ligand_ff::A;
      case autodock_ff::C: return autodock_ligand_ff::C;
      case autodock_ff::H: return autodock_ligand_ff::H;
      case autodock_ff::HD: return autodock_ligand_ff::HD;
      case autodock_ff::N: return autodock_ligand_ff::N;
      case autodock_ff::NA: return autodock_ligand_ff::NA;
      case autodock_ff::OA: return autodock_ligand_ff::OA;
      case autodock_ff::SA: return autodock_ligand_ff::SA;
      case autodock_ff::Cl: return autodock_ligand_ff::Cl;
      case autodock_ff::F: return autodock_ligand_ff::F;
      case autodock_ff::S: return autodock_ligand_ff::S;
      case autodock_ff::Br: return autodock_ligand_ff::Br;
      case autodock_ff::P: return autodock_ligand_ff::P;
      case autodock_ff::I: return autodock_ligand_ff::I;
      default: throw std::runtime_error("Missing Autodock ligand type");
    }
  }

  inline autodock_ff autodock_type_from_ligand(const autodock_ligand_ff& ligand_map) {
    switch (ligand_map) {
      case autodock_ligand_ff::A: return autodock_ff::A;
      case autodock_ligand_ff::C: return autodock_ff::C;
      case autodock_ligand_ff::H: return autodock_ff::H;
      case autodock_ligand_ff::HD: return autodock_ff::HD;
      case autodock_ligand_ff::N: return autodock_ff::N;
      case autodock_ligand_ff::NA: return autodock_ff::NA;
      case autodock_ligand_ff::OA: return autodock_ff::OA;
      case autodock_ligand_ff::SA: return autodock_ff::SA;
      case autodock_ligand_ff::Cl: return autodock_ff::Cl;
      case autodock_ligand_ff::F: return autodock_ff::F;
      case autodock_ligand_ff::S: return autodock_ff::S;
      case autodock_ligand_ff::Br: return autodock_ff::Br;
      case autodock_ligand_ff::P: return autodock_ff::P;
      case autodock_ligand_ff::I: return autodock_ff::I;
      default: throw std::runtime_error("Missing Autodock type");
    }
  }

  struct map_description {
    autodock_ligand_ff value;
    std::string_view name;
  };
  extern const std::array<map_description, num_autodock_ligand_types()> MAP_DICTIONARY;

  inline const map_description& get_description(const autodock_ligand_ff e) {
    assert(MAP_DICTIONARY[static_cast<int>(e)].value == e);
    return MAP_DICTIONARY[static_cast<int>(e)];
  }
  autodock_ligand_ff parse_map_symbol(const std::string_view symbol);
} // namespace mudock
