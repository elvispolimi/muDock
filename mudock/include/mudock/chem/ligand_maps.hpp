#pragma once

#include <mudock/chem/autodock_types.hpp>
#include <stdexcept>

namespace mudock {
  // TODO check from tests which one are needed
  enum class ligand_map_types : int { A = 0, C, H, HD, N, NA, OA, SA, Cl, F, S, Br, P, I };

  static constexpr int num_ligand_map_types() { return 14; }

  inline ligand_map_types map_from_autodock_type(const autodock_ff& autodock_type) {
    switch (autodock_type) {
      case autodock_ff::A: return ligand_map_types::A;
      case autodock_ff::C: return ligand_map_types::C;
      case autodock_ff::H: return ligand_map_types::H;
      case autodock_ff::HD: return ligand_map_types::HD;
      case autodock_ff::N: return ligand_map_types::N;
      case autodock_ff::NA: return ligand_map_types::NA;
      case autodock_ff::OA: return ligand_map_types::OA;
      case autodock_ff::SA: return ligand_map_types::SA;
      case autodock_ff::Cl: return ligand_map_types::Cl;
      case autodock_ff::F: return ligand_map_types::F;
      case autodock_ff::S: return ligand_map_types::S;
      case autodock_ff::Br: return ligand_map_types::Br;
      case autodock_ff::P: return ligand_map_types::P;
      case autodock_ff::I: return ligand_map_types::I;
      default: throw std::runtime_error("Missing map texture");
    }
  }

  inline autodock_ff autodock_type_from_map(const ligand_map_types& ligand_map) {
    switch (ligand_map) {
      case ligand_map_types::A: return autodock_ff::A;
      case ligand_map_types::C: return autodock_ff::C;
      case ligand_map_types::H: return autodock_ff::H;
      case ligand_map_types::HD: return autodock_ff::HD;
      case ligand_map_types::N: return autodock_ff::N;
      case ligand_map_types::NA: return autodock_ff::NA;
      case ligand_map_types::OA: return autodock_ff::OA;
      case ligand_map_types::SA: return autodock_ff::SA;
      case ligand_map_types::Cl: return autodock_ff::Cl;
      case ligand_map_types::F: return autodock_ff::F;
      case ligand_map_types::S: return autodock_ff::S;
      case ligand_map_types::Br: return autodock_ff::Br;
      case ligand_map_types::P: return autodock_ff::P;
      case ligand_map_types::I: return autodock_ff::I;
      default: throw std::runtime_error("Missing autodock type from texture!");
    }
  }

  struct map_description {
    ligand_map_types value;
    std::string_view name;
  };
  extern const std::array<map_description, num_ligand_map_types()> MAP_DICTIONARY;

  inline const map_description& get_description(const ligand_map_types e) {
    assert(MAP_DICTIONARY[static_cast<int>(e)].value == e);
    return MAP_DICTIONARY[static_cast<int>(e)];
  }
  ligand_map_types parse_map_symbol(const std::string_view symbol);
} // namespace mudock
