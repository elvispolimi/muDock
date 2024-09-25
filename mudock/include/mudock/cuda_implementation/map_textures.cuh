#pragma once

#include <array>
#include <mudock/chem/autodock_types.hpp>
#include <stdexcept>

namespace mudock {
  enum class cuda_map_textures : int { A = 0, C, H, HD, N, NA, OA, SA, Cl, F, S, Br, P, I };

  static constexpr int num_cuda_map_textures() { return 14; }

  inline cuda_map_textures map_from_autodock_type(const autodock_ff& autodock_type) {
    switch (autodock_type) {
      case autodock_ff::A: return cuda_map_textures::A;
      case autodock_ff::C: return cuda_map_textures::C;
      case autodock_ff::H: return cuda_map_textures::H;
      case autodock_ff::HD: return cuda_map_textures::HD;
      case autodock_ff::N: return cuda_map_textures::N;
      case autodock_ff::NA: return cuda_map_textures::NA;
      case autodock_ff::OA: return cuda_map_textures::OA;
      case autodock_ff::SA: return cuda_map_textures::SA;
      case autodock_ff::Cl: return cuda_map_textures::Cl;
      case autodock_ff::F: return cuda_map_textures::F;
      case autodock_ff::S: return cuda_map_textures::S;
      case autodock_ff::Br: return cuda_map_textures::Br;
      case autodock_ff::P: return cuda_map_textures::P;
      case autodock_ff::I: return cuda_map_textures::I;
      default:
        throw std::runtime_error("Missing map texture " + std::string{get_description(autodock_type).name});
    }
  }

  inline autodock_ff autodock_type_from_map(const cuda_map_textures& cuda_map) {
    switch (cuda_map) {
      case cuda_map_textures::A: return autodock_ff::A;
      case cuda_map_textures::C: return autodock_ff::C;
      case cuda_map_textures::H: return autodock_ff::H;
      case cuda_map_textures::HD: return autodock_ff::HD;
      case cuda_map_textures::N: return autodock_ff::N;
      case cuda_map_textures::NA: return autodock_ff::NA;
      case cuda_map_textures::OA: return autodock_ff::OA;
      case cuda_map_textures::SA: return autodock_ff::SA;
      case cuda_map_textures::Cl: return autodock_ff::Cl;
      case cuda_map_textures::F: return autodock_ff::F;
      case cuda_map_textures::S: return autodock_ff::S;
      case cuda_map_textures::Br: return autodock_ff::Br;
      case cuda_map_textures::P: return autodock_ff::P;
      case cuda_map_textures::I: return autodock_ff::I;
      default: throw std::runtime_error("Missing autodock type from texture!");
    }
  }
} // namespace mudock