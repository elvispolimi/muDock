#pragma once

#include <array>
#include <mudock/chem/autodock_types.hpp>
#include <stdexcept>

namespace mudock {
  enum class device_map_textures : int { A = 0, C, H, HD, N, NA, OA, SA, Cl, F, S, Br, P, I };

  static constexpr int num_device_map_textures() { return 14; }

  inline device_map_textures map_from_autodock_type(const autodock_ff& autodock_type) {
    switch (autodock_type) {
      case autodock_ff::A: return device_map_textures::A;
      case autodock_ff::C: return device_map_textures::C;
      case autodock_ff::H: return device_map_textures::H;
      case autodock_ff::HD: return device_map_textures::HD;
      case autodock_ff::N: return device_map_textures::N;
      case autodock_ff::NA: return device_map_textures::NA;
      case autodock_ff::OA: return device_map_textures::OA;
      case autodock_ff::SA: return device_map_textures::SA;
      case autodock_ff::Cl: return device_map_textures::Cl;
      case autodock_ff::F: return device_map_textures::F;
      case autodock_ff::S: return device_map_textures::S;
      case autodock_ff::Br: return device_map_textures::Br;
      case autodock_ff::P: return device_map_textures::P;
      case autodock_ff::I: return device_map_textures::I;
      default:
        throw std::runtime_error("Missing map texture " + std::string{get_description(autodock_type).name});
    }
  }

  inline autodock_ff autodock_type_from_map(const device_map_textures& cuda_map) {
    switch (cuda_map) {
      case device_map_textures::A: return autodock_ff::A;
      case device_map_textures::C: return autodock_ff::C;
      case device_map_textures::H: return autodock_ff::H;
      case device_map_textures::HD: return autodock_ff::HD;
      case device_map_textures::N: return autodock_ff::N;
      case device_map_textures::NA: return autodock_ff::NA;
      case device_map_textures::OA: return autodock_ff::OA;
      case device_map_textures::SA: return autodock_ff::SA;
      case device_map_textures::Cl: return autodock_ff::Cl;
      case device_map_textures::F: return autodock_ff::F;
      case device_map_textures::S: return autodock_ff::S;
      case device_map_textures::Br: return autodock_ff::Br;
      case device_map_textures::P: return autodock_ff::P;
      case device_map_textures::I: return autodock_ff::I;
      default: throw std::runtime_error("Missing autodock type from texture!");
    }
  }
} // namespace mudock