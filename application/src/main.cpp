#include <iostream>
#include <mudock/mudock.hpp>

auto test(const mudock::molecule auto& molecule) { return molecule.atoms.coordinates.x[7]; }

int main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[]) {
  std::cout << "Hello world" << std::endl;

  mudock::static_molecule ligand;
  test(ligand);

  return EXIT_SUCCESS;
}
