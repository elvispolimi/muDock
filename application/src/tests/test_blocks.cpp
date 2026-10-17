#include <array>
#include <cstddef>
#include <iostream>
#include <mudock/soa_template/soa_blocks.hpp>
#include <mudock/soa_template/soa_layout.hpp>
#include <vector>

GENERATE_SOA_LAYOUT(SoAPositionAtoms,
                    SOA_COLUMN(float, x),
                    SOA_COLUMN(float, y),
                    SOA_COLUMN(float, z))

GENERATE_SOA_LAYOUT(SoAIndexes,
                    SOA_COLUMN(int, indexes),
                    SOA_COLUMN(int, adjacents))

GENERATE_SOA_BLOCKS(SoABlocks,
                    SOA_BLOCK(positionAtoms, SoAPositionAtoms),
                    SOA_BLOCK(indexData, SoAIndexes))

int main() {
  using Blocks = SoABlocks<>;
  constexpr Blocks::size_type positionElements = 3;
  constexpr Blocks::size_type indexElements = 5;
  const std::array sizes{positionElements, indexElements};

  std::vector<std::byte> buffer(Blocks::computeDataSize(sizes));
  Blocks blocks(buffer.data(), sizes);
  Blocks::View view{blocks};

  for (Blocks::size_type index = 0; index < positionElements; ++index) {
    view.positionAtoms().x()[index] = static_cast<float>(index);
    view.positionAtoms().y()[index] = static_cast<float>(index + 10);
    view.positionAtoms().z()[index] = static_cast<float>(index + 20);
  }
  for (Blocks::size_type index = 0; index < indexElements; ++index) {
    view.indexData().indexes()[index] = index;
    view.indexData().adjacents()[index] = index + 1;
  }

  const auto viewSizes = view.metadata().size();
  std::cout << "positionAtoms (" << viewSizes[0] << " elements):\n";
  for (Blocks::size_type index = 0; index < viewSizes[0]; ++index)
    std::cout << view.positionAtoms().x(index) << ' ' << view.positionAtoms().y(index) << ' '
              << view.positionAtoms().z(index) << '\n';

  std::cout << "indexData (" << viewSizes[1] << " elements):\n";
  for (Blocks::size_type index = 0; index < viewSizes[1]; ++index)
    std::cout << view.indexData().indexes(index) << ' ' << view.indexData().adjacents(index) << '\n';
}