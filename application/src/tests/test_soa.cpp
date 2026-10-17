#include <cstddef>
#include <iostream>
#include <mudock/soa_template/soa_layout.hpp>
#include <vector>

GENERATE_SOA_LAYOUT(SoATest, SOA_COLUMN(int, x), SOA_COLUMN(double, z));

int main() {
  constexpr SoATest<>::size_type elements = 4;
  std::vector<std::byte> buffer(SoATest<>::computeDataSize(elements));

  SoATest<> soa(buffer.data(), elements);
  SoATest<>::View view{soa};
  for (SoATest<>::size_type index = 0; index < elements; ++index)
    {
      view.x()[index] = index;
      view.z()[index] = index * 2.0;
    }

  for (SoATest<>::size_type index = 0; index < view.metadata().size(); ++index)
    {
      std::cout << view.x(index) << '\n';
      std::cout << view.z(index) << '\n';
    }
}    