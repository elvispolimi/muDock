#include <alpaka/acc/Traits.hpp>
#include <alpaka/example/ExampleDefaultAcc.hpp>
#include <mudock/implementation_types.hpp>

#include <cstddef>
#include <cstdlib>
#include <iostream>

int main() {
#ifndef MUDOCK_USE_ALPAKA
  return EXIT_FAILURE;
#else
  const auto impl = mudock::get_impl_type("ALPAKA");
  if (impl != mudock::implementation_type::ALPAKA) {
    return EXIT_FAILURE;
  }

  using Dim = alpaka::DimInt<1u>;
  using Idx = std::size_t;
  using Acc = alpaka::ExampleDefaultAcc<Dim, Idx>;

  std::cout << "Alpaka default accelerator: " << alpaka::getAccName<Acc>() << std::endl;
  return EXIT_SUCCESS;
#endif
}
