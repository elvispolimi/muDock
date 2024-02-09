#ifndef MUDCOK_MOLECULE_HDR
#define MUDCOK_MOLECULE_HDR

#include <mudock/type_alias.hpp>
#include <span>

namespace mudock {

  struct molecule {
    std::span<coordinate_type> x;
    std::span<coordinate_type> y;
    std::span<coordinate_type> z;
  };

} // namespace mudock

#endif // MUDCOK_MOLECULE_HDR
