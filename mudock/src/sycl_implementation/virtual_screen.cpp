#include <mudock/sycl_implementation/virtual_screen.hpp>

namespace mudock {
  static constexpr std::size_t max_non_bonds{1024 * 3};

  virtual_screen_sycl::virtual_screen_sycl(const knobs k,
                                           std::shared_ptr<const grid_atom_mapper> &grid_atom_maps,
                                           std::shared_ptr<const grid_map> &electro_map,
                                           std::shared_ptr<const grid_map> &desolv_map)
      : configuration(k), index_maps(electro_map.get()->index), center_maps(electro_map.get()->center) {}

  void virtual_screen_sycl::operator()(batch &incoming_batch) {}
} // namespace mudock
