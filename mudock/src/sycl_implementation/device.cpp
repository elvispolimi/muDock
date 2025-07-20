#include <mudock/chem/autodock_protein.hpp>
#include <mudock/sycl_implementation/device.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {
  device::device(const sycl::device& dev, const autodock_protein& adt_protein)
      : center(adt_protein.get_center()),
        minimum(adt_protein.get_min()),
        maximum(adt_protein.get_max()),
        adt_protein(adt_protein),
        dev(dev),
        queue(get_queue()),
        atom_tex(queue, adt_protein) {
    queue.wait();
  }

  sycl::queue device::get_queue() const { return sycl::queue{dev, sycl::property::queue::in_order{}}; }

  int device::get_sub_group_size() const {
    return static_cast<int>(dev.get_info<sycl::info::device::sub_group_sizes>().at(0));
  }

} // namespace mudock
