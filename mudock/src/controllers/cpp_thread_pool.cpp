#include <cstdint>
#include <mudock/controllers/cpp_thread_pool.hpp>
#include <thread>

namespace mudock {
  template<>
  thread_pool<kernel_type::CPP>::thread_pool(std::shared_ptr<squeue<static_molecule>> i_queue,
                                             std::shared_ptr<squeue<static_molecule>> o_queue,
                                             const device_conf& dev_c) {
    if (dev_c.d_t != device_type::CPU) {
      throw std::runtime_error{"CPP version supports only CPU devices"};
    }
    for (const std::size_t id: dev_c.ids) {
      this->workers.emplace_back(worker<kernel_type::CPP>, i_queue, o_queue, id);
    }
  };

} // namespace mudock
