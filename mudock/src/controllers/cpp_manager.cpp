#include <mudock/controllers/cpp_manager.hpp>
#include <thread>

namespace mudock {
  template<>
  manager<language_types::CPP>::manager(std::shared_ptr<squeue<static_molecule>> i_queue,
                                        std::shared_ptr<squeue<static_molecule>> o_queue,
                                        const device_conf& dev_c) {
    if (dev_c.d_t != device_types::CPU) {
      throw std::runtime_error{"CPP version supports only CPU devices"};
    }
    for (index_type id: dev_c.ids) {
      this->workers.emplace_back(worker<language_types::CPP>, i_queue, o_queue, id);
    }
  };

} // namespace mudock
