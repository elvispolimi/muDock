#pragma once

#include <mudock/controllers/squeue.hpp>
#include <mudock/controllers/thread_pool.hpp>

namespace mudock {
  template<>
  thread_pool<kernel_type::CPP>::thread_pool(std::shared_ptr<squeue<static_molecule>> i_queue,
                                                std::shared_ptr<squeue<static_molecule>> o_queue,
                                                const device_conf& dev_c);
} // namespace mudock
