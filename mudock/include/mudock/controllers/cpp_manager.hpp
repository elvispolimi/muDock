#pragma once

#include <mudock/controllers/manager.hpp>
#include <mudock/controllers/squeue.hpp>

namespace mudock {
  template<>
  manager<language_types::CPP>::manager(std::shared_ptr<squeue<static_molecule>>,
                                        std::shared_ptr<squeue<static_molecule>>,
                                        const device_conf&);
} // namespace mudock
