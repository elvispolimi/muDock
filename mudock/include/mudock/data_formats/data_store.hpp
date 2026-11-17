#pragma once

#include <mudock/data_formats/data_store_soa.hpp>
#include <mudock/data_formats/data_layout_sets.hpp>

namespace mudock {

  template<typename queue_type, typename soa_type = all_data_layouts>
    requires std::derived_from<queue_type, queue>
  using data_store = data_store_soa<queue_type, soa_type>;

}  // namespace mudock
