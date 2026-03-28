#pragma once

#include <mudock/compute/precomputed_adt_score.hpp>
#include <mudock/compute/precomputed_adt_score_kernel.hpp>
#include <mudock/cpp_implementation/queue_cpp.hpp>

namespace mudock {
  
  
  template<>
  inline int get_precomputed_adt_score_batch<queue_cpp>(const int, std::shared_ptr<queue_cpp>) {
    // Lasciamo 10 per essere coerenti con l'originale.
    // Se la RAM esplode o la Cache fa i capricci, abbassiamo questo valore
    return 10;
  }


  template<>
  void precomputed_adt_score_kernel<queue_cpp>::operator()();
  
} // namespace mudock