#pragma once

#include <functional>
#include <mudock/chem/x_score_layer.hpp>

namespace mudock {

  /**
   * @brief This class assigns x-score's x_tool and x_logp types.
   * 
   * @tparam layer_type The specific type of the scoring layer. Must satisfy the is_x_score_layer concept.
   * @param layer The layer instance being configured.
   */

  // generic template for x_score layers
  template<typename layer_type>
    requires is_x_score_layer<layer_type>
  void assign_x_score_types(layer_type&);

  /**
   * @brief Ligand specialization (static) for assigning x-score types.
   *
   * @param layer The static ligand layer instance being configured.
   */
  template<>
  void assign_x_score_types(x_score_static_layer&);

  /**
   * @brief Protein specialization (dynamic) for assigning x-score types.
   *
   * @param layer The dynamic protein layer instance being configured.
   */
  template<>
  void assign_x_score_types(x_score_dynamic_layer&);

} // namespace mudock
