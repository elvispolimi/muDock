#pragma once

#include <mudock/chem/x_score_layer.hpp>
#include <mudock/chem/x_score_xtool_types.hpp>
#include <mudock/chem/x_score_xlogp_types.hpp>
#include <mudock/chem/x_score_residue_xtool_types.hpp>
#include <mudock/chem/residue.hpp>
#include <mudock/chem/assign_x_score_types.hpp> 
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>

namespace mudock { 
 struct x_score_protein: public x_score_dynamic_layer {

	x_score_protein(dynamic_molecule& _molecule) : x_score_dynamic_layer(_molecule) {

	}
 };

}