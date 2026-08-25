#pragma once

#include <mudock/chem/x_score_layer.hpp>
#include <mudock/chem/x_score_xtool_types.hpp>
#include <mudock/chem/x_score_xlogp_types.hpp>
#include <mudock/chem/x_score_residue_xtool_types.hpp>
#include <mudock/chem/assign_x_score_types.hpp> 
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>

namespace mudock { 
 struct x_score_ligand: public x_score_static_layer {

		x_score_ligand(static_molecule& _molecule) : x_score_static_layer(_molecule) {

		}
 };

}