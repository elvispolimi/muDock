#pragma once

#include <mudock/batch.hpp>
#include <mudock/molecule.hpp>
#include <mudock/compute/stage.hpp>
#include <mudock/compute/scratchpad.hpp>
#include <mudock/compute/scoring.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/format/writer.hpp>


namespace mudock {
  template<typename queue_t, template<typename> typename scoring_t>
    requires std::derived_from<queue_t, queue> && std::derived_from<scoring_t<queue_t>, scoring<queue_t>>
  struct local_search: public stage<queue_t> {
    local_search(std::shared_ptr<scratchpad<queue_t>> _scratch, 
                 std::shared_ptr<scoring_t<queue_t>> _score)
                 : stage<queue_t>(_scratch),
                   score_stage(_score) {
                   iterations = (*this->scratch).configuration.lsit;
                   };
    virtual void prepare(batch<static_molecule>&) = 0;
    virtual void operator()()                     = 0;

    virtual ~local_search() = default;
    
  protected:
    std::shared_ptr<scoring_t<queue_t>> score_stage;
    size_t iterations;
    bool standalone_local_search{false};
    std::optional<static_molecule> ligand_template;
    int counter = 0;
    

    // WARNING: ONLY FOR STANDALONE (ASSUMES 1 INDIVIDUAL)
    void dump_pose(const std::string complex_id,
                   const int ligand_index,
                   const fp_type* x_scratch_b,
                   const fp_type* y_scratch_b,
                   const fp_type* z_scratch_b,
                   const int batch_atoms) {
      assert(ligand_template.has_value() && "Ligand template was not stored before dump_pose");

      static_molecule pose = *ligand_template;
      const int num_atoms = pose.num_atoms();
      const int atom_stride = ligand_index * batch_atoms;

      const fp_type *__restrict__ scratch_x = x_scratch_b + atom_stride;
      const fp_type *__restrict__ scratch_y = y_scratch_b + atom_stride;
      const fp_type *__restrict__ scratch_z = z_scratch_b + atom_stride;

      std::memcpy(pose.x(), scratch_x, num_atoms * sizeof(fp_type));
      std::memcpy(pose.y(), scratch_y, num_atoms * sizeof(fp_type));
      std::memcpy(pose.z(), scratch_z, num_atoms * sizeof(fp_type));

      // TODO L if there is no (dump) directory, it doesn't save the pose. fix it.
      const std::string filename = "dump/pose_" + complex_id + "_" + std::to_string(ligand_index) + ".mol2";
      std::ofstream ofs(filename, std::ios::out);
      writer<supported_format::MOL2>(pose, ofs);
    }


    bool is_standalone_local_search(const knobs& conf) {
      return conf.population_number == 1 && conf.num_generations == 1 && conf.lsrate == 100;
    }
  };
} // namespace mudock