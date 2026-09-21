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
    

    // TODO L should i check if multiple ligands? i think so
    void dump_pose(std::size_t i) {
      assert(ligand_template.has_value() && "Ligand template was not stored before dump_pose");

      auto& x_scratch_b = (*this->scratch).template get<buffer_data_type::X_SCRATCH>();
      auto& y_scratch_b = (*this->scratch).template get<buffer_data_type::Y_SCRATCH>();
      auto& z_scratch_b = (*this->scratch).template get<buffer_data_type::Z_SCRATCH>();
      x_scratch_b.copy_device2host();
      y_scratch_b.copy_device2host();
      z_scratch_b.copy_device2host();
      (*this->scratch).get_queue()->synchronize();

      static_molecule pose = *ligand_template;
      const int num_atoms = pose.num_atoms();
      std::memcpy(pose.x(), x_scratch_b.host_pointer(), num_atoms * sizeof(fp_type));
      std::memcpy(pose.y(), y_scratch_b.host_pointer(), num_atoms * sizeof(fp_type));
      std::memcpy(pose.z(), z_scratch_b.host_pointer(), num_atoms * sizeof(fp_type));

      // TODO L if there is no (dump) directory, it doesn't save the pose. fix it.
      const std::string filename = "dump/pose_" + std::to_string(i) + ".mol2";
      std::ofstream ofs(filename, std::ios::out);
      writer<supported_format::MOL2>(pose, ofs);
    }

    void dump_best_pose_genetic(int ligand_index, int individual_index, int population_number, int batch_atoms) {
      auto& x_scratch_b = (*this->scratch).template get<buffer_data_type::X_SCRATCH>();
      auto& y_scratch_b = (*this->scratch).template get<buffer_data_type::Y_SCRATCH>();
      auto& z_scratch_b = (*this->scratch).template get<buffer_data_type::Z_SCRATCH>();
      
      x_scratch_b.copy_device2host();
      y_scratch_b.copy_device2host();
      z_scratch_b.copy_device2host();
      (*this->scratch).get_queue()->synchronize();

      fp_type *__restrict__ x_scratch_p = x_scratch_b.host_pointer();
      fp_type *__restrict__ y_scratch_p = y_scratch_b.host_pointer();
      fp_type *__restrict__ z_scratch_p = z_scratch_b.host_pointer();

      int atom_stride  = ligand_index * batch_atoms;

      fp_type *__restrict__ scratch_x = x_scratch_p + atom_stride * population_number;
      fp_type *__restrict__ scratch_y = y_scratch_p + atom_stride * population_number;
      fp_type *__restrict__ scratch_z = z_scratch_p + atom_stride * population_number;

      fp_type *__restrict__ scratch_x_l = scratch_x + individual_index * batch_atoms;
      fp_type *__restrict__ scratch_y_l = scratch_y + individual_index * batch_atoms;
      fp_type *__restrict__ scratch_z_l = scratch_z + individual_index * batch_atoms;
      
      static_molecule pose = *ligand_template;
      const int num_atoms = pose.num_atoms();
      std::memcpy(pose.x(), scratch_x_l, num_atoms * sizeof(fp_type));
      std::memcpy(pose.y(), scratch_y_l, num_atoms * sizeof(fp_type));
      std::memcpy(pose.z(), scratch_z_l, num_atoms * sizeof(fp_type));

      // TODO L if there is no (dump) directory, it doesn't save the pose. fix it.
      const std::string filename = "dump_best/pose_" + std::to_string(++counter) + ".mol2";
      std::ofstream ofs(filename, std::ios::out);
      writer<supported_format::MOL2>(pose, ofs);
    }

    bool is_standalone_local_search(const knobs& conf) {
      return conf.population_number == 1 && conf.num_generations == 1 && conf.lsrate == 100;
    }
  };
} // namespace mudock