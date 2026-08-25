#pragma once

#include <cstddef>
#include <cstring>
#include <mudock/batch.hpp>
#include <mudock/chem/x_score_hb.hpp>
#include <mudock/chem/x_score_hb_term.hpp>
#include <mudock/chem/x_score_ligand.hpp>
#include <mudock/chem/x_score_protein.hpp>
#include <mudock/chem/x_score_rt.hpp>
#include <mudock/chem/x_score_validity.hpp>
#include <mudock/compute/x_score_kernel.hpp>
#include <mudock/compute/x_score_terms.hpp>
#ifndef __CUDACC__
  #include <mudock/compute/buffer_utils.hpp>
  #include <mudock/compute/scoring.hpp>
  #include <mudock/compute/scratchpad.hpp>
#endif
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  template<typename queue_type>
  int get_x_score_batch(const int, std::shared_ptr<queue_type>, const size_t);

#ifndef __CUDACC__
  // XScore scoring stage: vdW, HB, HP and RT implemented; HM and HS not implemented
  // NB: XScore's protein define pocket filter seems unnecessary for the VDW term (8 Å cutoff is always within the 10 Å pocket) and for the HB term (5 Å cutoff is always within the 10 Å pocket)
  template<typename queue_type>
  struct x_score: public scoring<queue_type> {
    x_score(std::shared_ptr<scratchpad<queue_type>> _scratch,
            std::shared_ptr<scratchpad<queue_type>> _device_scratch,
            dynamic_molecule &protein)
        : scoring<queue_type>(_scratch),
          lig_x(_scratch->get_queue()),
          lig_y(_scratch->get_queue()),
          lig_z(_scratch->get_queue()),
          lig_vdw(_scratch->get_queue()),
          lig_scorable(_scratch->get_queue()),
          lig_hb(_scratch->get_queue()),
          lig_rt(_scratch->get_queue()),
          lig_hbt(_scratch->get_queue()),
          prot_x(_scratch->get_queue()),
          prot_y(_scratch->get_queue()),
          prot_z(_scratch->get_queue()),
          prot_vdw(_scratch->get_queue()),
          prot_scorable(_scratch->get_queue()),
          prot_hb(_scratch->get_queue()),
          terms(_scratch->get_queue()),
          device_scratch(_device_scratch) {
      x_score_protein xs_prot{protein};
      const auto &pmol = xs_prot.get_base_molecule();
      num_prot_atoms   = pmol.num_atoms();

      prot_x.alloc(num_prot_atoms);
      prot_y.alloc(num_prot_atoms);
      prot_z.alloc(num_prot_atoms);
      prot_vdw.alloc(num_prot_atoms);
      prot_scorable.alloc(num_prot_atoms);
      prot_hb.alloc(num_prot_atoms);

      for (int j = 0; j < num_prot_atoms; ++j) {
        prot_x()[j]   = pmol.x(j);
        prot_y()[j]   = pmol.y(j);
        prot_z()[j]   = pmol.z(j);
        prot_vdw()[j] = xs_prot.vdw_radius(j);
        prot_scorable()[j] =
            is_protein_scorable(xs_prot.valid(j), xs_prot.x_score_xtool_type(j)) ? 1 : 0;
        prot_hb()[j] = static_cast<int>(xs_prot.hb(j));
      }

      prot_x.copy_host2device();
      prot_y.copy_host2device();
      prot_z.copy_host2device();
      prot_vdw.copy_host2device();
      prot_scorable.copy_host2device();
      prot_hb.copy_host2device();

      // HB term: build the protein donor/acceptor atom list once, it's reused with every ligand
      prot_hb_atoms = build_protein_hb_atoms(xs_prot);
    }

    void prepare(batch<static_molecule> &batch) {
      batch_ligands                = batch.num_ligands;
      batch_atoms                  = batch.batch_max_atoms;
      const int tot_atoms_in_batch = batch_ligands * batch_atoms;
      auto q                       = (*this->scratch).get_queue();

      load_num_atoms(batch, this->scratch);

      const int scores_per_ligand =
          std::max(1, static_cast<int>((*this->scratch).configuration.population_number));

      auto &score_b = (*this->scratch).template get<buffer_data_type::SCORES>();
      if (!score_b.is_valid() ||
          score_b.num_elements() != static_cast<size_t>(batch_ligands * scores_per_ligand)) {
        score_b.alloc(batch_ligands * scores_per_ligand);
        score_b.set_valid();
      }

      lig_x.alloc(tot_atoms_in_batch);
      lig_y.alloc(tot_atoms_in_batch);
      lig_z.alloc(tot_atoms_in_batch);
      lig_vdw.alloc(tot_atoms_in_batch);
      lig_scorable.alloc(tot_atoms_in_batch);
      lig_hb.alloc(tot_atoms_in_batch);
      lig_rt.alloc(batch_ligands);
      lig_hbt.alloc(batch_ligands);
      terms.alloc(batch_ligands * static_cast<int>(x_term_count));

      for (int ligand_index = 0; ligand_index < batch_ligands; ++ligand_index) {
        auto &ligand           = *batch.molecules[ligand_index];
        const int stride_atoms = ligand_index * batch_atoms;
        const int num_atoms    = ligand.num_atoms();

        x_score_ligand xs_lig{ligand}; // assigns xtool types + vdw radius + validity

        const auto x = ligand.x(), y = ligand.y(), z = ligand.z();
        std::memcpy((void *) (lig_x() + stride_atoms), x, num_atoms * sizeof(fp_type));
        std::memcpy((void *) (lig_y() + stride_atoms), y, num_atoms * sizeof(fp_type));
        std::memcpy((void *) (lig_z() + stride_atoms), z, num_atoms * sizeof(fp_type));
        std::memcpy((void *) (lig_vdw() + stride_atoms), xs_lig.vdw_radius(), num_atoms * sizeof(fp_type));

        for (int i = 0; i < num_atoms; ++i) {
          lig_scorable()[stride_atoms + i] =
              is_ligand_scorable(xs_lig.valid(i), xs_lig.x_score_xtool_type(i)) ? 1 : 0;
          lig_hb()[stride_atoms + i] = static_cast<int>(xs_lig.hb(i));
        }

        // RT Term for the ligand
        lig_rt()[ligand_index] = compute_x_score_rt(xs_lig);

        // HB Term for the ligands
        // NB: Calculated once per ligand, but pose dependent
        x_score_hb_atoms lig_hb_atoms = build_ligand_hb_atoms(xs_lig);
        lig_hbt()[ligand_index] = compute_x_score_hb(lig_hb_atoms, prot_hb_atoms);
      }

      lig_x.copy_host2device();
      lig_y.copy_host2device();
      lig_z.copy_host2device();
      lig_vdw.copy_host2device();
      lig_scorable.copy_host2device();
      lig_hb.copy_host2device();
      lig_rt.copy_host2device();
      lig_hbt.copy_host2device();

      const int *num_atoms_b = (*this->scratch).template get<buffer_data_type::NUM_ATOMS>().dev_pointer();

      kernel = std::make_unique<x_score_kernel<queue_type>>(batch_ligands,
                                                            batch_atoms,
                                                            num_atoms_b,
                                                            lig_x.dev_pointer(),
                                                            lig_y.dev_pointer(),
                                                            lig_z.dev_pointer(),
                                                            lig_vdw.dev_pointer(),
                                                            lig_scorable.dev_pointer(),
                                                            lig_hb.dev_pointer(),
                                                            lig_rt.dev_pointer(),
                                                            lig_hbt.dev_pointer(),
                                                            num_prot_atoms,
                                                            prot_x.dev_pointer(),
                                                            prot_y.dev_pointer(),
                                                            prot_z.dev_pointer(),
                                                            prot_vdw.dev_pointer(),
                                                            prot_scorable.dev_pointer(),
                                                            prot_hb.dev_pointer(),
                                                            terms.dev_pointer(),
                                                            q);
    }

    void operator()() {
      assert(kernel && "Kernel method not yet prepared");
      (*kernel)();
    }

    static int get_ligand_mem(const int max_atoms, const knobs conf) {
      std::size_t mem{0};
      const int scores_per_ligand = std::max(1, static_cast<int>(conf.population_number));
      mem += sizeof(fp_type) * scores_per_ligand; // scores
      mem += sizeof(int);                         // num atoms
      mem += sizeof(fp_type) * max_atoms;         // lig x
      mem += sizeof(fp_type) * max_atoms;         // lig y
      mem += sizeof(fp_type) * max_atoms;         // lig z
      mem += sizeof(fp_type) * max_atoms;         // lig vdw radius
      mem += sizeof(int) * max_atoms;             // lig scorable flag
      mem += sizeof(int) * max_atoms;             // lig hb class
      mem += sizeof(fp_type);                     // lig rotor term
      mem += sizeof(fp_type);                     // lig HB term
      mem += sizeof(fp_type) * x_term_count;      // per-ligand XScore terms
      return static_cast<int>(mem);
    }

  private:
    int batch_ligands{0};
    int batch_atoms{0};
    int num_prot_atoms{0};

    // ligand atom data
    buffer_vector<fp_type, queue_type> lig_x;
    buffer_vector<fp_type, queue_type> lig_y;
    buffer_vector<fp_type, queue_type> lig_z;
    buffer_vector<fp_type, queue_type> lig_vdw;
    buffer_vector<int, queue_type> lig_scorable;
    buffer_vector<int, queue_type> lig_hb;
    buffer_vector<fp_type, queue_type> lig_rt;
    buffer_vector<fp_type, queue_type> lig_hbt; 

    // protein atom data
    buffer_vector<fp_type, queue_type> prot_x;
    buffer_vector<fp_type, queue_type> prot_y;
    buffer_vector<fp_type, queue_type> prot_z;
    buffer_vector<fp_type, queue_type> prot_vdw;
    buffer_vector<int, queue_type> prot_scorable;
    buffer_vector<int, queue_type> prot_hb;

    // protein HB donor/acceptor
    x_score_hb_atoms prot_hb_atoms;

    // per-ligand XScore terms (x_term_count values per ligand)
    buffer_vector<fp_type, queue_type> terms;

    std::shared_ptr<scratchpad<queue_type>> device_scratch;
    std::unique_ptr<x_score_kernel<queue_type>> kernel;

    static inline bool is_ligand_scorable(const x_score_validity v, const xtool_ff t) {
      if (v == x_score_validity::invalid)
        return false;
      return t != xtool_ff::H && t != xtool_ff::Hhb && t != xtool_ff::Hg;
    }

    static inline bool is_protein_scorable(const x_score_validity v, const xtool_ff t) {
      if (v == x_score_validity::invalid)
        return false;
      return t != xtool_ff::H && t != xtool_ff::Hhb && t != xtool_ff::Hg && t != xtool_ff::Ow;
    }

    void teardown_impl(batch<static_molecule> &batch) override {
      auto &scores_b               = (*this->scratch).template get<buffer_data_type::SCORES>();
      const auto scores_per_ligand = scores_b.num_elements() / batch.num_ligands;
      terms.copy_device2host();
      (*this->scratch).get_queue()->synchronize();
      for (int i = 0; i < batch.num_ligands; ++i) {
        auto &ligand      = *batch.molecules[i];
        const fp_type *t  = terms() + i * static_cast<int>(x_term_count);
        const fp_type vdw = t[x_term_vdw];
        const fp_type hb  = t[x_term_hb];
        const fp_type hp  = t[x_term_hp];
        const fp_type rt  = t[x_term_rt];
        const fp_type pkd = t[x_term_pkd];

        for (std::size_t s = 0; s < scores_per_ligand; ++s)
          scores_b()[i * scores_per_ligand + s] = pkd;

        ligand.properties.assign(property_type::SCORE,
                                 "VDW=" + std::to_string(vdw) + " HB=" + std::to_string(hb) +
                                     " HP=" + std::to_string(hp) + " RT=" + std::to_string(rt) +
                                     " PKD=" + std::to_string(pkd));
      }
    }
  };
#endif
} // namespace mudock