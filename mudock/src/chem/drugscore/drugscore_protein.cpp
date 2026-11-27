
#include <mudock/chem/drugscore/drugscore_protein.hpp>
#include <mudock/chem/drugscore/drugscore_model.hpp>
#include <mudock/chem/grid_const.hpp>

namespace mudock {

  void drugscore_precomp_protein::prepare() {
    const auto& molecule = this->get_base_molecule();

    // get the protein's atoms coordinate
    const auto prt_x = molecule.get_x();
    const auto prt_y = molecule.get_y();
    const auto prt_z = molecule.get_z();

    // check atoms amounts
    assert(prt_x.size() == prt_y.size());
    assert(prt_x.size() == prt_z.size());
    const auto num_atoms = prt_x.size();

    // bail out if there are no protein atoms
    if (num_atoms == std::size_t{0}) {
      throw std::runtime_error("Protein without atoms");
    }

    // find the protein boundaries
    auto min = point3D{std::ranges::min(prt_x), std::ranges::min(prt_y), std::ranges::min(prt_z)};
    std::transform(min.begin(), min.end(), min.begin(), [](const auto p) {
      return std::floor((p - (grid_spacing + maxdist())) * inv_spacing) / inv_spacing;
    });
    auto max = point3D{std::ranges::max(prt_x), std::ranges::max(prt_y), std::ranges::max(prt_z)};
    std::transform(max.begin(), max.end(), max.begin(), [](const auto p) {
      return std::ceil((p + (grid_spacing + maxdist())) * inv_spacing) / inv_spacing;
    });

    // declare the maps that will describe the protein
    const auto resolution = grid_spacing; // taking resolution from grid_const.hpp
    ds_grid               = {min, max, resolution};

    const auto ds_types = get_drugscore_type();

    const auto generic_map = get_atom_map(static_cast<ds_type>(0)); // just to get the map dimensions
    const auto size_x      = generic_map.size<0>();
    const auto size_y      = generic_map.size<1>();
    const auto size_z      = generic_map.size<2>();

    // pre-compute atom points3D
    std::vector<point3D> atom_points;
    atom_points.reserve(num_atoms);
    for (std::size_t atom_index = 0; atom_index < num_atoms; ++atom_index) {
      const auto atom_x = prt_x[atom_index];
      const auto atom_y = prt_y[atom_index];
      const auto atom_z = prt_z[atom_index];
      atom_points.emplace_back(atom_x, atom_y, atom_z);
    }

    // iterate over all grid points
#pragma omp parallel for collapse(3) schedule(static)
    for (size_t grid_x = 0; grid_x < size_x; ++grid_x) {
      for (size_t grid_y = 0; grid_y < size_y; ++grid_y) {
        for (size_t grid_z = 0; grid_z < size_z; ++grid_z) {
          std::array<fp_type, num_valued_ds_types()> voxel_values{};
          voxel_values.fill(fp_type{0});

          const auto grid_point = generic_map.to_coord(static_cast<fp_type>(grid_x),
                                                       static_cast<fp_type>(grid_y),
                                                       static_cast<fp_type>(grid_z));

          // iterate over all protein atoms
          for (std::size_t atom_index = 0; atom_index < num_atoms; ++atom_index) {
            // get protein atom drugscore type
            const auto atom_ds = ds_types[atom_index];
            if (atom_ds == ds_type::ND) {
              continue; // skip non-valued types
            }

            // compute distance bucket
            const auto dist                = grid_point.distance(atom_points[atom_index]);
            const unsigned int dist_bucket = bucketize(dist);
            const bool range_valid         = (dist_bucket >= min_bucket()) && (dist_bucket <= max_bucket());
            if (!range_valid) {
              continue; // skip out-of-range distances
            }

            const unsigned int protein_idx = ds2index(atom_ds);
            const unsigned int bucket      = dist_bucket - min_bucket();

            // for each possible ligand atom type, accumulate the score for this voxel
            for (unsigned int ligand_idx = 0; ligand_idx < num_valued_ds_types(); ++ligand_idx) {
              const auto score = DRUGSCORE_MODEL[ligand_idx][protein_idx][bucket];
              voxel_values[ligand_idx] += score;
            }
          }

          // assign voxel values to the maps
#pragma GCC unroll num_valued_ds_types()
          for (unsigned int i = 0; i < num_valued_ds_types(); ++i) {
            get_atom_map(static_cast<ds_type>(i)).get(grid_x, grid_y, grid_z) = voxel_values[i];
          }
        }
      }
    }
  }

} // namespace mudock
