#pragma once

#include <functional>
#include <mudock/chem/drugscore/drugscore_layer.hpp>
#include <mudock/grid.hpp>

#ifndef DRUGSCORE_PRECOMP
  #define DRUGSCORE_PRECOMP 1
#endif

namespace mudock {

  struct drugscore_protein: public drugscore_dynamic_layer {
    drugscore_protein(dynamic_molecule& _molecule, std::function<void(drugscore_dynamic_layer&)> f = {})
        : drugscore_dynamic_layer(_molecule, f) {}
  };

  struct drugscore_grid {
    drugscore_grid(const point3D min, const point3D max, const fp_type resolution)
        : index(min.difference(max)
                    .apply(std::abs<fp_type>)
                    .divide({resolution})
                    .apply(static_cast<fp_type (*)(fp_type)>(std::ceil))
                    .add({fp_type{1}})),
          data(index.size_x(), index.size_y(), index.size_z(), num_valued_ds_types()),
          _resolution(resolution),
          _min(min),
          _max(max),
          _center{max.difference(min).divide({fp_type{2}}).add(min)} {};
    drugscore_grid() {};

    inline auto get_atom_map(const ds_type type) {
      assert(static_cast<unsigned int>(type) < num_valued_ds_types());
      return space_grid_view<fp_type>{
          _min,
          _max,
          _center,
          _resolution,
          data.get_slice(md_index<4>{index.size_x(), index.size_y(), index.size_z(), static_cast<int>(type)},
                         index)};
    };

    inline auto get_atom_map(const ds_type type) const {
      return space_grid_view<const fp_type>{
          _min,
          _max,
          _center,
          _resolution,
          data.get_slice(md_index<4>{index.size_x(), index.size_y(), index.size_z(), static_cast<int>(type)},
                         index)};
    };

    [[nodiscard]] inline auto get_map_flat_size() const { return index.flat_size(); };
    [[nodiscard]] inline const auto* get_maps_pointer() const { return data.data(); };
    [[nodiscard]] inline auto* get_min_p() const { return _min.data(); };
    [[nodiscard]] inline auto* get_max_p() const { return _max.data(); };
    [[nodiscard]] inline auto* get_center_p() const { return _center.data(); };
    [[nodiscard]] inline auto get_center() const { return _center; };
    [[nodiscard]] inline auto get_min() const { return _min; };
    [[nodiscard]] inline auto get_max() const { return _max; };

    [[nodiscard]] inline auto get_size_x() const { return index.size_x(); };
    [[nodiscard]] inline auto get_size_xy() const { return index.size_xy(); };
    [[nodiscard]] inline auto get_size_xyz() const { return index.flat_size(); };

    md_index<3> index;
    md_container<std::vector<fp_type>, 4> data;
    fp_type _resolution = fp_type{0.5};
    point<fp_type, 3> _min, _max, _center;
  };

  struct drugscore_precomp_protein: public drugscore_dynamic_layer {
    drugscore_precomp_protein(dynamic_molecule& _molecule,
                              std::function<void(drugscore_dynamic_layer&)> f = {})
        : drugscore_dynamic_layer(_molecule, f) {
      drugscore_precomp_protein::prepare();
    }

    inline auto get_atom_map(const ds_type type) { return ds_grid.get_atom_map(type); };
    inline auto get_atom_map(const ds_type type) const { return ds_grid.get_atom_map(type); };
    [[nodiscard]] inline auto get_map_flat_size() const { return ds_grid.get_map_flat_size(); };
    [[nodiscard]] inline const auto* get_maps_pointer() const { return ds_grid.get_maps_pointer(); };

    [[nodiscard]] inline auto* get_min_p() const { return ds_grid.get_min_p(); };
    [[nodiscard]] inline auto* get_max_p() const { return ds_grid.get_max_p(); };
    [[nodiscard]] inline auto* get_center_p() const { return ds_grid.get_center_p(); };
    [[nodiscard]] inline auto get_center() const { return ds_grid.get_center(); };
    [[nodiscard]] inline auto get_min() const { return ds_grid.get_min(); };
    [[nodiscard]] inline auto get_max() const { return ds_grid.get_max(); };

    [[nodiscard]] inline auto get_size_x() const { return ds_grid.get_size_x(); };
    [[nodiscard]] inline auto get_size_xy() const { return ds_grid.get_size_xy(); };
    [[nodiscard]] inline auto get_size_xyz() const { return ds_grid.get_size_xyz(); };

  private:
    drugscore_grid ds_grid;
    void prepare();
  };

} // namespace mudock
