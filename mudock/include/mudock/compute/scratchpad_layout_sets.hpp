#pragma once

#include <concepts>
#include <tuple>
#include <type_traits>
#include <mudock/compute/scratchpad_data.hpp>

namespace mudock {

  template<typename... Types>
  struct type_list {};

  using all_scratch_layouts = type_list<ScratchInputCoordinates<>,
                                        ScratchCoordinates<>,
                                        ScratchBatchInfo<>,
                                        ScratchGenetic<>,
                                        ScratchScores<>,
                                        ScratchBest<>,
                                        ScratchProteinBounds<>,
                                        ScratchProteinDimensions<>,
                                        ScratchProteinMaps<>>;

  using scoring_scratch_layouts = type_list<ScratchCoordinates<>,
                                            ScratchBatchInfo<>,
                                            ScratchScores<>,
                                            ScratchProteinBounds<>,
                                            ScratchProteinDimensions<>,
                                            ScratchProteinMaps<>>;

  using genetic_scratch_layouts = all_scratch_layouts;

  template<typename Layout, typename LayoutList>
  struct layout_in_list;

  template<typename Layout, typename... Layouts>
  struct layout_in_list<Layout, type_list<Layouts...>>
      : std::bool_constant<(std::same_as<Layout, Layouts> || ...)> {};

}  // namespace mudock
