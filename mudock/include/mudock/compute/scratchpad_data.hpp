#pragma once

#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/soa_template/soa_blocks.hpp>
#include <mudock/soa_template/soa_layout.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  GENERATE_SOA_LAYOUT(ScratchInputCoordinates,
                      SOA_COLUMN(fp_type, x),
                      SOA_COLUMN(fp_type, y),
                      SOA_COLUMN(fp_type, z))

  GENERATE_SOA_LAYOUT(ScratchCoordinates,
                      SOA_COLUMN(fp_type, x),
                      SOA_COLUMN(fp_type, y),
                      SOA_COLUMN(fp_type, z))

  GENERATE_SOA_LAYOUT(ScratchBatchInfo,
                      SOA_COLUMN(int, num_atoms),
                      SOA_COLUMN(int, num_rotamers))

  GENERATE_SOA_LAYOUT(ScratchGenetic,
                      SOA_COLUMN(chromosome, chromosomes),
                      SOA_COLUMN(chromosome, next_chromosomes))

  GENERATE_SOA_LAYOUT(ScratchScores,
                      SOA_COLUMN(fp_type, scores))

  GENERATE_SOA_LAYOUT(ScratchBest,
                      SOA_COLUMN(chromosome, chromosomes),
                      SOA_COLUMN(fp_type, scores))

  GENERATE_SOA_LAYOUT(ScratchProteinBounds,
                      SOA_COLUMN(fp_type, minimum),
                      SOA_COLUMN(fp_type, maximum),
                      SOA_COLUMN(fp_type, center))

  GENERATE_SOA_LAYOUT(ScratchProteinDimensions,
                      SOA_COLUMN(int, size_x),
                      SOA_COLUMN(int, size_xy),
                      SOA_COLUMN(int, size_xyz))

  GENERATE_SOA_LAYOUT(ScratchProteinMaps,
                      SOA_COLUMN(fp_type, maps))

  GENERATE_SOA_LAYOUT(ScratchGeometryFragments,
                      SOA_COLUMN(int, values))

  GENERATE_SOA_BLOCKS(GeometryData,
                      SOA_BLOCK(fragments, ScratchGeometryFragments),
                      SOA_BLOCK(fragments_start, ScratchGeometryFragments),
                      SOA_BLOCK(fragments_start_atoms, ScratchGeometryFragments),
                      SOA_BLOCK(fragments_stop_atoms, ScratchGeometryFragments),
                      SOA_BLOCK(fragments_indices_start, ScratchGeometryFragments))

  GENERATE_SOA_BLOCKS(ScratchpadData,
                      SOA_BLOCK(coordinates, ScratchCoordinates),
                      SOA_BLOCK(scratch_coordinates, ScratchCoordinates),
                      SOA_BLOCK(batch_info, ScratchBatchInfo),
                      SOA_BLOCK(genetic, ScratchGenetic),
                      SOA_BLOCK(scores, ScratchScores),
                      SOA_BLOCK(protein_bounds, ScratchProteinBounds),
                      SOA_BLOCK(protein_dimensions, ScratchProteinDimensions),
                      SOA_BLOCK(protein_maps, ScratchProteinMaps))

}  // namespace mudock
