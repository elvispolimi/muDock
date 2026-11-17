#pragma once

#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/data_formats/batch_soa.hpp>
#include <mudock/data_formats/coordinates_soa.hpp>
#include <mudock/data_formats/protein_soa.hpp>
#include <mudock/soa_template/soa_blocks.hpp>
#include <mudock/soa_template/soa_layout.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {
  
  GENERATE_SOA_LAYOUT(Genetic,
                      SOA_COLUMN(chromosome, chromosomes),
                      SOA_COLUMN(chromosome, next_chromosomes))

  GENERATE_SOA_LAYOUT(Scores,
                      SOA_COLUMN(fp_type, scores))

  GENERATE_SOA_LAYOUT(Best,
                      SOA_COLUMN(chromosome, chromosomes),
                      SOA_COLUMN(fp_type, scores))

  GENERATE_SOA_BLOCKS(GeneticData,
                      SOA_BLOCK(coordinates, Coordinates),
                      SOA_BLOCK(scoring_coordinates, Coordinates),
                      SOA_BLOCK(batch_info, BatchInfo),
                      SOA_BLOCK(genetic, Genetic),
                      SOA_BLOCK(scores, Scores),
                      SOA_BLOCK(best, Best))

  GENERATE_SOA_BLOCKS(ScoringData,
                      SOA_BLOCK(scoring_coordinates, Coordinates),
                      SOA_BLOCK(batch_info, BatchInfo),
                      SOA_BLOCK(scores, Scores))

}
