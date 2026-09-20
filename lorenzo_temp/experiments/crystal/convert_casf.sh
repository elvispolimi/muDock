#!/usr/bin/env bash

MAX_JOBS=$(nproc)

BUILD=omp
DATA_DIR="./data/crystals"

for dir in ./data/crystals/*/; do
    PDBID=$(basename "$dir")

    INPUT="${DATA_DIR}/${PDBID}/${PDBID}_ligand.mol2"
    OUTPUT="${DATA_DIR}/${PDBID}/${PDBID}_ligand.adtmol2"

    ./builds/"$BUILD"/application/converter \
        --input "$INPUT" \
        --output "$OUTPUT"
done

printf "\nDone.\n"