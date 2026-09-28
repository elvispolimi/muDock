#!/usr/bin/env bash

MAX_JOBS=$(nproc)

BUILD=omp
DATA_DIR="./data/coreset_CASF_2016"

for dir in ./data/coreset_CASF_2016/*/; do
    PDBID=$(basename "$dir")

    INPUT="${DATA_DIR}/${PDBID}/${PDBID}_ligand.mol2"
    OUTPUT="${DATA_DIR}/${PDBID}/${PDBID}_ligand.adtmol2"

    ./builds/"$BUILD"/application/converter \
        --input "$INPUT" \
        --output "$OUTPUT"
done

printf "\nDone.\n"