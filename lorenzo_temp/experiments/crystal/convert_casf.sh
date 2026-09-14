#!/usr/bin/env bash

MAX_JOBS=$(nproc)

BUILD=omp

ids=("1a30" "1bcu" "1e66" "1f8b" "1f8c" "1f8d" "1gpk" "1h23" "1hfs" "1hnn")

for PDBID in "${ids[@]}"; do
    INPUT="./data/prova_crystals/${PDBID}/${PDBID}_ligand.mol2"
    OUTPUT="./data/prova_crystals/${PDBID}/${PDBID}_ligand.adtmol2"
                

    ./builds/"$BUILD"/application/converter \
        --input "$INPUT" \
        --output "$OUTPUT"
            
done

printf "\nDone.\n"