#!/usr/bin/env bash

BUILD=omp
DATA_DIR="./data/crystals"

dirs=($DATA_DIR/*/)
TOTAL_RUNS=${#dirs[@]}
COMPLETED_RUNS=0

echo "Running Local Search on crystals..."

for dir in "${dirs[@]}"; do
    PDBID=$(basename "$dir")

    PROTEIN="${DATA_DIR}/${PDBID}/${PDBID}_protein.pdb"
    LIGAND="${DATA_DIR}/${PDBID}/${PDBID}_ligand.adtmol2"

    ./builds/"$BUILD"/application/local_search/local_search \
        --protein "$PROTEIN" \
        --ligand "$LIGAND" \
        > /dev/null 2>&1
    
    COMPLETED_RUNS=$((COMPLETED_RUNS + 1))
    printf "\rProgress: [%d/%d] %3d%%" \
        "$COMPLETED_RUNS" \
        "$TOTAL_RUNS" \
        "$((COMPLETED_RUNS * 100 / TOTAL_RUNS))"
     
done

python3 lorenzo_temp/experiments/notify.py

printf "\nDone.\n"