#!/usr/bin/env bash

MAX_JOBS=$(nproc)

BUILD=omp

ids=("1fkb" "1hii" "2ya6" "3udd" "4few" "5cst" "5uez" "5wuk")

TOTAL_RUNS=${#ids[@]}
COMPLETED_RUNS=0

rm ./lorenzo_temp/experiments/test/*

echo "Docking..."
for PDBID in "${ids[@]}"; do
    PROTEIN="./data/${PDBID}/${PDBID}_protein.pdb"
    LIGAND="./data/${PDBID}/${PDBID}_ligand.adtmol2"
                
    OUT_PATH="./lorenzo_temp/experiments/test/${PDBID}_${LSRATE}_${LSIT}_${SEED}.txt"

    ./builds/"$BUILD"/application/local_search/local_search \
        --protein "$PROTEIN" \
        --ligand "$LIGAND" \
        > /dev/null 2>&1
    
    if [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]; then
        wait -n
        COMPLETED_RUNS=$((COMPLETED_RUNS + 1))
        printf "\rProgress: [%d/%d] %3d%%" \
            "$COMPLETED_RUNS" \
            "$TOTAL_RUNS" \
            "$((COMPLETED_RUNS * 100 / TOTAL_RUNS))"
    fi
            
done

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
    wait -n
    COMPLETED_RUNS=$((COMPLETED_RUNS + 1))
    printf "\rProgress: [%d/%d] %3d%%" \
        "$COMPLETED_RUNS" \
        "$TOTAL_RUNS" \
        "$((COMPLETED_RUNS * 100 / TOTAL_RUNS))"
done

python3 lorenzo_temp/experiments/notify.py

printf "\nDone.\n"