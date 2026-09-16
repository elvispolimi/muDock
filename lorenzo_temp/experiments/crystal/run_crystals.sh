#!/usr/bin/env bash

MAX_JOBS=4

BUILD=omp

ids=("1fkb" "1hii" "2ya6" "3udd" "4few" "5cst" "5uez" "5wuk")
lsits=("10" "25" "50" "100" "200")

TOTAL_RUNS=$((${#ids[@]} * ${#lsits[@]}))
COMPLETED_RUNS=0

rm ./lorenzo_temp/experiments/test/*

echo "Docking..."
for PDBID in "${ids[@]}"; do
    for LSIT in "${lsits[@]}"; do
        PROTEIN="./data/${PDBID}/${PDBID}_protein.pdb"
        LIGAND="./data/${PDBID}/${PDBID}_ligand.adtmol2"

        ./builds/"$BUILD"/application/local_search/local_search \
            --protein "$PROTEIN" \
            --ligand "$LIGAND" \
            --lsit "$LSIT" \
            > /dev/null 2>&1 \
            &
        
        if [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]; then
            wait -n
            COMPLETED_RUNS=$((COMPLETED_RUNS + 1))
            printf "\rProgress: [%d/%d] %3d%%" \
                "$COMPLETED_RUNS" \
                "$TOTAL_RUNS" \
                "$((COMPLETED_RUNS * 100 / TOTAL_RUNS))"
        fi
    done      
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