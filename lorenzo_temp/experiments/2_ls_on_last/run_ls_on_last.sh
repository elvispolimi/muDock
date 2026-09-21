#!/usr/bin/env bash

# ids=("1fkb" "1hii" "2ya6" "3udd" "4few" "5cst" "5uez" "5wuk")
# lsrates=("10" "25" "50" "75" "100")
# lsits=("50" "100" "150" "200" "250" "300")
# seeds=("0" "11111" "22222" "33333" "44444" "55555")

# for PDBID in "${ids[@]}"; do
# for LSRATE in "${lsrates[@]}"; do
# for LSIT in "${lsits[@]}"; do
# for SEED in "${seeds[@]}"; do

IGNORE_TOKEN="Experiment"
MAX_JOBS=$(nproc)
BUILD=omp
GENERATIONS=5
NUM_SEEDS=5

ids=("5uez")
lsrates=("10" "25" "50")
lsits=("10" "50" "100")
seeds=($(seq 0 $((NUM_SEEDS - 1))))

TOTAL_RUNS_GA=$((${#ids[@]} * ${#seeds[@]}))
TOTAL_RUNS_LGA=$((${#ids[@]} * ${#lsrates[@]} * ${#lsits[@]} * ${#seeds[@]}))

TOTAL_RUNS=$((${TOTAL_RUNS_GA} + ${TOTAL_RUNS_LGA}))
COMPLETED_RUNS=0

rm ./lorenzo_temp/experiments/test/*

echo "Docking..."

SEARCH=genetic
for PDBID in "${ids[@]}"; do
    PROTEIN="./data/${PDBID}/${PDBID}_protein.pdb"
    LIGAND="./data/${PDBID}/${PDBID}_ligand.adtmol2"
    for SEED in "${seeds[@]}"; do
        
        OUT_PATH="./lorenzo_temp/experiments/test/${PDBID}_${SEED}.txt"

        ./builds/"$BUILD"/application/muDock \
            --protein "$PROTEIN" \
            --ligand "$LIGAND" \
            --seed "$SEED" \
            --search "$SEARCH" \
            --generations "$GENERATIONS" \
            2>&1 | grep "$IGNORE_TOKEN" >> "$OUT_PATH" \
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


SEARCH=lga
for PDBID in "${ids[@]}"; do
    PROTEIN="./data/${PDBID}/${PDBID}_protein.pdb"
    LIGAND="./data/${PDBID}/${PDBID}_ligand.adtmol2"
    for LSRATE in "${lsrates[@]}"; do
        for LSIT in "${lsits[@]}"; do
            for SEED in "${seeds[@]}"; do
                
                OUT_PATH="./lorenzo_temp/experiments/test/${PDBID}_${LSRATE}_${LSIT}_${SEED}.txt"

                ./builds/"$BUILD"/application/muDock \
                    --protein "$PROTEIN" \
                    --ligand "$LIGAND" \
                    --seed "$SEED" \
                    --search "$SEARCH" \
                    --generations "$GENERATIONS" \
                    --lsrate "$LSRATE" \
                    --lsit "$LSIT" \
                    --ls_on_last 2 \
                    --ls_on_best 1 \
                    --ls_every 99999999 \
                    2>&1 | grep "$IGNORE_TOKEN" >> "$OUT_PATH" \
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



sed 's/^Experiment,//' ./lorenzo_temp/experiments/test/* > ./lorenzo_temp/experiments/results.csv
python3 lorenzo_temp/experiments/notify.py

printf "\nDone.\n"