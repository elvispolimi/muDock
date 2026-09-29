#!/usr/bin/env bash
IGNORE_TOKEN="Experiment"

MAX_JOBS=12

BUILD=build
DATA_DIR="./data/coreset_CASF_2016"
TEST_DIR="./lorenzo_temp/experiments/test"

ids=(3gy4 3ryj 3ehy 4de1 4lzs 2yki 3oe5 4ty7 3o9i 3uri)
crystal_scores=(-1.054222 -3.089948 -7.308385 -2.891916 -1.826324 -7.595546 -3.509636 -3.464769 -4.736212 -4.646878)

rhos=("0.80" "0.90" "0.95" )
epsilons=("0.01" "0.0001" )
NUM_SEEDS=20
seeds=($(seq 0 $((NUM_SEEDS - 1))))

TOTAL_RUNS=$((${#ids[@]} * ${#rhos[@]} * ${#epsilons[@]} * ${#seeds[@]}))
COMPLETED_RUNS=0

rm ${TEST_DIR}/*
echo "Docking..."
for SEED in "${seeds[@]}"; do
    for i in "${!ids[@]}"; do
        PDBID="${ids[$i]}"
        CRYSTAL_SCORE="${crystal_scores[$i]}"
        for RHO in "${rhos[@]}"; do
            for EPSILON in "${epsilons[@]}"; do

            PROTEIN="${DATA_DIR}/${PDBID}/${PDBID}_protein.pdb"
            LIGAND="${DATA_DIR}/${PDBID}/${PDBID}_ligand.adtmol2"
            
            OUT_PATH="${TEST_DIR}/${PDBID}_${RHO}_${EPSILON}_${SEED}.txt"

            ./builds/"$BUILD"/application/muDock \
                --protein "$PROTEIN" \
                --ligand "$LIGAND" \
                --seed "$SEED" \
                --search lga \
                --population 100 \
                --generations 750 \
                --lsrate 100 \
                --lsit 200 \
                --autostop 1 \
                --crystal_score "$CRYSTAL_SCORE" \
                --tolerance_window 50 \
                --best_score_diff_thld 0.0001 \
                --epsilon "$EPSILON" \
                --rho "$RHO" \
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

sed 's/^Experiment,//' ${TEST_DIR}/* > ./lorenzo_temp/experiments/results_rho_epsilon_max_ls.csv

python3 lorenzo_temp/experiments/notify.py

printf "\nDone.\n"