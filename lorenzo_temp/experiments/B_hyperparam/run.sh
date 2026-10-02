#!/usr/bin/env bash
IGNORE_TOKEN="Experiment"

MAX_JOBS=12

BUILD=build
DATA_DIR="./data/coreset_CASF_2016"
EXP_DIR="./lorenzo_temp/experiments"
TEST_DIR="${EXP_DIR}/test"

ids=(1mq6 1z95 2cet 3ary 3coz 3dx2 3fur 3gv9 3k5v 3lka 3nw9 3p5o 3qgy 3qqs 4bkt 4dld 4j21 4kz6 4ty7 5tmn)
crystal_scores=(-8.766284 -4.041984 -0.100647 -0.404427 -1.644277 -2.943025 -1.518824 -0.644809 -0.654886 -4.840909 -7.317592 -2.251232 -1.48791 -5.76047 -3.921186 -5.967893 -8.883362 -0.303014 -3.464769 -10.883055)

rhos=("0.80" "0.90" "0.95" )
epsilons=("0.01" "0.000001" )
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
                > "$OUT_PATH" 2>&1 &

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

grep -h "$IGNORE_TOKEN" "$TEST_DIR"/* |
    sed 's/^Experiment,//' \
    > ${EXP_DIR}/results_rho_epsilon_max_ls.csv

python3 ${EXP_DIR}/notify.py

printf "\nDone.\n"