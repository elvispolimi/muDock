#!/usr/bin/env bash
IGNORE_TOKEN="Experiment"

MAX_JOBS=16

BUILD=build
DATA_DIR="./data/coreset_CASF_2016"
EXP_DIR="./lorenzo_temp/experiments"
TEST_DIR="${EXP_DIR}/test"

ids=(1bcu 1oyt 2r9w 2zda 3ao4 3bv9 3ehy 3gr2 3gv9 3lka 3nx7 3tsk 3utu 3zso 3zsx 3zt2 4cig 4gr0 4jxs 4kz6)
crystal_scores=(-2.453103 -1.394959 -5.982522 -4.452805 -2.45466 55.23962 -7.889795 34.48156 -0.51968 -0.952467 -8.168454 -12.761602 -10.54346 -1.462153 -1.400095 -2.277798 -4.416257 -17.216919 -3.233058 -0.480943)

# COARSE
lsrates=("5" "50" "100")
lsits=("100" "200" "400" "800" "1600" "3200")

# FINE
# lsrates=("5" "25" "50" "75" "100")
# lsits=("10" "75" "150" "225" "300")

NUM_SEEDS=20
seeds=($(seq 0 $((NUM_SEEDS - 1))))

POPULATION=100
BUDGET=1000000

TOTAL_RUNS=$((${#ids[@]} * ${#lsrates[@]} * ${#lsits[@]} * ${#seeds[@]}))
COMPLETED_RUNS=0

rm ${TEST_DIR}/*
echo "Docking..."
for SEED in "${seeds[@]}"; do
    for i in "${!ids[@]}"; do
        PDBID="${ids[$i]}"
        CRYSTAL_SCORE="${crystal_scores[$i]}"
        for LSRATE in "${lsrates[@]}"; do
            for LSIT in "${lsits[@]}"; do

            POCKET="${DATA_DIR}/${PDBID}/${PDBID}_pocket.pdb"
            LIGAND="${DATA_DIR}/${PDBID}/${PDBID}_ligand.adtmol2"

            EVALS_PER_GENERATION=$((POPULATION + POPULATION * LSRATE * LSIT / 100))
            GENERATIONS=$((BUDGET / EVALS_PER_GENERATION))            
            
            OUT_PATH="${TEST_DIR}/${PDBID}_${LSRATE}_${LSIT}_${SEED}.txt"

            ./builds/"$BUILD"/application/muDock \
                --protein "$POCKET" \
                --ligand "$LIGAND" \
                --seed "$SEED" \
                --search lga \
                --population "$POPULATION" \
                --generations "$GENERATIONS" \
                --lsrate "$LSRATE" \
                --lsit "$LSIT" \
                --crystal_score "$CRYSTAL_SCORE" \
                --placement preserve \
                --max-search-box 5 \
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
    > ${EXP_DIR}/results_ls_config.csv

python3 ${EXP_DIR}/notify.py

printf "\nDone.\n"
