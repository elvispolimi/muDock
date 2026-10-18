#!/usr/bin/env bash
IGNORE_TOKEN="Experiment"

MAX_JOBS=12

BUILD=build
DATA_DIR="./data/coreset_CASF_2016"
EXP_DIR="./lorenzo_temp/experiments"
TEST_DIR="${EXP_DIR}/test"

ids=(3p5o 3u5j 4lzs 4ogj 4wiv 3ebp 3f3a 3f3c 3f3d 3f3e 4mme 4e6q 4f09 4gfm 4hge 3nw9 3oe4 3oe5 3ozs 3ozt)
crystal_scores=(-5.091403 -4.212269 -3.624627 -1.08664 -5.579877 -1.070766 -3.952088 -5.573053 -4.529698 -5.970263 -5.017744 -2.686873 -3.135863 -3.722687 -4.370185 -8.040822 -9.958447 -3.980942 -0.716404 -1.978494)

lsrates=("10" "50" "100")
lsits=("10" "150" "300")

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
