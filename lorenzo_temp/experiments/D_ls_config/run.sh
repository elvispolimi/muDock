#!/usr/bin/env bash
IGNORE_TOKEN="Experiment"

MAX_JOBS=12

BUILD=build
DATA_DIR="./data/coreset_CASF_2016"
EXP_DIR="./lorenzo_temp/experiments"
TEST_DIR="${EXP_DIR}/test"

ids=(3qgy 4m0z 3u8k 3u8n 3wtj 3zdg 4qac 2r9w 3gv9 4jxs 4kz6 3ao4 3zsx 3zt2 4cig 3ehy 3lka 3nx7 3tsk 4gr0)
crystal_scores=(-2.307344 -5.341486 -2.007125 -3.749675 -1.827195 -1.108771 -1.540321 -5.982524 -0.503303 -3.226412 -0.480906 -2.311208 -1.366008 -2.277802 -4.31475 -7.889799 -0.991583 -8.169148 -12.731228 -17.205564)

lsrates=("5" "50" "100")
lsits=("100" "200" "300")

NUM_SEEDS=20
seeds=($(seq 0 $((NUM_SEEDS - 1))))

POPULATION=100
BUDGET=500000

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
                --max-search-box 10 \
                --placement preserve \
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
