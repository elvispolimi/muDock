#!/usr/bin/env bash
IGNORE_TOKEN="Experiment"

MAX_JOBS=12

BUILD=build
DATA_DIR="./data"
TEST_DIR="./lorenzo_temp/experiments/test"

ids=(1fkb 1hii 2ya6 3udd 4few 5cst 5uez 5wuk)
# ids=(1fkb)

NUM_SEEDS=20
seeds=($(seq 0 $((NUM_SEEDS - 1))))
SEARCH=genetic
GENERATIONS=1000

TOLERANCE_WINDOW=30
BEST_SCORE_THLD=0.00001


echo "Docking..."
# TOTAL_RUNS=$(find "$DATA_DIR" -mindepth 1 -maxdepth 1 -type d | wc -l)
TOTAL_RUNS=$((${#ids[@]} * ${#seeds[@]}))
COMPLETED_RUNS=0

# for dir in ${DATA_DIR}/*/; do
    # PDBID=$(basename "$dir")
rm ${TEST_DIR}/AS/*
# rm ${TEST_DIR}/NOAS/*
AUTOSTOP=1
for i in "${!ids[@]}"; do
    PDBID="${ids[$i]}"
    for SEED in "${seeds[@]}"; do

        PROTEIN="${DATA_DIR}/${PDBID}/${PDBID}_protein.pdb"
        LIGAND="${DATA_DIR}/${PDBID}/${PDBID}_ligand.adtmol2"
        
        OUT_PATH="${TEST_DIR}/AS/${PDBID}_${SEED}_AS.txt"

        ./builds/"$BUILD"/application/muDock \
            --protein "$PROTEIN" \
            --ligand "$LIGAND" \
            --seed "$SEED" \
            --search "$SEARCH" \
            --generations "$GENERATIONS" \
            --autostop "$AUTOSTOP" \
            --tolerance_window "$TOLERANCE_WINDOW" \
            --best_score_diff_thld "$BEST_SCORE_THLD" \
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

# AUTOSTOP=0
# for i in "${!ids[@]}"; do
#     PDBID="${ids[$i]}"
#     for SEED in "${seeds[@]}"; do

#         PROTEIN="${DATA_DIR}/${PDBID}/${PDBID}_protein.pdb"
#         LIGAND="${DATA_DIR}/${PDBID}/${PDBID}_ligand.adtmol2"
        
#         OUT_PATH="${TEST_DIR}/NOAS/${PDBID}_${SEED}.txt"

#         ./builds/"$BUILD"/application/muDock \
#             --protein "$PROTEIN" \
#             --ligand "$LIGAND" \
#             --seed "$SEED" \
#             --search "$SEARCH" \
#             --generations "$GENERATIONS" \
#             --autostop "$AUTOSTOP" \
#             2>&1 | grep "$IGNORE_TOKEN" >> "$OUT_PATH" \
#             &

#         if [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]; then
#             wait -n
#             COMPLETED_RUNS=$((COMPLETED_RUNS + 1))
#             printf "\rProgress: [%d/%d] %3d%%" \
#                 "$COMPLETED_RUNS" \
#                 "$TOTAL_RUNS" \
#                 "$((COMPLETED_RUNS * 100 / TOTAL_RUNS))"
#         fi
#     done
# done

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
    wait -n
    COMPLETED_RUNS=$((COMPLETED_RUNS + 1))
    printf "\rProgress: [%d/%d] %3d%%" \
        "$COMPLETED_RUNS" \
        "$TOTAL_RUNS" \
        "$((COMPLETED_RUNS * 100 / TOTAL_RUNS))"
done

# cat ${TEST_DIR}/* > ./lorenzo_temp/experiments/A_hyperparam/results.csv
sed 's/^Experiment,//' ${TEST_DIR}/AS/* > ./lorenzo_temp/experiments/results_${TOLERANCE_WINDOW}_${BEST_SCORE_THLD}_AS.csv
# sed 's/^Experiment,//' ${TEST_DIR}/NOAS/* > ./lorenzo_temp/experiments/results_${TOLERANCE_WINDOW}_${BEST_SCORE_THLD}_NOAS.csv


python3 lorenzo_temp/experiments/notify.py

printf "\nDone.\n"