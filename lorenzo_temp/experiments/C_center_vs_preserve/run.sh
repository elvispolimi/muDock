#!/usr/bin/env bash
IGNORE_TOKEN="Experiment"

MAX_JOBS=12

BUILD=build
DATA_DIR="./data"
EXP_DIR="./lorenzo_temp/experiments"
TEST_DIR="${EXP_DIR}/test"

ids=("1fkb" "1hii" "2ya6" "3udd" "4few" "5cst" "5uez" "5wuk")
crystal_scores=(-15.797226 -15.021526 -10.792323 -18.935410 244.454712 -3.826046 -9.453270 -10.575790)

# ids=("1fkb")
# crystal_scores=(-15.797226)

placements=("preserve" "center_bbox")

NUM_SEEDS=10
seeds=($(seq 0 $((NUM_SEEDS - 1))))

TOTAL_RUNS=$((${#ids[@]} * ${#seeds[@]} * ${#placements[@]}))
COMPLETED_RUNS=0

rm -rf ${TEST_DIR}/*
echo "Docking..."
for PLACEMENT in "${placements[@]}"; do
    mkdir ${TEST_DIR}/${PLACEMENT}
    for i in "${!ids[@]}"; do
        PDBID="${ids[$i]}"
        CRYSTAL_SCORE="${crystal_scores[$i]}"
        for SEED in "${seeds[@]}"; do
            PROTEIN="${DATA_DIR}/${PDBID}/${PDBID}_protein.pdb"
            LIGAND="${DATA_DIR}/${PDBID}/${PDBID}_ligand.adtmol2"
            
            OUT_PATH="${TEST_DIR}/${PLACEMENT}/${PDBID}_${SEED}.txt"

            ./builds/"$BUILD"/application/muDock \
                --protein "$PROTEIN" \
                --ligand "$LIGAND" \
                --seed "$SEED" \
                --search lga \
                --population 100 \
                --generations 750 \
                --lsrate 50 \
                --lsit 100 \
                --autostop 1 \
                --crystal_score "$CRYSTAL_SCORE" \
                --tolerance_window 50 \
                --best_score_diff_thld 0.0001 \
                --max-search-box 10 \
                --placement "$PLACEMENT" \
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

    while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
        wait -n
        COMPLETED_RUNS=$((COMPLETED_RUNS + 1))
        printf "\rProgress: [%d/%d] %3d%%" \
            "$COMPLETED_RUNS" \
            "$TOTAL_RUNS" \
            "$((COMPLETED_RUNS * 100 / TOTAL_RUNS))"
    done

    grep -h "$IGNORE_TOKEN" "$TEST_DIR"/${PLACEMENT}/* |
        sed 's/^Experiment,//' \
        > ${EXP_DIR}/results_${PLACEMENT}.csv

done



python3 ${EXP_DIR}/notify.py

printf "\nDone.\n"