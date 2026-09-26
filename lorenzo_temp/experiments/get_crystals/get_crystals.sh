#!/usr/bin/env bash
IGNORE_TOKEN="Experiment"

MAX_JOBS=$(nproc)

BUILD=build

DATA_DIR="./data/crystals"

TEST_DIR="./lorenzo_temp/experiments/test"

echo "Docking..."
TOTAL_RUNS=$(find "$DATA_DIR" -mindepth 1 -maxdepth 1 -type d | wc -l)
COMPLETED_RUNS=0
rm ${TEST_DIR}/*

for dir in ${DATA_DIR}/*/; do
    PDBID=$(basename "$dir")

    PROTEIN="${DATA_DIR}/${PDBID}/${PDBID}_protein.pdb"
    LIGAND="${DATA_DIR}/${PDBID}/${PDBID}_ligand.adtmol2"
    
    OUT_PATH="${TEST_DIR}/${PDBID}_crystal.txt"

    ./builds/"$BUILD"/application/muDock \
        --protein "$PROTEIN" \
        --ligand "$LIGAND" \
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

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
    wait -n
    COMPLETED_RUNS=$((COMPLETED_RUNS + 1))
    printf "\rProgress: [%d/%d] %3d%%" \
        "$COMPLETED_RUNS" \
        "$TOTAL_RUNS" \
        "$((COMPLETED_RUNS * 100 / TOTAL_RUNS))"
done

cat ${TEST_DIR}/* > ./lorenzo_temp/experiments/get_crystals/crystals.csv

# sed 's/^Experiment,//' ${TEST_DIR}/* > ./lorenzo_temp/experiments/results.csv
python3 lorenzo_temp/experiments/notify.py

printf "\nDone.\n"