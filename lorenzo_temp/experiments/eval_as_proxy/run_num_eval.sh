#!/usr/bin/env bash

# ids=("1fkb" "1hii" "2ya6" "3udd" "4few" "5cst" "5uez" "5wuk")
# lsrates=("10" "25" "50" "75" "100")
# lsits=("50" "100" "150" "200" "250" "300")
# seeds=("0" "11111" "22222" "33333" "44444" "55555")

# for PDBID in "${ids[@]}"; do
# for LSRATE in "${lsrates[@]}"; do
# for LSIT in "${lsits[@]}"; do
# for SEED in "${seeds[@]}"; do

IGNORE_TOKEN="INFO All Done!"

MAX_JOBS=$(nproc)

BUILD=omp
POPULATION=100
GENERATIONS=1000
NUM_SEEDS=1
ids=("1fkb" "1hii" "2ya6" "3udd" "4few" "5cst" "5uez" "5wuk")
seeds=($(seq 0 $((NUM_SEEDS - 1))))

# budget in number of evaluations
NUM_EVALUATIONS_GA=$((POPULATION * GENERATIONS))
echo BUDGET: ${NUM_EVALUATIONS_GA} evaluations


TEST_DIR="./lorenzo_temp/experiments/test"


echo "Docking..."
TOTAL_RUNS=$((${#ids[@]} * ${#seeds[@]}))
COMPLETED_RUNS=0
for PDBID in "${ids[@]}"; do
    rm ${TEST_DIR}/*
    PROTEIN="./data/${PDBID}/${PDBID}_protein.pdb"
    LIGAND="./data/${PDBID}/${PDBID}_ligand.adtmol2"
    
    SEARCH=genetic
    for SEED in "${seeds[@]}"; do
        OUT_PATH="${TEST_DIR}/${PDBID}_GA_${SEED}.txt"

        ./builds/"$BUILD"/application/muDock \
            --protein "$PROTEIN" \
            --ligand "$LIGAND" \
            --seed "$SEED" \
            --search "$SEARCH" \
            --generations "$GENERATIONS" \
            --population "$POPULATION" \
            2>&1 | grep "$IGNORE_TOKEN" >> "$OUT_PATH" \
            # &
        
        if [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]; then
            wait -n
            COMPLETED_RUNS=$((COMPLETED_RUNS + 1))
            printf "\rProgress: [%d/%d] %3d%%" \
                "$COMPLETED_RUNS" \
                "$TOTAL_RUNS" \
                "$((COMPLETED_RUNS * 100 / TOTAL_RUNS))"
        fi
    done

    lsrates=("10" "50" "100")
    lsits=("10" "50" "100")
    SEARCH=lga
    for LSRATE in "${lsrates[@]}"; do
        for LSIT in "${lsits[@]}"; do
            # Compute the number of generations with the same evaluation budget as the GA
            GENERATIONS=$((NUM_EVALUATIONS_GA / (POPULATION * (1 + LSRATE * LSIT / 100))))
            for SEED in "${seeds[@]}"; do
                
                OUT_PATH="${TEST_DIR}/${PDBID}_${LSRATE}_${LSIT}_${SEED}.txt"

                ./builds/"$BUILD"/application/muDock \
                    --protein "$PROTEIN" \
                    --ligand "$LIGAND" \
                    --seed "$SEED" \
                    --search "$SEARCH" \
                    --generations "$GENERATIONS" \
                    --population "$POPULATION" \
                    --lsrate "$LSRATE" \
                    --lsit "$LSIT" \
                    2>&1 | grep "$IGNORE_TOKEN" >> "$OUT_PATH" \
                    # &
                
                if [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]; then
                    wait -n
                    COMPLETED_RUNS_LGA=$((COMPLETED_RUNS_LGA + 1))
                    printf "\rProgress: [%d/%d] %3d%%" \
                        "$COMPLETED_RUNS_LGA" \
                        "$TOTAL_RUNS_LGA" \
                        "$((COMPLETED_RUNS_LGA * 100 / TOTAL_RUNS_LGA))"
                fi
            
            done
        done
    done
    cat ${TEST_DIR}/* > ./lorenzo_temp/experiments/eval_as_proxy/${PDBID}.csv
done


TOTAL_RUNS_LGA=$((${#ids[@]} * ${#lsrates[@]} * ${#lsits[@]} * ${#seeds[@]}))


while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
    wait -n
    COMPLETED_RUNS_LGA=$((COMPLETED_RUNS_LGA + 1))
    printf "\rProgress: [%d/%d] %3d%%" \
        "$COMPLETED_RUNS_LGA" \
        "$TOTAL_RUNS_LGA" \
        "$((COMPLETED_RUNS_LGA * 100 / TOTAL_RUNS_LGA))"
done


# sed 's/^Experiment,//' ${TEST_DIR}/* > ./lorenzo_temp/experiments/results.csv
python3 lorenzo_temp/experiments/notify.py

printf "\nDone.\n"