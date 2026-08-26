#!/usr/bin/env bash

# ids=("1fkb" "1hii" "2ya6" "3udd" "4few" "5cst" "5uez" "5wuk")
# lsrates=("10" "25" "50" "75" "100")
# lsits=("50" "100" "150" "200" "250" "300")
# seeds=("0" "11111" "22222" "33333" "44444" "55555")

# for PDBID in "${ids[@]}"; do
# for LSRATE in "${lsrates[@]}"; do
# for LSIT in "${lsits[@]}"; do
# for SEED in "${seeds[@]}"; do

MAX_JOBS=8

POPULATION=100
GENERATIONS=100

ids=("5uez")
lsrates=("10" "50" "100")
lsits=("10" "150" "300")
seeds=("11111" "22222" "33333" "44444" "55555")

for PDBID in "${ids[@]}"; do
    PROTEIN="./data/${PDBID}/${PDBID}_protein.pdb"
    LIGAND="./data/${PDBID}/${PDBID}_ligand.adtmol2"
    for LSRATE in "${lsrates[@]}"; do
        for LSIT in "${lsits[@]}"; do
            for SEED in "${seeds[@]}"; do
                
                echo "=== Running ligand=${PDBID}, lsrate=${LSRATE}, lsit=${LSIT}, seed=${SEED} ==="

                OUT_PATH="./script/experiments/test/${PDBID}_${LSRATE}_${LSIT}_${SEED}.txt"

                ./builds/vanilla/application/muDock \
                    --protein "$PROTEIN" \
                    --ligand "$LIGAND" \
                    --seed "$SEED" \
                    --search lga \
                    --autostop 1 \
                    --generations "$GENERATIONS" \
                    --population "$POPULATION" \
                    --lsrate "$LSRATE" \
                    --lsit "$LSIT" \
                    --use CPP:CPU:0 \
                    2>&1 | grep Exp | sed "s/$/ $POPULATION $LSRATE $LSIT/" >> "$OUT_PATH" \
                    &
                
                if [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]; then
                    wait -n
                fi
            
            done
        done
    done
done

wait

cat ./script/experiments/test/* > ./script/experiments/results.txt

echo "Done."