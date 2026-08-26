#!/usr/bin/env bash

set -euo pipefail

ids=("1fkb" "1hii" "2ya6" "3udd" "4few" "5cst" "5uez" "5wuk")

for PDBID in "${ids[@]}"; do
    PROTEIN="./data/${PDBID}/${PDBID}_protein.pdb"
    LIGAND="./data/${PDBID}/${PDBID}_ligand.adtmol2"
    OUT_DIR="plots/lsit/${PDBID}"

    rm -rf "$OUT_DIR"

    mkdir -p "$OUT_DIR"

    for LSIT in 0 5 10 15 20 25 40 50 100 150 200 250 300; do
        echo "=== Running with lsit=${LSIT} ==="

        rm -f genetic.csv

        ./builds/omp/application/muDock \
            --protein "$PROTEIN" \
            --ligand "$LIGAND" \
            --seed 0 \
            --generations 100 \
            --population 100 \
            --search lga \
            --lsit "$LSIT"

        python3 plot_metric.py genetic.csv \
            --out-dir "$OUT_DIR" \
            -o "$LSIT"

        cp genetic.csv "${OUT_DIR}/${LSIT}_iter.csv"
    done

    rm -f genetic.csv
done

echo "Done."