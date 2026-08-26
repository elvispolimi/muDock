#!/usr/bin/env bash

set -euo pipefail

ids=("1fkb" "1hii" "2ya6" "3udd" "4few" "5cst" "5uez" "5wuk")

for PDBID in "${ids[@]}"; do
    PROTEIN="./data/${PDBID}/${PDBID}_protein.pdb"
    LIGAND="./data/${PDBID}/${PDBID}_ligand.adtmol2"
    OUT_DIR="plots/lsrate/${PDBID}"

    rm -rf "$OUT_DIR"

    mkdir -p "$OUT_DIR"

    for LSRATE in 0 1 2 5 10 25 50 75 100; do
        echo "=== Running with lsrate=${LSRATE} ==="

        rm -f genetic.csv

        ./builds/omp/application/muDock \
            --protein "$PROTEIN" \
            --ligand "$LIGAND" \
            --seed 0 \
            --generations 100 \
            --population 100 \
            --search lga \
            --lsrate "$LSRATE"

        python3 plot_metric.py genetic.csv \
            --out-dir "$OUT_DIR" \
            -o "$LSRATE"

        cp genetic.csv "${OUT_DIR}/${LSRATE}_rate.csv"
    done

    rm -f genetic.csv
done

echo "Done."