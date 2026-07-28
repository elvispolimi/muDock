#!/usr/bin/env bash

set -euo pipefail

ids=("1fkb" "1hii" "2ya6" "3udd" "4few" "5cst" "5uez" "5wuk")

for PDBID in "${ids[@]}"; do
    PROTEIN="./data/${PDBID}/${PDBID}_protein.pdb"
    LIGAND="./data/${PDBID}/${PDBID}_ligand.adtmol2"
    OUT_DIR="plots/local_search/${PDBID}"

    rm -rf "$OUT_DIR"

    mkdir -p "$OUT_DIR"

    echo "=== Running with id=${PDBID} ==="

    rm -f adadelta_scores.csv
    rm -f adadelta_com.csv

    ./builds/omp/application/local_search/local_search \
        --protein "$PROTEIN" \
        --ligand "$LIGAND" \

    python3 plot_metric.py adadelta_scores.csv \
        --out-dir "$OUT_DIR"

    python3 plot_metric.py adadelta_com.csv \
        --out-dir "$OUT_DIR"

    cp adadelta_scores.csv "${OUT_DIR}/adadelta_scores.csv"
    cp adadelta_com.csv "${OUT_DIR}/adadelta_com.csv"

done

echo "Done."