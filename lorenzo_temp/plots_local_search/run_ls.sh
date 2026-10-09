#!/usr/bin/env bash

set -euo pipefail

ids=("1fkb" "1hii" "2ya6" "3udd" "4few" "5cst" "5uez" "5wuk")
# ids=("1fkb")
# rhos=("0.80" "0.90" "0.99" )
rhos=("0.95" )
epsilons=("0.000001" "0.0001" "0.01")

for RHO in "${rhos[@]}"; do
    for EPSILON in "${epsilons[@]}"; do
        for PDBID in "${ids[@]}"; do
            PROTEIN="./data/${PDBID}/${PDBID}_protein.pdb"
            LIGAND="./data/${PDBID}/${PDBID}_ligand.adtmol2"
            OUT_DIR="./rho_${RHO}_epsilon_${EPSILON}/${PDBID}"

            rm -rf "$OUT_DIR"

            mkdir -p "$OUT_DIR"

            echo "=== Running with id=${PDBID} ==="

            rm -f adadelta_scores.csv
            rm -f adadelta_com.csv

            ./builds/omp/application/local_search/local_search \
                --protein "$PROTEIN" \
                --ligand "$LIGAND" \
                --placement preserve \
                --rho "$RHO" \
                --epsilon "$EPSILON"



            python3 ./lorenzo_temp/plots_local_search/plot_metric.py adadelta_scores.csv \
                --out-dir "$OUT_DIR"

            python3 ./lorenzo_temp/plots_local_search/plot_metric.py adadelta_com.csv \
                --out-dir "$OUT_DIR"

            cp adadelta_scores.csv "${OUT_DIR}/adadelta_scores.csv"
            cp adadelta_com.csv "${OUT_DIR}/adadelta_com.csv"

        done
    done
done

echo "Done."