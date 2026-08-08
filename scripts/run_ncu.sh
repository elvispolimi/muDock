#!/bin/bash
# Nsight Compute (ncu) Runtime Profiling Script for muDock
# Extrapolates HPC hardware metrics for all optimized kernels

REPO_DIR="/work/onedina/muDock_ON"
CUDA_BIN="$REPO_DIR/build/cuda/application/muDock"
ALPAKA_BIN="$REPO_DIR/build/alpaka-cuda/application/muDock"

# The regex matches any kernel whose name contains adt_score, genetic, or geom_transform
KERNEL_REGEX="adt_score|genetic|geom_transform"
METRICS="sm__warps_active.avg.pct_of_peak_sustained_active,sm__throughput.avg.pct_of_peak_sustained_elapsed,gpu__compute_memory_throughput.avg.pct_of_peak_sustained_elapsed,launch__registers_per_thread"

echo "=== NSIGHT COMPUTE RUNTIME PROFILER ==="
echo "Profiling CUDA Native..."
ncu --set default \
    --kernel-name regex:"$KERNEL_REGEX" \
    --metrics $METRICS \
    -o ncu_report_cuda \
    -f \
    $CUDA_BIN --use CUDA:GPU:0 --protein $REPO_DIR/data/1fkb/1fkb_pocket.pdbqt --ligand $REPO_DIR/data/1fkb/1fkb_ligand.adtmol2 --population 25 --generations 1 --seed 42 > ncu_cuda.txt

echo "Profiling Alpaka CUDA..."
ncu --set default \
    --kernel-name regex:"$KERNEL_REGEX" \
    --metrics $METRICS \
    -o ncu_report_alpaka \
    -f \
    $ALPAKA_BIN --use ALPAKA:GPU:0 --protein $REPO_DIR/data/1fkb/1fkb_pocket.pdbqt --ligand $REPO_DIR/data/1fkb/1fkb_ligand.adtmol2 --population 25 --generations 1 --seed 42 > ncu_alpaka.txt

echo "Done! Check ncu_cuda.txt and ncu_alpaka.txt for the text reports."
echo "You can also open the ncu_report_cuda.ncu-rep and ncu_report_alpaka.ncu-rep files in the Nsight Compute GUI on your local machine."
