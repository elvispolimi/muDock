#!/bin/bash
# run_all.sh
# Script to run all macro profilers (GPU and CPU) and generate plots.

set -e

echo "=========================================================="
echo "  1) Running GPU Profiler (macro_rapid)"
echo "=========================================================="
/home/onedina/.conda/envs/mudock_profiling/bin/python profiler_macro_rapid.py

echo "=========================================================="
echo "  2) Generating GPU Plots (macro_rapid)"
echo "=========================================================="
/home/onedina/.conda/envs/mudock_profiling/bin/python plotter_macro_rapid.py

echo "=========================================================="
echo "  3) Running CPU Profiler (macro_cpu)"
echo "=========================================================="
/home/onedina/.conda/envs/mudock_profiling/bin/python profiler_macro_cpu.py

echo "=========================================================="
echo "  4) Generating CPU Plots (macro_cpu)"
echo "=========================================================="
/home/onedina/.conda/envs/mudock_profiling/bin/python plotter_macro_cpu.py

echo "=========================================================="
echo "  All profiling suites completed successfully!"
echo "=========================================================="
