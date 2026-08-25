# Profiling Scripts

This folder contains the complete Python-based ecosystem used to profile and evaluate the muDock engine. I use these scripts to run benchmarks on the cluster and generate publication-ready plots.

- `plotter_advanced_cpu.py` -> Generates advanced analytical heatmaps for CPU performance penalties.
- `plotter_final_cpu.py` -> Creates standard latency and throughput charts comparing Native C++ against Alpaka on CPU.
- `plotter_final_macro.py` -> Aggregates results to create final paper-quality macro charts for GPU.
- `plotter_macro_cpu.py` -> Produces global CPU throughput and latency bar charts.
- `plotter_macro_global.py` -> Creates global throughput charts and speedup heatmaps for the large GPU datasets.
- `plotter_macro_rapid.py` -> Produces quick comparative plots for immediate feedback during development.
- `profiler_host_metrics.py` -> Monitors system-level resource consumption like RAM and Host CPU utilization.
- `profiler_macro_cpu.py` -> Executes the full CPU test matrix across different backends and parameters.
- `profiler_macro_global.py` -> Executes the heavy GPU benchmarking grid comparing Native CUDA against Alpaka CUDA.
- `profiler_macro_rapid.py` -> Runs a lightweight, single-configuration benchmark for rapid GPU validation.
- `profiler_micro_adt.py` -> Traces the execution time of the `calc_energy` kernel using Nsight Systems.
- `profiler_micro_genetic.py` -> Profiles the evolutionary algorithm kernels (`initialize`, `iterate`, `finalize`) using Nsight Systems.
- `profiler_micro_geom.py` -> Analyzes the `apply` kernel (geometric transformation) using Nsight Systems.
- `profiler_ncu.py` -> Parses Nsight Compute logs to extract hardware metrics like SM occupancy and memory bandwidth.
- `profiler_ptxas.py` -> Analyzes compiler output logs to extract static register limits and local memory spilling.
- `run_all.sh` -> Automates the execution of multiple profiling jobs sequentially.
- `run_ncu.sh` -> Wraps and executes the NVIDIA Nsight Compute hardware profiler.
