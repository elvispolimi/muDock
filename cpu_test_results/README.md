# CPU Profiling Results

I created this directory to store the final macro-profiling results for the CPU backends of muDock. Here, I compare the performance of the Native C++ implementation against the new Alpaka C++ version, focusing on both single-thread (Serial) and multi-thread (OpenMP) execution. The main goal of this data is to measure the total throughput in evaluations per second and to validate the performance improvements I made during the porting process.

### 2. Results in the `single/` Folder

The `single/` folder contains the definitive, most optimized results for the 1fkb dataset (a single ligand). These results include the final C++17 memory fix for the random number generator, which allowed the Alpaka Serial backend to become faster than the original Native C++ code.

To generate and analyze these results, I used the following scripts from the `scripts/` directory:

- **`profiler_macro_cpu.py`**: I ran this script on the cluster to execute the docking simulations and generate the raw data in `macro_results_cpu_single.csv`.
- **`plotter_final_cpu.py`**: I used this script to read the CSV data and create the standard PDF graphs that compare throughput and latency across all CPU backends.

### 3. Results in the `single_old_rng/` Folder

The `single_old_rng/` folder serves as a historical archive. It contains the baseline performance data for the same 1fkb dataset, but it was generated *before* I applied the final memory fix for the generator. In these old results, the Alpaka CPU version was about 28% slower than Native C++ because of a memory bottleneck caused by copying a 5KB state per thread.

I generated the data in this folder using the exact same scripts:

- **`profiler_macro_cpu.py`**: To run the simulations and generate the old CSV data.
- **`plotter_advanced_cpu.py`**: I used this advanced plotting script on this old data to generate specific heatmaps (`cpu_advanced_02_rng_gap_heatmap.pdf`) to study the 28% performance gap before I finally fixed it.
