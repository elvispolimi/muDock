# GPU Profiling Results

I created this directory to serve as the central archive for all the GPU profiling results of the muDock project. My main goal here was to analyze and compare the performance of the Native CUDA code against the Alpaka C++ GPU implementation. To keep everything organized, I divided the data into three separate subdirectories, each focusing on a different level of architectural analysis (macro, micro, and hardware).

### 2. Results in the `profiling_results_macro/` Folder

This folder contains the high-level performance data. I focused on evaluating the overall application throughput (Evaluations per second) and the total latency across various grid sizes (Populations x Generations).

To generate and analyze these results, I used the following scripts from the `scripts/` directory:

- **`profiler_macro_global.py`** and **`profiler_macro_rapid.py`**: I used these to run the docking simulations on the cluster and generate the raw CSV data.
- **`plotter_macro_global.py`** and **`plotter_final_macro.py`**: I used these to read the CSVs and generate the comparative PDFs and speedup heatmaps.

### 3. Results in the `profiling_results_micro/` Folder

This folder contains the results of the micro-profiling I performed. I used NVIDIA Nsight Systems (`nsys`) to isolate the execution time of individual GPU kernels to pinpoint architectural bottlenecks.

I generated the data in this folder using these scripts:

- **`profiler_micro_adt.py`**, **`profiler_micro_genetic.py`**, and **`profiler_micro_geom.py`**: I ran these to trace the specific kernels and extract the average micro-second latency per launch, producing the comparative bar charts.

### 4. Results in the `profiling_results_hardware/` Folder

This folder goes deeper than kernel latency, analyzing the exact hardware utilization and static compiler output. I looked at metrics like streaming multiprocessor occupancy, memory bandwidth, and register pressure.

To extract these low-level metrics, I used:

- **`profiler_ncu.py`**: To parse Nsight Compute logs for hardware metrics.
- **`profiler_ptxas.py`**: To analyze the compiled PTX logs and verify if loop unrolling was causing local memory spilling.
- **`profiler_host_metrics.py`**: To prove that Alpaka does not introduce host-side RAM overhead compared to Native CUDA.
