# Old Test Results Archive

I have created this directory to serve as a historical archive for the intermediate profiling data that I generated during the various optimization phases of the muDock project. The results stored here do not represent the final performance of the application. Instead, I have retained them strictly to scientifically document the performance delta and the impact of the specific patches I applied to the codebase and build system over time.

I organized this archive into several subdirectories to map the evolution of the project. In the pre-batch multiple folder, I kept the results generated before I implemented the occupancy optimization that batches multiple ligands per kernel launch. These baseline runs lacked the high-level parallelism required to fully saturate GPU warps. 

In the pre-cmake folder, I stored the results generated right before I fixed a critical bug in the CMake configuration. At that time, I had already applied the Alpaka C++ micro-optimizations, but due to the bug, the NVCC architectural fast-math flags were not being injected. This resulted in Alpaka performing worse than Native CUDA because mathematical operations were bottlenecked by emulated precision casts. 

Immediately following that, I created the post-cmake folder to document the results after I resolved the bug. This folder contains the historical proof of Alpaka's GPU success, demonstrating how the performance jumped from 6000 to over 7500 evaluations per second, completely closing the gap with Native CUDA. I later promoted the final, cleaned version of these results to the official GPU test results directory.

I also kept a pre-macro folder with results from before I removed the legacy C preprocessor macros, back when the codebase had not yet been refactored into cleaner Alpaka C++ templates. Finally, the old report folders contain the initial, unpolished iterations of my macro and micro-profiling scripts, mostly holding redundant raw data and obsolete graphs from the early stages of building my Python profiling suite.
