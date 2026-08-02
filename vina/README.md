# VINA

## Introduction

As detailed in the main [muDock](../README.md) documentation, muDock is a modular and high-performance molecular docking engine designed for flexibility, vectorization, and accelerator porting. Thanks to its stage-based architecture, muDock allows for the seamless integration and prototyping of alternative scoring models alongside its default energy engine.

This directory implements the **AutoDock Vina scoring function** as an alternative evaluation backend. Developed based on the empirical scoring principles established by Trott and Olson (2010), the Vina scoring model is widely recognized as an industry standard across modern molecular docking tools and protocols due to its accuracy and computational efficiency.

### Background & Academic Context

My work on muDock began a year ago as part of the *Multidisciplinary Project* course. The initial scope involved porting the AutoDock Vina scoring function from Python to C++ and developing a baseline, naive CUDA implementation. 

This foundation provided an ideal bridge into the *GPUs & Heterogeneous Systems* course, where the focus shifted from simple functional porting to deep architectural optimization. The primary objective of this second phase was to transform that initial prototype into a highly optimized, production-ready CUDA kernel, pushing execution performance, memory throughput, and hardware utilization to their limits.

## Scoring Algorithm

The total fitness score $E_{\text{total}}$ is obtained by combining the intermolecular energy ($\mathrm{score}_{\mathrm{inter}}$) and intramolecular energy ($\mathrm{score}_{\mathrm{intra}}$), normalized by the weighted number of active rotatable bonds ($N_{\text{rot}}$):

$$E_{\text{total}} = \frac{\text{score}_{\text{inter}} + \text{score}_{\text{intra}}}{1 + w_{\text{rot}} \cdot N_{\text{rot}}}$$

Where:
*   **$\text{score}_{\text{inter}}$**: Evaluates non-bonded interactions between receptor and ligand atoms.
*   **$\text{score}_{\text{intra}}$**: Evaluates internal steric clashes and interactions within the ligand itself.
*   **$w_{\text{rot}}$**: A weight coefficient penalizing conformational entropy lost upon binding.
*   **$N_{\text{rot}}$**: The number of active rotatable torsions in the ligand.

#### Pairwise Energy Function
Both $\mathrm{score}_{\mathrm{inter}}$ and $\mathrm{score}_{\mathrm{intra}}$ are calculated by summing the pairwise interaction energy $E_{\text{pair}}(d_{ij})$ over all valid **interacting pairs**:

$$\text{score} = \sum_{i,j \in \text{interacting pairs}} E_{\text{pair}}(d_{ij})$$

The pairwise energy $E_{\text{pair}}$ is computed as a function of the surface distance $d_{ij}$ between atom $i$ and atom $j$:

$$d_{ij} = r_{ij} - (R_{ti} + R_{tj})$$

Where $r_{ij}$ is the interatomic distance, and $R_{ti}, R_{tj}$ are the van der Waals radii of the respective atom types.

#### Interaction Terms (Gaussians and Repulsion)
The potential function $E_{\text{pair}}(d_{ij})$ is a weighted linear combination of spatial term components:

$$E_{\text{pair}}(d_{ij}) = w_1 \cdot \text{gauss}_1(d_{ij}) + w_2 \cdot \text{gauss}_2(d_{ij}) + w_{\text{rep}} \cdot \text{repulsion}(d_{ij}) + w_{\text{hbond}} \cdot e_{\text{hbond}}(d_{ij}) + w_{\text{hydrophobic}} \cdot e_{\text{hydrophobic}}(d_{ij})$$

*   **Gaussian 1 (Short-range attraction)**:
    $$\text{gauss}_1(d) = e^{-\left(\frac{d}{0.5\,\text{Å}}\right)^2}$$
*   **Gaussian 2 (Long-range attraction)**:
    $$\text{gauss}_2(d) = e^{-\left(\frac{d - 3\,\text{Å}}{2\,\text{Å}}\right)^2}$$
*   **Repulsion (Steric clash penalty)**:

$$
\text{repulsion}(d) = \begin{cases} d^2 & \text{if } d \lt 0 \\ 0 & \text{if } d \ge 0 \end{cases}
$$

*   **Hydrogen Bonding**: Distance and angle dependent attractive term applied specifically between designated donor and acceptor atom pairs.
*   **Hydrophobic Interactions**: Favorable energy term applied when both atoms $i$ and $j$ are flagged as hydrophobic.

#### Interacting Pairs Definition
An atom pair $(i, j)$ is classified as an **interacting pair** only if it fulfills two structural constraints:
1. The atoms are separated by more than 3 consecutive bonds (they are not 1-2, 1-3, or 1-4 topological neighbors).
2. The atoms belong to different **rigid fragments**.

A *rigid fragment* is a rigid block of atoms whose internal relative spatial coordinates remain invariant regardless of any rotatable bond rotations. Atoms within the same rigid fragment or connected by 3 or fewer bonds are structurally constrained, making pairwise non-bonded scoring physically meaningless for them.

## The files

The Vina scoring implementation is designed around a clean separation between high-level lifecycle management and backend-specific execution kernels. Below is an overview of the introduced files and their responsibilities:
- vina_score.hpp 
- vina_score.cpp 
- vina_score_kernel.hpp 
- vina_score_cpp.hpp 
- vina_score_cpp.cpp 
- vina_score_cuda.cuh 
- vina_score_cuda.cu 

### File Details

*   **`vina_score.hpp` / `vina_score.cpp`**: 
    Serves as the main high-level API for the Vina scoring engine.
    *   **Initialization**: Instantiated by passing the protein structure, setting up the device scratchpad containing all required static receptor parameters and lookup data.
    *   **Batch Preparation**: Features a `prepare()` method that accepts a batch of ligands and stores them into a specialized `vina_score_kernel` execution object.
    *   **Execution Driver**: Overloads `operator()` to delegate scoring execution to the underlying kernel object. It is templated to instantiate the correct device kernel based on the requested backend (e.g., `CPP:CPU` or `CUDA:GPU`).

*   **`vina_score_kernel.hpp`**: 
    Defines the base `vina_score_kernel` structure and interface. It acts as the data container holding the prepared ligand batch and the contract for backend kernel invocation (`operator()`).

*   **`vina_score_cpp.hpp` / `vina_score_cpp.cpp`**: 
    Contains the CPU C++ backend implementation. Provides the specialized host execution logic and helper functions for CPU execution paths.

*   **`vina_score_cuda.cuh` / `vina_score_cuda.cu`**: 
    Contains the GPU CUDA backend implementation. Defines the CUDA kernel entry points, GPU-side memory layout, and CUDA-specific helper utilities.

## Build and Configuration

Vina scoring does not require any additional build flags beyond the standard muDock configuration options, so you can refer to the main [muDock README](../README.md) for full instructions.

For convenience, here are the build commands for the **CUDA backend** (targeting e.g. SM80):

```bash
cmake -S /path/to/muDock -B /path/to/muDock/build \
  -DMUDOCK_ENABLE_CUDA=ON \
  -DMUDOCK_GPU_ARCHITECTURES=nvidia:sm_80 \
  -DCMAKE_BUILD_TYPE=Release
cmake --build /path/to/muDock/build
```

## How to run

To execute docking or scoring using the Vina energy model, specify `--score vina` when running the `muDock` application. You can control whether to perform a full docking search or a static evaluation using the `--search` flag:

*   `--search genetic`: Runs the full Genetic Algorithm search using Vina scoring.
*   `--search none`: Performs a single-point energy evaluation (scoring only) on the input ligand structures without search optimization.

### GPU (CUDA) Execution

```bash
./build/application/muDock \
  --protein <path_to_protein> \
  --ligand <path_to_ligand_dataset> \
  --use CUDA:GPU:0[:WORKERS][:DEVICE_MEMORY_BYTES] \
  --score vina \
  --search genetic
```

Some proteins can be found inside the `data` folder (es. data/1fkb/1fkb_protein.pdb)


## Performance Evaluation

To evaluate and benchmark the performance of the implemented Vina scoring engine, all tests were conducted using the **1FKB** protein target (`data/1fkb/1fkb_protein.pdbqt`). The input workload consisted of a dataset of **10,000 ligands** generated with an automated dataset pipeline.

To ensure an equitable distribution of workload and optimize resource utilization, ligands were grouped into **size-based buckets**, ensuring that molecules of similar atom counts and structural complexity were processed together within the same execution batch.

All profiling metrics and hardware counters were extracted using **NVIDIA Nsight Compute** with full metric collection enabled:

```bash
ncu --set full
```

Kernel execution parameters evolved across iterations:

- Early Implementations: Used warp-level optimizations, utilizing a block size of 32 threads (1 warp) and a bucket multiplier of 18.

- Intermediate Configurations: Where feasible, CUDA kernels were scaled to 256 threads per block along with a bucket multiplier of 36. This aggressive batching strategy pushed GPU VRAM utilization near saturation, particularly during active profiling with ncu-ui.

- Final Configurations: In the latest versions of muDock, the batch/bucket size is dynamically derived from the device memory limit passed as a CLI argument. For these versions, the kernels were kept at 256 threads per block. The `--use` flag was configured as follows:


```bash
--use CUDA:GPU:0:2:3000000000
```

All benchmark binaries were compiled using the following key optimization CMake flags:

```cmake
CMAKE_BUILD_TYPE            Release
MUDOCK_ENABLE_FAST          ON
MUDOCK_ATOM_CLUSTER_LEVEL   LARGE
# In older versions
MUDOCK_ENABLE_BUCKET        ON 
```
Each implementation was evaluated under these standardized build flags and dataset partitioning conditions to guarantee fair, highly reproducible throughput and execution time comparisons across hardware backends.

---

### Evolutionary Timeline & Optimization Journey

The optimization of the CUDA Vina kernel was an iterative process that evolved through several key architectural milestones. The main progression of the kernel's development can be traced through the following key commits:

1. **`e06a2a9` - Naive Baseline**: Integration of the initial naive CUDA kernel, into the newly updated muDock architecture.
2. **`3854c97` - Memory Footprint Reduction**: Eliminated unnecessary intermediate buffers to reduce memory overhead and latency.
3. **`0063f11` - Shared Memory & Grid Parallelism**: Moved ligand coordinates to shared memory (initial approach), replaced standard library math calls with hardware-accelerated intrinsics (e.g., `__expf`), and flipped the `inter_score` parallel execution model to iterate over receptor atoms instead of ligand limits.
4. **`3ebe289` - Occupancy & Shared Memory Tuning**: Increased thread block size to 256 with a bucket multiplier of 36 (tuned for hardware constraints, e.g., NVIDIA RTX 4050 6GB VRAM). Loaded full ligand structure into shared memory and introduced intrinsic micro optimizations (`fmaxf`, `fminf`).
5. **`412dd5c` - Shared Memory Layout Refactoring**: Consolidated 7 independent shared memory arrays into 2 compact `float4` arrays to mitigate bank conflicts and memory stalls. Applied `#pragma unroll` on ligand iterations within `inter_score` to hide memory access latencies.
6. **`current` - Architectural Integration & FMA**: Final state merged with the latest `dev` branch, taking advantage of updated muDock core pipeline features and leveraging `fmaf` intrinsics for fused multiply-add execution in a single clock cycle.

---

### Note on Development Complexity & Branch Merges

*(Context for code review)*

While the list above highlights the main performance milestones, the actual development history encompasses significantly more complexity than these discrete checkpoints suggest. 

The implementation was developed in parallel with major ongoing structural updates to the main muDock framework. Maintaining compatibility required frequent, complex merges with the `dev` branch. This constant evolution introduced several challenges:

*   **File Format & Parser Extensions**: To enable parallel parsing for tens of thousands of ligands, the native `adtmol2` format, converter, and parser had to be substantially extended to encode Vina-specific atomic metadata (van der Waals radii, hydrogen bond donors/acceptors, hydrophobicity, topological neighbor masks, etc.).
*   **External Dependencies & Build Breaking Changes**: Updates to third-party libraries, such as OpenBabel requiring Eigen3 as an explicit dependency in recent builds, mean that earlier commits may fail to compile directly without manual adjustments or backporting dependency fixes.
*   **Architectural Adapters & Non-Linear Progress**: Structural changes in the core library repeatedly invalidated prior assumptions. Optimizations that proved effective at one stage sometimes had to be redesigned or removed as the underlying execution pipeline evolved. Furthermore, edge-case bugs discovered later in development (originating both within the Vina scoring logic and deep inside core muDock components) make earlier intermediate commits unusable on certain input structures.

