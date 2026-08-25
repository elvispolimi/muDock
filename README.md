![muDock_icon](share/icon_200_186.png)

# muDock — Molecular Docking Microapp

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19384509.svg)](https://doi.org/10.5281/zenodo.19384509)

muDock is a compact, Autodock-style docking engine that uses a genetic algorithm and the Autodock 4.0 energy model. It was born as a benchmarking tool for Autodock-style workflows, and today it remains a small, focused codebase for experimenting with performance techniques (kernel porting, vectorization, accelerator backends, and approximation strategies) while staying usable as a docking tool.

The pipeline is intentionally split into clean stages (input parsing, scoring, search, and output), so individual pieces can be swapped or extended. That structure makes it practical to prototype new scoring functions, docking algorithms, or search strategies without rewriting the rest of the system.

## Repository layout

- `application` — CLI entry point and executable sources
- `mudock` — core library and domain logic
- `chem` — chemical knowledge (JSON) used for code generation
- `cmake` — CMake helpers and dependency setup
- `script` — utility scripts (formatting, code generation, etc.)
- `share` — icons and non-code assets
- `test` — tests (disabled by default)

## Dependencies

Required:
- CMake 3.25+
- A C++20-capable compiler (GCC/Clang/IntelLLVM)
- Boost (components: `program_options`, `graph`, `context`)
- oneTBB
- OpenBabel3

Required when using the bundled CMake presets:
- Ninja

Optional (enabled via build flags):
- MPI
- OpenMP (CPU parallelism)
- CUDA Toolkit (with `curand`) for CUDA backend
- HIP + hiprand for HIP backend (`rocrand` is also needed on AMD platforms)
- SYCL toolchain based on Intel LLVM / oneAPI
- Google Highway (`HWY`) and/or `xsimd` for CPU vectorization
- LIKWID for profiling

## Build and configuration

The project uses CMake and builds like a standard CMake package.

Basic build:

```bash
cmake -S /path/to/muDock -B /path/to/muDock/build -DCMAKE_INSTALL_PREFIX=/install/path
cmake --build /path/to/muDock/build
```

Convenience presets for CPU-only development and CI are also provided:

```bash
cmake --preset dev-cpu-debug
cmake --build --preset dev-cpu-debug
ctest --preset dev-cpu-debug
```

These presets set `"generator": "Ninja"` in `CMakePresets.json`, so install Ninja before using them. If you prefer another generator, use the manual `cmake -S ... -B ...` flow above instead.

Most useful configuration options:

- `CMAKE_BUILD_TYPE` — `Release` (default), `Debug`, or `RelWithDebInfo`
- `MUDOCK_ENABLE_TEST` — enable tests (not allowed in `Release`)
- `MUDOCK_GPU_ARCHITECTURES` — target accelerator in the form `platform:arch` (for example `nvidia:sm_86`, `amd:gfx90a`, `intel:gen12`)
- `MUDOCK_ENABLE_OMP` — OpenMP CPU parallelism
- `MUDOCK_ENABLE_MPI` — MPI frontend
- `MUDOCK_ENABLE_CUDA` — CUDA backend
- `MUDOCK_ENABLE_HIP` — HIP backend
- `MUDOCK_ENABLE_SYCL` — SYCL backend
- `MUDOCK_ENABLE_GH` — Google Highway vectorization
- `MUDOCK_ENABLE_XSIMD` — xsimd vectorization
- `MUDOCK_SYCL_WG_SIZE` — SYCL work-group size override
- `MUDOCK_ATOM_CLUSTER_LEVEL` — ligand atom-cluster granularity: `OFF`, `MEDIUM`, `LARGE`, `EXTREME`
- `MUDOCK_STAGE_BUCKET_POLICY` — stage bucket sizing policy: `DEVICE_ALIGNED`, `SM_ALIGNED`, `MAX_UTILIZATION`
- `MUDOCK_STAGE_BUCKET_OVERRIDE`, `MUDOCK_STAGE_BUCKET_MULTIPLE_OVERRIDE` — explicit stage bucket overrides
- `MUDOCK_DISABLE_UNROLL`, `MUDOCK_UNROLL_FACTOR` — control kernel loop unrolling
- `MUDOCK_ENABLE_STAGE_BUCKET_TRACE` — print stage bucket decisions at runtime
- `MUDOCK_ENABLE_LIKWID` — LIKWID profiling
- `MUDOCK_CPU_ARCHITECTURES`, `MUDOCK_CPU_TARGET`, `MUDOCK_CPU_TUNE` — fine-tune CPU code generation

Example: CUDA build targeting SM80:

```bash
cmake -S /path/to/muDock -B /path/to/muDock/build \
  -DMUDOCK_ENABLE_CUDA=ON \
  -DMUDOCK_GPU_ARCHITECTURES=nvidia:sm_80 \
  -DCMAKE_BUILD_TYPE=Release
cmake --build /path/to/muDock/build
```

Example: OpenMP CPU build:

```bash
cmake -S /path/to/muDock -B /path/to/muDock/build \
  -DMUDOCK_ENABLE_OMP=ON \
  -DCMAKE_BUILD_TYPE=Release
cmake --build /path/to/muDock/build
```

## Running

Main application (default uses AutoDock 4.0 scoring and Genetic Algorithm search):

```bash
./build/application/muDock --protein /path/to/protein.pdb --ligand /path/to/ligands.mol2 --use CPP:CPU:0
```

To use the X-Score scoring function instead (supported for scoring only, no search):

```bash
./build/application/muDock --protein /path/to/protein.pdb --ligand /path/to/ligands.mol2 --use CPP:CPU:0 --search none --score xscore
```

The `--use` flag maps implementations to devices using:

`IMPLEMENTATION:DEVICE:IDS[:WORKERS][:MEMORY_BYTES]`

For example:

- `CPP:CPU:0`
- `CUDA:GPU:0`
- `SYCL:GPU:0-1:2:1000000000`

For the TBB stream frontend, the most useful extra runtime flags are:

- `--tokens` — maximum number of in-flight TBB tokens
- `--bytes_per_token` — chunk size used by the stream reader
- `--queue_size` — bounded queue size between parsing and screening
- `--observer` — print periodic throughput and backlog information
- `--time_limit_sec` — stop admitting new work after a time limit and drain what is already in flight

Converter:

```bash
./build/application/converter --input input.mol2 --output output.pdbqt
```

Supported formats (by file extension):
- `mol2`
- `pdbqt`
- `pdb`
- `adtmol2`

Note on ligand parsing: `adtmol2` and `mol2` now use muDock's native parser and can be split and parsed in parallel. `adtmol2` still carries extra AutoDock-specific fields, while plain `mol2` remains a generic interchange format.

`mol2` writing is also available natively. The TBB stream frontend and the MPI frontend accept both `adtmol2` and `mol2` streams.

## Tests

Enable tests at configure time (Release is not allowed):

```bash
cmake -S /path/to/muDock -B /path/to/muDock/build -DMUDOCK_ENABLE_TEST=ON -DCMAKE_BUILD_TYPE=Debug
cmake --build /path/to/muDock/build
ctest --test-dir /path/to/muDock/build
```

The repository also ships CMake presets tailored for CPU-only testing:

- `dev-cpu-debug` for local development
- `ci-cpu-debug` for automation and reproducible CI runs

## CI

GitHub Actions CPU CI is configured in `.github/workflows/cpu-ci.yml`.

- Triggered on every `push` and `pull_request` (plus manual `workflow_dispatch`)
- Installs the required CPU-side dependencies on `ubuntu-24.04`
- Configures with the `ci-cpu-debug` preset
- Builds the project and runs `ctest --preset ci-cpu-debug`

The workflow intentionally uses a `Debug` build because muDock rejects `MUDOCK_ENABLE_TEST=ON` in `Release`.

## Devcontainer

A VS Code Dev Containers setup is provided in `.devcontainer/`.

To use it:

1. Open the repository in VS Code.
2. Run `Dev Containers: Reopen in Container`.
3. After the container is created, the project is configured automatically with `cmake --preset dev-cpu-debug`.

Common commands inside the devcontainer:

```bash
cmake --build --preset dev-cpu-debug
ctest --preset dev-cpu-debug
```

The devcontainer installs the same core CPU dependencies used by CI, plus common development tools such as `clangd`, `clang-format`, `gdb`, and `ripgrep`.

## Apptainer / Singularity

An Apptainer recipe is provided at `apptainer/mudock.def` for reproducible CPU-focused environments.

Build the image:

```bash
apptainer build mudock-dev.sif apptainer/mudock.def
```

If your system requires it, use `--fakeroot` or your site-specific Apptainer build flow instead.

Run an interactive shell with the repository mounted:

```bash
apptainer shell --bind "$PWD":/workspace mudock-dev.sif
cd /workspace
cmake --preset dev-cpu-debug
cmake --build --preset dev-cpu-debug
ctest --preset dev-cpu-debug
```

## Troubleshooting

- If CMake cannot find Boost through `BoostConfig.cmake`, muDock now falls back to CMake's standard `FindBoost` module automatically.
- If `ctest` reports no tests, confirm you configured with `MUDOCK_ENABLE_TEST=ON` or used the provided CPU debug presets.
- If OpenBabel is found at configure time but fails at runtime, verify that the runtime package is installed in addition to the development headers.
- If Apptainer image builds are blocked by permissions, retry with `apptainer build --fakeroot ...` or follow your cluster's container policy.

## References

- Gianmarco Accordi, Jens Domke, Theresa Pollinger, Davide Gadioli, Gianluca Palermo. "Towards High-Performance and Portable Molecular Docking on CPUs Through Vectorization." IEEE Cluster 2025. DOI: https://doi.org/10.1109/CLUSTER59342.2025.11186493
