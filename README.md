![muDock_icon](share/icon_200_186.png)

# muDock — Molecular Docking Microapp

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
- Boost (components: `program_options`, `graph`, `fiber`)
- OpenBabel3

Fetched automatically:
- Microsoft GSL (via CMake FetchContent)

Optional (enabled via build flags):
- OpenMP (CPU parallelism)
- CUDA Toolkit (with `curand`) for CUDA backend
- HIP + hiprand for HIP backend
- SYCL toolchain (oneAPI/LLVM)
- Google Highway (`HWY`) and/or `xsimd` for CPU vectorization
- LIKWID for profiling

## Build and configuration

The project uses CMake and builds like a standard CMake package.

Basic build:

```bash
cmake -S /path/to/muDock -B /path/to/muDock/build -DCMAKE_INSTALL_PREFIX=/install/path
cmake --build /path/to/muDock/build
```

Key configuration options:

- `CMAKE_BUILD_TYPE` — `Release` (default), `Debug`, or `RelWithDebInfo`
- `MUDOCK_ENABLE_FAST` — enables aggressive optimizations (auto-ON for Release)
- `MUDOCK_ENABLE_OMP` — OpenMP CPU parallelism
- `MUDOCK_ENABLE_CUDA` — CUDA backend
- `MUDOCK_ENABLE_HIP` — HIP backend
- `MUDOCK_ENABLE_SYCL` — SYCL backend
- `MUDOCK_ENABLE_GH` — Google Highway vectorization
- `MUDOCK_ENABLE_XSIMD` — xsimd vectorization
- `MUDOCK_ENABLE_LIKWID` — LIKWID profiling
- `MUDOCK_ENABLE_TEST` — enable tests (not allowed in Release)

GPU/accelerator target configuration:

- `MUDOCK_GPU_ARCHITECTURES` — format `platform:arch`
  - Examples: `cuda:sm_80`, `amd:gfx90a`, `intel:gen12`
- `MUDOCK_CPU_ARCHITECTURES`, `MUDOCK_CPU_TARGET`, `MUDOCK_CPU_TUNE` — fine-tune CPU codegen

Example: CUDA build targeting SM80:

```bash
cmake -S /path/to/muDock -B /path/to/muDock/build \
  -DMUDOCK_ENABLE_CUDA=ON \
  -DMUDOCK_GPU_ARCHITECTURES=cuda:sm_80 \
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
