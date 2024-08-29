![muDock_icon](share/icon_200_186.png)

# µDock - Molecular Docking microapp

Molecular docking is a well-known problem in scientific literature.
It aims to estimate the 3D displacement of a molecule with a small molecular weight, named _ligand_, when it interacts with the target protein.
The main issue is the number of degree of freedom involved in the operation.
Due to the difference in size between the protein and ligand, we need to find the correct position and orientation of the ligand inside a pocket cavity.
Moreover, the molecules are not rigid bodies, but we can rotate a subset of their bonds, changing their geometry, without changing their chemical and physical properties.
Thus, we have six degrees of freedom for rigid rototranslations in a 3D space and one additional degree of freedom for each of those bonds.
Arguably, the most popular tool that performs molecular docking is autodock ([https://autodock.scripps.edu/](https://autodock.scripps.edu/)).
Its implementation is composed of a feature rich ecosystem of tools, that might employ different optimization strategies to solve the exploration problem.

This repository provides a nimble Autodock-like implementation that uses a genetic algorithm as an optimization strategy and the energy model of Autodock 4.0.
The goal is to provide a microapp to use as a prototype for applying a wide range of optimizations, from kernel porting to a different language, to the application of approximate computation techniques.

## How to build from sources

We use CMake as a building system with a standard compilation procedure.
This example will compile and install µDock in the `/install/path` path, assuming that `/path/to/muDock` is the path of this repository root:

```bash
$ cmake -S /path/to/muDock -B /path/to/muDock/build -DCMAKE_INSTALL_PREFIX=/install/path
$ cmake --build /path/to/muDock/build
$ cmake --install /path/to/muDock/build
```

## How to use the application

This application reads from the standard input a ligand library in mol2 format. It will print on the standard output the best docked pose with the score as comment

USAGE: ./build/application/muDock --protein "protein.pdb" --use CPP:CPU:0 < "/path/to/ligands.mol2"

Available options:
--help print this help message
--protein arg (="protein.pdb") Path to the protein file (in PDB)
--use arg (=CPP:CPU:0) Devices configuration

### Devices Configuration

The device configuration options follows the following structure:

```
KERNEL:DEVICE:IDS
```

Where:

- KERNEL is one of the following
  - **CPP** which is the multi-threaded version of muDock
- IDS is like the printer dialogue, for example
  - **1,3,4,7-9,12** means that you are targeting devices with ids 1,3,4,7,8,9,12
  - **0-3,5,8-10** means that you are targeting devices with ids 0,1,2,3,5,8,9,10
- Device is one of the following
  - **CPU** use the CPU of the target machine, each IDS is a thread with the specify affinity

The full list of supported devices can be found in [here](./mudock/include/mudock/devices/device_types.hpp), while the kernel types list in [here](./mudock/include/mudock/devices/kernel_types.hpp)

## Repository structure

The repository is structured as follows:

- `application`: contains the source files of the µDock executable.

- `chem`: contains chemical knowledge in JSON format, typically used to automatically generate sources

- `cmake`: contains cmake helper function, e.g. to download third party dependencies

- `mudock`: contains the sources of the µDock library, which implements domain concerns

- `script`: contains utility scripts, e.g. beautifier or the script that generates sources

- `share`: contains graphical resources and other files that are not strictly related to the application.
