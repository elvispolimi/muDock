![muDock_icon](share/icon_200_186.png)
# µDock - Molecular Docking microapp

Molecular docking is a well-known problem in scientific literature.
It aims to estimate the 3D displacement of a molecule with a small molecular weight, named *ligand*, when it interacts with the target protein.
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

TBD

## Repository structure

The repository is structured as follows:

- `application`: this folder contains the source files of the µDock executable.

- `share`: this folder contains graphical resources and other files that are not strictly related to the application.
