# Feature-Aligned Parametrization in Penner Coordinates

<strong>Ryan Capouellez<sup>1</sup>, Rodrigo Singh<sup>1</sup>, Martin Heistermann<sup>2</sup>, David Bommes<sup>2</sup>, Denis Zorin<sup>1</sup></strong>

<small><sup>1</sup>New York University, <sup>2</sup>University of Bern</small>

An implementation of [Feature-Aligned Parametrization in Penner Coordinates](https://dl.acm.org/doi/10.1145/3731216).

![Challenging parametrizations](media/teaser.jpg)

### Overview

This method generates an approximately isometric seamless parameterization of an input `obj` mesh that is aligned to a provided feature set. Retriangulation is often necessary to satisfy these constraints, so the initial mesh is intrinsically refined to produce an output mesh with a compatible parameterization.

## Installation

To install this project on a Unix-based system, use the following standard CMake build procedure:

```bash
git clone --recurse-submodules https://github.com/rjc8237/feature-aligned-penner.git
cd feature-aligned-penner
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j 4
```

## Usage

The core parameterization method is `bin/parametrize_aligned`. This executable takes the following arguments:

|flag | description| default|
| --- | --- | --- |
|`--name` | name of the mesh (without `.obj` suffix) | `none`|
|`--input` | input directory with mesh | `./`|
|`--output` | output directory for parameterized mesh | `./`|
|`--show_parameterization` | open viewer to see parameterization | `false`|

The input mesh must be at the input path `<input>/<name>.obj`, and it must be a manifold surface with a single connected component.

### Library

Penner coordinates are global coordinates on the space of metrics on meshes with a fixed vertex set and topology, but varying connectivity, making it homeomorphic to the Euclidean space of dimension equal to the number of edges in the mesh, without any additional constraints imposed.

These coordinates underly the recent advances in parametrization with cone and holonomy angle constraints. To engender future work in this direction, we provide an independent library containing the data structures and methods for Penner coordinates at [geometryprocessing/penner-optimization](https://github.com/geometryprocessing/penner-optimization).

## Citation

```
@article{capouellez:2024:feature,
author = {Capouellez, Ryan and Singh, Rodrigo and Heistermann, Martin and Bommes, David and Zorin, Denis},
title = {Feature-Aligned Parametrization in Penner Coordinates},
year = {2025},
issue_date = {August 2025},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
volume = {44},
number = {4},
issn = {0730-0301},
url = {https://doi.org/10.1145/3731216},
doi = {10.1145/3731216},
journal = {ACM Trans. Graph.},
month = jul,
articleno = {110},
numpages = {21},
}
```
