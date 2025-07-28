# Feature-Aligned Parametrization in Penner Coordinates

<strong>Ryan Capouellez<sup>1</sup>, Rodrigo Singh<sup>1</sup>, Martin Heistermann<sup>2</sup>, David Bommes<sup>2</sup>, Denis Zorin<sup>1</sup></strong>

<small><sup>1</sup>New York University, <sup>2</sup>University of Bern</small>

An implementation of [Feature-Aligned Parametrization in Penner Coordinates](https://dl.acm.org/doi/10.1145/3731216).

![Challenging parametrizations](media/teaser.jpg)

### Overview

This method generates an approximately isometric seamless parameterization of an input `obj` mesh that is aligned to a provided feature set. TODO with parametric cone angle and holonomy signature constraints. Retriangulation is often necessary to satisfy these constraints, so the initial mesh is intrinsically refined to produce an output mesh with a compatible parameterization.

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
|`--mesh` | Mesh filepath| `none`|
|`--cones` | Target cone filepath| `none`|
|`--field` | Cross field rotation form filepath| `none`|
|`--output` | Output directory| `./`|
|`--remove_holonomy_constraints` | Only constrain cone angles (parametrization no longer seamless)| `false`|
|`--max_itr` | Maximum number of iterations| `500`| 
|`--error_eps` | Maximum allowed angle error| `1e-12`|
|`--max_triangle_quality` | Regularize initial metric until quality measure below value (`0` for no regularization) | `0`|
|`--use_delaunay` | Perform Newton in Delaunay connectivity | `true`|
|`--fit_field` | Fit new cross field instead of using cones and field from file| `false`|

The input mesh must be a manifold surface with a single connected component. The input cone and field files can be generated with `bin/generate_frame_field`, or they can generated algorithmically at runtime with the `--fit_field` flag. The output is a refined mesh with a parameterization and a file of metric coordinate values.

## Figure Reproduction

Scripts to generate the figures of the paper are included in `figures`.

![Some example figures](media/examples.jpg)

The models (with cone and field data) and cameras used in [Seamless Parametrization in Penner Coordinates](https://dl.acm.org/doi/10.1145/3658202) necessary for these scripts will be available soon;  `closed-myles`, `thingi10k-tetwild` and `cameras` must be copied to `data/closed-Myles`, `data/thingi10k-tetwild`, and `data/cameras` respectively.

A Conda environment must be activated (before compiling the code) with
```
conda env create -f environment.yml
conda activate curvature-metric
```
The figure bash scripts can then be run independently or in batch with
```
bash fig-all.sh
```

Note that most bash scripts generate an output directory with a JSON file specifying parameters for the parameterization and rendering pipeline python script `scripts/_pipeline.py`. Such JSON files can also be used for general batch parameterization and analysis.

### Library

Penner coordinates are global coordinates on the space of metrics on meshes with a fixed vertex set and topology, but varying connectivity, making it homeomorphic to the Euclidean space of dimension equal to the number of edges in the mesh, without any additional constraints imposed.

These coordinates underly the recent advances in parametrization with cone and holonomy angle constraints. To engender future work in this direction, we provide an independent library containing the data structures and methods for Penner coordinates at [geometryprocessing/penner-optimization](https://github.com/geometryprocessing/penner-optimization).

## Citation

```
@article{capouellez:2023:seamless,
author = {Capouellez, Ryan and Zorin, Denis},
title = {Seamless Parametrization in Penner Coordinates},
year = {2024},
issue_date = {July 2024},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
volume = {43},
number = {4},
issn = {0730-0301},
url = {https://doi.org/10.1145/3658202},
doi = {10.1145/3658202},
journal = {ACM Trans. Graph.},
month = {jul},
articleno = {61},
numpages = {13},
keywords = {parametrization, seamless, discrete metrics, cone metrics, conformal, intrinsic triangulation, penner coordinates}
}
```
