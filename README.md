# EQUIL

This repository contains the code used in

**G. Van Parys, J. P. Graves, A. Merle, C. Heiss, M. Bassanini, and T. Glauser,**  
*Investigation of finite aspect ratio effects in axisymmetric magnetic equilibria with toroidal rotation and pressure anisotropy*.

This branch is a frozen snapshot associated with the paper. The scripts in `ArticleResults/` reproduce the figures presented in the manuscript.

## Repository contents

- `src/`: core equilibrium solver and supporting routines
- `ArticleResults/`: scripts used to reproduce the figures of the paper
- `Mathematica/`: symbolic derivations and auxiliary notebook material

## Reproducing the paper figures

The figure scripts are self-contained and can be run directly from MATLAB after adding the repository to the path.

A typical session is:

```matlab
addpath(genpath('src'))
addpath(genpath('ArticleResults'))
```
Then run the desired script from `ArticleResults/`.

The main reproduction scripts are:

- `ArticleResults/numerical_convergence.m`: Fig. 1
- `ArticleResults/Solovev.m`: Fig. 2
- `ArticleResults/isotropic_1.m`: Figs. 3–5
- `ArticleResults/parallel_rotating_1.m`: Figs. 6–7
- `ArticleResults/perp_anisotropic.m`: Figs. 8–10

Each script contains a short header indicating which section and figures it reproduces.

## Citation

If you use this code, please cite the associated paper above.

## License

This code is released under the MIT License. See LICENSE.