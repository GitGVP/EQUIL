# EQUIL
EQUIL is an equilibrium code based on the geometrical inverse formulation: input profiles are assumed to be given as functions of a normalized flux coordinate $\hat r$ and the code finds the plasma diamagnetism and the geometry of the flux surfaces. For the anisotropic case, the magnetic field strength $B$ is also solved for self-consistently. The boundary conditions are imposed such that the aspect ratio $\epsilon:=a/R_0$ is fixed a priori, so the horizontal shift is obtained self-consistently by the code. The rest of the (up-down symmetric) shaping harmonics have Dirichlet outer BCs (fixed boundary).


## Using EQUIL

Start Matlab from `examples/` add to the path `src/` and sub-folders and then run a few examples to understand how to modify the inputs. The organization of the code is taken from LIUQE (MEQ).

The starting point should simply be to run `examples/simple_case.m` after having added the paths.

## Discretization
The Newton loop can be found in `src/equilY` which calls the different residuals and Jacobians. The most up to date pair is `residuals_noRepl` and `jacobian_noRepl` but for numerical convergence purposes I will provide a pair for which a full solution is available.
The discretization part is in `src/utils` and consists in `assemble_FE_matrices_bspline` the main routine constructing the discretization matrices. This function can be replaced by other discretizations (different FE, spectral method etc.) and the outputs are:
- Quadrature points $r_q$ (`r_q`), size $N_q$, where the residuals and Jacobian will be evaluated 
- Evaluation matrices $P_0,P_1,P_2$ (`P0_full`,`P1_full`,`P2_full`), size $N_q\times N_{\text{d.o.f.}}$, allows to go from the state vector $x$ to the profiles and their derivatives' values on the quadrature points $r_q$ (function, first derivative and second derivative respectively)
- Residual projection matrix $M$ (`M_profiles`), size $N_{\text{d.o.f.,tot}}\times N_{\text{prof}} N_q$, from the residual evaluated on $r_q$ to the residuals in terms of degrees of freedom, schematically $M_{i,q} = \Lambda_i w(r_q)$ where $w$ is the weight to evaluate the integral. This is the same matrix stacked $N_{\text{prof}}$ times for each residual.
- Jacobian assembly matrices (`M_extended`, `P_templates`, `P_extended`), matrices that are built from $M$ and $P_i$ above but with some operations to adjust their sizes etc. for the jacobain assembly. See jacobian equation in the article for more explanations, and also the last parts of `jacobian_noREpl.m` in `src/residuals_jacobians/` to understand how these are used, note that in the notation of the paper `jacTotal` is $\mathcal J_q^{\alpha,\gamma,(d)}$.


## Terms of Use

EQUIL is a scientific code under MIT license. If you use this code in a publication, please cite the associated paper:
Van Parys et al., "Investigation of finite aspect ratio effects in axisymmetric magnetic equilibria with toroidal rotation and pressure anisotropy" (submitted to PPCF).