# EQUIL
EQUIL is an equilibrium code based on the geometrical inverse formulation: input profiles are assumed to be given as functions of a normalized flux coordinate $\hat r$ and the code finds the plasma diamagnetism and the geometry of the flux surfaces. For the anisotropic case, the magnetic field strength $B$ is also solved for self-consistently. The boundary conditions are imposed such that the aspect ratio $\epsilon:=a/R_0$ is fixed a priori, so the horizontal shift is obtained self-consistently by the code. The rest of the (up-down symmetric) shaping harmonics have Dirichlet outer BCs (fixed boundary).


## Using EQUIL

Start Matlab from `examples/` add to the path `src/` and sub-folders and then run a few examples to understand how to modify the inputs. The organization of the code is taken from LIUQE (MEQ).

## Left to do
- [ ] Push to even higher $\beta_p$, and higher anisotropy. This may require to modify the discretization (from B-splines to spectral method ?) as it has been identified that when the code fails, the near-axis deviations from Dirichlet/Neumann BCs are to blame and not penalized enough with the current discretization.
- [ ] Make the residuals coming from the vectorial force balance as fast as the ones from $\hat G$. They are identical for isotropic and anisotropic plasmas, but the rotation adds some terms which make the current $\hat G$ formulation break.



