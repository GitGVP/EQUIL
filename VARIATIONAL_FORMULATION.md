# Variational EQUIL formulation
This document describes the normalized formulation implemented by `equilVariationalSol`. It is separate from the generated strong-form residuals used by `equilSol`. The main differences with the legacy code are: the magnetic field strength is solved for locally and not part of the state vector, the representation is generalized to the up-down asymmetric case, and the gauge-enforcing profile is adapted for the anisotropic case.
## Action, flux and geometrical variations
The section is organized as follows: we first introduce the coordinate-independent action, vary it with respect to the poloidal flux $\psi$ and the magnetic field strength $B$, and then pull back to general flux coordinates and applying at that point a geometrical variation. The specific representation and unknowns solved by the code are specified in the next section
The action from which the weak residuals are derived is:
\[
\mathcal A[\psi,B]=\int \mathrm{d}V\left[\frac{B^2}{2\mu_0}-P_\parallel(\psi,R,B)-\frac{T(\psi)}{\mu_0R}B_\phi(\nabla\psi,R,B)\right].
\]
with:
\[
  B_p^2[\nabla\psi, R]=\frac{|\nabla\psi|^2}{R^2},\qquad B_\phi[\nabla\psi, R, B]=\sqrt{B^2-B_p^2[\nabla\psi,R]}.
\]
where we put the explicit dependences of $B_\phi(\nabla\psi,R,B)$ as these will provide the kinetic term for the flux variation and contribute to the geometrical variations defined next. The anisotropy parameter is given by (subscripts indicate partial derivatives):
\[
\sigma_\parallel=\mu_0\frac{P_\parallel-P_\perp}{B^2}=\frac{\mu_0 P_{\parallel,B}}{B},
\]
and varying with respect to $B$ directly gives the non-linear algebraic equation for the local magnetic field strength ($\mathcal L$ is the integrand of the action):
\[
\left(\frac{\partial\mathcal L}{\partial B}\right)=\frac{B}{\mu_0}-P_{\parallel,B}-\frac{T}{\mu_0R}\frac{B}{B_\phi}=0\Rightarrow (1-\sigma_\parallel)RB_\phi=T(\psi)
\]
while varying with respect to $\psi \to \psi + \Phi$ at constant $B$ gives (after using the $B$ equation above):
\[
 \mathcal R_\psi[\Phi]= \int \mathrm{d}V\left[\frac{1-\sigma_\parallel}{\mu_0 R^2}\nabla\psi \cdot\nabla\Phi - \left(P_{\parallel,\psi}+\frac{TT_{,\psi}}{\mu_0R^2(1-\sigma_\parallel)}\right)\Phi\right]
\]
which is the weak form of the residuals used in e.g. the anisotropic-rotating version of the MEQ code. We now introduce a choice of independent coordinates $(r,\omega)$ for which $\psi=\psi(r)$ and importantly, use the pulled-back pressure and $T$ maps:
\[
 P_\parallel(r, R, B) = \tilde P_\parallel(\psi(r), R, B),\qquad T(r) = \tilde T(\psi(r)).
\]
and we drop the tildes next, their derivatives acquire a factor $T_{,\psi}= T_r \psi_r^{-1}$ so we are ready to give the weak residuals associated with the $\psi$ variation:
\[
 \mathcal R_\psi[\Phi]= 2\pi\iint \mathrm{d}r \mathrm d \omega\frac{J}{\psi_r}\left[(1-\sigma_\parallel)\frac{B_p^2}{\mu_0} \Phi_r- \left(P_{\parallel,r}+\frac{TT_{r}}{\mu_0R^2(1-\sigma_\parallel)}\right)\Phi\right].
\]
where $2\pi$ comes from the toroidal integration, and in these coordinates:
\[
 B_p^2 = \frac{g_{\omega\omega}}{J^2}\psi_r^2
\] 
This step is important as it introduces the $\psi_r^{-1}$ singularities and reminds us that the maps $T(r)$ and $P_\parallel(r,R,B)$ may not have arbitrary dependence in $r$. After specifying the functional form of $\psi(r)$ it will be clear that diverted cases for which $q\to\infty$ cannot be represented. Finally, geometrical variations give;
\[
\begin{aligned}
\delta_\mathrm{g} J={}&(R_rZ_\omega-R_\omega Z_r)\delta R+R\big(
Z_\omega\delta R_r+R_r\delta Z_\omega
-Z_r\delta R_\omega-R_\omega\delta Z_r\big),\\
\delta_\mathrm{g} g_{\omega\omega}={}&2R_\omega\delta R_\omega
+2Z_\omega\delta Z_\omega,\\
\delta_\mathrm{g} B_p^2&=\psi_r^2
\left(\frac{\delta_\mathrm{g} g_{\omega\omega}}{J^2}
-\frac{2g_{\omega\omega}\,\delta_\mathrm{g} J}{J^3}\right),\qquad\delta_\mathrm{g} B_\phi=-\frac{\delta_\mathrm{g} B_p^2}{2B_\phi}.
\end{aligned}
\]
Where $\delta_\mathrm{g}$ distinguishes between geometrical and flux variation, i.e. $\delta_\mathrm{g} B_p \neq \delta B_p$. Noting that $\mathrm d V = J \mathrm d r \mathrm d \omega$ we require:
\[
\delta_\mathrm{g}(J\mathcal L)=
\mathcal L\,\delta_\mathrm{g} J+J\left(-P_{\parallel,R}+\frac{TB_\phi}{\mu_0R^2}\right)\delta_\mathrm{g} R
-\frac{JT}{\mu_0R}\delta_\mathrm{g} B_\phi.
\]
to get the equations determining the functions parametrizing the flux surfaces.  We can compute the safety factor:
\[
q = \int_0^{2\pi} \mathrm d\omega \frac{\bm B\cdot\nabla\phi}{\bm B\cdot\nabla\omega} = \frac{T}{\psi_r}\int_0^{2\pi} \mathrm d\omega \frac{J}{R^2(1-\sigma_\parallel)}.
\]
Define the prescribed cylindrical-reference anisotropy and axis normalization by
\[
\sigma_{\rm cyl}(r)=\left.\sigma_\parallel(r,B)\right|_{B=1},
\qquad a_0=1-\sigma_{\rm cyl}(0).
\]
The poloidal flux is then given by
\[
\psi_r=\frac{r}{R_0}\frac{T}{q}\frac{1}{1-\sigma_{\rm cyl}(r)}
\]
solving $\mathcal R_\psi[\Phi]=0$ for the diamagnetism profile $T(r)$. Note that the $\psi_r^{-1}$ divergences now correspond to $q\to \infty$ (diverted flux surface) or $T=0$ i.e. toroidal field vanishing. In the next section we introduce our specific representation and normalized quantities.

## Up-down asymmetric $(\Delta, S_m, A_m, P)$ representation
We do not explicit the normalization which is the same as that of the EQUIL paper. The representation is generalized to allow for up-down asymmetry:
\[
\begin{aligned}
\hat R(\hat r,\omega)={}&1+\epsilon(\hat r+\epsilon^2 \hat P)\cos\omega-\epsilon^2\hat\Delta
 +\epsilon^2\sum_{m=2}^{N_s}\hat S_{m}\cos\bigl[(m-1)\omega\bigr]
 +\epsilon^2\sum_{m\in\mathcal M_A}\hat A_m\sin\bigl[(m-1)\omega\bigr],\\
\hat Z(\hat r,\omega)={}&\epsilon(\hat r+\epsilon^2 \hat P)\sin\omega+\epsilon^2\hat \Delta_Z
 -\epsilon^2\sum_{m=2}^{N_s}\hat S_{m}\sin\bigl[(m-1)\omega\bigr]
 +\epsilon^2\sum_{m\in\mathcal M_A}\hat A_m\cos\bigl[(m-1)\omega\bigr].
\end{aligned}
\]
The up-down asymmetry leads to a **vertical Shafranov shift** $\Delta_Z$. The unknowns related to the geometry are thus the functions $(\hat \Delta, \hat \Delta_Z, \hat S_m, \hat A_m, \hat P)$ to be found on $\hat r \in [0;1]$. Clearly $\hat P$ is a radial coordinate renormalization and since the action above is independent of the chosen independent coordinates, its equation is not derived from the action, but from the requirement that the Straight Field Line poloidal coordinate $\vartheta(\hat r,\omega)$ be a (periodic) angle:
\[
\frac{\partial \vartheta}{\partial \omega} = \frac{\hat J/\epsilon^2}{\hat r \hat R^2}\frac{1-\sigma_{\rm cyl}(\hat r)}{1-\sigma_\parallel}
\]
So that the requirement that $\vartheta(r,2\pi)-\vartheta(r,0)=2\pi$ means that:
\[
  \left\langle\frac{\hat J/\epsilon^2}{\hat R^2(1-\sigma_\parallel)}\right\rangle_\omega=\frac{\hat r}{1-\sigma_{\rm cyl}(\hat r)}
\]
which preserves the near axis SFL jacobian $J_\mathrm{SFL}\sim rR^2/R_0$ and gives the renormalization $\epsilon_\mathrm{eff} = \epsilon + \epsilon^3 \hat P(\hat r=1;\sigma_\parallel)$. The equations for $(\hat \Delta, \hat \Delta_Z, \hat S_m, \hat A_m)$ are obtained from the geometrical variation of the last section, and more explicitly:
\[
\begin{aligned}
\delta_\Delta R=-\Lambda,&\qquad\delta_\Delta Z=0, \\
\delta_{S_{m+1}} R=\Lambda\cos(m\omega),&\qquad \delta_{S_{m+1}} Z=-\Lambda\sin(m\omega), \\
\delta_{A_m} R=\Lambda\sin((m-1)\omega),&\qquad \delta_{A_m} Z=+\Lambda\cos((m-1)\omega), \\
\delta_{\Delta_Z} R=0,&\qquad \delta_{\Delta_Z} Z=\Lambda
\end{aligned}
\]
where $\Lambda(\hat r)$ is a test-function and e.g. $\delta_\Delta$ corresponds to the variation $\hat \Delta(\hat r)\to\hat\Delta(\hat r) + \Lambda(\hat r)$.

Postprocessing integrates the same exact constraint,
\[
\vartheta(r,\omega)=\int_0^\omega \mathrm d\omega'\,
\frac{[1-\sigma_{\rm cyl}(r)]J}
{\epsilon^2 rR^2(1-\sigma_\parallel)},
\]
and inverts it on a uniform periodic $\vartheta$ grid. The resulting fields are constructed independently of the large legacy `equil_SFL` routine.

The implementation is now isolated in `src/equilVariational_SFL.m`.  When `do_SFL=true`, the variational B-splines are also evaluated to second radial order and the Fourier geometry is differentiated analytically.  This gives

\[
R_s=R_r-\frac{\vartheta_r}{\vartheta_\omega}R_\omega,
\qquad
R_\vartheta=\frac{R_\omega}{\vartheta_\omega},
\]

and the analogous derivatives of $Z$, followed by $g_{11}$, $g_{12}$, $g_{22}$, $J_{\rm SFL}/s$, the contravariant magnetic field, pressure scalars, curvature, and the contravariant equilibrium current. The current's perpendicular part is evaluated from the solved anisotropic force balance; Ampere's law for $j^\vartheta$ fixes its field-aligned freedom and gives the exact isotropic MISHKA limit. The small set of tensor derivatives required by the Mishka coefficient files is also returned. The second-derivative B-spline
matrices are neither constructed nor evaluated when `do_SFL=false`.

Both $\boldsymbol\kappa=(\boldsymbol b\mathbin\cdot\nabla)\boldsymbol b$
and the static force represented by $\boldsymbol j\times\boldsymbol B$ are
analytically perpendicular to $\boldsymbol B$. In the SFL postprocessor,
their contravariant components are assembled from independently evaluated
first- and second-geometry derivatives. The resulting tiny parallel
component is therefore a reconstruction error, not another residual of the
variational equilibrium. In particular, its contribution to the pressure
force is multiplied by $P_\parallel-P_\perp$, so a fixed rejection threshold
spuriously rejects otherwise converged equilibria as the parallel anisotropy
is increased. The implementation projects the reconstructed curvature and
force perpendicular to $\boldsymbol B$ before constructing the current. It
retains the unprojected defects as `sfl_curvature_parallel_error_raw` and
`sfl_force_parallel_error_raw`, and the final three-component consistency as
`sfl_force_closure_error`; these are diagnostics and not convergence tests.

The legacy-compatible radial space does not impose polar Taylor regularity. Its first-span second derivatives can therefore be large and cancelling even when the sampled $R$ and $Z$ are smooth.  Only for that space, `equilVariational_SFL` replaces $\vartheta_\omega$ (and, in isotropy, $B$) inside the first span by a cubic continuation from its three outer resolved samples and the exact axis value.  The `axis_regular` space uses the native expressions without this continuation.

The variational postprocessor also returns the legacy-normalized global quantities `Wk`, `Wkpar`, `Wkperp`, `Wkrot`, `Wp`, `Ip`, `Ftt`, `Ft0`, `Ft`, `bp`, `bppar`, `bpperp`, `bprot`, `li`, and `bpli2`, together with the compatibility aliases `BBp2`, `BBt2`, and `BB2`. These are evaluated directly on the variational Gauss grid. The remaining legacy-only moments are deliberately not copied: `jphi` and the contravariant current require a separate weak reconstruction, while quantities such as `N0`, `Mm1`, and `Lm1` belong to the old analytical postprocessing rather than the equilibrium interface.



### Optional up-down asymmetric profiles

The $A_m$ profiles are the quadrature partners of $S_m$, rather than the $V_m$ family of the Fitzpatrick expansion. In complex form,
\[
w=R-iZ=\epsilon(r+\epsilon^2P)e^{-i\omega}
+\epsilon^2\sum_{m\ge2}(S_m-iA_m)e^{i(m-1)\omega}
-\epsilon^2(\Delta+i\Delta_Z).
\]
Their physical normal displacements are therefore the complete Fourier pair
\[
\boldsymbol e_r\cdot\delta_{S_m}\boldsymbol x=\cos(m\omega),
\qquad
\boldsymbol e_r\cdot\delta_{A_m}\boldsymbol x=\sin(m\omega).
\]
In particular, $(S_2,A_2)$ is an elliptic quadrupole with orientation
\[
\alpha_2=\frac12\operatorname{atan2}(A_2,S_2).
\]
Unequal upper/lower triangularity is not a single $A_3$. At small shaping, a smooth split-Miller boundary starts with linked $A_2$ and $A_4$ terms of opposite sign, followed by higher even modes. `A_modes` contains physical mode numbers $m\ge2$; `LX.Abc`, `LY.A`, and `LY.Ap` use that order. `Na` is the number of selected profiles and defaults to zero. It need not equal `Ns`: independence means that the cosine and sine coefficients are separate variations, not that every calculation must instantiate both. Thus `Na=0` recovers the symmetric problem, sparse `A_modes` are useful for measured boundaries, and `A_modes=2:(Ns+1)` gives the common paired truncation. The default axis class is
\[
A_m(r)=r^{m-1}(a_0+a_2r^2+\cdots).
\]
`A_leading_powers` is not needed to make the $S_m$ and $A_m$ variations independent. It is only an optional restriction on their magnetic-axis Taylor class. For example, `A_modes=[2 4]` with `A_leading_powers=[3 3]` represents unequal triangularity without allowing an $O(r)$ tilted-axis $A_2$ component. A genuinely tilted ellipse, and the standard MEQ benchmark, keep the default $A_2=O(r)$ class. For boundary fitting, `src/utils/get_Sbc.m` reconstructs the EQUIL angle by eliminating both real and imaginary negative Fourier modes with frequency two and above. It then returns
\[
S_m=\Re c_{m-1},\qquad A_m=-\Im c_{m-1},
\]
where $c_n$ is the positive-frequency coefficient of $w$. The angle solve accepts `angle_mode_limit` and `angle_initial_coefficients` and rejects requested modes above the conservative Nyquist cutoff.

The angle reconstruction is nonlinear. With a dense contour and many gauge harmonics it can possess several numerically accessible parameterizations, even though the physical curve is unchanged. Dense reconstructions should
therefore be seeded with the converged low-order angle coefficients from the native contour. `get_Sbc` accepts `angle_mode_limit` and `angle_initial_coefficients` for this purpose. It now rejects asymmetric modes
above the conservative Nyquist cutoff instead of failing later with an array index error.

`Analysis/get_fbt_boundary.m` can sample a surface on uniform geometrical rays from a spline or linear interpolation of FBT's existing Cartesian flux grid. This is postprocessing, not an increase of the FBT PDE resolution. It permits high-mode sensitivity tests, but only modes stable under point count, interpolation method, and angle-gauge truncation should be interpreted.

`Analysis/MEQ_related/asymmetric_variational_study.m` checks the Fourier gauge and the $A_m$ geometry derivatives, then exercises the FBT-to-EQUIL path for unequal upper/lower triangularity and tilted elongation.

For routine conversion rather than a detailed benchmark, `src/utils/meq_to_equil.m` accepts an FBT/MEQ equilibrium and returns its matched variational EQUIL state. `examples/match_meq.m` is the minimal user-facing example; detailed profile and flux-surface reports remain under `Analysis/MEQ_related`.

### MEQ benchmark

For the comparison with MEQ, the magnetic-axis position $(R_A,Z_A)$ is known and we take $R_0=R_A$. The retained MEQ boundary is first treated as a geometrical object independently of the EQUIL iteration. Defining
\[
W_\Gamma(\Theta)=\frac{R_\Gamma(\Theta)-R_A-i[Z_\Gamma(\Theta)-Z_A]}{R_0},
\]
where $\Theta$ is an arbitrary parametrization, which we then transform to the geometrical angle $\Theta\to\theta_\mathrm{geom}$, then `get_Sbc` reconstructs the angle $\omega$ for which
\[
W_\Gamma=\gamma_{-1}e^{-i\omega}+\gamma_0
+\sum_{k\geq1}\gamma_ke^{ik\omega},
\qquad \gamma_{-1}>0,
\]
and all harmonics with frequency smaller than $-1$ vanish. From these $\gamma_k$ we start EQUIL iterations to find $\epsilon$: at each trial value, the fixed boundary data are algebraically rescaled:
\[
S_m(1)=\frac{\Re\gamma_{m-1}}{\epsilon^2},\qquad
A_m(1)=-\frac{\Im\gamma_{m-1}}{\epsilon^2},\qquad
\Delta_Z(1)=-\frac{\Im\gamma_0}{\epsilon^2}.
\]
The profiles $q$ and $P_\parallel$'s push-back from $\hat\psi\to \hat r$ must be reversed, so the MEQ input profiles on the normalized flux are expressed as functions of $\hat r$ via the previous iteration's mapping of  $\hat \psi(\hat r)$. The axis value of the $q$ fit is taken from MEQ's separately evaluated `qA`, when available, rather than extrapolated from the first `iqQ` sample; this avoids selecting a horizontally displaced branch in strongly shaped cases. The benchmark updates $\epsilon$ by imposing:
\[
F_\epsilon=\epsilon[1+\epsilon^2\hat P(1)]-\gamma_{-1}=0.
\]
The horizontal shift is not prescribed. It remains a prediction of the force-balance problem and is checked through
\[
F_R=-\epsilon^2\hat \Delta(1)-\Re\gamma_0.
\]
Thus $F_\epsilon$ is an iteration residual whereas $F_R$ is a validation quantity.

### Magnetic field strength $B$

At each two-dimensional quadrature point $(\hat r_q, \omega_q)$ the scalar equation
\[
G(B)=1-\sigma(r,R,B)-\frac{T}{RB_\phi}=0
\]
is solved via newton steps, using the following jacobian:
\[
G_B=-\sigma_B+\frac{TB}{RB_\phi^3},\qquad
\sigma_{\parallel,B}=\frac{\Pi_{\parallel,BB}}{B}
-\frac{\Pi_{\parallel,B}}{B^2}.
\]
For isotropic cases \(\sigma_\parallel=0\) and the exact solution is used:
\[
B=\sqrt{B_p^2+\left(\frac{T}{R}\right)^2}.
\]
Anisotropic closures retain the safeguarded scalar Newton iteration and enforce \(B^2>B_p^2\), \(1-\sigma_\parallel>0\), and non-small \(|G_B|\) which physically correspond to a positive magnetic field norm, the firehose stability criteration and the mirror stability criterion.

For any global perturbation \(x_j\), local implicit differentiation gives
\[
B_{x_j}=-\frac{G_{x_j}}{G_B}.
\]
The implementation evaluates fixed-\(B\) directional derivatives and then forms the reduced response along \((x_j,B_{x_j})\). This is the scalar static-condensation or Schur-complement formula
\[
R^{\rm red}_x=R_x-R_BG_B^{-1}G_x.
\]

## Radial discretization

### Axis-regular B-spline space

Near-axis regularity in local Cartesian coordinates requires a polar harmonic with angular wave number $n$ to have the Taylor class
\[
f_n(\hat r)=\hat r^n(c_0+c_2\hat r^2+c_4\hat r^4+\cdots).
\]
The variational implementation uses open, uniform B-splines and constrains the polynomial on the first span to this class. In particular,
\[
t_2,\Delta=\mathcal O(\hat r^2),\qquad
P=p_1\hat r+p_3\hat r^3+\cdots,qquad
S_{n+1},A_{n+1}=\mathcal O(\hat r^n).
\]
For degree four, the allowed first-span powers are therefore $(\hat r^2,\hat r^4)$ for $t_2$ and $\Delta$, $(\hat r,\hat r^3)$ for $P$ and $S_2$, and only $\hat r^4$ for $S_5$. The constraint removes inadmissible axis coefficients; it is not a preconditioner for the equations away from the axis.

This is the default selected by `radial_discretization='axis_regular'`. The alternative `radial_discretization='legacy'` uses exactly the open, uniform B-spline construction of the default `equilSol` discretization, without the first-span Taylor-parity constraints. Every profile removes the first clamped B-spline and therefore vanishes on axis. The $\Delta$ space also removes the second coefficient, while each fixed-edge $S_m$, $A_m$, or $\Delta_Z$ space removes the final unknown and uses that B-spline as its prescribed boundary lift. Consequently all $S_m$ profiles have the same unconstrained near-axis spline space in this mode, as in legacy EQUIL.

For equal `(m,nq,spline_p)`, the legacy option has the same knot vector, radial Gauss points, and common profile matrices as `assemble_FE_matrices_bspline_neumann`: the $t_2$ and $P$ matrices equal legacy `P0` and `P1`; the $\Delta$ matrices equal their columns `2:end`; and the unknown $S_m$ matrices equal columns `1:end-1`, with `P0_end` as the lift. The historical `r_nodes` array with $2m+1$ entries is retained as metadata, although both assemblers use $m$ actual spline spans. Choosing the same `om_pts` also gives the identical angular grid.

If $d=m+p-1$, the legacy-compatible lengths are
\[
N_{t_2}=N_P=d,\qquad N_\Delta=d-1,\qquad
N_{S_m}=N_{A_m}=N_{\Delta_Z}=d-1.
\]
The variational formulation has no dummy isotropic $B_s$ block. Thus at `m=15`, `p=10`, and `Ns=4` it has 163 unknowns, while `equilSol` has the same 163 physical-profile unknowns plus 24 dummy $B_s$ coefficients, for 187 total.

The numerical parameters have distinct roles. `m` is the number of radial spans, `spline_p` is the polynomial degree, and `nq` is the number of Gauss points per span. Thus
\[
N_q=m\,n_q
\]
is the number of radial integration points, not the number of global radial unknowns. Increasing `nq` at fixed $(m,p)$ increases residual and local-$B$ work but leaves the global system size unchanged. This distinction is obscured in `equilSol` because its defaults set `spline_p=nq`; the two parameters remain mathematically independent there as well.

For the default axis-regular option, a leading power $n\le p$ gives
\[
N_{\rm free}(n)=m+\left\lfloor\frac{p-n}{2}\right\rfloor,
\qquad
N_{\rm edge}(n)=m-1+\left\lfloor\frac{p-n}{2}\right\rfloor
\]
coefficients for a free-edge and fixed-edge profile, respectively. With `Ns=4`, `spline_p=4` and no asymmetric profiles, the variational system consequently has $7m+1$ unknowns: 85 at `m=12`. The settings used by the legacy side of `examples/variational_vs_legacy.m`, namely `m=15`, `spline_p=nq=10`, `Ns=4` and `Nb=1`, give 187 unknowns. The difference is primarily radial degree and span count, not the angular or Gauss grids.

The default variational values remain `radial_discretization='axis_regular'`, `m=8`, `spline_p=4`, `nq=6`, and `om_pts=96`. In axis-regular mode the requested spline degree is raised when `Ns` or an asymmetric leading power cannot be represented on the first span. Legacy mode retains the requested degree, so the same `(m,nq,spline_p,om_pts)` can be passed unchanged to both solvers.

`examples/variational_analysis/matched_legacy_discretization_check.m` asserts equality of the two grids, projection matrices, profile lengths, boundary lift, and quadrature projection before solving the matched case-2 equilibrium.

## Analytical variational Jacobian

For the unmodified static isotropic equation of state,
\[
P_\parallel=P_\perp=P(\psi),\qquad \sigma_\parallel=0\Rightarrow RB_\phi=T(\psi)
\]
Substitution in the action gives the isotropic action (Lao, 1981)
\[
\mathcal A_{\rm iso}[\psi]
=\int\mathrm dV\left[
\frac{B_p^2}{2\mu_0}-P(\psi)-\frac{T(\psi)^2}{2\mu_0R^2}
\right].
\]
Geometrical and flux variations stay the same as in the anisotropic case of the first section. The geometrical variation residuals are given by:
\[
\mathcal R_{\rm g}[\delta R,\delta Z]
=2\pi\iint\mathrm dr\,\mathrm d\omega\left[
\mathcal L_{\rm iso}\,\delta J
+\frac{J}{2\mu_0}\delta B_p^2
+\frac{JT^2}{\mu_0R^3}\delta R
\right],
\]
For the analytical Jacobian, let a dot denote the trial direction generated by one global coefficient. Since $q(r)$ is fixed,
\[
\dot\psi_r=\frac{r}{R_0q}\dot T
\]
The required kinematic derivatives are
\[
\begin{aligned}
\dot J={}&(R_rZ_\omega-R_\omega Z_r)\dot R
+R\left(Z_\omega\dot R_r+R_r\dot Z_\omega
-Z_r\dot R_\omega-R_\omega\dot Z_r\right),\\
\dot g_{\omega\omega}={}&2R_\omega\dot R_\omega+2Z_\omega\dot Z_\omega,\\
\dot B_p^2={}&
2\psi_r\dot\psi_r\frac{g_{\omega\omega}}{J^2}
+\psi_r^2\frac{\dot g_{\omega\omega}}{J^2}
-2\psi_r^2\frac{g_{\omega\omega}\dot J}{J^3},\\
\dot{\mathcal L}_{\rm iso}={}&
\frac{\dot B_p^2}{2\mu_0}
-\frac{T\dot T}{\mu_0R^2}
+\frac{T^2\dot R}{\mu_0R^3}.
\end{aligned}
\]
In the code $T=a_0+\epsilon^2t_2$, and hence a $t_2$ trial function gives $\dot T=\epsilon^2\dot t_2$.
Holding the flux test function $\Phi$ fixed gives
\[
\begin{aligned}
\dot{\mathcal R}_\psi[\Phi]
={}&2\pi\iint\mathrm dr\,\mathrm d\omega\,
\frac{J}{\psi_r}\Bigg\{
\left(\frac{\dot J}{J}-\frac{\dot\psi_r}{\psi_r}\right)
\left[
\frac{B_p^2}{\mu_0}\Phi_r
-\left(P_r+\frac{TT_r}{\mu_0R^2}\right)\Phi
\right]\\
&\qquad+\frac{\dot B_p^2}{\mu_0}\Phi_r
-\left[
\frac{\dot T\,T_r+T\dot T_r}{\mu_0R^2}
-\frac{2TT_r}{\mu_0R^3}\dot R
\right]\Phi\Bigg\}.
\end{aligned}
\]
For a fixed geometrical test direction, the mixed variations are
\[
\begin{aligned}
\dot{\delta J}={}&
\delta R(\dot R_rZ_\omega+R_r\dot Z_\omega
-\dot R_\omega Z_r-R_\omega\dot Z_r)\\
&+\dot R(\delta R_rZ_\omega+R_r\delta Z_\omega
-\delta R_\omega Z_r-R_\omega\delta Z_r)\\
&+R(\delta R_r\dot Z_\omega+\dot R_r\delta Z_\omega
-\delta R_\omega\dot Z_r-\dot R_\omega\delta Z_r),\\
\dot{\delta g}_{\omega\omega}={}&
2\dot R_\omega\delta R_\omega+2\dot Z_\omega\delta Z_\omega.
\end{aligned}
\]
Equivalently, differentiating the compact coefficients above gives
\[
\begin{aligned}
\dot C_J={}&-\frac12\dot B_p^2-\frac{T\dot T}{R^2}
+\frac{U\dot R}{R^3},\\
\dot C_R={}&\frac{U\dot J}{R^3}
+\frac{2JT\dot T}{R^3}-\frac{3JU\dot R}{R^4},\\
\dot C_g={}&\epsilon^2\left(
\frac{\psi_r\dot\psi_r}{J}
-\frac{\psi_r^2\dot J}{2J^2}\right).
\end{aligned}
\]
Thus every geometrical Jacobian entry is assembled as
\[
\dot C_J\delta J+C_J\dot{\delta J}
+\dot C_R\delta R+C_R\dot{\delta R}
+\dot C_g\delta g_{\omega\omega}+C_g\dot{\delta g}_{\omega\omega},
\]
with the corresponding differentiated edge term. This is the compact normalized form of differentiating $\mathcal R_{\rm g}$ directly.

Finally, in isotropy the gauge residual and its tangent are
\[
\begin{aligned}
\mathcal R_P[\Lambda]
&=\int_0^1\mathrm dr\,\Lambda
\left[\left\langle\frac{J}{\epsilon^2R^2}\right\rangle_\omega-r\right],\\
\dot{\mathcal R}_P[\Lambda]
&=\int_0^1\mathrm dr\,\Lambda\left\langle
\frac{\dot J}{\epsilon^2R^2}
-\frac{2J\dot R}{\epsilon^2R^3}
\right\rangle_\omega.
\end{aligned}
\]
The $P$ equation is therefore unchanged; eliminating $B$ does not create the $\Delta$--$P$ coupling.

For a given trial profile, every differentiated physical field can be written in the separated form
\[
\dot f(r,\omega)=f_0(r,\omega)B_0(r)+f_1(r,\omega)B_1(r),
\]
where $B_0$ and $B_1$ contain all radial value and derivative trial functions for that profile. Consequently a typical value-test contribution is assembled as
\[
B_{0,{\rm test}}^T\operatorname{diag}\!\left(w_r\langle f_0\rangle_\omega\right)B_{0,{\rm trial}}
+B_{0,{\rm test}}^T\operatorname{diag}\!\left(w_r\langle f_1\rangle_\omega\right)B_{1,{\rm trial}},
\]
with $B_{1,{\rm test}}$ replacing $B_{0,{\rm test}}$ for derivative-test terms. The edge terms have the analogous outer-product form. Both analytical Jacobians therefore evaluate two-dimensional coefficient fields, perform the angular average, and only then apply the radial matrices. There is no loop over individual Jacobian columns and no three-dimensional quadrature-by-angle-by-coefficient tangent array. The result is returned sparse.

`equil_variational_finite_difference_jacobian` remains only as an explicitly selected centered-difference verification utility. Production dispatch selects the reduced residual and `equil_variational_isotropic_jacobian` for the native static isotropic closure. Every other native closure uses the general residual and `equil_variational_general_jacobian`; there is no EOS-name or `mach20` allowlist. A custom residual requires an explicitly supplied analytical Jacobian.

`equil_variational_general_jacobian` differentiates the same weak residual for the full $P_\parallel(r,R,B)$ closure without changing its normalization, vacuum subtraction, or natural edge terms. Static anisotropy is the special case in which all explicit-$R$ derivatives vanish. At fixed local $B$,

\[
\dot G|_B=-\sigma_R\dot R-\frac{\dot T}{RB_\phi}
+\frac{T\dot R}{R^2B_\phi}
-\frac{T\dot B_p^2}{2RB_\phi^3},
\qquad \sigma_R=\frac{P_{\parallel,BR}}{B}.
\]

The pointwise response and the remaining local tangents are

\[
\dot B=-\frac{\dot G|_B}{G_B},\qquad
\dot B_\phi=\frac{B\dot B-\dot B_p^2/2}{B_\phi},\qquad
\dot\sigma=\sigma_R\dot R+\sigma_B\dot B,
\]
\[
\dot P_\parallel=P_{\parallel,R}\dot R+P_{\parallel,B}\dot B,
\qquad
\dot P_{\parallel,r}=P_{\parallel,rR}\dot R+P_{\parallel,rB}\dot B,
\]
\[
\dot P_{\parallel,R}=P_{\parallel,RR}\dot R+P_{\parallel,RB}\dot B.
\]

The residual still contains the partial derivative $P_{\parallel,r}|_{R,B}$; no $R_rP_{\parallel,R}$ term is introduced. Its Newton tangent contains $P_{\parallel,rR}\dot R$ because the partial derivative itself changes with the trial geometry. The geometrical coefficient is

\[
C_R=J\left(-P_{\parallel,R}+\frac{TB_\phi}{R^2}\right),
\]

with the implemented vacuum subtraction applied algebraically before quadrature.

The normalized EOS adapter exposes `Pi`, `Pr`, `PR`, `PB`, `PBB`, `PrB`, `PRB`, `PrR`, and `PRR`. Native EOS functions return this named second-order interface when called with the `variational` request, so the weak Jacobian does not evaluate or require any legacy third derivatives. This includes the branch-consistent, $C^2$ low-$B$ `biMaxwellian` closure. Custom pressure maps must provide the same fields. The profile $\sigma_{\rm cyl}(r)$ and $a_0=1-\sigma_{\rm cyl}(0)$ remain prescribed during the global Newton solve.

Local condensation preserves the same separated representation used by the isotropic matrix assembler: every tangent is stored as two two-dimensional fields multiplying the radial trial value and derivative matrices. The implementation loops over physical profile blocks, averages in the angular direction, and then performs sparse radial matrix products. It does not loop over global coefficients or allocate an `Nq`-by-`Nomega`-by-`Ncoeff` tangent.

`examples/variational_analysis/general_jacobian_check.m` checks initial and converged matrices and random directions for rotating and static `de_Blank`, `thermal_and_deBlank`, and `biMaxwellian`, together with the isotropic limit and the low-$B$ EOS derivatives. A representative R2016b run gave full-matrix errors between $4.5\times10^{-10}$ and $1.5\times10^{-8}$, directional errors below $5.6\times10^{-9}$, and an isotropic-limit error of $1.2\times10^{-16}$. The rotating analytical and column-reference construction times were 0.0085 s and 0.101 s on the 25-degree-of-freedom check. Absolute timings remain machine dependent.

## Difference from standard EQUIL

The isotropic legacy path used in the comparison script solves a strong force-balance residual. It evaluates second radial derivatives, takes selected poloidal Fourier moments, and projects those moments with the radial mass matrix. Its `Bs` block is a dummy equation in `residuals_iso_static`: the coefficient is driven to zero while $RB_\phi=1+\epsilon^2t_2$ is used directly. The generated `jacobian_iso_static` differentiates this complete residual analytically.

The variational path instead evaluates the weak flux equation, the geometrical first variations of the action, and the weak SFL constraint. Only first radial derivatives occur. The global unknowns are $(t_2,\Delta,P,S_m)$, with optional $(A_m,\Delta_Z)$. In static isotropy $B$ is substituted in the reduced action before residual assembly; for other closures it is eliminated at every two-dimensional quadrature point. The public organization remains parallel: `equilVariationalSol` constructs the space, `equilVariationalX` supplies the prescribed profiles, `equilVariationalY` solves the system, and `equilVariationalPP` constructs `LY`.

For a smooth isotropic equilibrium, integration by parts and matching boundary conditions should relate the strong and weak force-balance statements. This does not make the discrete systems identical. They use different axis spaces, different test projections, and different treatments of the natural edge terms. The legacy space only eliminates selected endpoint B-spline coefficients; it does not impose the harmonic-by-harmonic Taylor parity of the variational space. Conversely, the variational equations contain the edge terms produced by the weak action and enforce the radial-coordinate constraint only in the $P$ test space. A residual that remains after both discretizations and nonlinear solves have independently converged is therefore possible, but it should not be interpreted until a stable variational $h$-sequence has been obtained.

There is also a postprocessing difference at the radial endpoints. `equilVariationalPP` evaluates every constrained spline at $\hat r=0$ and $1$. `equilPP` constructs `r_plt=[0;r_q;1]` but pads the free profiles $t_2$, $\Delta$, and $P$ at the edge with their value at the last interior Gauss point. The final legacy row is consequently not an evaluation at $\hat r=1$. This affects plotted edge values, the edge geometry, and `to_venus_cleaner`; interior comparisons should exclude both artificial endpoint rows and report the edge separately.

## Numerical behaviour of the variational solve

The previous efficiency deficit came from the Jacobian, not the axis prescription or the scalar $B$ solve. If the system has $N$ global coefficients, the reference `equil_variational_jacobian` perturbs every coefficient in both directions. Each column evaluates two global states, two fixed-$B$ constraints, the condensed response
\[
B_{x_j}=-G_B^{-1}G_{x_j},
\]
and two weak residual assemblies. It therefore costs approximately $2N$ residual assemblies and is retained as a verification reference, not as a production fallback.

The production analytical tangents no longer pay this cost. The reduced isotropic residual removes the $B$ constraint, while the general formulation condenses it pointwise; both evaluate all coefficients of each trial profile together and return sparse matrices. The first block-vectorized isotropic implementation was nevertheless slow at the matched legacy discretization: it repeatedly formed $N_q\times N_\omega\times n_c$ arrays, where $n_c$ is the number of coefficients in one trial-profile block. Profiling the case-2 initial state attributed about 52% of its 1.0 s Jacobian time to 1005 broadcast field multiplications; the linear solve itself took only about 1 ms. The excess per-iteration time was therefore memory traffic in the analytical assembly, not conditioning, the line search, the radial axis prescription, or the reduced isotropic equations.

The separated matrix assembly above removes those arrays. On the matched case-2 setting `m=15`, `nq=spline_p=10`, and `om_pts=300`, it agrees with the previous analytical matrix to a relative Frobenius error of $9.2\times10^{-16}$ and with centered finite differences to order $10^{-8}$. Representative warmed R2016b timings at the initial state are 0.010 s for the variational residual and 0.064 s for its Jacobian, compared with 0.016 s and 0.311 s for the generated legacy residual and Jacobian. The complete direct solves without SFL postprocessing took 0.57 s for seven variational Newton iterations and 3.44 s for eight legacy iterations on the same machine. These absolute timings are machine dependent, but they show that the matched variational iteration is no longer intrinsically slower; it is faster here despite retaining the weak residual, natural edge terms, and line search.

The analytical matrices also remove `jacobian_step` from normal solver operation. Centered differences remain useful as a verification oracle: checks at the comparison resolution give a relative Frobenius difference of order $10^{-8}$, consistent with truncation and roundoff in the finite-difference reference. Solving a small isotropic problem with the old general residual/Jacobian and with the reduced residual plus centered differences gives the same coefficient vector to approximately $2.5\times10^{-14}$.

### Why increasing `m` is difficult

The axis constraint removes singular Taylor components but does not remove the coordinate/gauge degeneracy of the action. The $P$ equation is the discrete condition that selects a radial parametrization, and it is coupled most strongly to the horizontal shift $\Delta$. With the exact isotropic Jacobian in case 2, the smallest right singular vector of the converged `m=16` matrix has 94% of its norm in the $\Delta$ block and 32% in $P$; the smallest left singular vector has 81% and 58%, respectively. The condition number grows from approximately $1.8\times10^6$ at `m=12` to $3.6\times10^7$ at `m=16`, and to $4.6\times10^8$ on the converged `m=20` branch. For comparison, the same-degree 115-variable legacy Jacobian has condition number about $2.7\times10^6$.

The same near-null direction is therefore present without finite-difference noise: it belongs to the discrete reduced equations. The old finite-difference tangent could mask or amplify it depending on `jacobian_step`, but it was not its source. In an exact-Jacobian transfer sequence `m=12,16,20`, the first two refinements converge, while the `m=20` solution moves to a state with $\min[J/(\epsilon^2r)]\simeq2.8\times10^{-2}$. Transfer to `m=24` then stalls near a residual of $3.5\times10^{-4}$ and a condition number of order $10^{10}$. Analytical differentiation consequently fixes the per-iteration work and tangent accuracy, but not the branch-selection problem.

The remaining difficulty is globalization. The line search accepts any residual-decreasing step for which $R>0$, $J>0$, and the local magnetic constraints remain admissible. It has no trust region, no lower margin on $J/(\epsilon^2r)$, no branch-distance criterion, and no pseudo-arclength parameter. In case 2 the converged `m=20` state reached
\[
\min\frac{J}{\epsilon^2r}=2.8\times10^{-2},
\]
compared with about $0.39$ at `m=12` and `m=16`; its profiles also moved farther from both the preceding mesh and the legacy result. Transfer to `m=24` then stagnated. Case 1 shows the same pattern more mildly: `m=16` converges from the transferred `m=12` state, while the converged `m=20` profiles move sharply away from both. Thus a small algebraic residual at one refined mesh is not, by itself, evidence of $h$-convergence.

The arbitrary sequence $m=12,16,20,\ldots$ also gives non-nested uniform knot sets, so a transferred state is a projection rather than an exact injection. This is a secondary effect: the projected residuals in the tests were already $O(10^{-4})$, after which the near-null Newton direction and residual-only line search determined the branch. The current evidence points to a coupled gauge/globalization problem in the discrete weak equations rather than a failure of the local axis Taylor prescription or an inaccurate isotropic tangent.

Radial and angular quadrature are not limiting these tests. Re-evaluating the converged case-2 `m=16`, `nq=8`, `om_pts=96` coefficients gives the same residual at `nq=8,10,12` and at 48, 96, or 192 angular points. Re-solving with `nq=12` and 192 angular points requires no Newton update. Increasing `nq` therefore increases cost without curing the observed `m` instability.

## Newton and beta continuation

When no `LX.x` is supplied, the present continuation solves beta zero, then
\[
\beta_{\rm low}=\operatorname{sign}(\beta)\min(|\beta|,1),
\]
and finally the target beta; failed jumps are bisected. If the finite-aspect-ratio beta-zero solve fails, a seed is obtained at $\epsilon=10^{-3}$ and retried at the requested $\epsilon$. Supplying `LX.x` bypasses this path.

This schedule remains unnecessarily expensive for the two mild $\beta=0.32$ comparison cases. With the analytical isotropic Jacobian, both direct target solves converge in 13 Newton iterations: representative times are $1.1$ s for case 1 and $0.9$ s for case 2. The default continuation takes approximately $3.6$ s and $2.9$ s, respectively. Case 1 still requires the small-$\epsilon$ beta-zero seed before accepting the finite-$\epsilon$ beta-zero and target stages; case 2 accepts beta zero directly but spends 33 iterations there before the 9-iteration target stage. The target is therefore easier than the nominal low-beta initialization in these cases, while `equilSol` also converges directly. The direct and continued variational profiles agree to better than $1.3\times10^{-6}$ in relative sampled norm, so the added stages do not select a measurably different equilibrium here.

A more appropriate strategy is to try the beta-aware target initial guess first and invoke continuation only after a diagnosed failure. Continuation should reuse the last accepted global state and, for anisotropic closures, local $B$; it should adapt its parameter to the failure mechanism and distinguish a pressure-path failure from a $\Delta$--$P$ geometry failure. A trust-region or pseudo-arclength method with a positive margin on $J/(\epsilon^2r)$ is more relevant to the latter than additional beta stages. Radial refinement should similarly be performed as an explicit mesh-continuation path with projection diagnostics. Reworking the schedule alone will reduce wasted solves, but it will not remove the exact near-null direction.

## The two comparison cases

The earlier unmatched comparison used the axis-regular variational `m=12`, `p=4` space and the legacy `m=15`, `p=10` space. On $0.01\le\hat r\le0.99$, case 1 then gave relative sampled $L^2$ differences of approximately 1.32% in $\Delta$, 0.14% in $t_2$, 0.59% in $P$, and 0.09% in $S_2$. Case 2 gave 8.33%, 1.52%, 3.77%, and 0.82%, respectively. These were not same-space consistency tests.

`examples/variational_vs_legacy.m` now selects the legacy-compatible variational space and passes `m=15`, `nq=spline_p=10`, and `om_pts=300` to both solvers, with beta continuation disabled and a common nonlinear tolerance. Both solvers then converge directly: case 1 takes 8 variational versus 9 legacy iterations, and case 2 takes 7 versus 8. On the common interior Gauss grid, the case-1 relative differences in $(\Delta,t_2,P,S_2)$ are approximately $(4.2\times10^{-5},7.7\times10^{-6},1.4\times10^{-5},6.0\times10^{-5})$. For case 2 they are $(8.4\times10^{-7},2.2\times10^{-6},1.6\times10^{-6},1.6\times10^{-5})$. The case-2 relative Frobenius difference in `theta_SFL` is $5.3\times10^{-6}$, with maximum differences $3.3\times10^{-4}$ in the angle, $3.7\times10^{-5}$ in `RR_sfl`, and $2.0\times10^{-4}$ in `ZZ_sfl`. Thus most of the former profile discrepancy was caused by the different radial trial/test spaces; the remaining same-space difference is the appropriate measure of the strong-versus-weak equations and their postprocessing conventions.

The case-2 discrepancy does not decrease monotonically along the separate axis-regular $h$-refinement sequence. At comparable global size, the converged 113-variable variational `m=16`, `p=4` state differs from the 115-variable legacy `m=12`, `p=4` state by about 58% in $\Delta$, 9.4% in $t_2$, 24% in $P$, and 4.8% in $S_2$ in the same sampled norm. The exact tangent confirms that this is a branch/refinement failure, not evidence that the lower-resolution `m=12` result merely needs more radial points. The matched legacy-space result above separates that axis-regular refinement issue from the much smaller strong-versus-weak residual.

`examples/variational_simple_case.m` remains useful as a smooth low-order comparison, while `examples/variational_analysis/radial_convergence.m` checks a milder `Ns=1`, $\epsilon=0.2$ problem. Neither test exercises the $\Delta$--$P$ refinement pathology exposed by the two `variational_vs_legacy` profiles. For this reason, the latter script is the relevant regression for solver robustness as well as for profile agreement.

## Straight-field-line geometry

For the isotropic comparison cases $a_0=1$ and $\sigma=0$, so the variational expression reduces to
\[
\frac{\partial\vartheta}{\partial\omega}
=\frac{J}{\epsilon^2\hat r R^2}.
\]
The legacy code evaluates exactly the same quantity as `JoverR/RR/epsilon^2/r`. Consequently, `theta_SFL` is not a new independent equilibrium discrepancy in these cases. It is determined by $(\Delta,P,S_m)$ and their first derivatives through $J$, and then affects `RR_sfl` and `ZZ_sfl` through inversion of $\vartheta(\omega)$. The $P$ discrepancy is especially relevant because $P$ is the coordinate-fixing profile.

There is one algorithmic difference after integration. The variational postprocessor records
\[
e_\vartheta=2\pi\left\langle
\frac{J}{\epsilon^2\hat rR^2}
\right\rangle_\omega-2\pi
\]
as `theta_SFL_periodicity_error` and rescales every non-axis row so that its integrated endpoint is exactly $2\pi$. The legacy postprocessor does not rescale. The weak $P$ equation makes the gauge error orthogonal to the discrete test space; it does not make its pointwise value identically zero. At `m=12`, the largest unrescaled period errors are about $2.5\times10^{-3}$ rad in case 1 and $2.0\times10^{-2}$ rad in case 2 for the variational solution, compared with $2.4\times10^{-3}$ and $3.2\times10^{-3}$ for the legacy solution.

After interpolation to the variational grids, the maximum differences are approximately
\[
\begin{array}{c|ccc}
 & |\delta R_{\rm SFL}|_{\max} & |\delta Z_{\rm SFL}|_{\max}
 & |\delta\vartheta|_{\max} \\
\hline
\text{case 1} & 1.29\times10^{-3} & 3.67\times10^{-3} & 9.06\times10^{-3} \\
\text{case 2} & 4.94\times10^{-3} & 5.63\times10^{-3} & 6.67\times10^{-2}
\end{array}
\]
in units normalized to $R_0$ for the geometry. The large $Z_{\rm SFL}$ maxima occur at the last legacy row and are amplified by its copied, rather than evaluated, free-profile edge values. Excluding both endpoint rows reduces the $Z_{\rm SFL}$ maxima to $7.1\times10^{-4}$ and $1.4\times10^{-3}$. Applying versus omitting only the variational $2\pi$ rescaling changes $Z_{\rm SFL}$ by about $6.3\times10^{-4}$ in both cases. Thus the SFL discrepancy has three separable sources: the equilibrium profiles, the weak periodicity error and rescaling convention, and the legacy edge padding.

## Fields entering the VENUS interfaces

The two writer functions do not implement the same schema or normalization. `equil_to_venus_variational` writes dimensionless geometry and profiles, their selected radial derivatives, and metadata. `to_venus_cleaner` inserts fixed dimensional values of $(R_0,B_0)$ and writes additional constant rotation and thermodynamic arrays. Those differences should be removed algebraically before comparing files.

The version-7 variational writer always exports the native total EOS partial
derivatives as `/anisotropy/Dparallel` and `/anisotropy/Dperp`. These are
$(\partial_rP_\parallel)_{R,B}$ and $(\partial_rP_\perp)_{R,B}$ and require no
thermal/species bookkeeping. They support the two-variable unsplit
perpendicular anisotropic VENUS model for every anisotropic EQUIL equation of
state.

Version 7 also exports the native SFL metric derivative
`/tensors/dg22dtheta`. It is required by the true perpendicular Galerkin
operator in VENUS and is assembled analytically as
\(2(R_\theta R_{\theta\theta}+Z_\theta Z_{\theta\theta})\), rather than
reconstructed by differentiating the HDF5 samples.

For the explicitly decomposed `thermal_and_deBlank` closure, the writer also
writes the collisional thermal pressure as the one-dimensional `/profiles/P`
field and subtracts its native derivative from the total EOS partial
derivatives before exporting `/anisotropy/Dparallel_A` and
`/anisotropy/Dperp_A`. Consequently these optional compatibility fields are
exactly $(\partial_r\Pi_{\parallel A})_{R,B}$ and
$(\partial_r\Pi_{\perp A})_{R,B}$; VENUS does not differentiate pressure
samples or follow the condensed-$B$ solution curve. The same file contains
the total pressure maps, contravariant curvature and current, and the native
SFL $(B^2)_r$. An anisotropic EOS without an explicit decomposition omits only
the species-A datasets and is rejected only by the older split stability
model.

The common equilibrium-dependent content is much smaller. Both writers use `RR_sfl` and `ZZ_sfl`; hence every discrepancy in $(\Delta,P,S_m)$ enters through the SFL geometry. For the anisotropic SFL representation define the auxiliary flux profile
\[
F_{\rm aux}=-\frac{a_0+\epsilon^2t_2}{1-\sigma_{\rm cyl}},
\qquad
g=\frac{\epsilon^2F_{\rm aux}}{q}.
\]
It gives the native field through $B^\vartheta=g/J_a$ and $B^\phi=qB^\vartheta$. In anisotropy, $R^2B^\phi=T/(1-\sigma)$ is
two-dimensional and is not the profile $F_{\rm aux}$; the two coincide only in isotropy. The writer exports native `Btheta` and `Bphi` as an independent cross-check and rejects a profile/native mismatch above $10^{-7}$. It also writes the analytical derivatives of $F_{\rm aux}$ and $g$. The safety factor $q$, its derivative, and the pressure are prescribed inputs and are not independent solved discrepancies when the same `LX` functions are used. The remaining density, gamma, temperature, and rotation arrays are interface choices rather than consequences of $(\Delta,t_2,P,S_m)$.

The uniform-$\vartheta$ remap interpolates primitive SFL fields separately.
After that interpolation, $J_a$ is reconstructed from the remapped $R$ and
$\sigma$.  Independently retaining the interpolated $B^\vartheta$ would
therefore evaluate the nonlinear product as
$J_{a,{\rm remap}}I(1/J_a)g$, which is not exactly $g$.  This produced a
maximum $B^\vartheta$ mismatch of about $6.7\times10^{-7}$ for an
`axis_regular`, `thermal_and_deBlank`, $N_s=10$, $\epsilon=0.4$ case even
though the identity held to roundoff before remapping.  The postprocessor now
reconstructs
\[
B^\vartheta=\frac{\epsilon^2T}
 {q(1-\sigma_{\rm cyl})J_a},\qquad
B^\phi=qB^\vartheta
\]
after reconstructing $J_a$.  This preserves both exported profile identities
and the anisotropic pointwise relation
$R^2B^\phi=T/(1-\sigma)$ without identifying that quantity with
$F_{\rm aux}$.

The surface-wise normalization of the weak SFL angle to a period of $2\pi$
can nevertheless leave a finite-resolution difference between the remapped
scalar $B^2$ and the norm reconstructed from the remapped metric and
contravariant field. This measures the pointwise defect of the weak gauge
representation; it is not a second residual of the converged equilibrium.
VENUS records the normalized difference as `magnetic_identity_error` and
warns when it exceeds $10^{-2}$, but does not reject the case. A large value
should still be treated as a resolution diagnostic, especially for a case
intended for quantitative use near marginality.

At the comparison resolution, converting the legacy $F$ and $g$ to the variational dimensionless normalization gives maximum differences of about $(6.2\times10^{-5},3.7\times10^{-6})$ in case 1 and $(4.7\times10^{-4},3.6\times10^{-5})$ in case 2. Geometry, particularly its SFL remapping and edge row, is therefore the larger interface discrepancy in these tests. A VENUS comparison should use a common physical normalization, discard schema-only fields, compare the open radial interval first, and then inspect the edge values separately.

## Diagnostic procedure

A refinement or formulation comparison should first fix $(q,\beta,\epsilon,N_s)$ and state `spline_p` independently of `nq`. It should report global degrees of freedom, quadrature sizes, direct and continued Newton histories, Jacobian construction time, the selected residual and Jacobian functions, the smallest singular direction by profile block, `min_J_over_eps2r`, `gauge_error`, and `theta_SFL_periodicity_error`. `jacobian_step` is relevant only to the explicitly selected finite-difference and column-reference utilities. A transferred mesh solve must be compared both with its preceding mesh and with the legacy solution; convergence of the algebraic residual alone is insufficient.

Profile errors should be evaluated on a common radial grid over the open interval, with exact-axis and exact-edge values reported separately. SFL comparisons should distinguish unrescaled $\vartheta$, the explicit $2\pi$ normalization, and the subsequent interpolation of $(R,Z)$. Finally, the two VENUS files should be reduced to common dimensionless $(R_{\rm SFL},Z_{\rm SFL},F,g,q,P)$ fields before any remaining discrepancy is attributed to the equilibrium formulation.
