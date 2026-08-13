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
and our approach is to assume that the poloidal flux is given by:
\[
\psi_r=\frac{r}{R_0}\frac{T}{q}\frac{1}{(1-\sigma_{\parallel,0})}
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
\frac{\partial \vartheta}{\partial \omega} = \frac{\hat J/\epsilon^2}{\hat r \hat R^2}\frac{1-\sigma_{\parallel,0}}{1-\sigma_\parallel}
\]
So that the requirement that $\vartheta(r,2\pi)-\vartheta(r,0)=2\pi$ means that:
\[
  \left\langle\frac{\hat J/\epsilon^2}{\hat R^2(1-\sigma_\parallel)}\right\rangle_\omega=\frac{\hat r}{1-\sigma_{\parallel,0}}
\]
which preserves the near axis SFL jacobian $J_\mathrm{SFL}\sim rR^2/R_0$ but leads to a renormalization of $\epsilon_\mathrm{eff} = \epsilon + \epsilon^3 \hat P(\hat r=1;\sigma_\parallel)$ which has the undesirable property that it is anisotropy dependent. The equations for $(\hat \Delta, \hat \Delta_Z, \hat S_m, \hat A_m)$ are obtained from the geometrical variation of the last section, and more explicitly:
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
\frac{(1-\sigma_{\parallel,0})J}
{\epsilon^2 rR^2(1-\sigma_\parallel)},
\]
and inverts it on a uniform periodic $\vartheta$ grid. The resulting fields are `theta_SFL`, `RR_sfl`, and `ZZ_sfl`. No stability-specific quantities from the legacy `equil_SFL` routine are evaluated.

The variational postprocessor also returns the legacy-normalized global quantities `Wk`, `Wkpar`, `Wkperp`, `Wkrot`, `Wp`, `Ip`, `Ftt`, `Ft0`, `Ft`, `bp`, `bppar`, `bpperp`, `bprot`, `li`, and `bpli2`, together with the compatibility aliases `BBp2`, `BBt2`, and `BB2`. These are evaluated directly on the variational Gauss grid. The remaining legacy-only fields are deliberately not copied: `Bs` and `Bsp` belonged to the old global Fourier representation of $B$; `jphi` and the contravariant current components require a separate weak reconstruction; higher radial derivatives such as `t2pp`, `deltapp`, `deltappp`, `Ppp`, `Pppp`, and `Spp` are only used by the old stability postprocessing. The stability-specific moments computed by `equil_SFL` are likewise outside the equilibrium interface requested here.



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

## Radial discretization: Axis-regular B-spline space

Near axis regularity in local Cartesian coordinates implies that $S_m,A_m\sim r^{|n|}(c_0+c_1r^2+c_2r^4+\dots)$. The implemented radial space therefore uses open, uniform B-splines in \(r\) with the first spline span constrained to have the required Taylor parity:
\[
t_2,\Delta=\mathcal O(\hat r^2),\qquad  \hat P=p_1\hat r+\mathcal O(\hat r^3),\qquad \hat S_{n+1},\hat A_{n+1}=\mathcal O \hat(r^n).
\]
For example, a degree-four first span contains \(r^2,r^4\) for \(t_2\) and \(\Delta\), \(r,r^3\) for \(\mathsf P\), and \(r,r^3\) for \(S_2\).
`m` is the number of radial spans, `spline_p` the B-spline degree, and `nq` the Gauss points per span. The requested default degree is $4$; internally it is raised to $\max(4,N_s,p_A)$ when the retained axis Taylor powers require it.

## Difference from standard EQUIL

| | `equilSol` | `equilVariationalSol` |
|---|---|---|
| Global fields | \(t_2,\Delta,\mathsf P,S_m,B_m\) | \(t_2,\Delta,\mathsf P,S_m\), optional \(A_m,\Delta_Z\) |
| Magnetic strength | Fourier harmonics in the radial FE basis | Local scalar solve at \((r_q,\omega_p)\) |
| Residual | Poloidal Fourier moments of strong force balance | Weak flux equation, action variations, SFL constraint |
| Radial derivatives | Strong residuals include second derivatives | First derivatives only |
| Jacobian source | Generated residual/Jacobian files | Direct condensed differentiation |
| Default radial basis | B-splines | B-splines with a parity-constrained first span |
| Initialization path | User-managed restarts | Low-beta seed and adaptive fixed-\(\epsilon\) continuation |
| Postprocessing | Strong-form currents, stability-specific SFL quantities, and integral diagnostics | Geometry, pressure, local \(B\), integral diagnostics, and the SFL angle/remapped geometry |

The public organization remains parallel: `equilVariationalSol` builds the
discretization, `equilVariationalX` supplies normalized inputs,
`equilVariationalY` solves, and `equilVariationalPP` constructs plotting
fields.

## Verification scripts

The focused checks are:

- `examples/variational_simple_case.m` solves the same isotropic beta-20
  case with both formulations, prints profile errors and compares the flux
  surfaces.
- `Analysis/test_variational_s2_axis_refinement.m` checks the enforced
  linear axis behavior of \(S_2\) under B-spline \(h\)-refinement.
- `Analysis/test_variational_beta_continuation.m` checks the accepted
  beta stages and rejected-attempt count.
- `examples/variational_analysis/radial_convergence.m` performs B-spline
  \(h\)-refinement for a smooth case and optionally for an axis-compatible
  pedestal.
- `examples/variational_analysis/epsilon_convergence.m` compares with the
  unchanged `equil_ana` LO and NLO profiles and reproduces the layout of
  the article's aspect-ratio plots.

The asymmetric and strong-anisotropy scripts remain examples, but they are
not part of the routine regression set.

## Numerical parameters and defaults

The routine defaults are deliberately moderate:

| Option | Default | Meaning |
|---|---:|---|
| `m` | 8 | uniform radial B-spline spans |
| `spline_p` | 4 | requested B-spline degree; raised automatically for axis regularity |
| `nq` | 6 | Gauss points per radial span |
| `om_pts` | 96 | periodic angular points |
| `Na` | 0 | number of asymmetric quadrature profiles |
| `A_modes` | `[]` | selected physical modes; defaults to `2:(Na+1)` |
| `A_leading_powers` | `[]` | optional higher axis powers for selected `A_m` profiles |
| `vertical_shift` | `false` | add the fixed-edge vertical-center profile |
| `epsilon_match` | `'boundary'` | MEQ benchmark scale update (`'volume'` selects legacy `Q2Q`) |
| `NLtol` | \(10^{-10}\) | global residual tolerance |
| `nk` | 30 | nominal Newton iterations |

Useful settings are:

| Task | `m` | `spline_p` | `nq` | `om_pts` | `NLtol` |
|---|---:|---:|---:|---:|---:|
| routine isotropic or rotating solve | 8 | 4 | 6 | 96 | \(10^{-10}\) |
| beta-20 legacy comparison | 8 | 4 | 6 | 96 | \(10^{-10}\) |
| aspect-ratio scan to \(10^{-3}\) | 8 | 4 | 8 | 128 | \(3\times10^{-13}\) |
| smooth radial convergence study | 3--14 | 4 | 6 | 64 | \(10^{-11}\) |

A smooth degree-\(p\) Galerkin problem is expected to approach
\(O(h^{p+1})\) in \(L^2\) and \(O(h^p)\) in \(H^1\) before nonlinear,
quadrature, or roundoff errors dominate. The radial-convergence script
prints the measured slopes instead of assuming that the asymptotic range
has been reached. With degree four and \(m=3,\ldots,10\), the smooth test
measured 4.81 for the \(t_2\) \(L^2\) error and 4.91 for its derivative
error. The \(\Delta\) sequence was still non-monotone and pre-asymptotic;
the script retains those curves as a diagnostic rather than presenting
their fitted slope as an achieved Galerkin order.

The beta-20 comparison at \(\epsilon=0.4\), zero shaping boundary values,
and the default radial space converges in the accepted beta stages
\([0,1,20]\). A representative final residual is below \(4\times10^{-11}\).
On \(0.01\le r\le0.99\), representative relative differences from the
legacy `iso_static` result are about 0.15% for \(t_2\), 0.39% for
\(\Delta\), 0.93% for \(\mathsf P\), 1.3% for \(S_2\), and 6.4% for
\(S_3\). At beta zero the corresponding absolute agreement is much
closer. Degree five gives nearly identical errors. The remaining relative
\(S_3\) difference corresponds to a maximum absolute difference below
\(9\times10^{-3}\).

## Newton and beta continuation

Continuation has no schedule or trial-budget options. If no initial `LX.x` is supplied and beta is nonzero, it uses:
1. beta zero at the requested \(\epsilon\);
2. a moderate nonzero seed
   \[
   \beta_{\rm low}=\operatorname{sign}(\beta)\min(|\beta|,1);
   \]
3. the requested beta.
A failed jump is bisected from the last accepted beta, and the target is retried immediately after every accepted midpoint. Thus routine beta-20 runs need three accepted states.

For strongly shaped boundaries the analytic beta-zero initial guess can fail directly at finite aspect ratio. Only in that case the solver obtains a beta-zero seed at \(\epsilon=10^{-3}\) and retries the requested \(\epsilon\) in one Newton solve. This is a fallback initialization, not an epsilon continuation scan. The attempted beta and epsilon values are returned in `LY.beta_continuation_attempts` and `LY.beta_continuation_attempt_epsilons`.

A continuation attempt is rejected early only when eight successive residuals reduce by less than 3% while remaining more than \(10^4\) times the tolerance. Monotonically decreasing histories are allowed to extend beyond `nk` when their logarithmic slope predicts convergence, up to the hard limit \(3\,nk\). This separates a residual stuck near \(10^{-1}\) from a slow attempt that is already close to the requested tolerance. A failed line search is rejected immediately.

The continuation regression deliberately sets `nk=10`: the beta-20 final stage needs 13 iterations and is extended successfully, while both beta 3 and beta 20 still use only three accepted pressure states and no rejected jumps.

The continuation trials always disable `do_ana` and `do_shift_NLO`. If requested, the analytical calculation is performed once after the final accepted equilibrium. Supplying `LX.x` bypasses beta continuation and is the preferred mechanism for a physics scan.

## Aspect-ratio convergence

`examples/variational_analysis/epsilon_convergence.m` evaluates `equil_ana`, including the optional NLO shift solve, exactly once. It then warm-starts the numerical states from \(\epsilon=0.4\) down to
\(10^{-3}\), carrying both the global coefficients and the local magnetic field.

The NLO shift is the most radially sensitive quantity in this scan. The example therefore uses 16 radial spans and spline degree six. With only eight spans, the fixed-space error in the numerical $\epsilon\to0$ coefficient is about four percent; subtracting the analytic first-order term then hides the expected $O(\epsilon^2)$ remainder even though the nonlinear solve is converged.

The expected pre-plateau behavior is

\[
\Delta-\Delta_0=O(\epsilon),\qquad
\Delta-\Delta_0-\epsilon\Delta_1=O(\epsilon^2),
\]

with the same orders for \(S_2\), while the first omitted isotropic correction to \(t_{2,0}\) is \(O(\epsilon^2)\). A reduced six-point check with the documented scan settings gave LO slopes 1.06 for \(\Delta\) and
0.98 for \(S_2\), with global residuals between \(7\times10^{-16}\) and \(1.3\times10^{-13}\). The script also forms an independent numerical \(\epsilon\to0\) reference so that an `equil_ana` integration floor can
be distinguished from a different NLO target.

The convergence plots intentionally use the Euclidean norm of the sampled
radial values, matching the article scripts and their vertical scale.
Quadrature-weighted norms change the ordinate but not the asymptotic order.
