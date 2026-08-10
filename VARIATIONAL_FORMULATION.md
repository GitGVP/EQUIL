# Variational EQUIL formulation

This document describes the normalized formulation implemented by
`equilVariationalSol`. It is separate from the generated strong-form
residuals used by `equilSol`.

## Normalized fields, geometry, and action

The global radial state is

\[
x=\left(t_2,\Delta,\mathsf P,S_2,\ldots,S_{N_s+1},
V_2,\ldots,V_{N_h+1}\right).
\]

Here \(\mathsf P\) is the straight-field-line coordinate gauge, not a
thermodynamic pressure. The optional \(V_n\) profiles describe up-down
asymmetry; `Nh=0` omits those blocks entirely. Magnetic-field strength is
not a global unknown. It is eliminated independently at each
\((r,\omega)\) quadrature point.

With \(\hat r\in[0,1]\), the normalized mapping is

\[
\begin{aligned}
\hat R={}&1+\epsilon\hat r\cos\omega-\epsilon^2\hat\Delta
 +\epsilon^2\sum_{n=1}^{N_s}\hat S_{n+1}\cos(n\omega)
 +\epsilon^2\sum_{n=1}^{N_h}\hat V_{n+1}\sin(n\omega)
 +\epsilon^3\hat{\mathsf P}\cos\omega,\\
\hat Z={}&\epsilon\hat r\sin\omega
 -\epsilon^2\sum_{n=1}^{N_s}\hat S_{n+1}\sin(n\omega)
 -\epsilon^2\sum_{n=1}^{N_h}\hat V_{n+1}\cos(n\omega)
 +\epsilon^3\hat{\mathsf P}\sin\omega.
\end{aligned}
\]

With \(\epsilon=a/R_0\). Hereafter we drop the hats for convenience.
Define

\[
J=R(R_rZ_\omega-R_\omega Z_r),\qquad
g_{\omega\omega}=R_\omega^2+Z_\omega^2.
\]

The normalized legacy equations of state return quantities on the beta
scale. Internally the action uses

\[
\Pi_\parallel=\frac{\mu_0p_\parallel}{B_0^2}
=\epsilon^2\beta_\parallel.
\]

Consequently every pressure value and pressure derivative supplied by a
legacy equation of state is multiplied by \(\epsilon^2\). A native
`pressure_map` instead supplies \(\Pi_\parallel\) and its derivatives
directly.

The anisotropy and on-axis normalization are

\[
\sigma=\frac{\Pi_{\parallel,B}}{B},\qquad
\sigma_0=\Pi_{\parallel,B}(0,1,1),
\]

and

\[
T=1-\sigma_0+\epsilon^2t_2.
\]

For the prescribed safety factor,

\[
\boxed{\psi_r=\frac{\epsilon rT}{q(r)(1-\sigma_0)}}.
\]

In the code, `a0 = 1-sigma0` is used.

Because derivatives in the code are taken with respect to normalized
\(r\), the poloidal-field term is evaluated as

\[
B_p^2=\frac{(\epsilon\psi_r)^2g_{\omega\omega}}{J^2},
\qquad
B_\phi=\sqrt{B^2-B_p^2}.
\]

The normalized action density and action are

\[
\boxed{\mathcal L=\frac{B^2}{2}-\Pi_\parallel-\frac{T}{R}B_\phi},
\qquad
\mathcal A=\int dr\,d\omega\,J\mathcal L.
\]

## Local magnetic-field elimination

At each two-dimensional quadrature point the scalar equation

\[
G(B)=1-\sigma(r,R,B)-\frac{T}{RB_\phi}=0
\]

is solved locally. Its Newton derivative is

\[
G_B=-\sigma_B+\frac{TB}{RB_\phi^3},\qquad
\sigma_B=\frac{\Pi_{\parallel,BB}}{B}
-\frac{\Pi_{\parallel,B}}{B^2}.
\]

For `isotropic` and `isotropic_rotating`, \(\sigma=0\) and the exact
solution is used:

\[
B=\sqrt{B_p^2+\left(\frac{T}{R}\right)^2}.
\]

Thus an isotropic calculation cannot be limited by a numerical local-
\(B\) residual. Anisotropic closures retain the safeguarded scalar Newton
iteration and enforce \(B^2>B_p^2\), \(1-\sigma>0\), and non-small
\(|G_B|\).

## Weak residuals

Write \(J_s=J/\epsilon\). After dropping common constants, the weak flux
equation for a radial test function \(v\) is

\[
\begin{aligned}
R_{t_2}[v]=\int dr\,d\omega\,\bigg[
&\frac{(1-\sigma)g_{\omega\omega}}{J_s}\,\psi_rv_r\\
&-\left(
\frac{J_s\Pi_{\parallel,r}}{\psi_r}
+\frac{J_sTT_r}{R^2(1-\sigma)\psi_r}
\right)v\bigg].
\end{aligned}
\]

The partial derivative \(\Pi_{\parallel,r}\) is taken at fixed \(R,B\).
Only \(t_2\) and \(t_{2,r}\) appear. There is no derivative of the
prescribed factor \(\epsilon r/[q(1-\sigma_0)]\): \(\psi_r\) occurs
algebraically in this weak equation. A derivative of that factor would be
introduced only by returning to the strong radial equation and expanding
its derivative. In the first term the dependence on \(T\) is carried by
\(\psi_r\). The explicit \(T\) multiplying \(T_r\) in the final term is
retained because it comes from the action term \(-TB_\phi/R\), not from
re-expanding \(\psi_r\).

For a geometrical variation \((\delta R,\delta Z)\),

\[
\begin{aligned}
\delta J={}&(R_rZ_\omega-R_\omega Z_r)\delta R+R\big(
Z_\omega\delta R_r+R_r\delta Z_\omega
-Z_r\delta R_\omega-R_\omega\delta Z_r\big),\\
\delta g_{\omega\omega}={}&2R_\omega\delta R_\omega
+2Z_\omega\delta Z_\omega,\\
\delta B_p^2={}&(\epsilon\psi_r)^2
\left(\frac{\delta g_{\omega\omega}}{J^2}
-\frac{2g_{\omega\omega}\,\delta J}{J^3}\right),\\
\delta B_\phi={}&-\frac{\delta B_p^2}{2B_\phi}.
\end{aligned}
\]

At fixed independent \(B\), stationarity of the locally reduced action
gives

\[
\delta(J\mathcal L)=
\mathcal L\,\delta J
+J\left(-\Pi_{\parallel,R}+\frac{TB_\phi}{R^2}\right)\delta R
-\frac{JT}{R}\delta B_\phi.
\]

The shift variation uses

\[
\delta R=-v,\qquad\delta Z=0,
\]

and harmonic \(S_{n+1}\) uses

\[
\delta R=v\cos(n\omega),\qquad
\delta Z=-v\sin(n\omega).
\]

The up-down asymmetric harmonic \(V_{n+1}\) uses

\[
\delta R=v\sin(n\omega),\qquad
\delta Z=-v\cos(n\omega).
\]

These compact variations are assembled directly. No expanded Mathematica
residual is used, and no \(\Delta_{rr}\), \(S_{m,rr}\), or
\(\mathsf P_{rr}\) enters. Natural outer-edge terms are retained for
profiles whose test space does not vanish at \(r=1\). The constant
toroidal-vacuum contribution is removed analytically before quadrature to
avoid cancellation when \(\epsilon\ll1\).

The gauge profile is determined by exact SFL periodicity:

\[
C_{\rm SFL}(r)=
\left\langle\frac{J}{R^2(1-\sigma)}\right\rangle_\omega
-\frac{\epsilon^2r}{1-\sigma_0}=0.
\]

The code divides this equation by \(\epsilon^2\), so its assembled kernel
is

\[
\frac{1}{\epsilon^2}
\left\langle\frac{J}{R^2(1-\sigma)}\right\rangle_\omega
-\frac{r}{1-\sigma_0}.
\]

All residual blocks are additionally scaled by known powers of \(\epsilon\).
This changes conditioning but not the root.

## Optional up-down asymmetric profiles

The generalized mapping follows equations (19) and (20) of Fitzpatrick,
*Physics of Plasmas* **31**, 082505 (2024). His \(\hat H_n\) profiles map
to the symmetric \(S_n\) profiles above, his relabeling function
\(\hat L\) maps to \(\mathsf P\), and his \(\hat V_n\) are the new
asymmetric profiles. The existing EQUIL mapping uses different
poloidal-angle and profile-sign conventions; after translating those
conventions, the code defines each \(V_n\) by the sine/cosine pair shown
in the generalized mapping.

At leading order in aspect ratio, equation (34) of that paper gives

\[
V_n''=-\left(2\frac{f_1'}{f_1}-\frac1r\right)V_n'
      +(n^2-1)\frac{V_n}{r^2},\qquad n>1,
\]

with the axis condition from its equation (40),

\[
V_n(r)=V_{nc}r^{n-1}+o(r^{n-1}).
\]

The \(n=1\) solution is excluded because a nonzero constant would move
the magnetic axis. In the finite-\(\epsilon\) implementation the ODE is
not inserted as a separate strong residual. The full weak equation is
instead obtained from the same action variation as every other geometry
profile, using the \(\delta R,\delta Z\) pair above. Thus it includes the
finite-aspect-ratio, flow, and anisotropic terms automatically and still
contains only first radial derivatives.

The code option `Nh` counts \(V_2,\ldots,V_{Nh+1}\), while `LX.Vbc`
prescribes their values at \(r=1\). Each radial space is

\[
V_{n+1}(r)=r^n(1-r)v(r)+V_{n+1}(1)r^n.
\]

Consequently the axis regularity and outer boundary value are exact, and
the unknown test functions vanish at the edge. No natural edge term is
therefore added to a \(V_n\) row. `LY.V` and `LY.Vp` contain the solved
profiles. `examples/variational_asymmetric_case.m` demonstrates two small
nonzero asymmetric boundary harmonics, ramped from the symmetric state.

## Condensed Jacobian

For any global perturbation \(x_j\), local implicit differentiation gives

\[
B_{x_j}=-\frac{G_{x_j}}{G_B}.
\]

The implementation evaluates fixed-\(B\) directional derivatives and then
forms the reduced response along \((x_j,B_{x_j})\). This is the scalar
static-condensation or Schur-complement formula

\[
R^{\rm red}_x=R_x-R_BG_B^{-1}G_x.
\]

There are therefore no global \(B_m\) rows or columns. The differentiation
step is reduced near geometry or mirror constraints, but the nonlinear
solver intentionally remains a short proof-of-principle Newton/backtracking
method.

## Axis-regular B-spline space

Smoothness in local Cartesian coordinates implies the magnetic-axis
Taylor class

\[
F_n(r)=r^{|n|}\left(c_0+c_1r^2+c_2r^4+\cdots\right).
\]

This is a local regularity statement, not a claim that an equilibrium
profile should be represented by one global polynomial in \(r^2\). The
implemented radial space therefore uses open, uniform B-splines in \(r\).
Only the first spline span is constrained to have the required Taylor
parity:

\[
t_2,\Delta=O(r^2),\qquad
\mathsf P=p_1r+O(r^3),\qquad
S_{n+1},V_{n+1}=O(r^n).
\]

For example, a degree-four first span contains \(r^2,r^4\) for \(t_2\)
and \(\Delta\), \(r,r^3\) for \(\mathsf P\), and \(r,r^3\) for \(S_2\).
The \(p_1r\) term is the radial-normalization gauge associated with the
globally imposed SFL constraint. The physical first-harmonic deformation
left after fixing that normalization begins at \(r^3\). Omitting \(p_1\)
forces the other axis coefficients to compensate and is the main cause of
the previously observed finite-beta disagreement with legacy EQUIL.
Beyond the first span the basis is an ordinary B-spline space. The same
constrained functions are used as Galerkin test functions, and the outer
values of every prescribed shaping profile are imposed by eliminating the
last spline coefficient.

There is one variational radial implementation; the former global
polynomial and parity switches are not solver options. `m` is the number
of radial spans, `spline_p` the B-spline degree, and `nq` the Gauss points
per span. Thus changing `m` is genuine \(h\)-refinement. The default
degree is four because it is the cheapest space that resolves the
high-beta legacy comparison; degree five gives essentially the same branch
and profile errors with four additional unknowns.

An axis-compatible pedestal must likewise be prescribed as a smooth
function of \(s=r^2\), for example

\[
\beta(r)=f(r^2),\qquad \beta_r=2r f'(r^2).
\]

## Difference from standard EQUIL

| | `equilSol` | `equilVariationalSol` |
|---|---|---|
| Global fields | \(t_2,\Delta,\mathsf P,S_m,B_m\) | \(t_2,\Delta,\mathsf P,S_m\), optional \(V_n\) |
| Magnetic strength | Fourier harmonics in the radial FE basis | Local scalar solve at \((r_q,\omega_p)\) |
| Residual | Poloidal Fourier moments of strong force balance | Weak flux equation, action variations, SFL constraint |
| Radial derivatives | Strong residuals include second derivatives | First derivatives only |
| Jacobian source | Generated residual/Jacobian files | Direct condensed differentiation |
| Default radial basis | B-splines | B-splines with a parity-constrained first span |
| Initialization path | User-managed restarts | Low-beta seed and adaptive fixed-\(\epsilon\) continuation |
| Postprocessing | Strong-form currents and integral diagnostics | Plot-ready geometry, pressure, local \(B\), flux and shape fields |

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
| `spline_p` | 4 | B-spline polynomial degree |
| `nq` | 6 | Gauss points per radial span |
| `om_pts` | 96 | periodic angular points |
| `NLtol` | \(10^{-10}\) | global residual tolerance |
| `nk` | 30 | nominal Newton iterations |

With `Ns=3` this gives 50 global radial unknowns. The first-span parity
constraints make \(S_2=O(r)\), \(S_3=O(r^2)\), and \(S_4=O(r^3)\)
independently of quadrature. Increasing `nq` or `om_pts` therefore does
not repair an axis-space error; `m` and `spline_p` control radial
approximation, `nq` controls radial integration, and `om_pts` controls
angular integration.

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

The matching axis-refinement test gives:

| spans | total DOFs | fitted near-axis power of \(S_2\) | relative \(S_2\) error from legacy |
|---:|---:|---:|---:|
| 4 | 26 | 0.9983 | 5.3% |
| 6 | 38 | 0.9988 | 2.1% |
| 8 | 50 | 0.9973 | 1.3% |

Thus the correct linear axis trace is already present in the smallest
space; refinement improves the finite-radius profile rather than repairing
an incorrect axis exponent.

## Newton and beta continuation

Continuation has no schedule or trial-budget options. If no initial
`LX.x` is supplied and beta is nonzero, it uses:

1. beta zero at the requested \(\epsilon\);
2. a moderate nonzero seed
   \[
   \beta_{\rm low}=\operatorname{sign}(\beta)\min(|\beta|,1);
   \]
3. the requested beta.

A failed jump is bisected from the last accepted beta, and the target is
retried immediately after every accepted midpoint. Thus routine beta-20
runs need three accepted states, not a long prescribed ramp.

For strongly shaped boundaries the analytic beta-zero initial guess can
fail directly at finite aspect ratio. Only in that case the solver obtains
a beta-zero seed at \(\epsilon=10^{-3}\) and retries the requested
\(\epsilon\) in one Newton solve. This is a fallback initialization, not an
epsilon continuation scan. The attempted beta and epsilon values are
returned in `LY.beta_continuation_attempts` and
`LY.beta_continuation_attempt_epsilons`.

A continuation attempt is rejected early only when eight successive
residuals reduce by less than 3% while remaining more than \(10^4\) times the
tolerance. Monotonically decreasing histories are allowed to extend beyond
`nk` when their logarithmic slope predicts convergence, up to the hard
limit \(3\,nk\). This separates a residual stuck near \(10^{-1}\) from a
slow attempt that is already close to the requested tolerance. A failed
line search is rejected immediately.

The continuation regression deliberately sets `nk=10`: the beta-20 final
stage needs 13 iterations and is extended successfully, while both beta 3
and beta 20 still use only three accepted pressure states and no rejected
jumps.

The continuation trials always disable `do_ana` and `do_shift_NLO`.
If requested, the analytical calculation is performed once after the final
accepted equilibrium. Supplying `LX.x` bypasses beta continuation and is
the preferred mechanism for a physics scan.

## Aspect-ratio convergence

`examples/variational_analysis/epsilon_convergence.m` evaluates
`equil_ana`, including the optional NLO shift solve, exactly once. It
then warm-starts the numerical states from \(\epsilon=0.4\) down to
\(10^{-3}\), carrying both the global coefficients and the local magnetic
field.

The expected pre-plateau behavior is

\[
\Delta-\Delta_0=O(\epsilon),\qquad
\Delta-\Delta_0-\epsilon\Delta_1=O(\epsilon^2),
\]

with the same orders for \(S_2\), while the first omitted isotropic
correction to \(t_{2,0}\) is \(O(\epsilon^2)\). A reduced six-point check
with the documented scan settings gave LO slopes 1.06 for \(\Delta\) and
0.98 for \(S_2\), with global residuals between \(7\times10^{-16}\) and
\(1.3\times10^{-13}\). The script also forms an independent numerical
\(\epsilon\to0\) reference so that an `equil_ana` integration floor can
be distinguished from a different NLO target.

The convergence plots intentionally use the Euclidean norm of the sampled
radial values, matching the article scripts and their vertical scale.
Quadrature-weighted norms change the ordinate but not the asymptotic order.
