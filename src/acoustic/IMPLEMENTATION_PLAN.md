# Two-dimensional mixed-grid acoustic solver: implementation plan

## Objective

Implement the node-centered finite-volume formulation of `textbook/Chap5-4-acoustic.md` for linear triangles and quadrilaterals in the same mesh. The solver will solve

\[
\frac{1}{\rho_0c_0^2}p'_{tt}
=\nabla\cdot\left(\frac{1}{\rho_0}\nabla p'\right)
-\frac{2}{\rho_0c_0}\frac{\partial(\mathbf M\cdot\nabla p')}{\partial t}
-\frac{1}{\rho_0}\mathbf M\cdot\nabla(\mathbf M\cdot\nabla p')+\mathcal S,
\]

for a steady, spatially uniform, subsonic background in the initial implementation. Cell-wise storage will be used for density, sound speed, and Mach number so the assembly follows the chapter and can later admit non-uniform background fields.

## 1. Mesh and element formulation

- Store acoustic pressure and its time derivatives at mesh nodes; retain the existing cell-centered flow fields unchanged.
- Construct each node's median-dual control area from the node, adjacent edge midpoints, and cell centroid. Sum the exact sub-area of every triangle and quadrilateral and verify that element sub-areas sum to element area.
- Evaluate triangle gradients exactly. Evaluate bilinear quadrilateral gradients at flux-integration points through the isoparametric Jacobian, including inverse mapping from physical to reference coordinates.
- Check `sum(grad N) = 0` at every integration point and write per-cell geometry diagnostics.
- Form the internal dual-face flux with the full tensor `A = I - M ⊗ M`; do not drop mixed derivatives when mean flow is present.

## 2. Discrete operators and physics

Assemble the global system

\[
\mathsf M\ddot p+\mathsf C\dot p+\mathsf Kp=F(t).
\]

- `M`: diagonal lumped acoustic mass `S_k,e / (rho_e c_e^2)`.
- `C`: the centered mixed space-time mean-flow term from Eq. (5.235), plus absorbing-boundary damping.
- `K`: minus the divergence form of `(I - M ⊗ M) grad(p) / rho`, with one-point integration on each straight internal dual-face segment.
- `F`: optional spatially normalized Gaussian volume source with harmonic time dependence.

The initial implementation will support uniform `density`, `sound_speed`, and `mean_velocity`; reject `|M| >= 1` and invalid material properties.

## 3. Boundary and initial conditions

Implement patch-level conditions in `case.in`:

- `rigid`: zero normal acoustic flux (hard wall).
- `absorbing`: the outgoing-wave relation `p_t / (c + U_n) + p_n = 0`, with its normal derivative contributing implicit damping and the full mean-flow tensor retaining tangential flux terms.
- `fixedValue`: prescribed constant pressure by degree-of-freedom elimination.
- `planeWave`: prescribed pressure from the same analytical plane wave used by the initial condition.

Provide zero, rigid-cavity-mode, plane-wave, and directed Gaussian-pulse initial fields. Initialize pressure velocity consistently and compute initial acceleration from the discrete balance, including boundary acceleration for prescribed pressure.

## 4. Time integration and linear algebra

- Implement Newmark with configurable `beta` and `gamma`, defaulting to average acceleration (`1/4`, `1/2`), using Eqs. (5.249)–(5.257).
- Implement the centered three-level scheme of Eq. (5.248), with a second-order Taylor startup and a fixed time step.
- Reuse PETSc for non-diagonal systems. Use a direct nodal update when the centered scheme has only diagonal mass/damping (static media with rigid or absorbing boundaries).
- Eliminate prescribed-pressure columns and adjust the right-hand side before pinning boundary rows.
- Report Courant number, pressure extrema, discrete energy, linear residual, and iteration count. Select `dt` automatically from `cfl` and the minimum element edge when it is not specified.

## 5. Program integration and output

- Add `physics = acoustic` selection in `src/main.cpp`; retain `physics = flow` as the default for existing cases.
- Add `src/acoustic/*.cpp` to the Makefile and provide `validate-acoustic` and acoustic mesh-generation targets.
- Write nodal pressure, `p_t`, `p_tt`, and dual area in mixed Tecplot triangle/quadrilateral zones.
- Write per-step summary and probe CSV files, plus element geometry checks.

## 6. Reproducible validation

Generate conforming triangle, quadrilateral, and half-triangle/half-quadrilateral meshes with a standard-library Python script. Validate at several resolutions:

1. **Rigid rectangular cavity:** initialize the `(2,1)` cosine mode; compare pressure contours, nodal area-weighted errors, and frequency from probe zero crossings with `omega = pi c sqrt((m/Lx)^2 + (n/Ly)^2)`. Check convergence on coarse/medium/fine mixed grids and cross-check both time schemes and all three element compositions.
2. **Plane wave without flow:** compare the full field and centerline with `p = A sin(k·x - omega t + phi)`, `omega = c |k|`, after an integer number of periods.
3. **Plane wave with uniform mean flow:** repeat with `omega = c |k| + U·k`; verify downstream phase speed and mean-flow matrix assembly against the analytical field.
4. **Directed Gaussian pulse and absorbing boundary:** compare the propagating pulse with `f(n·x - (c + U_n)t)` before it reaches the outlet. Compare late-time probe histories for absorbing and rigid outlets and measure reflected/incident peak ratio.

The validation driver will run the executable, require finite results and prescribed error/reflection gates, and write comparison CSVs and SVG figures without third-party Python packages. Figures will include numerical/analytical contours, centerline overlays, convergence, and boundary-reflection histories.

## 7. Acceptance and documentation

- Build cleanly in the supported WSL/PETSc environment and keep all existing flow regressions passing.
- Require the mixed-grid cavity and plane-wave errors to decrease under refinement, with observed rates reported rather than assumed.
- Require absorbing-boundary reflection to be substantially below rigid-wall reflection; report the measured ratio and the limitations of first-order absorbing conditions.
- Write `validation/acoustic_validation.md` with governing equations, case parameters, mesh statistics, quantitative metrics, contour/comparison figures, sources, and explicit limitations (serial execution, linear elements, uniform background, no PML, no fluid–structure coupling).

## References

- Project formulation: `textbook/Chap5-4-acoustic.md`, Eqs. (5.227)–(5.267).
- COMSOL, [Gaussian Pulse in 2D Uniform Flow: Convected Wave Equation and Absorbing Layers](https://doc.comsol.com/6.4/doc/com.comsol.help.models.aco.gaussian_pulse_absorbing_layers/gaussian_pulse_absorbing_layers.html), reviewed as a classical convected-wave/open-boundary benchmark.
- Fabbri, Bruzzone, and Rosso, [pyCAFE: A Finite Element Framework for Solving Acoustic Problems in Python](https://openresearchsoftware.metajnl.com/articles/10.5334/jors.677), including analytical rigid-rectangular-cavity validation.
