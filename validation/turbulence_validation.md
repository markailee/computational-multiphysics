# Validation of standard \(k\)-\(\varepsilon\) and SST \(k\)-\(\omega\) turbulence models

## 1. Scope and conclusion

This report documents the implementation and validation of two incompressible RANS closures in the two-dimensional unstructured finite-volume solver:

- standard high-Reynolds-number \(k\)-\(\varepsilon\) with logarithmic-law wall functions;
- Menter SST-2003 \(k\)-\(\omega\) with strain-based eddy-viscosity limiting and low-\(y^+\) wall integration.

Both models were exercised through complete SIMPLE/RANS calculations on three mesh/near-wall-spacing levels for a zero-pressure-gradient flat plate and a 90° elbow. The principal findings are:

- All cases satisfy the coupled convergence, finite-value, positivity, and mass-conservation checks.
- Every flat-plate case meets the prescribed 10% normalized-RMSE target for \(C_f(Re_\theta)\) and the Coles velocity profile.
- The converged solutions have **zero cells at the viscosity-ratio cap**, so the reported flat-plate accuracy is not produced by the numerical safeguard.
- The 75° elbow comparison gives normalized RMSE of approximately 18% for both models. This is a planar symmetry-plane approximation to a three-dimensional pipe experiment, not a full 3D elbow validation.

The baseline revision before these changes was `23b8f04`. The calculations used WSL Ubuntu, PETSc 3.21.0, MPICH 4.2.0, and the C++17 release build with `-O3`.

## 2. Implemented formulation

### 2.1 Coupled RANS system

The momentum equation uses the Boussinesq closure

\[
-\rho\overline{u_i'u_j'}=2\mu_t\left(S_{ij}-\frac13S_{kk}\delta_{ij}\right)-\frac23\rho k\delta_{ij}.
\]

The implementation includes face-interpolated \(\mu+\mu_t\), the variable-eddy-viscosity transpose-stress correction, and the frozen \(-\nabla(2\rho k/3)\) term. Thus `p` is static pressure rather than modified pressure. Turbulence equations use conservative mass-flux convection, an implicit first-order upwind matrix with second-order upwind deferred reconstruction, variable diffusion, and the same steady/Euler/BDF2 history mechanism as momentum. Most supplied cases enable the Barth–Jespersen limiter for convective reconstruction; the fine SST elbow uses the default unlimited reconstruction.

Each outer iteration performs momentum prediction, pressure/velocity correction, \(k\) and scale-equation solves, and an eddy-viscosity update. Convergence requires velocity, pressure, turbulence, and cell-continuity criteria simultaneously.

### 2.2 Standard \(k\)-\(\varepsilon\)

The constants are \(C_\mu=0.09\), \(C_{\varepsilon1}=1.44\), \(C_{\varepsilon2}=1.92\), \(\sigma_k=1.0\), and \(\sigma_\varepsilon=1.3\), with \(\mu_t=C_\mu\rho k^2/\varepsilon\).

At no-slip walls, momentum uses the standard logarithmic wall function with \(\kappa=0.41\), \(E=9.8\), and a viscous/log-law switch at \(y^+=11.53\). The near-wall equations use zero-normal-gradient \(k\), wall-equilibrium production, and

\[
\varepsilon_P=\frac{C_\mu^{3/4}k_P^{3/2}}{\kappa y_P}.
\]

The nonzero wall-cell constraint is eliminated while retaining its contribution to neighboring right-hand sides. A `maximum_eddy_viscosity_ratio = 1000` startup safeguard is implemented, but its use is reported separately from the model equations.

### 2.3 Menter SST-2003

The SST implementation includes \(F_1/F_2\), cross-diffusion, blended coefficients, limited production in both equations, and

\[
\mu_t=\frac{\rho a_1k}{\max(a_1\omega,SF_2)},\qquad a_1=0.31.
\]

The constants are \(\beta^*=0.09\), \(\gamma_1=5/9\), \(\gamma_2=0.44\), \(\sigma_{k1}=0.85\), \(\sigma_{k2}=1.0\), \(\sigma_{\omega1}=0.5\), \(\sigma_{\omega2}=0.856\), \(\beta_1=0.075\), and \(\beta_2=0.0828\). Wall distances are exact minimum point-to-segment distances. The no-slip conditions are \(k=0\) and \(\omega_w=6\nu/(\beta_1y_P^2)\), without the standard \(k\)-\(\varepsilon\) wall function.

Implementation references are the project teaching chapter, the [NASA TMR SST definition](https://tmbwg.github.io/turbmodels/sst.html), and Menter, Kuntz, and Langtry, *Ten Years of Industrial Experience with the SST Turbulence Model* (2003).

## 3. Cases and reference data

### 3.1 Flat plate

The domain follows the NASA 2DZPG layout: symmetry run-up to \(x=0\), no-slip plate \(0\le x\le2\), unit domain height, \(U_\infty=1\), \(\rho=1\), and \(\mu=2\times10^{-7}\). Thus \(Re_{x=1}=5\times10^6\). The local freestream velocity is used when integrating each momentum-thickness profile because the slip farfield permits a small acceleration.

The comparison data are unmodified NASA files:

- [Kármán–Schoenherr \(C_f(Re_\theta)\)](https://tmbwg.github.io/turbmodels/FlatPlate_validation/cf_K-S.dat), used for \(4000<Re_\theta<13000\);
- [Coles \(u^+(y^+)\) at \(Re_\theta=10000\)](https://tmbwg.github.io/turbmodels/FlatPlate_validation/u+y+.dat);
- [typical \(Re_\theta(x)\)](https://tmbwg.github.io/turbmodels/FlatPlate_validation/retheta_variation_typical.dat), used only as a setup cross-check.

These are theoretical correlations calibrated to experiments, not pointwise raw measurements. See `reference/flat_plate/README.md` for provenance.

### 3.2 Elbow

The planar elbow has \(R_c/D=2.8\), \(Re_D=43000\), 5D straight inlet, and 10D straight outlet. Inlet turbulence uses intensity 5% and length scale \(0.07D\). The primary comparison is the 75° radial velocity-magnitude profile from Enayet et al. (1982), distributed publicly in [Ansys VMFL030](https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/fbu_vm/Hlp_VMFL030.html).

The CSV was manually digitized from the experimental symbols. Estimated reading uncertainties are ±0.005 in \(y/d\) and ±0.02 in normalized velocity. The figure, CSV, and extraction notes are in `reference/elbow_enayet_1982/`. As documented there, \(y/d=0\) is mapped to the inner wall and \(y/d=1\) to the outer wall without fitting or reversal. This mapping is an approximation, not a reconstruction of the full 3D measurement plane.

## 4. Numerical integrity

All final cases have `converged = 1`, finite velocity/pressure/turbulence values, positive \(k\) and scale fields, maximum cell mass imbalance below `1e-8`, and absolute net boundary flux below `1e-8`. Across the matrix, cell imbalance and net boundary flux are at round-off levels (order `1e-15` or smaller). The final solutions have no cells at the viscosity-ratio cap.

The first-cell \(y^+\) values reported below are minima/maxima over the relevant no-slip patch, excluding the first 0.1 m of the plate. Low-shear elbow regions can fall below the nominal logarithmic-layer target; the wall function switches to its viscous relation below \(y^+=11.53\).

## 5. Flat-plate results

### 5.1 Standard \(k\)-\(\varepsilon\)

| Level | Cells | Iterations | \(y^+\) min–max | \(C_f\) RMSE | \(C_f\) max relative error | \(u^+\) RMSE | Drag coefficient |
|---|---:|---:|---:|---:|---:|---:|---:|
| Coarse | 4,500 | 409 | 69.81–89.65 | 6.836% | 10.349% | 1.684% | 0.00144102 |
| Medium | 9,380 | 722 | 43.60–54.88 | 6.850% | 10.465% | 1.839% | 0.00144115 |
| Fine | 20,000 | 1,389 | 26.17–33.04 | 6.812% | 10.623% | 2.052% | 0.00144132 |

The \(C_f\) and integrated-drag results are nearly insensitive to this grid/\(y^+\) sequence. The Coles-profile error remains near 2%. This is a grid and near-wall-spacing sensitivity study, not a formal order-of-accuracy measurement.

### 5.2 SST \(k\)-\(\omega\)

| Level | Cells | Iterations | \(y^+\) min–max | \(C_f\) RMSE | \(C_f\) max relative error | \(u^+\) RMSE | Drag coefficient |
|---|---:|---:|---:|---:|---:|---:|---:|
| Coarse | 4,500 | 1,247 | 1.803–2.338 | 2.775% | 4.593% | 6.812% | 0.00157507 |
| Medium | 9,380 | 4,015 | 0.880–1.111 | 4.785% | 8.283% | 2.395% | 0.00147422 |
| Fine | 20,000 | 3,358 | 0.437–0.551 | 6.451% | 10.198% | 1.398% | 0.00144721 |

SST improves the Coles profile as the near-wall mesh is refined, while the \(C_f\) RMSE increases from its especially good coarse-grid value. Fine-grid drag is close to the \(k\)-\(\varepsilon\) value. The coarse SST grid is already in a low-\(y^+\) regime, so these trends should not be interpreted as a monotone formal convergence rate.

### 5.3 Comparison plots

For a representative point at \(x=0.9743\) on the medium grids, the reference is interpolated at each model's computed \(Re_\theta\):

| Model | Computed \(Re_\theta\) | Numerical \(C_f\) | Reference \(C_f\) | Relative error |
|---|---:|---:|---:|---:|
| \(k\)-\(\varepsilon\) | 5367.3 | 0.0027122 | 0.0029651 | -8.53% |
| SST | 5513.7 | 0.0027679 | 0.0029495 | -6.16% |

- [Skin friction versus momentum-thickness Reynolds number](figures/plate_cf_comparison.svg)
- [Coles velocity profile comparison](figures/plate_uplus_comparison.svg)
- Per-case numerical/reference rows: `validation/data/plate_*_cf.csv` and `validation/data/plate_*_uplus.csv`.

## 6. Elbow results

| Model | Level | Cells | Iterations | \(y^+\) min–max | 75° normalized RMSE | Maximum relative error | Pressure drop | Loss coefficient |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| \(k\)-\(\varepsilon\) | Coarse | 3,960 | 765 | 23.24–37.44 | 17.749% | 41.339% | 0.116391 | 0.232782 |
| \(k\)-\(\varepsilon\) | Medium | 8,500 | 1,313 | 16.73–27.67 | 17.923% | 41.674% | 0.117090 | 0.234181 |
| \(k\)-\(\varepsilon\) | Fine | 16,800 | 2,295 | 11.93–20.73 | 18.069% | 41.747% | 0.117174 | 0.234348 |
| SST | Coarse | 3,960 | 741 | 1.524–4.047 | 17.474% | 41.196% | 0.138490 | 0.276980 |
| SST | Medium | 8,500 | 1,263 | 0.688–2.263 | 17.580% | 41.590% | 0.123978 | 0.247955 |
| SST | Fine | 16,800 | 2,136 | 0.329–1.257 | 17.861% | 41.148% | 0.120446 | 0.240892 |

The pressure drop is the static-pressure difference between uniformly sampled section averages 0.5D inside the straight inlet and 0.5D inside the straight outlet. The loss coefficient uses \(0.5\rho U_b^2\) and includes the complete straight and bent domain; it is not an elbow-only loss coefficient and is not compared directly with a correlation.

The medium-grid minimum streamwise velocity at 75° is 0.621 for \(k\)-\(\varepsilon\) and 0.297 for SST. Both models reproduce the broad velocity redistribution, but neither can represent circumferential Dean vortices or the measured three-dimensional secondary flow. The experimental digitization uncertainty is also much smaller than the 2D/3D modeling discrepancy, so the reported errors must not be used to rank the closures as if this were full 3D validation.

[75° numerical/experimental comparison](figures/elbow_75deg_comparison.svg)

Representative medium-grid values, normalized by unit bulk speed, are:

| \(y/d\) | Enayet experiment | \(k\)-\(\varepsilon\) | SST |
|---:|---:|---:|---:|
| 0.050 | 0.880 | 1.007 | 0.993 |
| 0.210 | 1.050 | 1.120 | 1.153 |
| 0.410 | 1.120 | 1.075 | 1.088 |
| 0.550 | 1.150 | 1.037 | 1.041 |
| 0.700 | 0.700 | 0.994 | 0.992 |
| 0.835 | 0.820 | 0.908 | 0.904 |
| 0.930 | 0.640 | 0.781 | 0.758 |

## 7. Reproduction

From the project directory in WSL:

```bash
export PETSC_DIR=/home/lwa/program/petsc-3.21.0
export PATH=/home/lwa/program/mpich-4.2.0/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin
export LD_LIBRARY_PATH=/home/lwa/program/petsc-3.21.0/lib:/home/lwa/program/lapack-3.12.1/lib:/home/lwa/program/mpich-4.2.0/lib:${LD_LIBRARY_PATH:-}

make -j
python3 cases/generate_meshes.py
make validate-turbulence
```

`make validate-turbulence` executes all twelve cases and regenerates `validation/metrics.csv`, the per-case CSVs, and SVG plots. `python3 cases/validate_turbulence.py` repeats only the post-processing. `make validate` additionally runs the original cavity, backstep, and cylinder regressions.

## 8. Limitations

- The solver is two-dimensional. The elbow comparison cannot validate three-dimensional secondary flow, swirl, or cross-sectional vortex structure.
- Both closures are linear eddy-viscosity models without transition, compressibility, buoyancy, rotation/curvature correction, DES, or Reynolds-stress transport.
- All supplied cases use second-order upwind reconstruction with an implicit first-order upwind matrix. The Barth–Jespersen limiter is enabled except for the fine SST elbow, which uses the default unlimited reconstruction. Grid and first-cell spacing are varied, but no formal spatial convergence order is claimed.
- The flat-plate references are correlations, while the elbow reference contains figure-digitization and planar-mapping uncertainty.
- Inlet turbulence and farfield conditions affect RANS predictions; the documented values are part of the case definition and should not be changed to tune agreement.
- The meshes are modest teaching/validation grids, not production grid-independent certifications.

## 9. Final assessment

The standard \(k\)-\(\varepsilon\) with wall functions and Menter SST-2003 implementations are operational in the current solver and pass the specified flat-plate reference checks at every tested level. Their elbow solutions are converged and physically meaningful within the planar approximation, with approximately 18% profile RMSE. Full elbow certification requires a three-dimensional solver and matching experimental boundary conditions.
