# 2D Multiphysics Finite-Volume Solver

A C++17 finite-volume project for two-dimensional multiphysics calculations on unstructured meshes. The executable implements **incompressible, Boussinesq thermal-dilation, and ideal-gas compressible Navier–Stokes solvers** using cell-centered SIMPLE pressure–velocity coupling and a **linear acoustic solver** using the node-centered dual-control-volume formulation of `textbook/Chap5-4-acoustic.md`. The mesh, discretization, linear-algebra, and output components are shared foundations for further physics modules.

The primary purpose of this solver is to support **course teaching for computational multiphysics**. It is also intended for numerical-method development. It emphasizes conservative face fluxes, explicit mesh connectivity, readable operator implementations, and reproducible benchmark cases.

**Copyright:** © Wanai Li, Sino-French Institute of Nuclear Engineering and Technology, Sun Yat-Sen University. All rights reserved.

> **Current scope:** serial two-dimensional incompressible/Boussinesq flow, single-component laminar ideal-gas compressible flow, passive inert multispecies transport, and linear acoustics are implemented as separate physics modules, with one-way Curle-type flow-to-acoustic source coupling for the cylinder aeolian-tone case. Flow and species use cell-centered fields; acoustics uses nodal fields on median-dual control areas. Composition-driven variable-density flow, reactions, real-gas models, compressible turbulence, elasticity, two-way fluid–structure/acoustic coupling, non-uniform acoustic media, MPI-distributed calculations, and the VOF/DPM/Euler–Euler solve routines remain future work.

## Contents

- [Implemented features](#implemented-features)
- [Numerical formulation](#numerical-formulation)
- [Project layout](#project-layout)
- [Prerequisites](#prerequisites)
- [Build and run in WSL](#build-and-run-in-wsl)
- [Case configuration](#case-configuration)
- [Boundary conditions](#boundary-conditions)
- [Material and multispecies configuration](#material-and-multispecies-configuration)
- [Acoustic configuration](#acoustic-configuration)
- [Gmsh meshes](#gmsh-meshes)
- [Output files](#output-files)
- [Validation cases](#validation-cases)
- [Adding a new case](#adding-a-new-case)
- [Adding another physics solver](#adding-another-physics-solver)
- [Current limitations](#current-limitations)
- [To-do work](#to-do-work)
- [Troubleshooting](#troubleshooting)

## Implemented features

### Unstructured mesh support

- Two-dimensional cell-centered control volumes.
- Linear Gmsh triangles and quadrilaterals.
- ASCII Gmsh 2.2 and 4.1 input.
- Physical boundary names read from Gmsh line elements.
- Owner/neighbour face connectivity and cell-to-face connectivity.
- Cell volumes, centroids, face areas, outward normals, interpolation weights, and non-orthogonal geometry data.
- Conservative internal-face fluxes: one flux is stored and used with opposite signs for the two adjacent cells.

### General finite-volume infrastructure

The reusable code under `src/fvm/` includes:

- Scalar and vector cell fields with boundary values.
- Green–Gauss gradient reconstruction with reconstructed-face corrections.
- Weighted least-squares gradient reconstruction.
- Steady, first-order backward-Euler, and constant-step BDF2 time terms. BDF2 uses backward Euler for its first step.
- Second-order upwind face reconstruction with an implicit upwind/deferred-correction convection operator and an optional Barth–Jespersen limiter.
- Over-relaxed orthogonal diffusion with explicit non-orthogonal correction.
- PETSc sparse-matrix, vector, KSP, and preconditioner wrappers.
- Tecplot ASCII output for unstructured triangles and quadrilaterals.
- A `PhysicsSolver` interface for future physics modules.

### Pressure-based Navier–Stokes solvers

The code under `src/flow/` implements:

- The SIMPLE algorithm on collocated unstructured grids.
- A shared cell-centered velocity and pressure arrangement.
- Rhie–Chow momentum-weighted face interpolation to suppress odd–even pressure decoupling.
- Conservative mass-flux updates in the pressure-correction equation.
- Separate velocity and pressure under-relaxation.
- Configurable PETSc linear solvers and preconditioners.
- Fixed-value, no-slip, zero-gradient, slip/symmetry, and parabolic velocity boundaries.
- Fixed-value and zero-gradient pressure boundaries.
- Automatic pressure-reference fixing when all pressure boundaries are Neumann conditions.
- Optional steady or transient calculations and a spatially uniform body force.
- Boussinesq thermal dilation with reference-density transport, hydrostatically modified pressure, integrated buoyancy, and coupled temperature transport.
- Single-component ideal-gas compressible flow with absolute pressure, cell-variable density, density storage and pressure convection in pressure correction, and conservative total-enthalpy energy with pressure work, kinetic-energy diffusion correction, and viscous-work fluxes.
- Pressure-based compressible SIMPLE for subsonic, transonic, and supersonic flow, with upwind pressure convection and consistent extrapolated-outflow density correction.
- Fixed-temperature, adiabatic/zero-gradient, and inward-positive fixed-heat-flux thermal boundaries, plus uniform volumetric heating.
- Standard high-Reynolds-number \(k\)-\(\varepsilon\) turbulence with logarithmic-law wall functions.
- Menter SST-2003 \(k\)-\(\omega\) turbulence with strain-based eddy-viscosity limiting and low-\(y^+\) wall integration.
- Conservative transport of turbulence variables with implicit upwind matrices and second-order upwind deferred reconstruction.
- Turbulence fields, wall shear, friction velocity, and \(y^+\) diagnostics in Tecplot and CSV output.

### Material models and multispecies transport

The independent `src/material/` layer provides constant-property fluids, a single-component ideal gas for compressible flow, ideal-gas mixtures, and volume-fraction mixture property rules. Ideal-gas mixtures support mass/mole conversion, Wilke viscosity, arithmetic/harmonic thermal conductivity, mass-weighted heat capacity, and constant or mixture-averaged diffusivities. All inputs use SI units.

The optional `species_model = transport` module adds two or more inert species to the existing cell-centered SIMPLE solver. It reuses steady, backward-Euler, and BDF2 discretizations, supports fixed-value and zero-gradient composition boundaries, includes turbulent diffusivity through a configurable Schmidt number, and solves the first `N-1` species with a conservative correction that enforces zero summed diffusion flux, while the last species is computed from the unit-sum constraint and its flux is closed by the residual. Species convergence, boundedness, composition closure, and mixture properties are included in output. The current flow coupling is explicitly **passive and constant-property**; composition-dependent low-Mach pressure coupling is not yet implemented.

`src/flow/multiphase/` defines independent VOF, algebraic-slip, Euler–Euler, and DPM interfaces and their coupling context. Their solve routines are intentionally unavailable at this stage and throw a clear exception if selected.

### Linear acoustic solver

The code under `src/acoustic/` implements the formulation of `textbook/Chap5-4-acoustic.md`:

- Nodal acoustic pressure on conforming mixed linear-triangle/quadrilateral meshes.
- Median-dual control areas and diagonal lumped acoustic mass.
- Full `(I - M ⊗ M)` fluxes and the centered mixed time-space term for a uniform subsonic mean flow.
- Rigid, first-order absorbing, fixed-pressure, and analytical-plane-wave boundary conditions.
- Zero, rigid-cavity-mode, plane-wave, and directed Gaussian-pulse initial conditions.
- Optional spatially normalized harmonic Gaussian source.
- Newmark and centered three-level time integration with consistent initial acceleration.
- PETSc solves for non-diagonal systems and direct nodal updates for diagonal centered-difference systems.
- Nodal Tecplot pressure/output rates, per-step summaries, probe histories, and geometry checks.

## Numerical formulation

### Governing equations

With the default `flow_model = incompressible`, the solver uses constant density \(\rho\) and constant dynamic viscosity \(\mu\). The incompressible momentum and continuity equations are

$$
\frac{\partial(\rho\mathbf{U})}{\partial t}
+\nabla\cdot(\rho\mathbf{U}\otimes\mathbf{U})
=-\nabla p+\mu\nabla^2\mathbf{U}+\rho\mathbf{f},
$$

$$
\nabla\cdot\mathbf{U}=0.
$$

The `viscosity` setting is the **dynamic viscosity** \(\mu\), not the kinematic viscosity \(\nu=\mu/\rho\). The supplied cases use \(\rho=1\), so their numerical dynamic and kinematic viscosities have the same value.

Integrating the momentum equation over control volume \(P\) gives a balance of time accumulation, convection, diffusion, pressure, and body-force terms. Each internal face has a single owner-oriented normal and mass flux, which preserves discrete conservation.

With `flow_model = boussinesq`, density remains at its reference value except for the buoyancy source \(-\rho_0\beta(T-T_0)\mathbf gV_P\), and `p` is hydrostatically modified pressure. Boussinesq flow solves the conservative temperature/static-enthalpy equation. With `flow_model = compressible`, the energy unknown is total enthalpy \(H=c_pT+|\mathbf U|^2/2\); temperature is recovered from \(H\) and the latest velocity, and `material_model = ideal_gas` evaluates \(\rho=p/(RT)\) from absolute pressure and recovered temperature. The pressure equation additionally contains density storage and upwind pressure convection. The same pressure-based SIMPLE formulation is used through Mach one and for supersonic flow; extrapolated outflows include the outgoing density correction in both the pressure matrix and mass-flux update. The full equations and implementation notes are in `textbook/Chap5-1-SIMPLE-scheme-basic-flow.md`, sections 5.7.7 and 5.7.8.

### Gradient reconstruction

Two cell-gradient options are available:

- `gradient = gauss`: Green–Gauss integration using face-interpolated values, followed by reconstructed-face correction sweeps.
- `gradient = least_squares`: a weighted least-squares fit to neighbouring cell differences, with Dirichlet boundary constraints where applicable. This is the default.

The same gradient interface is used for scalar fields and both components of vector fields. Set `gradient_limiter = barth_jespersen` to limit these gradients for convective reconstruction; the default is `none`. The Barth–Jespersen coefficient is computed independently for every scalar or vector component from the cell and its immediate neighbours. Face-centroid reconstructions are constrained to the local minimum and maximum, and the whole gradient of that component is scaled by one coefficient in `[0,1]`. Pressure correction, Rhie–Chow interpolation, and diffusive gradients remain unlimited.

### Convection

Convection uses second-order upwind reconstruction by default. For upwind cell \(U\), the face value is

$$
\phi_f^{\mathrm{SOU}}=\phi_U+\alpha_U\mathbf g_U\cdot(\mathbf x_f-\mathbf x_U),
$$

where \(\alpha_U=1\) without a limiter and is the Barth–Jespersen coefficient when enabled. The implicit matrix retains first-order upwind convection for robust diagonal properties. The deferred correction

$$
F_f^{\mathrm{SOU}}=F_f^{\mathrm{upwind}}
+\dot m_f\left(\phi_f^{\mathrm{SOU}}-\phi_U\right)
$$

is evaluated with the current outer-iteration field. Momentum, energy, turbulence, and species transport all use this operator; there is no central/blended convection option.

### Diffusion

The face-normal diffusive derivative uses an over-relaxed orthogonal decomposition. The implicit contribution is based on the centre-to-centre or centre-to-boundary distance. The remaining non-orthogonal contribution is evaluated explicitly from reconstructed cell gradients.

### Time derivative

The `time_scheme` setting accepts:

- `steady`: no time-derivative contribution.
- `euler` or `backward_euler`: first-order implicit backward Euler.
- `bdf2`: second-order backward differentiation for constant time steps, with backward-Euler startup.

For transient cases, SIMPLE outer iterations are repeated within each time step until the requested update and continuity tolerances are reached.

### SIMPLE and Rhie–Chow coupling

Each outer iteration performs the following operations:

1. Reconstruct velocity and pressure gradients.
2. Assemble and solve the two momentum components using the current pressure and face mass flux.
3. Recompute face mass fluxes from the predicted velocity using Rhie–Chow interpolation.
4. Assemble and solve the pressure-correction equation so that the corrected face fluxes satisfy discrete continuity; compressible flow also includes density storage and pressure convection.
5. Correct cell velocities, pressure, and face mass fluxes, then update density from the equation of state in compressible mode.
6. Solve energy when enabled, update temperature-dependent density and properties, and check velocity, pressure, temperature, density, and cell mass imbalances.

The Rhie–Chow normal velocity contains the difference between the face-normal pressure derivative and the interpolated cell pressure gradient. Its momentum weighting uses the relaxed momentum diagonal.

If no fixed-pressure boundary exists, the pressure-correction equation has a constant nullspace. The code removes it by fixing the pressure correction in cell zero. With a fixed-pressure outlet, no additional reference cell is needed.

### Linear systems

`src/fvm/PetscLinearSystem.hpp` provides the sparse-matrix and PETSc interface. The present implementation uses sequential `MATSEQAIJ` matrices and sequential vectors. Momentum and pressure-correction systems have separate persistent matrix/vector/KSP objects, but their nonzero values are refreshed during outer iterations.

Default linear-solver choices are:

| Equation | KSP | Preconditioner |
|---|---|---|
| Momentum | `bcgs` | `ilu` |
| Pressure correction | `cg` | `icc` |

Compressible pressure correction defaults to `bcgs` with `ilu` because density convection makes its matrix nonsymmetric. Energy also defaults to `bcgs`/`ilu`. These choices can be changed in `case.in`. Standard PETSc command-line options are passed through solver initialization. Turbulence equations default to `bcgs` with `ilu`.

### RANS equations and eddy viscosity

With `turbulence_model` different from `laminar`, the solver uses the Reynolds-averaged Navier–Stokes equations. The Boussinesq closure is

$$
-\rho\overline{u_i'u_j'}=2\mu_t\left(S_{ij}-\frac13S_{kk}\delta_{ij}\right)-\frac23\rho k\delta_{ij}.
$$

Momentum assembly uses \(\mu+\mu_t\) in the normal diffusion operator and includes the variable-\(\mu_t\) transpose-stress correction. The isotropic \(-2\rho k/3\) gradient is assembled explicitly with the current turbulent kinetic energy, so `p` remains static pressure rather than modified pressure. Production uses the full strain expression; no curvature or rotation correction is applied.

Each SIMPLE outer iteration solves momentum, pressure correction, \(k\), and the selected turbulence scale equation sequentially. Turbulence source terms and eddy viscosity are frozen within each linear solve and updated between outer iterations. Turbulence equations support the same steady, backward-Euler, and BDF2 time histories as momentum.

### Standard \(k\)-\(\varepsilon\) with wall functions

The standard model solves

$$
\mu_t=C_\mu\rho k^2/\varepsilon,\qquad C_\mu=0.09,\quad C_{\varepsilon1}=1.44,\quad
C_{\varepsilon2}=1.92,\quad \sigma_k=1.0,\quad \sigma_\varepsilon=1.3.
$$

At a no-slip wall, tangential momentum diffusion is replaced by wall-function shear. With \(u_\tau=C_\mu^{1/4}\sqrt{k_P}\), the implementation uses the molecular relation below \(y^+=11.53\) and

$$
U^+=\frac1\kappa\ln(Ey^+),\qquad \kappa=0.41,\quad E=9.8,
$$

above it. The near-wall kinetic-energy equation has zero normal gradient, consistent wall production, and the equilibrium value

$$
\varepsilon_P=\frac{C_\mu^{3/4}k_P^{3/2}}{\kappa y_P}.
$$

The nonzero wall-cell constraint is eliminated with its value retained in neighboring right-hand sides. The first-cell center should lie in the logarithmic layer, normally \(30\lesssim y^+\lesssim200\); the implementation does not make a low-Reynolds-number \(k\)-\(\varepsilon\) model out of a fine near-wall mesh.

A configurable `maximum_eddy_viscosity_ratio = 1000` is available as a numerical realizability safeguard against singular startup values. It is not a model constant. Cap usage is reported by the validation driver and should be negligible in a converged, physically meaningful solution.

### Menter SST-2003 \(k\)-\(\omega\)

The SST implementation uses the 2003 strain-invariant form, \(F_1/F_2\) blending, cross-diffusion, limited production in both equations, and

$$
\mu_t=\frac{\rho a_1k}{\max(a_1\omega,SF_2)},\qquad a_1=0.31.
$$

Constants are \(\beta^*=0.09\), \(\gamma_1=5/9\), \(\gamma_2=0.44\), \(\sigma_{k1}=0.85\), \(\sigma_{k2}=1.0\), \(\sigma_{\omega1}=0.5\), \(\sigma_{\omega2}=0.856\), \(\beta_1=0.075\), and \(\beta_2=0.0828\). Wall distance is the exact minimum point-to-segment distance to a no-slip wall. At the wall, \(k=0\) and

$$
\omega_w=\frac{6\nu}{\beta_1y_P^2}.
$$

SST is integrated to the wall and should use a wall-resolved mesh, normally \(y^+\approx1\) or smaller. The standard \(k\)-\(\varepsilon\) wall function is not used for SST.

### Turbulence boundary and initial settings

Select the model and initialize positive fields, for example:

```ini
turbulence_model = k_epsilon       # laminar, k_epsilon, or sst_k_omega
initial_velocity = 1 0
initial_k = 1.5e-4
initial_epsilon = 1.125           # use initial_omega for SST
turbulence_relaxation = 0.3
gradient_limiter = barth_jespersen
turbulence_tolerance = 2.0e-6
maximum_eddy_viscosity_ratio = 1000
```

Every patch of a turbulent case needs `k` and the appropriate scale entry in addition to `velocity` and `pressure`:

```ini
[boundary.inlet]
velocity = fixedValue 1 0
pressure = zeroGradient
k = fixedValue 1.5e-4
epsilon = fixedValue 1.125

[boundary.outlet]
velocity = zeroGradient
pressure = fixedValue 0
k = zeroGradient
epsilon = zeroGradient

[boundary.walls]
velocity = noSlip
pressure = zeroGradient
k = zeroGradient                 # SST: fixedValue 0
epsilon = wallFunction           # SST: omega = wallFunction
```

For an inlet with turbulence intensity \(I\) and eddy-viscosity ratio \(r=\mu_t/\mu\),

$$
k=\frac32(IU)^2,\qquad
\varepsilon=\frac{C_\mu\rho k^2}{r\mu},\qquad
\omega=\frac{\rho k}{r\mu}.
$$

Alternatively, \(\varepsilon=C_\mu^{3/4}k^{3/2}/\ell\) and \(\omega=\varepsilon/(\beta^*k)\) can be used with a chosen turbulence length scale \(\ell\). The supplied plate cases use \(I=1\%\) and \(r=0.009\) for \(k\)-\(\varepsilon\), and \(I=0.039\%\) and \(r=0.009\) for SST; the elbow cases use \(I=5\%\) and \(\ell=0.07D\).

Turbulence lower bounds default to `minimum_k = 1e-12`, `minimum_epsilon = 1e-12`, and `minimum_omega = 1e-8`. `turbulence_ksp` and `turbulence_pc` select the scalar linear solver and preconditioner. The summary's `turbulence_change` is the larger infinity-norm update of \(k\) and its scale, normalized by the corresponding current field maximum.

## Project layout

```text
.
├── Makefile                    Build, mesh-generation, and validation targets
├── README.md                   This document
├── src/
│   ├── main.cpp                Program entry and PETSc lifetime management
│   ├── mesh/
│   │   ├── Vector2.hpp         Two-dimensional vector arithmetic
│   │   ├── Mesh.hpp/.cpp       Geometry and unstructured connectivity
│   │   └── GmshReader.hpp/.cpp Gmsh 2.2/4.1 ASCII reader
│   ├── fvm/
│   │   ├── Field.hpp           Scalar/vector fields and boundary values
│   │   ├── Gradient.hpp/.cpp   Gauss and least-squares reconstruction
│   │   ├── Operators.hpp/.cpp  Interpolation, convection, and diffusion data
│   │   ├── ScalarTransport.*   Conservative variable-diffusivity scalar equations
│   │   ├── TimeDerivative.*    Steady, Euler, and BDF2 coefficients
│   │   ├── PetscLinearSystem.* Sparse matrices and PETSc KSP wrapper
│   │   ├── TecplotWriter.*     Unstructured Tecplot ASCII output
│   │   ├── Config.hpp/.cpp     Case-file parser
│   │   └── PhysicsSolver.hpp   Common interface for physics modules
│   ├── material/
│   │   ├── MaterialModel.*     Material interface, constant fluid, and factory
│   │   ├── IdealGasMixture.*   Ideal gas, Wilke viscosity, and mixture diffusivities
│   │   └── VolumeFractionMixture.* Phase-fraction property mixing
│   ├── flow/
│   │   ├── BoundaryConditions.* Flow boundary-condition handling
│   │   ├── EnergyModel.*       Temperature boundaries, total-enthalpy transport, and pressure/viscous work
│   │   ├── TurbulenceModel.*   Standard k-epsilon and SST-2003 implementations
│   │   ├── MultiSpeciesModel.* Species fields, transport, flux closure, and diagnostics
│   │   ├── multiphase/        VOF, ASM, Euler–Euler, DPM interfaces (not implemented)
│   │   └── SimpleSolver.*      SIMPLE, Rhie–Chow, RANS, and optional module scheduling
│   └── acoustic/
│       ├── AcousticBoundaryConditions.* Acoustic patch parsing
│       ├── AcousticSolver.*    Nodal dual-volume operators and time integration
│       └── IMPLEMENTATION_PLAN.md Numerical implementation and verification plan
├── cases/
│   ├── generate_meshes.py      Reproducible supplied Gmsh meshes
│   ├── validate.py             Laminar conservation and benchmark checks
│   ├── validate_turbulence.py  Turbulence runs, reference errors, and plots
│   ├── generate_acoustic_meshes.py Mixed acoustic mesh generation
│   ├── validate_acoustic.py    Analytical acoustic runs, metrics, and SVG plots
│   ├── generate_multispecies_meshes.py Structured species-verification meshes
│   ├── validate_multispecies.py Species analytical/property checks and grid refinement
│   ├── generate_thermal_meshes.py Structured thermal-cavity, nozzle, and bump meshes
│   ├── validate_thermal_compressible.py Boussinesq/compressible analytical validation
│   ├── boussinesq_cavity/      de Vahl Davis natural-convection grids
│   ├── compressible_conduction/ Zero-flow ideal-gas conduction
│   ├── compressible_closed_heating/ Transient constant-mass heating
│   ├── compressible_nozzle/    Subsonic isentropic nozzle grids
│   ├── compressible_bump/      Subsonic/transonic/supersonic compressible-SIMPLE bump cases
│   ├── species_stefan_tube/    Binary analytical Stefan diffusion at three resolutions
│   ├── species_flux_correction/ Unequal-diffusivity analytical flux-closure checks
│   ├── species_diffusion/      Three-species transient analytical diffusion
│   ├── species_air_properties/ Ideal-gas and Wilke property verification
│   ├── acoustic_*/             Cavity, plane-wave, convection, and absorption cases
│   ├── cylinder_aeolian/       Transient cylinder flow to Curle-dipole acoustic propagation
│   ├── run_cylinder_aeolian.py Hybrid flow/source/acoustic workflow and plots
│   ├── cavity/                 Lid-driven cavity configuration and mesh
│   ├── backstep/               Backward-facing step configuration and mesh
│   ├── cylinder/               Flow-around-cylinder configuration and mesh
│   ├── plate_*/                Flat-plate cases for both turbulence models
│   └── elbow_*/                90-degree elbow cases for both models
├── reference/                 NASA and Enayet/Ansys comparison data
└── validation/               Turbulence report, metrics, data, and figures
```

Additional teaching material is under `textbook/` and `texbook-ppt/`. The elastic directory is reserved for future development.

## Prerequisites

The supported build environment is WSL Ubuntu with:

- A C++17 compiler.
- GNU Make.
- MPICH providing `mpic++`.
- PETSc 3.21.0, including its KSP library and external-link configuration.
- Python 3 for mesh generation and validation.

The configured installation paths are:

```text
PETSc: /home/lwa/program/petsc-3.21.0
MPICH: /home/lwa/program/mpich-4.2.0
```

Python uses only the standard library. Gmsh is **not required** to run the supplied cases because their `.msh` files are already included. Gmsh is only needed if you want to generate independent meshes from your own geometry.

## Build and run in WSL

### 1. Configure the environment

From Git Bash, WSL commands may need `MSYS_NO_PATHCONV=1` so that Linux paths are not rewritten. Inside WSL, set:

```bash
export PETSC_DIR=/home/lwa/program/petsc-3.21.0
export PATH=/home/lwa/program/mpich-4.2.0/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin
export LD_LIBRARY_PATH=/home/lwa/program/petsc-3.21.0/lib:/home/lwa/program/lapack-3.12.1/lib:/home/lwa/program/mpich-4.2.0/lib:${LD_LIBRARY_PATH:-}
```

Then enter the project directory. On the current machine it is:

```bash
cd "/mnt/d/LessonTeaching/计算多物理场/课件"
```

A non-interactive Bash invocation does not necessarily source the profile that normally defines these paths.

### 2. Build

```bash
make -j
```

The executable is created as:

```text
multiphysics
```

Object and dependency files are placed under `build/`. To use a different PETSc prefix:

```bash
make -j PETSC_DIR=/path/to/petsc
```

For a clean rebuild:

```bash
make clean
make -j
```

### 3. Run a case

Pass either a case directory or the full path to a `case.in` file:

```bash
./multiphysics cases/cavity
./multiphysics cases/backstep/case.in
./multiphysics cases/cylinder
./multiphysics cases/plate_k_epsilon
./multiphysics cases/plate_sst_k_omega
./multiphysics cases/elbow_k_epsilon
./multiphysics cases/elbow_sst_k_omega
./multiphysics cases/species_stefan_tube/medium
./multiphysics cases/species_flux_correction/medium
./multiphysics cases/species_diffusion/medium
./multiphysics cases/species_air_properties
./multiphysics cases/boussinesq_cavity/fine
./multiphysics cases/compressible_conduction
./multiphysics cases/compressible_closed_heating
./multiphysics cases/compressible_nozzle/fine
./multiphysics cases/compressible_bump/subsonic
./multiphysics cases/compressible_bump/transonic
./multiphysics cases/compressible_bump/supersonic
./multiphysics cases/acoustic_cavity
./multiphysics cases/acoustic_plane_wave
./multiphysics cases/acoustic_convected
./multiphysics cases/acoustic_absorbing
```

Mesh and output paths in `case.in` are resolved relative to that case directory, not relative to the shell's working directory.

Run the solver directly, without `mpirun`; the current implementation intentionally rejects multi-rank execution rather than allowing independent ranks to overwrite the same output.

### Make targets

| Command | Purpose |
|---|---|
| `make` or `make -j` | Build `multiphysics`. |
| `make meshes` | Regenerate all supplied Gmsh meshes. This overwrites the generated `.msh` files. |
| `make validate` | Run all laminar regressions and the complete plate/elbow turbulence matrix. |
| `make validate-turbulence` | Run both models at all plate/elbow resolutions, check conservation and accuracy, and generate comparison data and SVG plots. |
| `make validate-acoustic` | Run all analytical acoustic cases, check errors/reflection, and generate comparison data and SVG plots. |
| `make validate-multispecies` | Run ten species/property cases, check analytical errors, grid order, conservation, and generate metrics and SVG plots. |
| `make validate-thermal` | Run the thermal and nozzle cases plus subsonic, transonic, and supersonic bump flows. |
| `make validate-bump` | Run the three compressible-SIMPLE circular-arc bump regimes and reference comparisons. |
| `make cylinder-aeolian` | Run transient `Re = 150` cylinder flow, construct a Curle dipole from its force history, propagate the sound, and plot the results. |
| `make clean` | Remove the executable and `build/`. It does not remove case results or `validation/`. |

## Case configuration

A case is defined by a `case.in` file. Settings use `key = value`; `#` starts a comment. Boundary sections use the Gmsh physical name:

```ini
[boundary.top]
velocity = fixedValue 1 0
pressure = zeroGradient
```

### Main settings

| Setting | Default | Meaning |
|---|---:|---|
| `physics` | `flow` | `flow` for cell-centered SIMPLE or `acoustic` for the nodal acoustic solver. |
| `flow_model` | `incompressible` | `incompressible`, `boussinesq`, or `compressible` (pressure-based ideal-gas SIMPLE). |
| `mesh` | required | Gmsh file, relative to the case directory unless absolute. |
| `output` | `results/flow.plt` | Tecplot output path. |
| `time_scheme` | `steady` | Flow: `steady`, `euler`, `backward_euler`, or `bdf2`. See the acoustic section for its schemes. |
| `density` | `1.0` | Constant density \(\rho\). Must be positive. |
| `viscosity` | `0.01` | Constant dynamic viscosity \(\mu\). |
| `initial_velocity` | `0 0` | Initial cell velocity, two components. |
| `initial_pressure` | `0.0` (`101325` compressible) | Initial pressure; compressible cases use positive absolute pressure. |
| `initial_temperature` | `300` | Initial absolute temperature for thermal models. |
| `turbulence_model` | `laminar` | `laminar`, `k_epsilon`, or `sst_k_omega`. |
| `initial_k` | `1e-6` | Initial turbulent kinetic energy; must exceed its floor. |
| `initial_epsilon` | `1e-6` | Initial dissipation for `k_epsilon`. |
| `initial_omega` | `1.0` | Initial specific dissipation for SST. |
| `turbulence_relaxation` | `0.5` | Under-relaxation for turbulence equations in `(0,1]`. |
| `turbulence_tolerance` | value of `tolerance` | Maximum normalized turbulence update. |
| `turbulence_ksp` | `bcgs` | PETSc KSP type for turbulence transport. |
| `turbulence_pc` | `ilu` | PETSc preconditioner for turbulence transport. |
| `maximum_eddy_viscosity_ratio` | `1000` | Numerical \(\mu_t/\mu\) safeguard for standard \(k\)-\(\varepsilon\). |
| `body_force` | `0 0` | Body force per unit mass, two components. |
| `gradient` | `least_squares` | `least_squares`/`lsq` or `gauss`/`green_gauss`. |
| `gradient_limiter` | `none` | `none` or `barth_jespersen` (`barth-jespersen` and `bj` are aliases), applied to convective reconstruction. |
| `velocity_relaxation` | `0.7` | Momentum under-relaxation in `(0,1]`. |
| `pressure_relaxation` | `0.3` | Pressure under-relaxation in `(0,1]`. |
| `maximum_iterations` | `2000` steady, `100` transient | Maximum SIMPLE iterations per steady solve or time step. |
| `tolerance` | `1e-6` | Maximum absolute velocity and pressure-update tolerance. |
| `continuity_tolerance` | `1e-8` | Maximum absolute cell mass-imbalance tolerance. |
| `linear_tolerance` | `1e-9` | PETSc relative linear residual tolerance. |
| `linear_absolute_tolerance` | `1e-12` | PETSc absolute linear residual tolerance. |
| `maximum_linear_iterations` | `1000` | Maximum iterations for each PETSc solve. |
| `momentum_ksp` | `bcgs` | PETSc KSP type for both momentum components. |
| `momentum_pc` | `ilu` | PETSc preconditioner for momentum. |
| `pressure_ksp` | `cg` | PETSc KSP type for pressure correction. |
| `pressure_pc` | `icc` | PETSc preconditioner for pressure correction. |
| `log_interval` | `10` | SIMPLE console-output interval. |

The outer `tolerance` and `continuity_tolerance` are absolute quantities. They should be reconsidered if the mesh, density, velocity scale, or units change substantially.

### Thermal and compressible settings

| Setting | Default | Meaning |
|---|---:|---|
| `material_model` | `constant` | Use `constant` for Boussinesq or `ideal_gas` for compressible flow. |
| `thermal_expansion` | `0` | Boussinesq \(\beta\) in 1/K; required to be positive. |
| `reference_temperature` | `initial_temperature` | \(T_0\) in the Boussinesq buoyancy source. |
| `gravity` | `0 0` | Gravity vector. Boussinesq pressure is hydrostatically modified; compressible flow applies \(\rho\mathbf g\) directly. |
| `thermal_conductivity` | `0` | Molecular conductivity \(k\) in W/(m K). |
| `specific_heat` | `0` | Constant \(c_p\) in J/(kg K); must be positive for energy transport. |
| `gas_constant` | `287` | Specific gas constant \(R\) for the single-component ideal gas. |
| `heat_source` | `0` | Uniform volumetric heat source in W/m³. |
| `viscous_dissipation` | `true` compressible, `false` Boussinesq | Include viscous dissipation in energy. |
| `energy_relaxation` | `0.7` | Energy-equation under-relaxation in `(0,1]`. |
| `turbulent_prandtl` | `0.9` | Turbulent Prandtl number for Boussinesq/RANS heat transport. |
| `energy_tolerance` | `tolerance` | Maximum absolute temperature update in K. |
| `density_tolerance` | `tolerance` | Maximum absolute density update in kg/m³ for compressible convergence. |
| `energy_ksp`, `energy_pc` | `bcgs`, `ilu` | Linear solver and preconditioner for energy. |

Every patch in a thermal case needs a `temperature` entry:

| Syntax | Behaviour |
|---|---|
| `temperature = fixedValue T` | Fixed positive absolute temperature in K. |
| `temperature = zeroGradient` | Zero normal temperature gradient; `adiabatic` is an alias. |
| `temperature = fixedHeatFlux q` | Prescribed heat flux in W/m², positive into the fluid. |

For compressible flow, `initial_pressure` and every `pressure = fixedValue ...` are absolute pressures in Pa. A steady compressible case requires at least one fixed-pressure boundary. Compressible mode currently supports laminar single-component ideal gas only; Boussinesq mode uses the existing laminar/RANS models with constant reference density.

### Transient settings

| Setting | Default | Meaning |
|---|---:|---|
| `dt` | `0.01` | Time-step size. |
| `end_time` | `1.0` | Physical end time. |
| `maximum_steps` | `1000000` | Safety limit on the number of time steps. |
| `allow_incomplete_time_steps` | `false` | Allow limited-iteration intermediate steps for transient-to-steady calculations. |
| `final_maximum_iterations` | `maximum_iterations` | Required fully converged iteration budget for the final step. |
| `restart` | disabled | Optional Tecplot file providing cell-centered `U`, `p`, and `T`, relative to the case directory. |
| `output_interval` | `10` | Write a Tecplot file every this many steps. |
| `force_patch` | disabled | Physical patch on which to integrate pressure and viscous force. |
| `force_output` | `results/forces.csv` | Force-history CSV path when `force_patch` is set. |
| `force_reference_velocity` | `1.0` | Reference velocity for `Cd` and `Cl`. |
| `force_reference_length` | `1.0` | Reference length/span for `Cd` and `Cl`. |
| `output_vorticity` | `false` | Add cell-centered `dv/dx - du/dy` to Tecplot output. |

The final step is shortened if necessary to finish exactly at `end_time`. Constant-step BDF2 is intended for a fixed `dt`; variable-step BDF2 coefficients are not yet implemented. Force histories include total, pressure, and viscous components. Coefficients use `0.5 * density * force_reference_velocity^2 * force_reference_length` as their normalization.

### Cavity example

```ini
mesh = mesh/cavity.msh
output = results/cavity.plt

time_scheme = steady
density = 1.0
viscosity = 0.01
gradient = least_squares

velocity_relaxation = 0.7
pressure_relaxation = 0.3
maximum_iterations = 3000
tolerance = 1.0e-6
continuity_tolerance = 1.0e-8
linear_tolerance = 1.0e-9
linear_absolute_tolerance = 1.0e-12
maximum_linear_iterations = 1000
log_interval = 25

[boundary.top]
velocity = fixedValue 1 0
pressure = zeroGradient

[boundary.walls]
velocity = noSlip
pressure = zeroGradient
```

## Boundary conditions

For flow cases, every Gmsh physical boundary patch must have both a `velocity` and a `pressure` entry. Acoustic patches require only their acoustic `pressure` entry, as described below. Patch names in `case.in` must exactly match the mesh's physical names, including case.

### Velocity conditions

| Syntax | Behaviour |
|---|---|
| `velocity = fixedValue ux uy` | Prescribed boundary velocity vector. |
| `velocity = noSlip` | Zero velocity at the boundary. |
| `velocity = zeroGradient` | Extrapolated velocity; normally paired with fixed pressure at an outlet. |
| `velocity = slip` | Free-slip/symmetry treatment with zero face-normal mass flux. |
| `velocity = symmetry` | Alias for `slip`. |
| `velocity = parabolic meanVelocity` | X-directed parabolic profile on a vertical inlet. |
| `velocity = parabolicX meanVelocity` | Explicit alias for `parabolic`. |

For a vertical inlet with \(y\) normalized between its endpoints as \(\eta\), the parabolic condition is

$$
u=6U_{\mathrm{mean}}\eta(1-\eta),\qquad v=0.
$$

Thus `meanVelocity` is the profile's mean value, not its centreline maximum.

### Pressure conditions

| Syntax | Behaviour |
|---|---|
| `pressure = fixedValue p` | Prescribed pressure and zero boundary pressure correction. |
| `pressure = zeroGradient` | Zero normal pressure-gradient boundary condition. |

A common open-boundary combination is:

```ini
[boundary.outlet]
velocity = zeroGradient
pressure = fixedValue 0
```

For a closed flow such as the cavity, all pressure boundaries can be `zeroGradient`; the solver then supplies the required discrete pressure reference.

## Material and multispecies configuration

For single-component compressible flow, use `flow_model = compressible` and `material_model = ideal_gas` without a `species` list; `gas_constant`, `specific_heat`, `thermal_conductivity`, and `viscosity` define the fluid. For inert multispecies transport, set `material_model = ideal_gas`, provide a `species` list, and set `species_model = transport`. In the latter passive module, the legacy `density` and `viscosity` settings remain the constant carrier properties used by momentum, continuity, and species convection. `flow_coupling = passive` is the only currently implemented species coupling. `temperature` and `thermo_pressure` are the passive material temperature and absolute thermodynamic pressure; `thermo_pressure` is not the incompressible SIMPLE gauge pressure.

| Setting | Default | Meaning |
|---|---:|---|
| `material_model` | `constant` | `constant`, `ideal_gas`, or `volume_fraction_mixture`. The last currently provides a property API only. |
| `species_model` | `none` | `none` or `transport`. |
| `flow_coupling` | `passive` | Passive constant-density/viscosity coupling; other values are rejected. |
| `species` | required for mixtures | Space-separated unique names; use letters, digits, and underscores, starting with a letter or underscore. |
| `diffusion_model` | `constant` | `constant` per-species diffusivities or `mixture_averaged` binary-pair mixing. |
| `temperature` | `300` | Constant material temperature in K; no energy equation is solved. |
| `thermo_pressure` | `101325` | Absolute pressure in Pa for the ideal-gas equation of state. |
| `initial_mass_fractions` | required | One non-negative value per species, summing to one within `1e-8`. |
| `species_relaxation` | `0.7` | Species linear-equation under-relaxation in `(0,1]`. |
| `turbulent_schmidt` | `0.7` | Turbulent Schmidt number; adds `mu_t/Sc_t` to each species diffusion coefficient. |
| `species_tolerance` | `tolerance` | Maximum absolute mass-fraction update. |
| `species_composition_tolerance` | `1e-8` | Allowed sum-closure error and bound tolerance at outer convergence. |
| `species_ksp`, `species_pc` | `bcgs`, `ilu` | Linear solver and preconditioner shared by the species solves. |
| `multiphase_model` | `none` | `vof`, `asm`, `euler_euler`, or `dpm` select interfaces that currently throw an explicit not-implemented exception. |

Each `[species.NAME]` section requires SI-valued `molar_mass` (kg/mol), `viscosity` (Pa s), `thermal_conductivity` (W/(m K)), `specific_heat` (J/(kg K)), and `diffusivity` (m²/s). For `diffusion_model = mixture_averaged`, every upper-triangular binary pair is also required. Pair settings use the first species as a section name:

```ini
material_model = ideal_gas
species_model = transport
flow_coupling = passive
species = N2 O2
diffusion_model = mixture_averaged
temperature = 300
thermo_pressure = 101325
density = 1.1719579191113527
viscosity = 1.8188965989347983e-5
initial_mass_fractions = 0.7670824885424203 0.23291751145757966
species_relaxation = 1.0
turbulent_schmidt = 0.7

[species.N2]
molar_mass = 0.0280134
viscosity = 1.76e-5
thermal_conductivity = 0.0259
specific_heat = 1040
diffusivity = 2.0e-5

[species.O2]
molar_mass = 0.0319988
viscosity = 2.04e-5
thermal_conductivity = 0.0266
specific_heat = 918
diffusivity = 2.0e-5

[binary_diffusion.N2]
O2 = 2.0e-5
```

Every flow boundary patch additionally needs one composition condition:

```ini
[boundary.inlet]
velocity = fixedValue 0.001 0
pressure = zeroGradient
mass_fractions = fixedValue 0.7670824885424203 0.23291751145757966

[boundary.outlet]
velocity = zeroGradient
pressure = fixedValue 0
mass_fractions = zeroGradient
```

All species on a patch share either `fixedValue` or `zeroGradient`; fixed values must sum to one. Zero-gradient outlets with substantial reverse flow need a more complete backflow treatment and are not a general outflow boundary model.

`VolumeFractionMixture` independently implements volume-weighted density, viscosity and thermal conductivity, and mass-weighted heat capacity from `phases` and `[phase.NAME]` data. It is available to future multiphase implementations; selecting it in the current executable without an active multiphase model is rejected. The full mathematical formulation and module ownership are documented in `textbook/Chap5-2-SIMPLE-scheme-multiphase-multispecies.md`, sections 5.11.7–5.11.9.

## Acoustic configuration

Select the nodal mixed-grid solver with `physics = acoustic`. It uses the convected pressure-wave equation and median-dual control areas from `textbook/Chap5-4-acoustic.md`, not the cell-centered SIMPLE discretization. The implementation and verification plan is `src/acoustic/IMPLEMENTATION_PLAN.md`.

### Main acoustic settings

| Setting | Default | Meaning |
|---|---:|---|
| `density` | `1.0` | Uniform background density; must be positive. |
| `sound_speed` | `1.0` | Uniform background sound speed \(c_0\); must be positive. |
| `mean_velocity` | `0 0` | Uniform background velocity; `|u0| < c0` is required. |
| `time_scheme` | `newmark` | `newmark` or `central_difference` (`central`/`explicit` aliases). |
| `newmark_beta` | `0.25` | Newmark displacement coefficient; must be positive. |
| `newmark_gamma` | `0.5` | Newmark velocity coefficient; must be positive. |
| `dt` | automatic | Target maximum time step. If zero/negative, computed from `cfl`. |
| `cfl` | `0.4` | Automatic-step Courant factor using minimum element edge and `c0 + |u0|`. |
| `end_time` | `1.0` | Physical end time. The solver uses equal steps `end_time/ceil(end_time/dt)`. |
| `output_interval` | final step | Write numbered Tecplot files every this many steps and at the final step. |
| `probe_position` | disabled | Two coordinates; the nearest node is written to `<stem>_probe.csv`. |
| `probe_positions` | disabled | Space-separated `x y` pairs; all nearest-node histories are written to `<stem>_probes.csv`. |
| `acoustic_ksp` | `bcgs` | PETSc KSP for non-diagonal acoustic systems. |
| `acoustic_pc` | `ilu` | PETSc preconditioner for non-diagonal acoustic systems. |
| `linear_tolerance` | `1e-10` | PETSc relative residual tolerance. |
| `linear_absolute_tolerance` | `1e-13` | PETSc absolute residual tolerance. |
| `maximum_linear_iterations` | `1000` | Maximum PETSc iterations per time step. |

### Initial conditions and sources

`initial_condition` accepts `zero`, `cavity_mode`, `plane_wave`, or `gaussian_pulse`. Relevant settings are:

| Setting | Default | Meaning |
|---|---:|---|
| `initial_amplitude` | `1.0` | Pressure amplitude for the selected analytical field. |
| `initial_pressure` | `0.0` | Uniform pressure for `zero`. |
| `initial_pressure_rate` | `0.0` | Uniform initial pressure rate; also used by `cavity_mode`. |
| `mode_indices` | `2 1` | Non-negative integer cavity-mode indices, not both zero. Domain extents are inferred from the mesh. |
| `wave_direction` | `1 0` | Plane-wave or pulse propagation direction; normalized internally. |
| `wave_number` | `2*pi` | Plane-wave wavenumber magnitude. `omega = c0 |k| + u0·k`. |
| `wave_phase` | `0.0` | Plane-wave initial phase in radians. |
| `gaussian_center` | domain center | Two coordinates for the directed pulse center. |
| `gaussian_width` | `0.1` | Gaussian standard deviation; must be positive. |
| `source_type` | `none`/`harmonic` | `none`, `harmonic` monopole, or `flow_dipole` from a force-history CSV. Defaults to `harmonic` when `source_amplitude` is nonzero. |
| `source_amplitude` | `0.0` | Integrated harmonic source amplitude; zero disables the source. |
| `source_file` | required for `flow_dipole` | CSV with increasing `time`, `force_x`, and `force_y` columns, relative to the case directory. |
| `source_scale` | `1.0` | Multiplier applied to the flow-force history. |
| `source_ramp_time` | `0.0` | Cosine startup-ramp duration; zero disables ramping. |
| `source_position` | domain center | Gaussian source center. |
| `source_width` | `0.1` | Gaussian source standard deviation. |
| `source_frequency` | `1.0` | Harmonic source frequency in cycles per unit time. |
| `source_phase` | `0.0` | Harmonic source phase in radians. |

The Gaussian source is normalized over nodal dual areas so its discrete weights sum to one. The directed Gaussian pulse uses `p_t = -(c0 + u0·n) n·grad(p)` initially. For `flow_dipole`, the implementation evaluates `div(F g / rho)`, where `F` is the force of the fluid on the body and the force of the body on the fluid is `-F`. `cases/run_cylinder_aeolian.py` supplies this convention after mean removal. If `end_time` is omitted for a flow dipole, the last source timestamp determines it.

### Acoustic boundaries

Every physical patch needs one `pressure` entry:

```ini
[boundary.left]
pressure = planeWave

[boundary.right]
pressure = absorbing

[boundary.bottom]
pressure = rigid

[boundary.top]
pressure = fixedValue 0
```

- `rigid` (`hardWall`, `wall`, and `zeroGradient` aliases): zero normal acoustic flux.
- `absorbing` (`nonreflecting`/`impedance` aliases): local outgoing characteristic relation. It is exact only in the normal-incidence cases for which it is derived; oblique and near-field reflection must be assessed separately.
- `fixedValue p`: prescribed constant pressure with degree-of-freedom elimination.
- `planeWave`: prescribed analytical pressure using the global wave settings; also provides consistent initial boundary rates/accelerations.

### Cavity example

```ini
physics = acoustic
mesh = mesh/medium_mixed.msh
output = results/cavity.plt
density = 1.0
sound_speed = 1.0
mean_velocity = 0 0
time_scheme = central_difference
dt = 0.005
end_time = 1.4142135623730951
initial_condition = cavity_mode
mode_indices = 2 1
probe_position = 0.25 0.25

[boundary.left]
pressure = rigid
[boundary.right]
pressure = rigid
[boundary.bottom]
pressure = rigid
[boundary.top]
pressure = rigid
```

Acoustic Tecplot zones contain nodal `p`, `p_t`, `p_tt`, and `dual_area`. Each run also writes `<stem>_summary.csv` with pressure extrema, quadratic energy, linear diagnostics, and Courant number, plus `<stem>_geometry.csv` with element area closure and `sum(grad N)` checks.

## Gmsh meshes

### Supported input

- ASCII Gmsh format 2.2 or 4.1.
- Two-dimensional linear triangles and quadrilaterals.
- Linear line elements for boundary faces.
- Physical names for robust boundary assignment.
- Z coordinates are read but ignored by the two-dimensional geometry.

Curved/high-order elements, binary Gmsh files, three-dimensional elements, and mixed 2D/3D models are not supported. Every geometric boundary edge must have a corresponding Gmsh line element; otherwise mesh construction stops with an error rather than silently omitting a boundary condition.

### Supplied meshes

Regenerate the flow and acoustic meshes with:

```bash
make meshes
```

Alternatively, run `python3 cases/generate_meshes.py` for flow or `python3 cases/generate_acoustic_meshes.py` for acoustics. The generators use only Python's standard library and write ASCII Gmsh 2.2 files.

| Case | Geometry and mesh |
|---|---|
| Cavity | Unit square, 3,200 triangles. |
| Backstep | Inlet region `[-1,0] × [0.5,1]` and expanded region `[0,8] × [0,1]`, 4,896 triangles. |
| Cylinder | Diameter-one cylinder at the origin, outer domain `[-4,16] × [-4,4]`, 2,304 body-fitted quadrilaterals. |
| Flat plate | Symmetry run-up and plate `0 ≤ x ≤ 2`; 4,500 / 9,380 / 20,000 quadrilaterals by level. |
| Elbow | Planar 90° annulus with `Rc/D = 2.8`, 5D upstream and 10D downstream; 3,960 / 8,500 / 16,800 quadrilaterals. |

The cylinder generator classifies inner-wall edges using both endpoint radii. This avoids mistakenly assigning a chord's midpoint to a far-field patch. Plate and elbow meshes use different first-cell heights for wall-function and wall-resolved models.

## Output files

### Tecplot fields

The default output is ASCII Tecplot. Nodal coordinates and cell-centered fields use finite-element triangle or quadrilateral zones:

| Variable | Location | Meaning |
|---|---|---|
| `X`, `Y` | Nodes | Mesh coordinates. |
| `U_x`, `U_y` | Cell centres | Velocity components. |
| `p` | Cell centres | Static pressure; hydrostatically modified for Boussinesq and absolute for compressible flow. |
| `T` | Cell centres | Absolute temperature for thermal models; recovered from total enthalpy in compressible flow. |
| `H0` | Cell centres | Total enthalpy \(c_pT+|\mathbf U|^2/2\) for compressible flow. |
| `rho` | Cell centres | Reference or equation-of-state density for thermal models. |
| `Mach` | Cell centres | Local Mach number for compressible models. |
| `lambda`, `cp` | Cell centres | Molecular conductivity and specific heat for thermal models. |
| `mass_imbalance` | Cell centres | Signed integrated cell mass imbalance, including transient density accumulation for compressible flow. |
| `vorticity` | Cell centres | `dv/dx - du/dy`, enabled by `output_vorticity = true`. |
| `k` | Cell centres | Turbulent kinetic energy, for turbulent cases. |
| `epsilon` | Cell centres | Dissipation, for `k_epsilon`. |
| `omega` | Cell centres | Specific dissipation, for SST. |
| `mu_t` | Cell centres | Turbulent dynamic eddy viscosity. |
| `F1`, `F2` | Cell centres | SST blending functions. |
| `Y_NAME` | Cell centres | Mass fraction of species `NAME`, for multispecies cases. |
| `D_NAME` | Cell centres | Molecular mixture-averaged diffusivity of species `NAME`. |
| `rho_mix`, `mu_mix` | Cell centres | Material density and dynamic viscosity; diagnostic in the current passive coupling. |
| `lambda_mix`, `cp_mix` | Cell centres | Material thermal conductivity and mass-specific heat capacity. |

Mixed meshes are written as separate triangle and quadrilateral zones. A steady case writes the configured file directly. A transient case appends the time-step number to the stem, for example `flow_100.plt`.

Turbulent cases also write `<output-stem>_wall.csv` with patch name, face-center coordinates, wall distance, shear-stress vector, friction velocity, and \(y^+\) for every no-slip face.

### Convergence summary

Each run also writes `<output-stem>_summary.csv` with:

- `time`
- `iterations`
- `velocity_change`
- `pressure_change`
- `max_cell_imbalance`
- `net_boundary_flux`
- `turbulence_change`
- `turbulence_linear_iterations`
- `species_change`
- `species_composition_error`
- `minimum_mass_fraction`, `maximum_mass_fraction`
- `species_linear_iterations`
- `temperature_change`, `density_change`, `energy_linear_iterations`
- `converged`

The console log uses `dU` for the maximum cell-velocity update, `dp` for the maximum relaxed pressure update, `dTurb` for the normalized turbulence update, `continuity` for the maximum absolute cell mass imbalance, and `linear` for the total PETSc iterations used by that outer iteration's momentum, pressure, energy, turbulence, and species solves. Thermal cases additionally print `dTemp` and `drho`; multispecies cases print `dY` and `Ysum`. Single-component laminar output retains zero turbulence diagnostics.

## Validation cases

Run the complete workflow:

```bash
make validate
```

This runs the three laminar regressions, all twelve plate/elbow turbulence cases, the analytical acoustic matrix, the ten multispecies/material cases, and the eight thermal/compressible cases. `cases/validate.py` checks the laminar cases; `cases/validate_turbulence.py --run` executes the turbulence matrix and writes `validation/metrics.csv`, comparison CSVs, and SVG plots. The detailed model-validation reports are under `validation/`.

### Current results

| Case | Cells | SIMPLE iterations | Maximum cell imbalance | Selected result |
|---|---:|---:|---:|---|
| Cavity | 3,200 | 1,173 | `4.998e-16` | Interpolated `u(0.5,0.5) = -0.24359` |
| Backstep | 4,896 | 1,488 | `1.390e-15` | Near-wall minimum `u = -0.07173`; estimated reattachment `x = 2.35117` |
| Cylinder | 2,304 | 333 | `2.320e-14` | Wake minimum `u = -0.10540`; estimated reattachment `x = 2.16189` |

The absolute net boundary mass fluxes are zero for the closed cavity to printed precision, `1.051e-14` for backstep, and `9.548e-15` for cylinder.

### Checks performed

For every case, the script requires:

- SIMPLE convergence.
- Finite velocity and pressure values.
- Maximum cell mass imbalance below `1e-8`.
- Absolute net boundary mass flux below `1e-8`.

Case-specific checks are:

- **Cavity:** inverse-distance interpolation estimates the horizontal velocity at `(0.5,0.5)`. It is compared with the Ghia-type Re = 100 value `-0.20581`, using an absolute error allowance of `0.08` on this coarse mesh. The current absolute difference is approximately `0.03778`.
- **Backstep:** a negative near-wall velocity and a subsequent reattachment inside the downstream domain are required.
- **Cylinder:** a negative centreline wake velocity and subsequent reattachment are required.

These are useful regression and physical-sanity checks, but they are **not a grid-convergence study**. The backstep and cylinder checks are currently qualitative; they are not yet quantitative drag, lift, or reattachment-length benchmark certifications.

The supplied physical parameters are:

| Case | Parameters |
|---|---|
| Cavity | Lid speed `1`, `rho = 1`, `mu = 0.01`, giving Re = 100 for unit length. |
| Backstep | Mean inlet speed `1`, inlet height `0.5`, `rho = 1`, `mu = 0.005`, giving Re = 100 based on inlet height. |
| Cylinder | Freestream speed `1`, diameter `1`, `rho = 1`, `mu = 0.025`, giving Re = 40. |

### Turbulence validation

The flat-plate study compares both models with NASA's Kármán–Schoenherr \(C_f(Re_\theta)\) relation and Coles \(u^+(y^+)\) profile at \(Re_\theta=10{,}000\), using three grid/near-wall-spacing levels. The elbow study compares the 75° streamwise section with Enayet's experimental points, also at three levels. All turbulence cases require coupled convergence, finite positive turbulence fields, local/global mass conservation, and appropriate near-wall spacing. Flat-plate normalized RMSE above 10% is a validation failure; elbow errors are reported without a hard accuracy gate because the solver is planar while the experiment is three-dimensional.

The reference files, digitization uncertainty, and source URLs are under `reference/`. See `validation/turbulence_validation.md` for final values, grid sensitivity, cap usage, and the limitations of comparing a 2D elbow with a 3D pipe experiment.

### Acoustic analytical validation

Run `make validate-acoustic` for 14 cavity, plane-wave, mean-flow, and pulse/outlet runs on triangle, quadrilateral, and mixed grids. The driver checks geometry closure, analytical field and frequency errors, grid convergence, energy drift, and outlet reflection. It writes `validation/acoustic_metrics.csv`, comparison CSVs under `validation/data/`, and SVG contours/profiles under `validation/figures/`.

| Finest-grid or refined case | Selected result |
|---|---|
| Rigid cavity `(2,1)` mode | Relative L2 error `0.06281%`; observed grid orders `1.954`, `2.239` |
| Stationary plane wave | Relative L2 error `0.10500%`; phase speed `0.999696` versus exact `1` |
| Plane wave with `M = 0.3` | Relative L2 error `0.24392%`; phase speed `1.299116` versus exact `1.3` |
| Gaussian pulse at `t = 1` | Relative L2 error `2.05364%` |
| Absorbing versus rigid outlet | Reflected/incident peak ratios `0.234999%` versus `99.558540%` |

See `validation/acoustic_validation.md` for the complete metrics, contours, comparison plots, acceptance criteria, references, and limitations. `make validate` runs the laminar, turbulence, acoustic, and multispecies checks.

### Multispecies analytical and material validation

`make validate-multispecies` runs three binary Stefan-tube grids, three unequal-diffusivity flux-correction grids, three transient three-species diffusion grids, and one uniform N₂/O₂ property case. The driver checks SIMPLE convergence, analytical mass-fraction profiles, grid-error reduction, boundedness, composition closure, continuity, ideal-gas density, Wilke viscosity, thermal conductivity, heat capacity, and diffusivity. It writes `validation/multispecies_metrics.csv`, comparison CSVs, and SVG profiles/convergence plots.

The finest-grid relative L2 errors are approximately `3.41e-5` for Stefan diffusion, `4.14e-5` for unequal-diffusivity flux correction, and `2.28e-4` for transient three-species diffusion, with observed grid orders near two. The uniform air case gives `rho = 1.1719579191 kg/m³` and `mu = 1.8188965989e-5 Pa s`. These are analytical and material-rule verifications, not full Maxwell–Stefan, reacting-flow, or composition-driven low-Mach validation. See `validation/multispecies_validation.md` for equations, measured metrics, acceptance criteria, references, and limitations.

### Boussinesq and compressible-flow validation

`make validate-thermal` runs three de Vahl Davis natural-convection grids (`Ra = 10^4`, `Pr = 0.71`), zero-flow ideal-gas conduction, transient closed-volume heating, and three subsonic isentropic-nozzle grids. The driver checks coupled convergence, mass conservation, analytical temperature/pressure/density profiles, Nusselt number, mid-plane velocity maxima, and grid refinement. The finest cavity gives `Nu = 2.27680` versus `2.243` (`1.507%` error). The finest nozzle has pressure, density, temperature, and velocity relative L2 errors of `0.969%`, `1.003%`, `0.0741%`, and `0.695%`, respectively. See [`validation/thermal_compressible_validation.md`](validation/thermal_compressible_validation.md) for equations, all eleven runs, acceptance criteria, source links, and model limitations.

### Bump-flow validation

`make validate-bump` runs the 10% circular-arc bump at inlet Mach `0.5` (subsonic) and `0.675` (transonic), and the 4% bump at Mach `1.4` (supersonic) using the pressure-based compressible SIMPLE solver. All three use `90 x 24` body-fitted quadrilateral grids and second-order upwind convection with `gradient_limiter = barth_jespersen`. The subsonic case checks fore/aft symmetry and isentropic exit Mach number; the transonic case compares exit Mach numbers `0.670` and `0.615` from Favini et al.; the supersonic case checks shock topology and approximate wall Mach profiles. The clean-run results and source links are documented in `validation/thermal_compressible_validation.md` and `cases/compressible_bump/README.md`.

### Cylinder vortex shedding and aeolian tone

`make cylinder-aeolian` runs the hybrid case under `cases/cylinder_aeolian/`: transient `Re = 150` flow, force-history extraction over three measured shedding periods, and acoustic propagation at `M = 0.2`. The workflow compares `St`, mean `Cd`, lift amplitude, and drag second harmonic with Inoue and Hatakeyama's published values and checks transverse-dipole directivity. Case equations, source conventions, commands, and limitations are in `cases/cylinder_aeolian/README.md`; final metrics and figures are described in `validation/cylinder_aeolian.md`.

## Adding a new case

1. Create a directory under `cases/`.
2. Add a supported ASCII Gmsh mesh with named physical line boundaries.
3. Create `case.in` with at least `mesh`, material properties, and one boundary section per physical patch.
4. Choose the gradient, convection, relaxation, and convergence settings.
5. Run:

   ```bash
   ./multiphysics cases/your-case
   ```

6. Inspect the console convergence, `<case>_summary.csv`, and the Tecplot result.
7. Check both local cell mass imbalance and global boundary mass conservation before interpreting the flow field.

For difficult cases, enable `gradient_limiter = barth_jespersen`, lower the SIMPLE relaxation factors, or improve the mesh. Tightening linear-solver tolerances does not by itself guarantee nonlinear SIMPLE convergence.

## Adding another physics solver

The general FVM layer is intended to be shared by flow, elastic, and acoustic modules. A new solver can reuse:

- `Mesh`, `Cell`, `Face`, and `BoundaryPatch` for geometry and connectivity.
- `ScalarField` and `VectorField` (or `Field<double>` and `Field<Vector2>`) for unknowns and boundary data.
- `GradientScheme` for scalar/vector reconstruction.
- `TimeDerivative` for time-term coefficients.
- Interpolation and flux utilities in `Operators.hpp`.
- `SparseMatrix` and `PetscLinearSystem` for algebraic systems.
- `TecplotWriter` for output.
- `PhysicsSolver` as the basic lifecycle interface:

```cpp
class PhysicsSolver {
public:
    virtual ~PhysicsSolver() = default;
    virtual std::string name() const = 0;
    virtual void initialize() = 0;
    virtual void advance(double dt) = 0;
    virtual void write(const std::filesystem::path& file) const = 0;
};
```

`src/main.cpp` selects `SimpleSolver` or `AcousticSolver` from `physics` in `case.in`. Inter-solver field exchange and coupled scheduling still need to be added before this becomes a coupled multiphysics executable.

## Current limitations

- Serial execution only; there is no domain decomposition or distributed PETSc matrix assembly.
- Two-dimensional linear triangles and quadrilaterals only.
- Incompressible/Boussinesq transport uses reference density; ideal-gas compressible flow updates density from pressure and temperature. Passive species and mixture properties do not feed variable properties back into SIMPLE.
- Compressible flow is limited to single-component laminar ideal gas with constant heat capacity and transport properties; real-gas, temperature-dependent-property, and compressible-turbulence models are unavailable. All compressible cases use pressure-based SIMPLE with second-order upwind convection. The optional Barth–Jespersen limiter controls reconstruction overshoots, but shock thickness and peak values remain grid- and dissipation-dependent.
- Multispecies transport is inert and isothermal, with constant or mixture-averaged Fick diffusion; no Maxwell–Stefan matrix, Soret/Dufour effects, reactions, or species-energy coupling.
- VOF, algebraic-slip, Euler–Euler, and DPM modules are interfaces only; no multiphase field or particle solve is implemented.
- Two-dimensional RANS only: three-dimensional secondary flow, swirl, and cross-sectional Dean vortices cannot be represented.
- Only standard \(k\)-\(\varepsilon\) and SST-2003; no transition, compressibility correction, buoyancy production, curvature/rotation correction, DES, or Reynolds-stress model.
- No elastic field solver, fluid–structure coupling, or two-way fluid–acoustic coupling yet. The cylinder case is one-way Curle-type coupling.
- Acoustic media and mean flow are currently uniform and steady; the implementation requires `|M| < 1`.
- Acoustic absorption is a local first-order condition, not a PML; general frequency-dependent impedance and frequency-domain Helmholtz solves are not implemented.
- No convection limiter, so central reconstruction can produce oscillations in strongly convection-dominated or poor-quality-mesh cases.
- No curved/high-order Gmsh element support.
- No binary Gmsh input or binary Tecplot output.
- Restart is limited to cell-centered `U`, `p`, and `T` from an existing Tecplot result; no portable checkpoint metadata or parallel restart support.
- No variable-time-step BDF2 coefficients.
- Grid-refinement studies are case-specific; there is no complete Navier–Stokes manufactured-solution or temporal-order study yet.
- The supplied validation meshes are deliberately modest in size and should not be treated as grid-independent production results.

## To-do work

### Physics modules and multiphysics coupling

- [x] Add independent material models and passive multispecies transport with conservative flux closure.
- [x] Define VOF, algebraic-slip, Euler–Euler, and DPM module interfaces.
- [ ] Implement the multiphase solve routines and composition-driven low-Mach pressure coupling.
- [ ] Add Maxwell–Stefan diffusion, species energy, and reaction-source interfaces.
- [ ] Implement the two-dimensional elastic solid solver, including displacement, stress, traction, and fixed-displacement boundaries.
- [x] Implement the two-dimensional acoustic solver with nodal pressure and dual-control-volume fluxes.
- [x] Add configuration-driven physics selection in `main.cpp`.
- [ ] Design shared field registration and solver scheduling for coupled calculations.
- [ ] Add conservative field/force interpolation between different physics meshes.
- [ ] Implement fluid–structure interface conditions and an explicit or partitioned fluid–structure coupling algorithm.
- [ ] Add PML, frequency-dependent impedance, vibrating-wall input, and frequency-domain acoustic solutions.
- [ ] Investigate strong/monolithic coupling where explicit partitioned coupling is unstable.

### Mesh and geometry

- [ ] Add curved Gmsh elements and high-order geometry representation.
- [ ] Support binary Gmsh files if input size becomes a limitation.
- [ ] Add mesh-quality reporting for aspect ratio, skewness, non-orthogonality, and minimum face delta.
- [ ] Add mesh-refinement and grid-convergence workflows.
- [ ] Consider general polygonal control volumes and hanging-node meshes.
- [ ] Add moving-mesh support for elastic boundaries and fluid–structure interaction.

### Finite-volume numerics

- [x] Add bounded second-order upwind convection with an optional Barth–Jespersen gradient limiter.
- [ ] Add additional convection schemes such as QUICK.
- [ ] Improve robustness of non-orthogonal corrections on highly skewed meshes.
- [x] Add ideal-gas density and cell-property fields to SIMPLE.
- [ ] Add temperature-dependent viscosity, conductivity, and heat capacity where appropriate.
- [ ] Add variable-step BDF2 and adaptive time-step control.
- [ ] Add residual normalization based on problem scales in addition to absolute update tolerances.
- [ ] Review and expand outlet/backflow treatments for strongly recirculating outflows.
- [ ] Add Robin and more general mixed boundary conditions to the shared discretization layer.

### Linear algebra and parallel computing

- [ ] Replace sequential matrices with distributed PETSc `MPIAIJ` matrices.
- [ ] Add mesh partitioning, local/global numbering, ghost cells, and halo exchange.
- [ ] Support distributed gradient reconstruction and conservative face-flux assembly.
- [ ] Add configurable parallel preconditioners and solver prefixes for individual equations.
- [ ] Study linear and nonlinear convergence scaling with mesh size and rank count.
- [ ] Add parallel-aware output without concurrent writes to the same files.

### Verification and validation

- [ ] Add unit tests for polygon geometry, face orientation, Gmsh 2.2/4.1 reading, and boundary connectivity.
- [ ] Verify that Gauss and least-squares gradients reproduce exact linear fields on disturbed meshes.
- [ ] Add manufactured-solution tests for convection, diffusion, time integration, and complete Navier–Stokes residuals.
- [ ] Measure spatial convergence on at least three systematically refined meshes.
- [ ] Verify backward-Euler and BDF2 temporal convergence rates.
- [ ] Compare complete cavity centreline profiles with benchmark data.
- [ ] Add quantitative backstep reattachment comparisons over a Reynolds-number range.
- [ ] Add cylinder drag, lift, wake length, and pressure coefficients.
- [ ] Add transient cylinder cases with lift history and Strouhal-number validation above the steady-wake regime.
- [ ] Add automated regression runs in a clean WSL environment.

### Usability and productivity

- [x] Add checkpoint/restart support for steady and transient calculations.
- [ ] Add optional binary or parallel Tecplot output for large meshes.
- [ ] Add command-line overrides for frequently changed case parameters.
- [ ] Add more diagnostic fields such as vorticity, Courant number, forces, and boundary flux reports.
- [ ] Add example transient and coupled-physics cases as those features become available.

## Troubleshooting

### `mpic++` is not found

Set the MPICH binary directory explicitly before running Make:

```bash
export PATH=/home/lwa/program/mpich-4.2.0/bin:$PATH
```

### PETSc configuration or libraries cannot be found

Confirm that this file exists:

```text
/home/lwa/program/petsc-3.21.0/lib/petsc/conf/variables
```

If PETSc is installed elsewhere, set `PETSC_DIR` or pass it to Make:

```bash
make PETSC_DIR=/your/petsc/prefix
```

For shared-library loading errors, include PETSc, LAPACK, and MPICH library directories in `LD_LIBRARY_PATH` as shown in the build section.

### `Only ASCII Gmsh files are supported`

Export the mesh as ASCII. The reader supports Gmsh 2.2 and 4.1, but not binary `.msh` files.

### `Only linear Gmsh lines, triangles, and quadrilaterals are supported`

Generate first-order two-dimensional elements, or add high-order support before using curved elements.

### `A geometric boundary face has no Gmsh line element`

Every exterior and internal physical boundary edge must belong to a Gmsh line element. Define the corresponding curves and physical groups in the geometry/mesh generator.

### `Missing required case setting: boundary...`

Check that every physical curve name has a matching section such as `[boundary.walls]`, and that both `velocity` and `pressure` are specified. Names are case-sensitive.

### SIMPLE converges slowly or diverges

- Check mesh quality and boundary-condition consistency.
- Reduce `velocity_relaxation` and `pressure_relaxation`.
- Enable `gradient_limiter = barth_jespersen` if reconstructed convection produces local overshoots.
- Confirm that at least one pressure reference exists, either explicitly or through the automatic closed-domain reference.
- Inspect whether the PETSc KSP itself is diverging; nonlinear relaxation cannot repair an invalid mesh or an inconsistent linear system.

### The residual is small but the result looks unphysical

Check the actual mesh boundary patches, local and global mass conservation, Reynolds-number definition, and whether the grid is adequate. Algebraic convergence alone does not establish physical accuracy.
