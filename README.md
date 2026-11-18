# 2D Multiphysics Finite-Volume Solver

A C++17, cell-centered finite-volume project for two-dimensional multiphysics calculations on unstructured meshes. The current executable implements an **incompressible Navier–Stokes solver** using the SIMPLE pressure–velocity coupling algorithm. The mesh, field, discretization, linear-algebra, and output components are reusable foundations for the planned elastic and acoustic solvers.

The project is intended for both numerical-method development and teaching. It emphasizes conservative face fluxes, explicit mesh connectivity, readable operator implementations, and reproducible benchmark cases.

> **Current scope:** the working physics solver is serial, two-dimensional, constant-property, incompressible flow. Elasticity, acoustics, fluid–structure interaction, and MPI-distributed calculations are planned work; they are not yet implemented as coupled solvers.

## Contents

- [Implemented features](#implemented-features)
- [Numerical formulation](#numerical-formulation)
- [Project layout](#project-layout)
- [Prerequisites](#prerequisites)
- [Build and run in WSL](#build-and-run-in-wsl)
- [Case configuration](#case-configuration)
- [Boundary conditions](#boundary-conditions)
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
- Central face reconstruction with an implicit upwind/deferred-correction convection operator.
- Over-relaxed orthogonal diffusion with explicit non-orthogonal correction.
- PETSc sparse-matrix, vector, KSP, and preconditioner wrappers.
- Tecplot ASCII output for unstructured triangles and quadrilaterals.
- A `PhysicsSolver` interface for future physics modules.

### Incompressible Navier–Stokes solver

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
- Standard high-Reynolds-number \(k\)-\(\varepsilon\) turbulence with logarithmic-law wall functions.
- Menter SST-2003 \(k\)-\(\omega\) turbulence with strain-based eddy-viscosity limiting and low-\(y^+\) wall integration.
- Conservative transport of turbulence variables with bounded implicit upwind convection and optional deferred correction.
- Turbulence fields, wall shear, friction velocity, and \(y^+\) diagnostics in Tecplot and CSV output.

## Numerical formulation

### Governing equations

The solver uses constant density \(\rho\) and constant dynamic viscosity \(\mu\). The incompressible momentum and continuity equations are

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

### Gradient reconstruction

Two cell-gradient options are available:

- `gradient = gauss`: Green–Gauss integration using face-interpolated values, followed by reconstructed-face correction sweeps.
- `gradient = least_squares`: a weighted least-squares fit to neighbouring cell differences, with Dirichlet boundary constraints where applicable. This is the default.

The same gradient interface is used for scalar fields and both components of vector fields.

### Convection

The implicit matrix uses first-order upwind convection because of its robust diagonal properties. A deferred correction restores the selected higher-order face value during the SIMPLE iterations:

$$
F_f^{\mathrm{HO}}=F_f^{\mathrm{upwind}}
+\beta\left(F_f^{\mathrm{central}}-F_f^{\mathrm{upwind}}\right).
$$

The `convection_blending` parameter is \(\beta\):

- `1.0`: reconstructed central convection; the intended second-order smooth-flow option.
- `0.0`: first-order upwind convection.
- Values between zero and one: a blend of the two.

The correction is evaluated with the current outer-iteration field, while the conservative upwind part is assembled into the PETSc matrix. No flux limiter is currently applied.

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
4. Assemble and solve the pressure-correction equation so that the corrected face fluxes satisfy discrete continuity.
5. Correct cell velocities, pressure, and face mass fluxes.
6. Check velocity updates, pressure updates, and cell mass imbalances.

The Rhie–Chow normal velocity contains the difference between the face-normal pressure derivative and the interpolated cell pressure gradient. Its momentum weighting uses the relaxed momentum diagonal.

If no fixed-pressure boundary exists, the pressure-correction equation has a constant nullspace. The code removes it by fixing the pressure correction in cell zero. With a fixed-pressure outlet, no additional reference cell is needed.

### Linear systems

`src/fvm/PetscLinearSystem.hpp` provides the sparse-matrix and PETSc interface. The present implementation uses sequential `MATSEQAIJ` matrices and sequential vectors. Momentum and pressure-correction systems have separate persistent matrix/vector/KSP objects, but their nonzero values are refreshed during outer iterations.

Default linear-solver choices are:

| Equation | KSP | Preconditioner |
|---|---|---|
| Momentum | `bcgs` | `ilu` |
| Pressure correction | `cg` | `icc` |

They can be changed in `case.in`. Standard PETSc command-line options are also passed through the solver initialization. Turbulence equations default to `bcgs` with `ilu`.

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
turbulence_convection_blending = 0.0
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
│   └── flow/
│       ├── BoundaryConditions.* Flow boundary-condition handling
│       ├── TurbulenceModel.*   Standard k-epsilon and SST-2003 implementations
│       └── SimpleSolver.*      SIMPLE, Rhie–Chow, and RANS assembly
├── cases/
│   ├── generate_meshes.py      Reproducible supplied Gmsh meshes
│   ├── validate.py             Laminar conservation and benchmark checks
│   ├── validate_turbulence.py  Turbulence runs, reference errors, and plots
│   ├── cavity/                 Lid-driven cavity configuration and mesh
│   ├── backstep/               Backward-facing step configuration and mesh
│   ├── cylinder/               Flow-around-cylinder configuration and mesh
│   ├── plate_*/                Flat-plate cases for both turbulence models
│   └── elbow_*/                90-degree elbow cases for both models
├── reference/                 NASA and Enayet/Ansys comparison data
└── validation/               Turbulence report, metrics, data, and figures
```

Additional teaching material is under `textbook/` and `texbook-ppt/`. The empty physics directories are reserved for future development.

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
| `mesh` | required | Gmsh file, relative to the case directory unless absolute. |
| `output` | `results/flow.plt` | Tecplot output path. |
| `time_scheme` | `steady` | `steady`, `euler`, `backward_euler`, or `bdf2`. |
| `density` | `1.0` | Constant density \(\rho\). Must be positive. |
| `viscosity` | `0.01` | Constant dynamic viscosity \(\mu\). |
| `initial_velocity` | `0 0` | Initial cell velocity, two components. |
| `initial_pressure` | `0.0` | Initial cell pressure. |
| `turbulence_model` | `laminar` | `laminar`, `k_epsilon`, or `sst_k_omega`. |
| `initial_k` | `1e-6` | Initial turbulent kinetic energy; must exceed its floor. |
| `initial_epsilon` | `1e-6` | Initial dissipation for `k_epsilon`. |
| `initial_omega` | `1.0` | Initial specific dissipation for SST. |
| `turbulence_relaxation` | `0.5` | Under-relaxation for turbulence equations in `(0,1]`. |
| `turbulence_convection_blending` | `0.0` | Zero for upwind, one for deferred central correction. |
| `turbulence_tolerance` | value of `tolerance` | Maximum normalized turbulence update. |
| `turbulence_ksp` | `bcgs` | PETSc KSP type for turbulence transport. |
| `turbulence_pc` | `ilu` | PETSc preconditioner for turbulence transport. |
| `maximum_eddy_viscosity_ratio` | `1000` | Numerical \(\mu_t/\mu\) safeguard for standard \(k\)-\(\varepsilon\). |
| `body_force` | `0 0` | Body force per unit mass, two components. |
| `gradient` | `least_squares` | `least_squares`/`lsq` or `gauss`/`green_gauss`. |
| `convection_blending` | `1.0` | Zero for upwind, one for central, or an intermediate blend. |
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

### Transient settings

| Setting | Default | Meaning |
|---|---:|---|
| `dt` | `0.01` | Time-step size. |
| `end_time` | `1.0` | Physical end time. |
| `maximum_steps` | `1000000` | Safety limit on the number of time steps. |
| `output_interval` | `10` | Write a Tecplot file every this many steps. |

The final step is shortened if necessary to finish exactly at `end_time`. Constant-step BDF2 is intended for a fixed `dt`; variable-step BDF2 coefficients are not yet implemented.

### Cavity example

```ini
mesh = mesh/cavity.msh
output = results/cavity.plt

time_scheme = steady
density = 1.0
viscosity = 0.01
gradient = least_squares
convection_blending = 1.0

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

Every Gmsh physical boundary patch must have both a `velocity` and a `pressure` entry. Patch names in `case.in` must exactly match the mesh's physical names, including case.

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

## Gmsh meshes

### Supported input

- ASCII Gmsh format 2.2 or 4.1.
- Two-dimensional linear triangles and quadrilaterals.
- Linear line elements for boundary faces.
- Physical names for robust boundary assignment.
- Z coordinates are read but ignored by the two-dimensional geometry.

Curved/high-order elements, binary Gmsh files, three-dimensional elements, and mixed 2D/3D models are not supported. Every geometric boundary edge must have a corresponding Gmsh line element; otherwise mesh construction stops with an error rather than silently omitting a boundary condition.

### Supplied meshes

Regenerate all three meshes with:

```bash
python3 cases/generate_meshes.py
# or
make meshes
```

The generator uses only Python's standard library and writes ASCII Gmsh 2.2 files.

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
| `p` | Cell centres | Static pressure. |
| `mass_imbalance` | Cell centres | Signed integrated cell mass imbalance in the current mass-flux units. |
| `k` | Cell centres | Turbulent kinetic energy, for turbulent cases. |
| `epsilon` | Cell centres | Dissipation, for `k_epsilon`. |
| `omega` | Cell centres | Specific dissipation, for SST. |
| `mu_t` | Cell centres | Turbulent dynamic eddy viscosity. |
| `F1`, `F2` | Cell centres | SST blending functions. |

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
- `converged`

The console log uses `dU` for the maximum cell-velocity update, `dp` for the maximum relaxed pressure update, `dT` for the normalized turbulence update, `continuity` for the maximum absolute cell mass imbalance, and `linear` for the total PETSc iterations used by that outer iteration's momentum, pressure, and turbulence solves. Laminar output retains zero turbulence diagnostics.

## Validation cases

Run the complete workflow:

```bash
make validate
```

This runs the three laminar regressions and all twelve plate/elbow turbulence cases. `cases/validate.py` checks the laminar cases; `cases/validate_turbulence.py --run` executes the turbulence matrix, checks conservation and accuracy, and writes `validation/metrics.csv`, comparison CSVs, and SVG plots. The detailed model-validation report is [`validation/turbulence_validation.md`](validation/turbulence_validation.md).

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

For difficult cases, first reduce `convection_blending`, lower the SIMPLE relaxation factors, or improve the mesh. Tightening linear-solver tolerances does not by itself guarantee nonlinear SIMPLE convergence.

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

At present, `src/main.cpp` directly constructs `SimpleSolver`. A general physics-selection mechanism and inter-solver field exchange still need to be added before this becomes a coupled multiphysics executable.

## Current limitations

- Serial execution only; there is no domain decomposition or distributed PETSc matrix assembly.
- Two-dimensional linear triangles and quadrilaterals only.
- Constant density and dynamic viscosity.
- Two-dimensional RANS only: three-dimensional secondary flow, swirl, and cross-sectional Dean vortices cannot be represented.
- Only standard \(k\)-\(\varepsilon\) and SST-2003; no transition, compressibility, buoyancy, curvature/rotation correction, DES, or Reynolds-stress model.
- No elastic or acoustic field solver yet.
- No fluid–structure or fluid–acoustic coupling yet.
- No convection limiter, so central reconstruction can produce oscillations in strongly convection-dominated or poor-quality-mesh cases.
- No curved/high-order Gmsh element support.
- No binary Gmsh input or binary Tecplot output.
- No restart/checkpoint facility.
- No variable-time-step BDF2 coefficients.
- No automated grid-refinement, temporal-order, or manufactured-solution study yet.
- The supplied validation meshes are deliberately modest in size and should not be treated as grid-independent production results.

## To-do work

### Physics modules and multiphysics coupling

- [ ] Implement the two-dimensional elastic solid solver, including displacement, stress, traction, and fixed-displacement boundaries.
- [ ] Implement the two-dimensional acoustic solver with appropriate pressure/velocity or displacement-potential formulation.
- [ ] Add configuration-driven physics selection instead of directly constructing `SimpleSolver` in `main.cpp`.
- [ ] Design shared field registration and solver scheduling for coupled calculations.
- [ ] Add conservative field/force interpolation between different physics meshes.
- [ ] Implement fluid–structure interface conditions and an explicit or partitioned fluid–structure coupling algorithm.
- [ ] Add acoustic radiation, impedance, and non-reflecting boundary conditions.
- [ ] Investigate strong/monolithic coupling where explicit partitioned coupling is unstable.

### Mesh and geometry

- [ ] Add curved Gmsh elements and high-order geometry representation.
- [ ] Support binary Gmsh files if input size becomes a limitation.
- [ ] Add mesh-quality reporting for aspect ratio, skewness, non-orthogonality, and minimum face delta.
- [ ] Add mesh-refinement and grid-convergence workflows.
- [ ] Consider general polygonal control volumes and hanging-node meshes.
- [ ] Add moving-mesh support for elastic boundaries and fluid–structure interaction.

### Finite-volume numerics

- [ ] Add bounded second-order convection with flux or gradient limiters.
- [ ] Add additional convection schemes such as second-order upwind and QUICK.
- [ ] Improve robustness of non-orthogonal corrections on highly skewed meshes.
- [ ] Add variable-property density and viscosity fields where appropriate.
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

- [ ] Add checkpoint/restart support for steady and transient calculations.
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
- Start with `convection_blending = 0`, then increase it for higher-order accuracy.
- Confirm that at least one pressure reference exists, either explicitly or through the automatic closed-domain reference.
- Inspect whether the PETSc KSP itself is diverging; nonlinear relaxation cannot repair an invalid mesh or an inconsistent linear system.

### The residual is small but the result looks unphysical

Check the actual mesh boundary patches, local and global mass conservation, Reynolds-number definition, and whether the grid is adequate. Algebraic convergence alone does not establish physical accuracy.
