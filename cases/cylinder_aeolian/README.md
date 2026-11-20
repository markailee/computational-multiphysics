# Cylinder vortex shedding and aeolian-tone propagation

This case couples the existing incompressible flow solver to the nodal acoustic solver for flow-induced sound from a circular cylinder.

## Benchmark

The reference is O. Inoue and N. Hatakeyama, [Sound generation by a two-dimensional circular cylinder in a uniform flow](https://metodologiauntref.wordpress.com/wp-content/uploads/2015/09/inoue2002.pdf), *Journal of Fluid Mechanics* 471, 285–314, 2002, DOI `10.1017/S0022112002002124`.

The selected parameters are:

- diameter `D = 1`, density `rho = 1`, and freestream `U = 1`;
- `Re = U D / nu = 150`, so `nu = 1/150`;
- acoustic sound speed `c = 5`, giving `M = U/c = 0.2`;
- published `St = f D/U = 0.183`, mean drag coefficient about `1.32`, and lift amplitude about `0.52`;
- a dominant lift dipole at the shedding frequency and a weaker drag dipole at twice that frequency.

The supplied `cases/cylinder/` configuration is steady at `Re = 40` and does not shed vortices. This case instead runs a transient BDF2 calculation with a small initial transverse perturbation.

## Workflow

1. `flow.in` solves the near field on the existing body-fitted cylinder mesh and writes `results/flow_forces.csv` and vorticity snapshots.
2. `cases/run_cylinder_aeolian.py` detects the shedding period from positive-going lift zero crossings. It takes the last three periods, subtracts the mean force, and writes `source/flow_dipole.csv` starting at acoustic time zero.
3. `acoustic.in` reads and linearly interpolates that force history on a mixed acoustic mesh covering `[-20,20] x [-15,15]`. The source duration determines the acoustic end time automatically.
4. The driver measures force coefficients, acoustic probe amplitudes, and lift-dipole directivity, and writes metrics and SVG figures under `validation/`.

Run the complete workflow in WSL:

```bash
make -j2
python3 cases/run_cylinder_aeolian.py --run
# or
make cylinder-aeolian
```

To reuse the existing flow result and rerun only the acoustic propagation:

```bash
python3 cases/run_cylinder_aeolian.py --acoustic-only
```

To regenerate the source and plots without rerunning either solver:

```bash
python3 cases/run_cylinder_aeolian.py
```

## Curle-type source

The flow force output uses the force **of the fluid on the cylinder**, per unit span:

\[
C_D=\frac{F_x}{\tfrac12\rho U^2D},\qquad
C_L=\frac{F_y}{\tfrac12\rho U^2D}.
\]

The force of the cylinder on the fluid is its negative. With a normalized compact Gaussian `g`, the acoustic right-hand side is therefore

\[
-\nabla\cdot\left(-\frac{\mathbf Fg}{\rho}\right)
=\nabla\cdot\left(\frac{\mathbf Fg}{\rho}\right).
\]

The acoustic solver evaluates this as a dipole using the Gaussian gradient. It does not differentiate the force explicitly; retardation through the wave equation produces the radiating time variation. A one-time-unit cosine ramp suppresses the artificial startup from abruptly switching on the source.

The acoustic medium retains `U = (1,0)` and `c = 5`, so the convected wave operator includes the low-Mach Doppler effect. All exterior acoustic patches use the available first-order absorbing condition.

## Outputs

- `results/flow_forces.csv`: total, pressure, and viscous forces plus `Cd` and `Cl`.
- `results/flow_*.plt`: near-field pressure, velocity, vorticity, and mass imbalance.
- `source/flow_dipole.csv`: mean-subtracted, time-shifted force input to the acoustic solver.
- `results/acoustic_*.plt`: propagating acoustic pressure and pressure rates.
- `results/acoustic_probes.csv`: histories at axial and transverse points.
- `validation/cylinder_aeolian_metrics.csv`: measured benchmark and directivity metrics.
- `validation/figures/cylinder_aeolian_*.svg`: wake, force, pressure, and directivity figures.

## Limitations

This is a one-way hybrid Curle-type model, not a compressible DNS. The coarse near-field mesh, bounded flow domain, initial transient, compact Gaussian regularization, low-order acoustic mesh, and first-order absorbing boundaries affect accuracy. The acoustic domain is intended to demonstrate propagation and directivity, not to reproduce the reference paper's very large far-field domain. Agreement with the published force coefficients is checked quantitatively, but small differences should not be interpreted as a grid-independent DNS certification.
