# Neutron transport validation

The cases exercise conservative SN upwind transport, exact axis reflection, multigroup diffusion, fixed-source and eigenvalue iteration, Dirichlet and Marshak boundaries, and triangle/quadrilateral meshes. Eigenvalue references use the matching discrete-grid buckling; continuous references are also listed.

| Case | k_eff | Reference k | k error | Distribution L2 | Balance | Fission production | Absorption | Net leakage |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| attenuation_16 | — | — | — | 2.506e-02 | 6.059e-16 | 0.000000000e+00 | 1.777003524e+00 | -1.777003524e+00 |
| attenuation_64 | — | — | — | 6.506e-03 | 7.791e-16 | 0.000000000e+00 | 1.802617511e+00 | -1.802617511e+00 |
| attenuation_triangles | — | — | — | 7.161e-03 | 1.991e-15 | 0.000000000e+00 | 1.802617511e+00 | -1.802617511e+00 |
| sn_critical | 1.2000000000 | 1.2000000000 | 2.987e-12 | 1.433e-13 | 2.490e-12 | 1.200000000e+00 | 1.000000000e+00 | 0.000000000e+00 |
| diffusion_16 | 1.0640421394 | 1.0640421394 | 7.501e-13 | 1.006e-13 | 5.382e-13 | 2.200000000e+01 | 1.000000000e+00 | 1.967587287e+01 |
| diffusion_32 | 1.0616039328 | 1.0616039328 | 1.013e-12 | 1.097e-13 | 7.936e-13 | 2.200000000e+01 | 1.000000000e+00 | 1.972335955e+01 |
| diffusion_two_group | 1.0061060454 | 1.0061060454 | 1.109e-12 | 1.871e-13 | 9.943e-13 | 1.804982963e+01 | 1.314057176e+00 | 1.662622825e+01 |
| diffusion_marshak | — | — | — | 1.154e-04 | 1.203e-12 | 0.000000000e+00 | 6.788985886e-01 | -6.788985886e-01 |

One-group continuous-eigenvalue refinement order: **2.0019** (expected approximately 2).
Attenuation L2 refinement order: **1.9458**.

The two-group reference solves $A_{11}\phi_1=\chi_1R$ and $A_{22}\phi_2-\Sigma_{s,1\to2}\phi_1=\chi_2R$, with removal $\Sigma_{r,g}=\Sigma_{t,g}-\Sigma_{s,g\to g}$. The balance check independently evaluates $P_f/k=A+L$ for eigenvalue problems and $P_f+Q=A+L$ for fixed-source problems.

## Distribution plots

- [SN attenuation](figures/neutron_attenuation.svg)
- [One-group diffusion eigenmode](figures/neutron_diffusion.svg)
- [Two-group diffusion eigenmode](figures/neutron_two_group.svg)
- [Marshak fixed-source diffusion](figures/neutron_diffusion_marshak.svg)

**All neutron validation checks passed.**
