# 90-degree elbow reference data

## Primary experiment

M. M. Enayet, M. M. Gibson, A. M. K. P. Taylor, and M. Yianneskis, “Laser-Doppler measurements of laminar and turbulent flow in a pipe bend,” *International Journal of Heat and Fluid Flow* 3, 213–219 (1982).

The turbulent case has diameter-based bulk Reynolds number 43,000, bend mean-radius/diameter ratio 2.8, and bulk speed approximately 0.92 m/s in the physical experiment. The Ansys VMFL030 benchmark uses a dynamically similar scaled configuration.

## Digitized 75-degree profile

`velocity_75deg_experimental.csv` was manually digitized from the experimental point symbols in the public [Ansys VMFL030](https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/fbu_vm/Hlp_VMFL030.html) figure “Comparison of Velocity Magnitude (m/s) at 75° Along the Bend.” The original image is retained as `ansys_vmfl030_75deg.jpg`; the domain schematic is `ansys_vmfl030_domain.svg`.

- `y_over_d` is the figure’s horizontal coordinate, from 0 at one wall to 1 at the opposite wall.
- The velocity columns retain the plotted values. Normalization by unit bulk speed is used for the nondimensional solver comparison.
- Estimated digitization uncertainty is ±0.005 in `y_over_d` and ±0.02 in velocity. The wall values are exact boundary values shown on the plot axes, not measured finite-speed wall data.
- For the 2D calculation, `y_over_d=0` is mapped to the inner-radius wall and `y_over_d=1` to the outer-radius wall. This is an explicit planar-mapping assumption; the public figure alone does not reconstruct the full 3D measurement-plane orientation. No reversal or fitting is performed to improve agreement.

The comparison is not a full three-dimensional validation: the experiment contains circumferential secondary flow that a planar solver cannot represent.

## Independent cross-check

C. E. Alvarado-Rodríguez, L. Di G. Sigalotti, J. Klapp, C. Fierro, F. Aragón, and A. R. Uribe-Ramírez, “SPH simulations of turbulent flow in curved pipes with different geometries. A comparison with experiments,” *Journal of Fluids Engineering* 143(9), 091503 (2021), [arXiv:2107.13729](https://arxiv.org/pdf/2107.13729), CC BY 4.0.

`alvarado_fig5_maps.png` contains independent numerical velocity maps at 30°, 60°, 75°, and downstream stations. `alvarado_fig6_upstream.png` contains the authors’ comparison with Enayet’s experimental upstream velocity profiles. These images are cross-checks only; the CSV error metric uses the experimental Ansys/Enayet points above.
