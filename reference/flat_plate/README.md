# Flat-plate reference data

These files are unmodified downloads from the NASA Langley Turbulence Modeling Resource (TMR), [2D Zero Pressure Gradient Flat Plate Validation Case](https://tmbwg.github.io/turbmodels/flatplate_val.html).

- `cf_karman_schoenherr.dat`: skin-friction coefficient versus momentum-thickness Reynolds number from the Kármán–Schoenherr relation. Original URL: <https://tmbwg.github.io/turbmodels/FlatPlate_validation/cf_K-S.dat>.
- `uplus_coles_retheta_10000.dat`: Coles mean-velocity law with van Driest near-wall damping at \(Re_\theta=10{,}000\). Original URL: <https://tmbwg.github.io/turbmodels/FlatPlate_validation/u+y+.dat>.
- `retheta_variation_typical.dat`: typical \(Re_\theta(x)\) for the TMR case, used only as a setup cross-check and not as solver input. Original URL: <https://tmbwg.github.io/turbmodels/FlatPlate_validation/retheta_variation_typical.dat>.

The validation case has \(Re=5\times10^6\) per unit length and a no-slip plate from \(x=0\) to \(x=2\). The present incompressible setup uses \(U_\infty=1\), \(\rho=1\), and \(\mu=2\times10^{-7}\). The Kármán–Schoenherr and Coles curves are theoretical correlations calibrated to experimental data, not pointwise measurements; their own uncertainty is discussed in the validation report.

The SST model definition and constants were cross-checked against the [TMR SST page](https://tmbwg.github.io/turbmodels/sst.html) and Menter, Kuntz, and Langtry (2003).
