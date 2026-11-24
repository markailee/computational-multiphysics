# Circular-arc bump references

## Geometry and flow conditions

The implemented cases use a channel of length `3` and height `1`, with a circular-arc bump centred on the lower wall between `x = 1` and `x = 2`. For chord length `c = 1` and height `h`, the circular radius is

`R = (h^2 + (c/2)^2)/(2h)`, and `y_lower(x) = h - R + sqrt(R^2 - (x-1.5)^2)` on the bump.

- Subsonic: `h/c = 0.10`, inlet Mach `0.5`.
- Transonic: `h/c = 0.10`, inlet Mach `0.675`.
- Supersonic: `h/c = 0.04`, inlet Mach `1.4`.

Useful descriptions and comparison results:

1. [SU2, Inviscid Bump in a Channel](https://su2code.github.io/tutorials/Inviscid_Bump/). Geometry, slip-wall treatment, subsonic inlet/outlet conditions, and the transonic total-condition settings. The tutorial's transonic values are `T0 = 316.224 K`, `p0 = 140513.23 Pa`, and outlet static pressure `101300 Pa`.
2. [Ferziger, Peric, and Street, Computational Methods for Fluid Dynamics (2020)](https://elmoukrie.com/wp-content/uploads/2022/05/joel-h.-ferziger-milovan-peric-robert-l.-street-computational-methods-for-fluid-dynamics-springer-international-publishing-2020.pdf). Chapter 11 describes subsonic, transonic, and supersonic circular-bump flows, including the nearly symmetric subsonic solution.
3. [Kasturi Rangan and Ghosh (2022), A face-based immersed boundary method for compressible flows using a uniform interpolation stencil](https://www.frontiersin.org/journals/mechanical-engineering/articles/10.3389/fmech.2022.903492/full). Section 4.1 and Table 2 give transonic exit Mach numbers from Favini et al. (1996): **top 0.670, bottom 0.615**. The paper also compares shock position and surface Mach profiles.
4. [Djavareshkian and Farnak (2008), Application of UNO Scheme for Steady and Transient Compressible Flow in Pressure-Based Method](https://scialert.net/fulltext/?doi=jas.2008.3361.3370). Transonic `M = 0.675`, 10% bump and supersonic `M = 1.4`, 4% bump cases, including wall Mach and pressure distributions. `uno_supersonic.gif` is Figure 4 from this page.
5. [Ma, Emad, and Chen (2008), A Local Meshless Method for Solving Compressible Euler Equations](https://scialert.net/fulltext/?doi=srj.2008.1.16). Pressure contours and wall pressure coefficients for the same transonic/supersonic geometries. `lmm_wall_pressure.gif` is Figure 8 from this page.
6. [A parallel pressure implicit splitting of operators algorithm applied to flows at all speeds](https://eprints.soton.ac.uk/21853/1/bres_01.pdf). Documents the symmetric `M = 0.5` subsonic solution and grid sensitivity of the wall Mach profiles.

The SU2 configuration used for the transonic case is available in the [tutorial repository](https://github.com/su2code/Tutorials/blob/master/compressible_flow/Inviscid_Bump/inv_channel.cfg).

## Digitized data

`transonic_exit_mach.csv` contains the direct values from Table 2 of Kasturi Rangan and Ghosh (2022), citing Favini et al. (1996).

`supersonic_mach_approx.csv` contains approximate points read from Figure 4(a) of Djavareshkian and Farnak (2008). These points are **plot digitization, not tabulated author data**, and should be treated as having approximately `0.05` Mach-unit reading uncertainty in addition to uncertainty in identifying the plotted wall. They are used for qualitative/profile comparison, not as a high-precision acceptance target.
