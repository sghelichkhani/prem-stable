# A gravitationally stable PREM for GIA models

This note explains the rationale and methods for producing a gravitationally stable PREM profile for use in
glacial isostatic adjustment (GIA) modelling studies. We provide both the original profile `PREM.csv` and
a preferred modified version `PREM-stable.csv` It has been developed as part of the GIAMIP modelling efforts
and, in addition to density, all other parameters ($g$, $V_P$, $V_S$, $\kappa$, and $\mu$) should also be
read in from the revised profile.

## Why is there an issue?

Self-gravitating GIA models require the viscoelastic properties and density of the solid Earth to be prescribed,
and 1D elastic and density structure have commonly been taken from PREM (Dziewonski & Anderson, 1981), although
several alternative published radial profiles are also used. PREM's elastic moduli and density are obtained from
fitting body wave and free-oscillation data. However, there is no requirement in its construction for the density
profile to be gravitationally stable (and indeed, the mantle convects vigorously over geological timescales).
This fact can therefore result in gravitational instabilities occurring in GIA models that are undesirable 
(Plag & Jüttner, 1995; Vermeersen & Mitrovica, 2000). In particular, when we consider a parcel of rock during
unloading that moves upwards compared to it surroundings, if it's density difference (i.e., buoyancy) with respect
to the new surroundings increases (for example, where PREM prescribes a negative density gradient with depth),
then it will want to rise further, and so on, leading to a runaway overturn. In practice, this means that if we
perturb such an Earth model during GIA loading cycles, it will not return to it's original structure.

### The stability condition

The reference Earth in GIA is considered to be in hydrostatic equilibrium, such that the change in hydrostatic
pressure with radius is,

$$\frac{dP}{dr} = -\rho g,$$

where $r$ is the radius, $P$ the pressure, $\rho$ the density, and $g$ the gravitational acceleration. When a small
parcel of rock moves upwards, the pressure on it decreases causing it to
expand. If the parcel exchanges no heat with its surroundings (i.e., an adiabatic process), the bulk
modulus $\kappa$ controls this expansion, $d\rho/\rho = dP/\kappa$. Referring to the
hydrostatic equilibrium equation, the density of the parcel therefore changes with radius as

$$\left(\frac{d\rho}{dr}\right)_{\mathrm{ad}} = -\frac{\rho^2 g}{\kappa},$$

which is the Adams–Williamson equation (Williamson and Adams, 1923). By comparing this adiabatic density
gradient with the density gradient of PREM at the same depth, we see that when the density of the Earth model
decreases upwards faster than the adiabatic density gradient, the parcel becomes denser than its
new surroundings and will sink back down (i.e., the layer is stable). However, if the density of the
Earth model decreases more slowly (or increases upward) compared to the adiabatic density gradient, the parcel
becomes more buoyant and continues to rise (i.e., the layer is unstable). A metric for is stability can be defined as

$$e = 1 + \frac{\kappa}{\rho^2 g}\ \frac{d\rho}{dr},$$

where $e$ is one minus the Bullen parameter. A layer with $e < 0$ is stable, a layer
with $e = 0$ is adiabatic, and a layer with $e > 0$ is unstable.

In a rigid, elastic layer, shear strength holds the parcel in place, such that even when $e > 0$, the gravitational
instability has no effect on the GIA model. In a layer that flows, however, the parcel will start to move at a rate
controlled by the viscosity and displacement then grows exponentially with time. This phenomenon is a Rayleigh–Taylor
instability, and its growth rate increases as the viscosity decreases. In a Love-number
calculation, it appears as a mode with a positive growth rate. Even in the absence of any loads, it never relaxes such
that the response does not relax to its fluid limit (Plag & Jüttner, 1995; Vermeersen & Mitrovica, 2000). The same
growing modes also occur in numerical GIA models (e.g., Huang et al., 2023).

It is important to note that, for an incompressible GIA model where $\kappa \to \infty$, the adiabatic density gradient is
zero and the condition becomes simple: density in the Earth model must never increase upwards. Incompressible models are therefore
only affected where PREM's density increases upward. Compressible models, however, can have $e > 0$ even in places where PREM's
density reduces upwards, so more care is required to ensure that the profile is gravitationally stable.

### Where is PREM gravitationally unstable?

| Depth (km) | $e$ in PREM | Density excess over adiabatic (kg/m³) |
|---|---|---|
| 0–24.4 (crust) | 1 (constant density in each layer) | ~20 |
| 24.4–220 | 1.13 (density inverted) | ~190 |
| 220–400 | 0.17 to 0.22 | 26 |
| 400–600 | −0.98 to −0.73 (fully stable) | 0 |
| 600–670 | 0.63 | 28 |
| 670–2891 | −0.01 to 0.03 | <9 |

The 24.4–220 km layer is the most significant problem, as it is unstable in both
compressible and incompressible models and includes the asthenosphere, which is often given the lowest viscosities in GIA models
(so instabilities develop quickly). The 400–600 km layer is always stable, while the other layers are unstable in compressible models only.

Since GIAMIP does not prescribe compressibility, groups with compressible models are likely to see more of these instabilities than 
groups with incompressible models and these difference will appear in the model spread. Moreover, since some of the experiments leave the
viscosity open, then any gravitational instabilities will develop at different rates between models. For example, in our Love number code,
the PREM response at low degrees misses its fluid limit by up to 3.4%, while a compressible G-ADOPT finite element model with a low viscosity
mantle developed a flow that continued to accelerate.

## The ANU-GADOPT team's approach to fixing this problem

A self-gravitating GIA model generally requires inputs of density $\rho$, bulk modulus $\kappa$, shear modulus $\mu$, and viscosity. Gravitational
acceleration $g$ can be computed self-consistently from the density. We note that some codes instead read in $V_P$ and $V_S$ and internally compute 
the bulk and shear moduli using $\mu = \rho V_S^2$ and $\kappa = \rho V_P^2 - \tfrac{4}{3}\mu$ (PREM actually provides $V_P$ and $V_S$ rather than
the moduli).

In our modification, we have chosen to leave $\kappa$ and $\mu$ consistent with those in the original PREM profile (i.e., the elastic
stiffness of the model is PREM's) and only change $\rho$ (which therefore results in small associated changes in $V_P$ and $V_S$ of
approximately $-\tfrac{1}{2}\, \Delta\rho/\rho$). The $V_P$ and $V_S$ in `PREM-stable.csv` correspond to these revised wavespeeds. The resulting
profile is, therefore, purely for the purposes of GIA rather than seismic modelling.

Our changes have modified density in the crust and throughout the whole mantle. While a purely elastic lithosphere cannot flow and so carries no 
instabilities, the lithosphere can be very thin in some regions and so the crust and shallowest mantle structure can still be susceptible to them.
In models that give this lid a high but finite viscosity, the lid flows slowly and a long run can show the instability. For this resason,
we therefore include the crust and lithospheric mantle in our modification (our revised profile changes crustal density by about 10 kg/m³ and
maintains the large density jumps at 15~km and 24.4 km). Deeper than 670 km, PREM is already close to adiabatic
($e < 0.03$), such that the required changes in density in the lower mantle are only a few kg/m³. Instabilities at these depths in the original PREM
are generally small and take a long time to occur for typical lower mantle viscosities (Vermeersen & Mitrovica, 2000). Nevertheless, we have decided
to still fix them as long-duration GIA models with a low viscosity lower mantle can still produce them. The properties of the core remain the same 
as those in the original PREM and, following our corrections, all layers of the crust and mantle are gravitationally stable.

### Optimisation framework

Our aim is to keep density as close to the original PREM as possible and we infer density at each node between the surface and the core-mantle boundary
as the solution of the least-squares problem

$$\min \sum_i w_i \left(\rho_i - \rho_i^{\mathrm{PREM}}\right)^2,$$

where $i$ represents each node and $w_i$ is the mass of the shell at that node. Three conditions are applied:

1. Every interval between two nodes is either stable or adiabatic:

   $$\frac{\rho_{i+1} - \rho_i}{r_{i+1} - r_i} + \frac{\bar\rho^2\,\bar g}{\bar\kappa} \le 0,$$

   where the bar is the mean value across the two nodes. This is the discrete form of $e \le 0$.
2. At every discontinuity, density in the layer below cannot be lower than that above.
3. The total mass of the crust and the mantle is equal to those values in PREM. Therefore, Earth's total mass and gravitational acceleration at the
   surface and within the core remain consistent. It is, however, possible for mass to move between the upper and the lower mantle.

Condition 1 contains $g$, which itself depends on density, so we solve the problem iteratively
in a loop. $g$ is computed based upon the current density profile, the least-squares problem
with these linear conditions results in a new density profile, from which a $g$ is revised and the loop repeats. Five
iterations is sufficient to bring the density profile within 2e-5 kg/m³ of its final values.

## The result

![Density of PREM and PREM-stable](prem_density.png)

*Density of PREM as originally published (blue), in our revised PREM-stable (orange), and in VM5i (grey; Table 2 of Spada and Melini, 2019). a. The crust
and whole mantle, 0 to 2891 km. b. Zoom of the upper 700 km, where the modifications are largest. VM5i averages PREM into 11 layers of constant density. 
In a compressible model, it is therefore not stable.*

| Depth (km) | Range of density changes (kg/m³) |
|---|---|
| 0–24.4 (crust) | −9 to +10 |
| 24.4–220 | −109 at 24.4 km, +81 at 220 km |
| 220–400 | +4 to +32 |
| 400–600 | -12 to +0.3 |
| 600–670 | −12 to +16 |
| 670–771 | −7 to −5 |
| 771–2741 | −5 to +3 |
| 2741–2891 (D'') | +2 to +2.5 |

The largest $e$ remaining anywhere in the crust and mantle is 2e-9.
The density increase across the 220 km discontinuity becomes zero, the jump at 400 km reduces from
+180 to +149 kg/m³, and the jump at 670 km reduces from +389 to +366 kg/m³. The lithospheric mantle
layer between 24.4 and 80 km becomes 1.7–3.2% lighter. $g$ changes by a maximum of 0.0055 m/s², which occurs at about 140 km
depth. $V_P$ and $V_S$ change by between −1.2 and +1.7%.

## How to use the files

`PREM.csv` is almost exactly the same as the original published version except that the 2 km thick ocean layer has been replaced with upper crust.
`PREM-stable.csv` is our revised, gravitationally stable model. The two files contain the same 339 depth nodes, so can be compared row by row.

| Column | Quantity | Unit |
|---|---|---|
| `radius_m` | radius | m |
| `depth_km` | depth | km |
| `region` | index of the PREM region in `prem1981.py` | |
| `rho_kg_m3` | density | kg/m³ |
| `vp_m_s` | P-wave speed | m/s |
| `vs_m_s` | S-wave speed | m/s |
| `kappa_Pa` | bulk modulus | Pa |
| `mu_Pa` | shear modulus | Pa |
| `g_m_s2` | gravitational acceleration | m/s² |

- Rows go from Earth's centre to its surface.
- Every boundary between original PREM regions is represented with two rows with the same radius. The
  first row is the value below the boundary and the second row is above.
  Some of these boundaries have no density jump in our modified profile, resulting in the two rows being
  equal.
- Nodes are spaced at most 10 km apart between the surface and core-mantle
  boundary, and at most 100 km apart within the core.
- Between the surface and core-mantle boundary, the density should be interpolated
  linearly between nodes (the gravitational acceleration column is computed assuming this density interpolation). If your code computes
  $g$ from density using the trapezoid rule, the result will differ from our column by at most 3e-5 of its value.
- Use pur $\kappa$ and $\mu$ (consistent with original PREM). If you read $V_P$ and $V_S$, read them
  from the our file too (different to original PREM). Do not combine the density of
  `PREM-stable.csv` with the $V_P$ and $V_S$ of another PREM profile as that
  combination gives elastic moduli that are not consistent with the original PREM.
- Our profiles do not include viscosity and attenuation.


## References

Dziewonski, A. M. and Anderson, D. L. (1981). Preliminary reference Earth
model. *Physics of the Earth and Planetary Interiors* 25, 297–356.

Huang, P., Steffen, R., Steffen, H., Klemann, V., Wu, P., van der Wal, W.,
Martinec, Z. and Tanaka, Y. (2023). A commercial finite element approach to
modelling Glacial Isostatic Adjustment on spherical self-gravitating
compressible earth models. *Geophysical Journal International* 235, 2231–2256.

Plag, H.-P. and Jüttner, H.-U. (1995). Rayleigh–Taylor instabilities of a
self-gravitating Earth. *Journal of Geodynamics* 20, 267–288.

Spada, G. and Melini, D. (2019). SELEN4 (SELEN version 4.0): a Fortran program for
solving the gravitationally and topographically self-consistent sea-level equation
in glacial isostatic adjustment modeling. Geoscientific Model Development 12, 5055–5075.

Vermeersen, L. L. A. and Mitrovica, J. X. (2000). Gravitational stability of
spherical self-gravitating relaxation models. *Geophysical Journal
International* 142, 351–360.

Williamson, E. D. and Adams, L. H. (1923). Density distribution in the Earth.
*Journal of the Washington Academy of Sciences* 13, 413–428.
