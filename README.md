# A gravitationally stable PREM for GIA models

This repository provides two radial Earth models for use in GIA simulations.
`PREM.csv` is the orignal version of PREM and `PREM-stable.csv` is a version of PREM
that has been modified so that it is gravitationally stable in every layer of the crust and the mantle.
`METHOD.md` explains why PREM needs this change, how the change has been made, and
how to use the files.

## Source of PREM

PREM is sourced from Table I of Dziewonski & Anderson (1981). `prem1981.py` holds
the polynomial coefficients taken directly from that paper. Two modifications have been made:

- Between 24.4 and 220 km depth, PREM is transversely isotropic. Our files use
  the isotropic wavespeeds from the footnote of Table I,
  V_P = 4.1875 + 3.9382x and V_S = 2.1519 + 2.3481x, with x = r / 6371 km.
- The 3 km ocean layer is replaced by a continuation of upper crust (2600 kg/m³, 5.8 km/s,
  3.2 km/s) to the surface. A GIA model applies the ocean as a surface load on top of
  the topography/bathymetry, so the Earth model includes no such fluid layer at the top.

The values are computed for the reference period of 1 s. The files do not contain
attenuation or viscosity.

## Files

| File | Content |
|---|---|
| `PREM.csv` | PREM from the Table I polynomials |
| `PREM-stable.csv` | The stable model, on the same nodes as `PREM.csv` |
| `METHOD.md` | The reasons, method, result, and column formats |
| `VM5i_Spada_Melini_2019.csv` | The 11-layer model VM5i, for comparison only |
| `prem_density.png` | Figure showing density of PREM, PREM-stable and VM5i |
| `prem1981.py` | The Table I coefficients of D&A81 and a function that evaluates them |
| `build_prem_stable.py` | Calculates the stable density and writes both CSV files |
| `plot_density.py` | Plots `prem_density.png` |

## Build the files

The scripts need Python 3, NumPy, SciPy and Matplotlib. Run them from this
directory:

```
python build_prem_stable.py
python plot_density.py
```

`build_prem_stable.py` should take a few seconds. It prints the density change, 
stability measure, total mass, surface gravitational acceleration, and density jumps across discontinuities.
It writes `PREM.csv` and `PREM-stable.csv`. `plot_density.py` writes
`prem_density.png`.

## Comparison model

`VM5i_Spada_Melini_2019.csv` is the Earth model VM5i from Table 2 of Spada & Melini
(2019). It consists of 11 layers of constant density, rigidity, and viscosity,
and it is the test model of the sea-level code SELEN4. The values were typed
from Table 2 and checked against the published PDF. The paper is published
under the Creative Commons Attribution 4.0 License.

The rows go from the surface down. Radius and depth are in km and density is
in kg/m³. Rigidity is in units of 1e11 Pa and viscosity in units of 1e21 Pa s.
The viscosity of the elastic lithosphere is written `inf`. The depth columns
are 6371 km minus the radii.

`prem_density.png` shows VM5i for comparison. No file in this repository is
built from it.

## References

Dziewonski, A. M. and Anderson, D. L. (1981). Preliminary reference Earth
model. *Physics of the Earth and Planetary Interiors* 25, 297–356.

Spada, G. and Melini, D. (2019). SELEN4 (SELEN version 4.0): a Fortran program
for solving the gravitationally and topographically self-consistent sea-level
equation in glacial isostatic adjustment modeling. *Geoscientific Model
Development* 12, 5055–5075.
