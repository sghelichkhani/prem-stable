# A gravitationally stable PREM for GIA models

This repository gives two Earth models for the GIA Model Intercomparison
Project (GIAMIP). `PREM.csv` is PREM. `PREM-stable.csv` is PREM with a density
that is gravitationally stable in every layer of the crust and the mantle.
`METHOD.md` explains why PREM needs this change, how the change is made, and
how to use the files.

## Source of PREM

The source is Table I of Dziewonski and Anderson (1981). `prem1981.py` holds
the polynomial coefficients of that table, typed from the printed paper. No
other PREM table is used. Two choices differ from the table as printed:

- Between 24.4 and 220 km depth, PREM is transversely isotropic. The files use
  the isotropic speeds from the footnote of Table I,
  V_P = 4.1875 + 3.9382x and V_S = 2.1519 + 2.3481x, with x = r / 6371 km.
- The 3 km ocean layer is replaced by the upper crust (2600 kg/m³, 5.8 km/s,
  3.2 km/s) up to the surface. A GIA model applies the ocean as a load, so the
  Earth model has no fluid layer at the top.

The values are for the reference period of 1 s. The files contain no
attenuation and no viscosity.

## Files

| File | Content |
|---|---|
| `PREM.csv` | PREM from the Table I polynomials |
| `PREM-stable.csv` | The stable model, on the same nodes as `PREM.csv` |
| `METHOD.md` | The reasons, the method, the result and the column format |
| `VM5i_Spada_Melini_2019.csv` | The 11-layer model VM5i, for comparison only |
| `prem_density.png` | The density of PREM, PREM-stable and VM5i |
| `prem1981.py` | The Table I coefficients and a function that evaluates them |
| `build_prem_stable.py` | Calculates the stable density and writes both CSV files |
| `plot_density.py` | Plots `prem_density.png` |

## Design decisions

- Only the density changes. The bulk and shear moduli stay as in PREM, so the
  elastic stiffness of the model is PREM's. V_P and V_S are calculated again
  from the moduli and the new density. The table is for GIA and not for
  seismology.
- The change covers the crust and the whole mantle, from the surface to the
  core-mantle boundary. A table that is stable from the surface down is stable
  for every lithosphere thickness and every lithosphere viscosity. The core
  stays PREM.
- The new density is the density closest to PREM, in a mass-weighted
  least-squares sense, that has no layer steeper than the adiabatic gradient
  and no density jump that increases upward. The mass of the crust and the
  mantle stays PREM's, so the total mass and the surface gravity do not change.
- The nodes are at most 10 km apart from the surface to the core-mantle
  boundary and at most 100 km apart in the core. Each boundary between PREM
  regions has two rows at the same radius. The fit is done on these nodes, so
  the stability condition holds for the file as it is written.

## Comparison model

`VM5i_Spada_Melini_2019.csv` is the Earth model VM5i of Spada and Melini
(2019), Table 2. It has 11 layers of constant density, rigidity and viscosity,
and it is the test model of the sea-level code SELEN4. The values were typed
from Table 2 and checked against the published PDF. The paper is published
under the Creative Commons Attribution 4.0 License.

The rows go from the surface down. Radius and depth are in km and density is
in kg/m³. Rigidity is in units of 1e11 Pa and viscosity in units of 1e21 Pa s.
The viscosity of the elastic lithosphere is written `inf`. The depth columns
are 6371 km minus the radii.

`prem_density.png` shows VM5i for comparison. No file in this repository is
built from it.

## Build the files

The scripts need Python 3, NumPy, SciPy and Matplotlib. Run them from this
directory:

```
python build_prem_stable.py
python plot_density.py
```

`build_prem_stable.py` takes a few seconds. It prints the density change, the
stability measure, the total mass, the surface gravity and the density jumps,
and it writes `PREM.csv` and `PREM-stable.csv`. `plot_density.py` writes
`prem_density.png`.

## References

Dziewonski, A. M. and Anderson, D. L. (1981). Preliminary reference Earth
model. *Physics of the Earth and Planetary Interiors* 25, 297–356.

Spada, G. and Melini, D. (2019). SELEN4 (SELEN version 4.0): a Fortran program
for solving the gravitationally and topographically self-consistent sea-level
equation in glacial isostatic adjustment modeling. *Geoscientific Model
Development* 12, 5055–5075. https://doi.org/10.5194/gmd-12-5055-2019
