# A gravitationally stable PREM for GIA models

This directory builds a version of PREM whose density is stable in the flowing
mantle. The target is the GIA Model Intercomparison Project (GIAMIP), which
prescribes PREM as the elastic structure but gives no PREM file.

## The problem

In PREM as published, density increases upward between the Moho and 220 km
depth. PREM is also steeper than the adiabatic gradient between 220 and 400 km
and between 600 and 670 km. In a compressible viscoelastic model, a layer like
that is gravitationally unstable where it flows. It carries Rayleigh-Taylor
modes that grow with time. An incompressible model is unstable only where the
density increases upward, which in PREM is the layer between 24.4 and 220 km.
`METHOD.md` explains the condition, the correction and the files.

The evidence:

- Vermeersen and Mitrovica (2000), `sources/vermeersen-mitrovica-2000.pdf`.
- Huang et al. (2023), `sources/huang-et-al-2023.pdf`, Section 5.
- lovejx, `~/Workplace/lovejx/NOTES/STRATIFICATION.md`. It measures the growing
  mode of PREM in the Love numbers and explains the 1 to 4 per cent low-degree
  error of the legacy Fortran.
- G-ADOPT runs by Will Scott, February 2026. A compressible model with a
  layered density and a weak mantle developed a growing flow. The same model
  with a density that increases smoothly with depth relaxed to equilibrium.

## The files

| File | Purpose |
|---|---|
| `prem1981.py` | PREM from the polynomials of Dziewonski and Anderson (1981), Table I |
| `check_prem1981.py` | Checks the transcription against ObsPy's `prem.nd` and the lovejx deck |
| `build_prem_stable.py` | Fits the stable density and writes `PREM.csv` and `PREM-stable.csv` |
| `PREM.csv` | PREM on the output node table, centre outward, SI units |
| `PREM-stable.csv` | The stable PREM on the same nodes |
| `plot_fit.py` | Plots e and the density adjustment of the two CSV files, `prem_stable_fit.png` |
| `plot_prem_stable.py` | Plots the older lovejx PREM-stable against `prem.nd`, `prem_vs_prem_stable.png` |
| `METHOD.md` | Why and how PREM is changed, and how to use the files. For the GIAMIP groups |
| `DECISIONS.md` | Every decision, with the date, the reason and the numbers |
| `sources/` | The scan of Table I and the two papers. Not tracked by git |
| `rejected/` | Two rules we tested and did not keep |
| `email/` | The message to Holly Han about PREM |

## How to rebuild

Use the project Python, `~/Workplace/python3.12/bin/python3.12`. Run the
commands from this directory.

1. Run `check_prem1981.py`. Every density misfit against the lovejx deck must
   be below 1e-5 g/cm^3.
2. Run `build_prem_stable.py`. It takes about 2 seconds. It prints the density
   change, the stability measure, the mass and the jumps.
3. Run `plot_fit.py`.

`check_prem1981.py` needs ObsPy and a checkout of lovejx at
`~/Workplace/lovejx`. The build needs only NumPy and SciPy.

## Status

PREM-stable.csv is stable from the surface to the core-mantle boundary. The
core is PREM. The columns are radius_m, depth_km, region, rho_kg_m3, vp_m_s, vs_m_s,
kappa_Pa, mu_Pa and g_m_s2. Every region boundary has two rows at one radius,
the lower side first. `DECISIONS.md` lists the open items.
