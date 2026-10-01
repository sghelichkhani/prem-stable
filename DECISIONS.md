# Decisions

Each entry gives the date, the decision, the reason and the numbers behind it.
Sia Ghelichkhani made the decisions. The numbers come from the scripts in this
directory unless an entry names another source.

## 2026-09-29: GIAMIP gives no PREM file

GIAMIP protocol v2.4 names PREM as the elastic structure of Exp01, Exp02 and
Exp05 to Exp07. It cites Dziewonski and Anderson (1981) only. The shared input
folder holds ice histories and topographies only. The model report form has
"PREM" filled in and no version. Each group therefore uses its own PREM.

The differences between PREM versions are small, except for density. PREM's
density is unstable in the flowing upper mantle, and the protocol leaves
compressibility open. A compressible group sees the instability and an
incompressible group does not. The difference appears in the model spread.

## 2026-09-30: Build from the published polynomials

**Decision.** Build the stable PREM from Table I of Dziewonski and Anderson
(1981), not from a table that someone else derived.

**Reason.** A derived table can carry changes that nobody recorded. The check
below found one: the wave speeds of the lovejx PREM deck between 24.4 and
220 km match no form in the paper. The paper is the one source that everyone
can check.

**Check.** `check_prem1981.py` compares the transcription with two tables that
were made independently of it.

- ObsPy's `prem.nd`. The misfit is at most 1.3e-4 in any quantity. In the LID
  and the LVZ, `prem.nd` uses the Voigt average of the transversely isotropic
  coefficients. Our transcription of those coefficients reproduces it to 9e-5.
- The lovejx legacy deck. The density misfit is at most 6.5e-6 g/cm^3 in every
  region. This is fine enough to find a wrong last digit in any density
  coefficient.

## 2026-09-30: Isotropic speeds from the footnote of Table I

**Decision.** Between 24.4 and 220 km, use the isotropic speeds of the Table I
footnote, V_P = 4.1875 + 3.9382x and V_S = 2.1519 + 2.3481x.

**Reason.** The footnote is the isotropic form that the paper itself gives.
The Voigt average that `prem.nd` uses differs from it by about 1e-3 km/s.

## 2026-09-30: Replace the ocean with crust

**Decision.** Replace the 3 km ocean layer by the upper crust (2.6 g/cm^3,
5.8 and 3.2 km/s) up to the surface.

**Reason.** A GIA model loads the ocean explicitly, so a fluid layer in the
Earth model is not wanted. `prem.nd` makes the same change.

## 2026-09-30: Adjust the density between 60 and 670 km

Replaced on 2026-10-01 by the entry "Extend the adjustment to the surface".

**Decision.** Change the density only between 60 and 670 km. Keep PREM exactly
above 60 km and below 670 km. Include the 600 to 670 km layer.

**Reason.** 60 km is the base of the VM5a lid, which is the thinnest lid of the
GIAMIP viscosity models. Material above 60 km does not flow in any GIAMIP
experiment, so its stratification carries no unstable mode. The stability
measure of PREM, from the polynomials, is:

| Depth (km) | e = 1 - Bullen parameter | Density excess over the adiabat |
|---|---|---|
| 24.4 to 220 | 1.13 | about 190 kg/m^3 |
| 220 to 400 | 0.17 to 0.22 | 26 kg/m^3 |
| 400 to 600 | -0.98 to -0.73 (stable) | 0 |
| 600 to 670 | 0.63 | 28 kg/m^3 |
| 670 to 2000 | 0.005 to 0.029 | 9 kg/m^3 |
| 2000 to 2891 | -0.011 to 0.008 | under 1 kg/m^3 |

A value of e = 0 is the adiabatic gradient. A positive value is unstable. The
lovejx notes measure the unstable mode of the 600 to 670 km layer as small. The
layer is included because a table for other groups must not contain a known
unstable layer.

## 2026-09-30: The smallest change that makes PREM stable

The rule in this entry still applies. The range is now 0 to 670 km, see the
entry of 2026-10-01. The numbers below are for the range 60 to 670 km.

**Decision.** Between 60 and 670 km, the density is the one closest to PREM in
the mass-weighted least-squares sense, subject to three conditions:

1. In every interval, the density gradient is at most the adiabatic gradient.
   Bulk modulus is PREM's. Gravity comes from the new density.
2. At 220, 400 and 670 km, the density below the discontinuity is at least the
   density above it.
3. The mass between 60 and 670 km is PREM's.

Bulk and shear moduli stay PREM's everywhere. The jump at 60 km is free,
because the lid above it does not flow. This is the rule of lovejx PREM-stable,
applied to the published polynomials with the lid at 60 km.

**Reason.** Three simpler rules fail.

| Rule | Density change | Failure |
|---|---|---|
| Each unstable layer keeps its own PREM mass | -76 to +80 kg/m^3 | Inverted jumps between flowing layers: 18 kg/m^3 at 220 km, 14 kg/m^3 at 600 km |
| 60 to 400 km keeps its mass, PREM jump at 220 km kept | -123 to +60 kg/m^3 | Large changes |
| 60 to 670 km keeps its mass, PREM jumps at 220 and 400 km kept | -146 to +61 kg/m^3 | Large changes, whole transition zone moved by +32 kg/m^3 |

The first rule was a test in an earlier version of `build_prem_stable.py`. The
other two are in `rejected/`. The size of the 60 to 220 km layer is the cause.
PREM's density at 220 km is about 150 kg/m^3 below the value that an adiabat
from 60 km gives. A rule that keeps the mass must spread that difference up and
down.

**Result.** From `build_prem_stable.py`:

| Depth (km) | Density change (kg/m^3) |
|---|---|
| 60 to 220 | -85 at 60 km, +71 at 220 km |
| 220 to 400 | -5.5 to +22 |
| 400 to 600 | 0, except near 600 km |
| 600 to 670 | -12 to +16 |

- The largest e between 60 and 670 km is -9e-14. PREM's is +1.13.
- The total mass changes by 1e-15 of itself. Surface gravity does not change.

**Consequences.**

- The density jump at 220 km becomes zero. PREM has +76 kg/m^3. The speed
  jumps at 220 km stay, because the moduli are PREM's.
- At 60 km there is a density step of 85 kg/m^3, with the lid denser than the
  mantle below it. It carries no unstable mode because the lid does not flow.
  The table is stable only for lids of 60 km or thicker.
- With fixed moduli, V_P and V_S change by about half the relative density
  change, up to about 1.3 per cent. The table is for GIA, not for seismology.

## 2026-10-01: Extend the adjustment to the surface

**Decision.** Apply the rule of 2026-09-30 between the surface and 670 km,
crust included. Conserve the mass of 0 to 670 km. Keep the lower mantle as
PREM. Keep the bulk and shear moduli as PREM.

The lower limit of 670 km in this entry was replaced later on 2026-10-01 by
the entry "Extend the adjustment to the whole mantle". The numbers below are
for the range 0 to 670 km.

**Reason.** The 60 km table is stable only for models whose lid is at least
60 km thick and does not flow. That is not true for every model that will use
the table. In a 3D model the lithosphere is thinner than 60 km in some places,
and there the 24.4 to 60 km layer flows with e = 1.13. G-ADOPT gives the lid a
high but finite viscosity, so the lid flows slowly in a long run. A table that
is stable from the surface down is stable for every lid. The extension also
removes the 85 kg/m^3 step at 60 km, where a denser lid sat on a lighter
mantle.

**Cost.** lovejx measured the cost of a lighter lithosphere on its 171-node
deck (`~/Workplace/lovejx/NOTES/STRATIFICATION.md`, Section 2 and the note of
2026-09-05). A lithosphere 2.3 per cent lighter gave half of the 0.5 to 0.9
per cent change of the degree-32 response over 120 kyr, through the elastic k.
lovejx then kept PREM's lithosphere, because for a rigid lid the change has no
benefit. That argument does not hold for a table that other groups use with
their own lids. The crust is two layers of constant density. That is
super-adiabatic (e = 1), and the fit gives each layer a small gradient.

**Node table.** Nodes are evenly spaced in each region of Table I, at most
10 km apart from the surface to the CMB and at most 100 km apart in the core.
Every region boundary has two rows at one radius. The table has 339 nodes. The
fit is done on this table, so the discrete stability measure of lovejx passes
on the file as it is written. Inside 0 to 670 km the density is linear between
nodes. Mass and gravity are integrated exactly for that density, and for the
PREM polynomials elsewhere.

**Result.** From `build_prem_stable.py`:

| Depth (km) | Density change (kg/m^3) |
|---|---|
| 0 to 15 | -9.1 at the surface, +10.1 at 15 km |
| 15 to 24.4 | -4.7 to +5.7 |
| 24.4 to 220 | -109 at 24.4 km, +81 at 220 km |
| 220 to 400 | +4.6 to +32 |
| 400 to 600 | +0.5, except -5 to -12 between 590 and 600 km |
| 600 to 670 | -12 to +16 |

| Discontinuity (km) | Density jump, below minus above (kg/m^3) | PREM |
|---|---|---|
| 15 | 285 | 300 |
| 24.4 | 366 | 481 |
| 220 | 0 | 76 |
| 400 | 149 | 180 |
| 670 | 372 | 389 |

- The largest e between 0 and 670 km is 5e-10. PREM's is +1.13.
- The total mass changes by 1e-15 of itself. Surface gravity does not change.
  Gravity changes by at most 5.4e-3 m/s^2, at about 140 km.
- V_P and V_S change by -1.2 to +1.7 per cent.
- `check_stratification` of lovejx, with every solid node marked as flowing,
  gives "stable" for PREM-stable.csv and "unstable" for PREM.csv. The largest
  judged e of PREM-stable is 0.029, in the lower mantle, under the 0.05
  threshold.

**Consequences.**

- The lid between 24.4 and 80 km is 57 to 109 kg/m^3 lighter than PREM, 1.7 to
  3.2 per cent. This changes the elastic response at high degree more than the
  60 km table did. The size of this change is not measured yet.
- No density step remains at 60 km, so the table makes no assumption about the
  lid thickness.

**Output.** `PREM.csv` and `PREM-stable.csv`, on the same nodes, centre
outward, in SI units. Columns: radius_m, depth_km, region (index into
`prem1981.REGIONS`), rho_kg_m3, vp_m_s, vs_m_s, kappa_Pa, mu_Pa, g_m_s2. The
attenuation (Q) of Table I is not in `prem1981.py` and is not in the files.
GIA codes do not use it.

## 2026-10-01: Extend the adjustment to the whole mantle

**Decision.** Apply the rule of 2026-09-30 between the surface and the
core-mantle boundary (2891 km), D'' included. Conserve the mass of the crust
and the mantle with one constraint, so mass can move between the upper and the
lower mantle. Keep the core as PREM. Keep the bulk and shear moduli as PREM.

**Reason.** Below 670 km PREM is slightly super-adiabatic (e up to 0.029,
excess up to 9 kg/m^3). The lovejx threshold of e = 0.05 accepts this, but that
threshold was set for lower-mantle viscosities of 5e20 to 5e21 Pa s. A lower
mantle with a lower viscosity can show the instability over a long run. With
the change, the table is stable in every layer of the crust and the mantle.
The core is excluded: the outer core is a fluid without shear strength, and
GIA codes treat it separately.

**Alternatives.** Two mass constraints, one for 0 to 670 km and one for 670 to
2891 km, give densities that differ from the single constraint by less than
0.5 kg/m^3 everywhere. A D'' kept as PREM was rejected, so that the whole
mantle satisfies one rule.

**Result.** From `build_prem_stable.py`. Above 670 km the density moves by
less than 0.5 kg/m^3 from the 0 to 670 km fit.

| Depth (km) | Density change (kg/m^3) |
|---|---|
| 0 to 24.4 | -9.3 to +9.9 |
| 24.4 to 220 | -109 at 24.4 km, +81 at 220 km |
| 220 to 400 | +4.4 to +32 |
| 400 to 600 | +0.3, except down to -12 between 590 and 600 km |
| 600 to 670 | -12 to +16 |
| 670 to 771 | -6.7 to -5.2 |
| 771 to 2741 | -5.2 to +3.2 |
| 2741 to 2891 | +2.3 to +2.5 |

- The largest e in the crust and the mantle is 1.9e-9. PREM's is +1.13.
- The jump at 670 km is 366 kg/m^3 (PREM 389). The jump at the core-mantle
  boundary is 4335 kg/m^3 (PREM 4337).
- The total mass changes by 1e-15 of itself. Surface gravity and gravity in
  the core do not change. Gravity changes by at most 5.5e-3 m/s^2, at 140 km,
  and by 1.4e-4 m/s^2 at 670 km.
- V_P and V_S change by -1.2 to +1.7 per cent.
- `check_stratification` of lovejx, with every solid node marked as flowing,
  gives "stable" for PREM-stable.csv, with the largest judged e 1.9e-9, and no
  inverted jumps.
- Between about 2350 and 2891 km PREM is stable (e down to -0.011). The fit
  gives the adiabatic gradient there (e = 0). The reason is the objective:
  after the change above, the adiabat is the stable profile closest to PREM.

**Figure.** `plot_fit.py` now shows two quantities only, e and the density
adjustment, in two depth ranges: 0 to 1000 km and 1000 to 2891 km. Each range
has its own horizontal scale, so the lower-mantle change of a few kg/m^3 is
visible.

## Open items

- **The cost in the Love numbers.** Run lovejx on PREM.csv and PREM-stable.csv
  with the GIAMIP viscosity profiles. Compare the elastic k and h at high degree
  and the response over a glacial cycle. The draft email states "well under a
  per cent" and needs this number.
- **A lovejx deck.** The CSV files have no viscosity. A lovejx legacy deck needs
  a node at the base of the lid of each viscosity profile.
- **Where the table lives.** Either `giamip/resources/models/` in gia-mip or a
  g-drift dataset. When it is final, the `gadopt` blocks of the gia-mip case
  files can name it in place of `prem`.
- **The email to Holly Han.** The draft in `email/` says that we have a stable
  table. PREM-stable.csv is that table.
