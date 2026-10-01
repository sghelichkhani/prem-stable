"""Check the transcribed PREM polynomials of prem1981.py against two tables.

1. ObsPy's prem.nd, the TauP table. It gives values to 5 decimals. In the LID
   and the LVZ (24.4 to 220 km) it uses the Voigt average of the transversely
   isotropic coefficients, not the footnote formulas, so there it is compared
   with `prem1981.voigt_isotropic`. Elsewhere it agrees with the polynomials to
   about 1e-4, which is larger than its rounding: prem.nd was not computed
   exactly from the printed 4-decimal coefficients.
2. The legacy 171-node PREM deck that lovejx ships. It gives density to about
   seven significant figures, so it can detect a mistyped last digit of any
   density coefficient. Its wave speeds in the LID and the LVZ match no form in
   the paper and are not compared.

Run from this directory with the project Python.
"""

import pathlib

import numpy as np

import prem1981 as P

PREM_ND = pathlib.Path.home() / "Workplace/python3.12/lib/python3.12/site-packages/obspy/taup/data/prem.nd"
DECK = pathlib.Path.home() / "Workplace/lovejx/lovejx/models/model.prem.l96.ump5.lm5"
LID_LVZ = ("LVZ, 220 to 80 km", "LID, 80 to 24.4 km")


def side_at(r, i, surface_first):
    """Which side of a discontinuity node i belongs to, from the node order."""
    nxt = i + 1 < len(r) and r[i + 1] == r[i]
    prv = i > 0 and r[i - 1] == r[i]
    if not (nxt or prv):
        return "below"
    first = nxt
    # Surface-first tables list the upper side first; centre-first tables the lower.
    return ("above" if first else "below") if surface_first else ("below" if first else "above")


def report(title, table):
    print(title)
    for name, errs in table.items():
        e = np.array(errs)
        print(f"  {name:34s} rho {e[:, 0].max():.1e}  vp {e[:, 1].max():.1e}  vs {e[:, 2].max():.1e}  n={len(e)}")
    print()


# 1. prem.nd, surface first, km and g/cm^3.
rows = [l.split() for l in open(PREM_ND) if l.strip() and l.split()[0][0].isdigit()]
depth, vp, vs, rho, _, _ = np.array(rows, dtype=float).T
r = P.A_KM - depth
table = {}
for i, ri in enumerate(r):
    k = P.region_of(ri, side_at(r, i, surface_first=True))
    name = P.REGIONS[k][2]
    ours = list(P.evaluate(k, ri))
    if name in LID_LVZ:
        ours[1], ours[2] = P.voigt_isotropic(ri)
    table.setdefault(name, []).append(np.abs(np.array(ours) - [rho[i], vp[i], vs[i]]))
report("prem.nd (LID and LVZ speeds against the Voigt average):", table)

# 2. lovejx legacy deck, centre first, SI units.
t = np.loadtxt(DECK)
r = t[:, 0] / 1e3
table = {}
for i, ri in enumerate(r):
    if ri > 6368.0:            # the deck's ocean-free top; same crust as ours
        continue
    k = P.region_of(ri, side_at(r, i, surface_first=False))
    name = P.REGIONS[k][2]
    ours = np.array(P.evaluate(k, ri))
    theirs = t[i, 1:4] / 1e3
    err = np.abs(ours - theirs)
    if name in LID_LVZ:
        err[1:] = np.nan       # speeds of unknown origin, not compared
    table.setdefault(name, []).append(err)
print("lovejx deck (LID and LVZ speeds not compared):")
for name, errs in table.items():
    e = np.array(errs)
    print(f"  {name:34s} rho {np.nanmax(e[:, 0]):.1e}  vp {np.nanmax(e[:, 1]) if not np.all(np.isnan(e[:, 1])) else float('nan'):.1e}  "
          f"vs {np.nanmax(e[:, 2]) if not np.all(np.isnan(e[:, 2])) else float('nan'):.1e}  n={len(e)}")
