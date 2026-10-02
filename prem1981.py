"""PREM evaluated from the polynomials of Dziewonski and Anderson (1981), Table I.

Source: Dziewonski, A. M. and Anderson, D. L. (1981). Preliminary reference
Earth model. Physics of the Earth and Planetary Interiors 25(4), 297-356,
Table I, transcribed from the printed table.

Every coefficient below is typed from that table and nothing else. To check a
coefficient, compare it with Table I of the paper.

Conventions of the table:

- x = r / a with a = 6371 km.
- Density in g/cm^3, wave speeds in km/s, at a reference period of 1 s.
- Radii in km. Each region is closed at both ends; at a shared radius the two
  regions give the two sides of a discontinuity.

Two choices are made here, both stated in the table or its footnote:

- Isotropic wave speeds. Between 24.4 and 220 km depth (the LID and the LVZ,
  radii 6151 to 6346.6 km) PREM is transversely isotropic. The footnote of
  Table I gives the effective isotropic speeds used here,
  V_P = 4.1875 + 3.9382x and V_S = 2.1519 + 2.3481x. Density in that range is
  isotropic and needs no approximation.
- No ocean. The 3 km ocean (6368 to 6371 km) is replaced by the upper crust,
  2.6 g/cm^3, 5.8 and 3.2 km/s, extended to the surface. A GIA model loads the
  ocean explicitly, so a fluid layer in the Earth model is not wanted.
"""

from __future__ import annotations

import numpy as np

A_KM = 6371.0

# One entry per region: (r_bottom_km, r_top_km, name, rho, vp, vs), where each
# of rho, vp, vs is the coefficient list (c0, c1, c2, c3) of c0 + c1 x + c2 x^2
# + c3 x^3. Coefficients absent from the table are zero.
REGIONS = [
    (0.0, 1221.5, "inner core",
     (13.0885, 0.0, -8.8381, 0.0),
     (11.2622, 0.0, -6.3640, 0.0),
     (3.6678, 0.0, -4.4475, 0.0)),
    (1221.5, 3480.0, "outer core",
     (12.5815, -1.2638, -3.6426, -5.5281),
     (11.0487, -4.0362, 4.8023, -13.5732),
     (0.0, 0.0, 0.0, 0.0)),
    (3480.0, 3630.0, "lower mantle, D''",
     (7.9565, -6.4761, 5.5283, -3.0807),
     (15.3891, -5.3181, 5.5242, -2.5514),
     (6.9254, 1.4672, -2.0834, 0.9783)),
    (3630.0, 5600.0, "lower mantle",
     (7.9565, -6.4761, 5.5283, -3.0807),
     (24.9520, -40.4673, 51.4832, -26.6419),
     (11.1671, -13.7818, 17.4575, -9.2777)),
    (5600.0, 5701.0, "lower mantle, top",
     (7.9565, -6.4761, 5.5283, -3.0807),
     (29.2766, -23.6027, 5.5242, -2.5514),
     (22.3459, -17.2473, -2.0834, 0.9783)),
    (5701.0, 5771.0, "transition zone, 670 to 600 km",
     (5.3197, -1.4836, 0.0, 0.0),
     (19.0957, -9.8672, 0.0, 0.0),
     (9.9839, -4.9324, 0.0, 0.0)),
    (5771.0, 5971.0, "transition zone, 600 to 400 km",
     (11.2494, -8.0298, 0.0, 0.0),
     (39.7027, -32.6166, 0.0, 0.0),
     (22.3512, -18.5856, 0.0, 0.0)),
    (5971.0, 6151.0, "upper mantle, 400 to 220 km",
     (7.1089, -3.8045, 0.0, 0.0),
     (20.3926, -12.2569, 0.0, 0.0),
     (8.9496, -4.4597, 0.0, 0.0)),
    (6151.0, 6291.0, "LVZ, 220 to 80 km",
     (2.6910, 0.6924, 0.0, 0.0),
     (4.1875, 3.9382, 0.0, 0.0),      # isotropic, footnote of Table I
     (2.1519, 2.3481, 0.0, 0.0)),     # isotropic, footnote of Table I
    (6291.0, 6346.6, "LID, 80 to 24.4 km",
     (2.6910, 0.6924, 0.0, 0.0),
     (4.1875, 3.9382, 0.0, 0.0),      # isotropic, footnote of Table I
     (2.1519, 2.3481, 0.0, 0.0)),     # isotropic, footnote of Table I
    (6346.6, 6356.0, "lower crust",
     (2.900, 0.0, 0.0, 0.0),
     (6.800, 0.0, 0.0, 0.0),
     (3.900, 0.0, 0.0, 0.0)),
    # The upper crust, 6356 to 6368 km in the table, extended to the surface in
    # place of the ocean (6368 to 6371 km: 1.020 g/cm^3, 1.45 km/s, V_S = 0).
    (6356.0, 6371.0, "upper crust, ocean replaced",
     (2.600, 0.0, 0.0, 0.0),
     (5.800, 0.0, 0.0, 0.0),
     (3.200, 0.0, 0.0, 0.0)),
]


def _poly(c, x):
    """c0 + c1 x + c2 x^2 + c3 x^3."""
    return c[0] + x * (c[1] + x * (c[2] + x * c[3]))


def evaluate(region: int, r_km):
    """Density (g/cm^3), V_P and V_S (km/s) of one region at radii r_km."""
    lo, hi, _, rho, vp, vs = REGIONS[region]
    r = np.asarray(r_km, dtype=float)
    if np.any(r < lo - 1e-9) or np.any(r > hi + 1e-9):
        raise ValueError(f"radius outside region {REGIONS[region][2]}")
    x = r / A_KM
    return _poly(rho, x), _poly(vp, x), _poly(vs, x)


def region_of(r_km: float, side: str = "below") -> int:
    """The region that holds radius r_km.

    At a discontinuity two regions share the radius. `side` picks the region
    below it ("below") or above it ("above").
    """
    for i, (lo, hi, *_rest) in enumerate(REGIONS):
        if side == "below" and lo < r_km <= hi:
            return i
        if side == "above" and lo <= r_km < hi:
            return i
    if r_km == 0.0:
        return 0
    if r_km == A_KM:
        return len(REGIONS) - 1
    raise ValueError(r_km)
