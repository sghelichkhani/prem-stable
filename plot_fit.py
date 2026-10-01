"""Plot the stability measure and the density adjustment of PREM-stable.csv.

Both files come from `build_prem_stable.py` and share one node table, so the
two models are compared row by row without interpolation. The figure has two
columns and two rows:

- left column: e = 1 - Bullen parameter per interval, for PREM and PREM-stable;
- right column: the density adjustment, PREM-stable minus PREM;
- top row: 0 to 1000 km depth, where PREM is far from stable and the
  adjustment reaches about 100 kg/m^3;
- bottom row: 1000 km to the core-mantle boundary, where PREM is close to
  adiabatic and the adjustment is a few kg/m^3. This row has its own x-axis
  range so that the small values are visible.

The moduli are equal in the two files and are not plotted. The changes of V_P,
V_S and gravity follow from the density and are listed in DECISIONS.md.

Run from this directory with the project Python. Writes prem_stable_fit.png.
"""

import pathlib

import matplotlib.pyplot as plt
import numpy as np

from plot_prem_stable import style, C_PREM, C_STABLE, INK, MUTED
from build_prem_stable import stability_e

HERE = pathlib.Path(__file__).resolve().parent
SPLIT_KM = 1000.0       # depth that separates the two rows
CMB_KM = 2891.0         # core-mantle boundary, bottom of the adjusted range


def read(name):
    """Read one of the CSV models as a structured array, centre outward."""
    return np.genfromtxt(HERE / name, delimiter=",", names=True)


def staircase(depth, i, e):
    """Draw a per-interval value as a constant over each interval.

    `i` indexes the intervals (nodes i and i+1) and `e` holds one value per
    interval. Returns x and y arrays for `plot`: every interval becomes a
    vertical segment from its upper to its lower node, so the value is shown
    over the depth it applies to, and not at an interpolated point.
    """
    x = np.repeat(e, 2)
    y = np.column_stack([depth[i], depth[i + 1]]).ravel()
    return x, y


def main():
    p, s = read("PREM.csv"), read("PREM-stable.csv")
    depth, r = p["depth_km"], p["radius_m"]

    # The moduli must be the same in both files; the comparison relies on it.
    assert np.allclose(p["kappa_Pa"], s["kappa_Pa"], rtol=1e-9)
    assert np.allclose(p["mu_Pa"], s["mu_Pa"], rtol=1e-9)

    # Stability measure on every interval of the crust and the mantle. The
    # fluid outer core (mu = 0) and the inner core are not part of the plot.
    mantle = (p["mu_Pa"] > 0) & (depth <= CMB_KM + 1e-9)
    i_p, e_p = stability_e(r, p["rho_kg_m3"], p["kappa_Pa"], p["g_m_s2"], mantle)
    i_s, e_s = stability_e(r, s["rho_kg_m3"], s["kappa_Pa"], s["g_m_s2"], mantle)
    d_rho = s["rho_kg_m3"] - p["rho_kg_m3"]                     # kg/m^3

    fig, ax = plt.subplots(2, 2, figsize=(9, 9.5), constrained_layout=True,
                           gridspec_kw=dict(height_ratios=[1, 1.9]))

    # Each row: depth range, x-range of e, x-range of the density adjustment.
    rows = [(0.0, SPLIT_KM, (-1.3, 1.3), (-120, 100)),
            (SPLIT_KM, CMB_KM, (-0.05, 0.05), (-8, 8))]
    for row, (top, bottom, e_lim, d_lim) in enumerate(rows):
        a_e, a_d = ax[row]

        # e per interval. The shaded half is unstable where the rock flows.
        a_e.axvspan(0, e_lim[1], color="#f6e3dc", zorder=0)
        a_e.axvline(0, color=MUTED, lw=1)
        a_e.plot(*staircase(depth, i_p, e_p), color=C_PREM, lw=2, label="PREM")
        a_e.plot(*staircase(depth, i_s, e_s), color=C_STABLE, lw=2, label="PREM-stable")
        a_e.set(xlim=e_lim, ylim=(bottom, top), ylabel="Depth (km)")

        # Density adjustment. Only PREM-stable changes, so one series. Only
        # the mantle and crust nodes are drawn: the core is unchanged, and a
        # line to the core node would draw a false step at the CMB.
        a_d.axvline(0, color=MUTED, lw=1)
        a_d.plot(d_rho[mantle], depth[mantle], color=C_STABLE, lw=2)
        a_d.set(xlim=d_lim, ylim=(bottom, top))

        for a in (a_e, a_d):
            style(a)

    # Labels: titles on the top row, x-axis labels on the bottom row.
    ax[0, 0].set_title("a. Stratification, 0 to 1000 km", loc="left", color=INK)
    ax[0, 1].set_title("b. Density adjustment, 0 to 1000 km", loc="left", color=INK)
    ax[1, 0].set_title("c. Stratification, 1000 to 2891 km", loc="left", color=INK)
    ax[1, 1].set_title("d. Density adjustment, 1000 to 2891 km", loc="left", color=INK)
    ax[1, 0].set_xlabel("e = 1 − Bullen parameter")
    ax[1, 1].set_xlabel("PREM-stable − PREM (kg/m³)")
    ax[0, 0].text(0.65, 980, "unstable", ha="center", va="bottom", color=INK, fontsize=9)
    ax[0, 0].text(-0.65, 980, "stable", ha="center", va="bottom", color=INK, fontsize=9)
    # The legend goes in the empty lower right part of panel a, between the
    # PREM line near e = 0 and the "unstable" label.
    ax[0, 0].legend(frameon=False, loc="lower left", bbox_to_anchor=(0.6, 0.1), labelcolor=INK)

    out = HERE / "prem_stable_fit.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out.name}")

    # Short numerical summary next to the figure.
    print(f"density adjustment: {d_rho.min():+.1f} to {d_rho.max():+.1f} kg/m^3")
    print(f"largest e in the crust and mantle: PREM {e_p.max():+.3f}, "
          f"PREM-stable {e_s.max():+.1e}")


if __name__ == "__main__":
    main()
