"""Plot PREM-stable.csv against PREM.csv: every quantity that the fit changes.

Both files come from `build_prem_stable.py` and share one node table, so the
two models are compared row by row without interpolation. The bulk and shear
moduli are equal in the two files and are not plotted. The quantities that
change are density, the stability measure e that follows from it, the wave
speeds and gravity.

With the moduli fixed, V_P = sqrt((kappa + 4 mu / 3) / rho) and
V_S = sqrt(mu / rho) both scale as rho^(-1/2). Their relative changes are
therefore identical, (rho_PREM / rho_stable)^(1/2) - 1, and one curve shows
both.

Run from this directory with the project Python. Writes prem_stable_fit.png.
"""

import pathlib

import matplotlib.pyplot as plt
import numpy as np

from plot_prem_stable import style, C_PREM, C_STABLE, INK, MUTED
from build_prem_stable import stability_e

HERE = pathlib.Path(__file__).resolve().parent
DEPTH_MAX_KM = 800.0    # plotted depth range; below 670 km the models are equal
ADJ_BOTTOM_KM = 670.0   # bottom of the adjusted range


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

    # Stability measure on every solid interval below the surface. The fluid
    # outer core is excluded (mu = 0), as it is in the lovejx check.
    solid = p["mu_Pa"] > 0
    i_p, e_p = stability_e(r, p["rho_kg_m3"], p["kappa_Pa"], p["g_m_s2"], solid)
    i_s, e_s = stability_e(r, s["rho_kg_m3"], s["kappa_Pa"], s["g_m_s2"], solid)

    d_rho = s["rho_kg_m3"] - p["rho_kg_m3"]                     # kg/m^3
    d_v = 100.0 * (s["vp_m_s"] / p["vp_m_s"] - 1.0)             # per cent
    d_g = 1e3 * (s["g_m_s2"] - p["g_m_s2"])                     # mm/s^2

    fig, ax = plt.subplots(1, 5, figsize=(16, 6.4), sharey=True, constrained_layout=True)

    # a. Density of both models.
    ax[0].plot(p["rho_kg_m3"], depth, color=C_PREM, lw=2, label="PREM")
    ax[0].plot(s["rho_kg_m3"], depth, color=C_STABLE, lw=2, label="PREM-stable")
    ax[0].set(xlim=(2500, 4500), xlabel="Density (kg/m³)", ylabel="Depth (km)")
    ax[0].set_title("a. Density", loc="left", color=INK)
    ax[0].legend(frameon=False, loc="lower left", labelcolor=INK)

    # b. Stability measure. e > 0 is denser than the adiabat would allow.
    ax[1].axvspan(0, 1.3, color="#f6e3dc", zorder=0)
    ax[1].axvline(0, color=MUTED, lw=1)
    ax[1].plot(*staircase(depth, i_p, e_p), color=C_PREM, lw=2)
    ax[1].plot(*staircase(depth, i_s, e_s), color=C_STABLE, lw=2)
    ax[1].set(xlim=(-1.3, 1.3), xlabel="e = 1 − Bullen parameter")
    ax[1].set_title("b. Stratification", loc="left", color=INK)
    ax[1].text(0.65, 785, "unstable", ha="center", va="bottom", color=INK, fontsize=9)
    ax[1].text(-0.65, 785, "stable", ha="center", va="bottom", color=INK, fontsize=9)

    # c. to e. The differences. One series each, so no legend: the title names it.
    for a, x, label, title in [
            (ax[2], d_rho, "PREM-stable − PREM (kg/m³)", "c. Density change"),
            (ax[3], d_v, "Relative change (%)", "d. V_P and V_S change"),
            (ax[4], d_g, "PREM-stable − PREM (mm/s²)", "e. Gravity change")]:
        a.axvline(0, color=MUTED, lw=1)
        a.plot(x, depth, color=C_STABLE, lw=2)
        a.set(xlabel=label)
        a.set_title(title, loc="left", color=INK)

    for a in ax:
        style(a)
        a.set_ylim(DEPTH_MAX_KM, 0)
        # Bottom of the adjusted range. Below it both models are PREM.
        a.axhline(ADJ_BOTTOM_KM, color=MUTED, lw=1, ls=(0, (4, 3)))
    ax[2].text(ax[2].get_xlim()[0], ADJ_BOTTOM_KM + 8, " PREM below 670 km",
               ha="left", va="top", color=MUTED, fontsize=9)

    out = HERE / "prem_stable_fit.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out.name}")

    # Short numerical summary next to the figure.
    adj = depth <= ADJ_BOTTOM_KM
    print(f"density change: {d_rho.min():+.1f} to {d_rho.max():+.1f} kg/m^3")
    print(f"wave-speed change: {d_v.min():+.2f} to {d_v.max():+.2f} per cent")
    print(f"gravity change: {d_g.min():+.2f} to {d_g.max():+.2f} mm/s^2")
    print(f"largest e above 670 km: PREM {e_p[adj[i_p]].max():+.3f}, "
          f"PREM-stable {e_s[adj[i_s]].max():+.1e}")


if __name__ == "__main__":
    main()
