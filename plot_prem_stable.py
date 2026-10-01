"""Compare PREM as published with the stabilised PREM of lovejx.

Two sources:

- PREM as published: the TauP table that ships with ObsPy,
  `obspy/taup/data/prem.nd`. Columns are depth (km), V_P, V_S (km/s),
  density (g/cm^3), Q_P and Q_S. A discontinuity is two lines at the same
  depth. The top 15 km is crust (2.6 g/cm^3); the table has no ocean layer.
- PREM-stable: `lovejx/models/model.prem-stable.l96.ump5.lm5`, the 171-node
  legacy PREM deck with the density gradient inside every flowing mantle region
  clipped at the adiabatic gradient. Bulk and shear moduli are unchanged; only
  density, and through it the wave speeds and gravity, differ. The construction
  is recorded in lovejx `NOTES/STRATIFICATION.md`.

The stability measure per interval between two nodes is

    A = d rho/dr + rho^2 g / kappa                    [kg/m^4]
    e = A kappa / (rho^2 g) = 1 - eta_Bullen          [dimensionless]

with r the radius, so d rho/dr < 0 is density increasing with depth. A parcel
moved upward decompresses along -rho^2 g / kappa. If the surrounding density
falls faster than that with height, the parcel arrives denser and sinks back
(A < 0, stable). A > 0 is super-adiabatic and gravitationally unstable where the
material flows. e = 0 is the Adams-Williamson (adiabatic) gradient; lovejx
flags a flowing interval as unstable above e = 0.05.

Run from anywhere with the project Python; writes a PNG next to this file.
"""

import pathlib

import matplotlib.pyplot as plt
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
PREM_ND = pathlib.Path.home() / "Workplace/python3.12/lib/python3.12/site-packages/obspy/taup/data/prem.nd"
DECKS = pathlib.Path.home() / "Workplace/lovejx/lovejx/models"
G = 6.674e-11           # m^3 kg^-1 s^-2
R_EARTH_KM = 6371.0

# Colours: the first two categorical slots of the dataviz reference palette.
C_PREM = "#2a78d6"      # blue, PREM as published
C_STABLE = "#eb6834"    # orange, PREM-stable
INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3dc"


def read_nd(path):
    """Read a TauP .nd table and convert to SI.

    Returns depth (km), radius (m), density (kg/m^3), kappa and mu (Pa), all
    ordered from the surface down, with discontinuities kept as repeated depths.
    """
    rows = [l.split() for l in open(path) if l.strip() and l.split()[0][0].isdigit()]
    depth, vp, vs, rho, _, _ = np.array(rows, dtype=float).T
    rho = rho * 1e3                     # g/cm^3 -> kg/m^3
    vp, vs = vp * 1e3, vs * 1e3         # km/s -> m/s
    mu = rho * vs**2
    kappa = rho * vp**2 - 4.0 * mu / 3.0
    return depth, (R_EARTH_KM - depth) * 1e3, rho, kappa, mu


def gravity_from_density(radius, rho):
    """Integrate g(r) = G M(r) / r^2 from the centre for a surface-down table.

    Density is taken as linear between nodes, which is what the table implies;
    PREM's own polynomials are at most cubic and the nodes are dense enough in
    the mantle for this to matter only at the 1e-3 level.
    """
    r, d = radius[::-1], rho[::-1]                  # centre outward
    shell = np.pi * 4.0 * np.diff(r) * 0.5 * (d[1:] * r[1:]**2 + d[:-1] * r[:-1]**2)
    mass = np.concatenate([[0.0], np.cumsum(shell)])
    g = np.where(r > 0, G * mass / np.where(r > 0, r, 1.0)**2, 0.0)
    return g[::-1]


def read_deck(name):
    """Read a lovejx legacy deck (radius, rho, vp, vs, g, eta), surface first."""
    r, rho, vp, vs, g, eta = np.loadtxt(DECKS / name).T[:, ::-1]
    mu = rho * vs**2
    kappa = rho * vp**2 - 4.0 * mu / 3.0
    return R_EARTH_KM - r / 1e3, r, rho, kappa, g, eta


def stability(depth, radius, rho, kappa, g):
    """The measure e per interval of finite width, on interval means.

    Returns the mid-depth of each interval and e. Zero-width intervals are the
    discontinuities; they are skipped here because a jump is judged on its own
    sign, not on a gradient.
    """
    dr = np.diff(radius)
    keep = np.abs(dr) > 1.0
    drho_dr = np.diff(rho)[keep] / dr[keep]
    rm = 0.5 * (rho[1:] + rho[:-1])[keep]
    km = 0.5 * (kappa[1:] + kappa[:-1])[keep]
    gm = 0.5 * (g[1:] + g[:-1])[keep]
    # Only the solid mantle and crust matter here; the fluid core has mu = 0
    # but kappa > 0, so the formula still holds and is simply not plotted.
    a = drho_dr + rm**2 * gm / km
    e = a * km / (rm**2 * gm)
    return 0.5 * (depth[1:] + depth[:-1])[keep], e


def style(ax):
    """Recessive axes: light grid, muted spines, depth increasing downward."""
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelcolor=INK)


def main():
    nd_depth, nd_r, nd_rho, nd_kappa, _ = read_nd(PREM_ND)
    nd_g = gravity_from_density(nd_r, nd_rho)
    p_depth, p_r, p_rho, p_kappa, p_g, _ = read_deck("model.prem.l96.ump5.lm5")
    s_depth, s_r, s_rho, s_kappa, s_g, s_eta = read_deck("model.prem-stable.l96.ump5.lm5")

    # The lovejx PREM deck and the TauP table are the same model; check that
    # before using one as the other.
    interp = np.interp(p_depth, nd_depth, nd_rho)   # both tables run surface down
    # Skip nodes on a discontinuity: there np.interp cannot tell which side of
    # the jump a node belongs to.
    jumps = nd_depth[:-1][np.diff(nd_depth) == 0]
    off_jump = np.min(np.abs(p_depth[:, None] - jumps[None, :]), axis=1) > 0.5
    mantle = (p_depth > 25) & (p_depth < 2885) & off_jump
    print("max |deck - prem.nd| density in the mantle: %.1f kg/m^3"
          % np.max(np.abs(p_rho - interp)[mantle]))
    print("max |stable - PREM| density: %.1f kg/m^3" % np.max(np.abs(s_rho - p_rho)))
    print("kappa unchanged to %.1e relative"
          % np.max(np.abs(s_kappa - p_kappa) / p_kappa))

    # Depth of the elastic lid in the deck (viscosity 1e42 and above).
    lid_km = s_depth[s_eta >= 1e40].max()

    nd_mid, nd_e = stability(nd_depth, nd_r, nd_rho, nd_kappa, nd_g)
    s_mid, s_e = stability(s_depth, s_r, s_rho, s_kappa, s_g)

    fig, axes = plt.subplots(1, 3, figsize=(12, 6.2), constrained_layout=True)
    zoom = 800.0

    # (a) Density in the upper mantle.
    ax = axes[0]
    ax.plot(nd_rho, nd_depth, color=C_PREM, lw=2, label="PREM as published")
    ax.plot(s_rho, s_depth, color=C_STABLE, lw=2, label="PREM-stable")
    ax.set_ylim(zoom, 0)
    ax.set_xlim(3250, 4050)
    ax.set_xlabel("Density (kg/m³)")
    ax.set_ylabel("Depth (km)")
    ax.set_title("a. Density", loc="left", color=INK)
    ax.legend(frameon=False, loc="lower left")

    # (b) Stability measure. Right of zero is unstable.
    ax = axes[1]
    ax.axvspan(0.0, 1.3, color="#f6e3dc", zorder=0)
    ax.axvline(0.0, color=MUTED, lw=1)
    ax.axvline(0.05, color=MUTED, lw=1, ls=":")
    ax.step(nd_e, nd_mid, where="mid", color=C_PREM, lw=2)
    ax.step(s_e, s_mid, where="mid", color=C_STABLE, lw=2)
    ax.set_ylim(zoom, 0)
    ax.set_xlim(-1.3, 1.3)
    ax.set_xlabel("e = 1 − Bullen parameter")
    ax.set_title("b. Stratification", loc="left", color=INK)
    ax.text(0.65, 20, "unstable\n(super-adiabatic)", ha="center", va="top",
            color=INK, fontsize=9)
    ax.text(-0.65, 20, "stable", ha="center", va="top", color=INK, fontsize=9)

    # (c) Density change over the whole mantle.
    ax = axes[2]
    ax.axvline(0.0, color=MUTED, lw=1)
    ax.plot(s_rho - p_rho, s_depth, color=C_STABLE, lw=2)
    ax.set_ylim(2891, 0)
    ax.set_xlabel("PREM-stable − PREM (kg/m³)")
    ax.set_title("c. Density change, whole mantle", loc="left", color=INK)

    for ax in axes:
        style(ax)
        # The lid does not flow, so its stratification carries no mode.
        ax.axhspan(0, lid_km, color="#ecebe4", zorder=0)
    axes[0].text(3640, lid_km / 2, f"elastic lid (0–{lid_km:.0f} km)",
                 va="center", color=MUTED, fontsize=8.5)

    out = HERE / "prem_vs_prem_stable.png"
    fig.savefig(out, dpi=160)
    print("wrote", out)


if __name__ == "__main__":
    main()
