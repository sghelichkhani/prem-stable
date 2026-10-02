"""Plot the density of PREM, of PREM-stable and of VM5i, the simplest view of the change.

PREM is the density as published. PREM-stable is the density that the fit in
`build_prem_stable.py` gives: the profile closest to PREM that is stable
everywhere in the crust and the mantle. VM5i is the 11-layer model of Spada
and Melini (2019, Table 2), a volume average of PREM in uniform layers that is
used in incompressible GIA codes. It is drawn for comparison only. Two panels:

- a. the whole crust and mantle, 0 to 2891 km. At this scale the two profiles
  differ visibly only above 400 km, which shows how small the change is
  overall;
- b. the top 700 km, where the change is largest: the density inversion of
  PREM between 24.4 and 220 km becomes a density that increases with depth.

The core is not plotted. It is the same in PREM and PREM-stable, and its
density is far outside the plotted range.

Run from this directory with the project Python. Writes prem_density.png.
"""

import pathlib

import matplotlib.pyplot as plt
import numpy as np


HERE = pathlib.Path(__file__).resolve().parent

# Colours: the first two categorical slots of the dataviz reference palette,
# checked for colour-blind separation and contrast against a light background.
C_PREM = "#2a78d6"      # blue, PREM as published
C_STABLE = "#eb6834"    # orange, PREM-stable
# VM5i is a reference and not one of the two models of this repository, so it
# gets a neutral grey and not a third categorical colour. The grey has a
# contrast of about 5:1 against the white background.
C_VM5I = "#6b6a64"      # grey, VM5i of Spada and Melini (2019)
INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3dc"
CMB_KM = 2891.0         # core-mantle boundary
ZOOM_KM = 700.0         # depth range of the right panel


def read(name):
    """Read one of the CSV models as a structured array, centre outward."""
    return np.genfromtxt(HERE / name, delimiter=",", names=True)


def read_vm5i():
    """Read the VM5i layers and return them as a step profile for `plot`.

    The file holds one row per layer, from the surface down, with the depth of
    the top and the bottom of each layer in km and one density in kg/m^3. The
    density is constant inside a layer, so each layer becomes a vertical
    segment, and consecutive segments join with a horizontal step at the
    layer boundary. The core row is left out.

    Returns
    -------
    (rho, depth)
        Arrays for `plot(rho, depth)`, two points per layer.
    """
    t = np.genfromtxt(HERE / "VM5i_Spada_Melini_2019.csv", delimiter=",",
                      names=True, dtype=None, encoding="utf-8")
    t = t[t["layer"] != "Core"]
    rho = np.repeat(t["density_kg_m3"], 2)
    depth = np.column_stack([t["depth_top_km"], t["depth_bottom_km"]]).ravel()
    return rho, depth


def style(ax):
    """Light grid, muted spines, so the two density profiles stand out."""
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelcolor=INK)


def main():
    p, s = read("PREM.csv"), read("PREM-stable.csv")
    depth = p["depth_km"]
    # Crust and mantle nodes only. The outer-core node at the core-mantle
    # boundary has the same depth as the base of the mantle, so the shear
    # modulus (zero in the fluid outer core) separates the two. The nodes run
    # from the centre outward, so this selection is one contiguous block.
    mantle = (depth <= CMB_KM + 1e-9) & (p["mu_Pa"] > 0)
    z = depth[mantle]
    rho_p = p["rho_kg_m3"][mantle]      # kg/m^3
    rho_s = s["rho_kg_m3"][mantle]      # kg/m^3
    rho_v, z_v = read_vm5i()            # kg/m^3 and km, step profile

    fig, ax = plt.subplots(1, 2, figsize=(9, 6.5), constrained_layout=True)
    panels = [(ax[0], CMB_KM, (2500, 5700), "a. Crust and mantle"),
              (ax[1], ZOOM_KM, (2500, 4500), "b. Top 700 km")]
    for a, bottom, xlim, title in panels:
        # PREM is drawn wider and underneath, so that it shows on both sides
        # of PREM-stable where the two lines coincide.
        a.plot(rho_p, z, color=C_PREM, lw=3.5, label="PREM, as published")
        a.plot(rho_s, z, color=C_STABLE, lw=1.8, label="PREM-stable")
        # VM5i on top, thin, so its steps stay visible where they cross the
        # two smooth profiles.
        a.plot(rho_v, z_v, color=C_VM5I, lw=1.4,
               label="VM5i (Spada and Melini, 2019)")
        a.set(xlim=xlim, ylim=(bottom, 0), xlabel="Density (kg/m³)")
        a.set_title(title, loc="left", color=INK)
        style(a)
    ax[0].set_ylabel("Depth (km)")
    # The lower left of panel a is empty: below 1000 km all three profiles
    # are denser than 4600 kg/m^3, on the right of the panel.
    ax[0].legend(frameon=False, loc="lower left", labelcolor=INK)

    out = HERE / "prem_density.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out.name}")


if __name__ == "__main__":
    main()
