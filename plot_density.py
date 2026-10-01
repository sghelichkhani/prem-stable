"""Plot the density of PREM and of PREM-stable, the simplest view of the change.

PREM is the density as published. PREM-stable is the density that the fit in
`build_prem_stable.py` gives: the profile closest to PREM that is stable
everywhere in the crust and the mantle. Two panels:

- a. the whole crust and mantle, 0 to 2891 km. At this scale the two profiles
  differ visibly only above 400 km, which shows how small the change is
  overall;
- b. the top 700 km, where the change is largest: the density inversion of
  PREM between 24.4 and 220 km becomes a density that increases with depth.

The core is not plotted, because it is the same in both models.

Run from this directory with the project Python. Writes prem_density.png.
"""

import pathlib

import matplotlib.pyplot as plt
import numpy as np

from plot_prem_stable import style, C_PREM, C_STABLE, INK

HERE = pathlib.Path(__file__).resolve().parent
CMB_KM = 2891.0         # core-mantle boundary
ZOOM_KM = 700.0         # depth range of the right panel


def read(name):
    """Read one of the CSV models as a structured array, centre outward."""
    return np.genfromtxt(HERE / name, delimiter=",", names=True)


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

    fig, ax = plt.subplots(1, 2, figsize=(9, 6.5), constrained_layout=True)
    panels = [(ax[0], CMB_KM, (2500, 5700), "a. Crust and mantle"),
              (ax[1], ZOOM_KM, (2500, 4500), "b. Top 700 km")]
    for a, bottom, xlim, title in panels:
        # PREM is drawn wider and underneath, so that it shows on both sides
        # of PREM-stable where the two lines coincide.
        a.plot(rho_p, z, color=C_PREM, lw=3.5, label="PREM, as published")
        a.plot(rho_s, z, color=C_STABLE, lw=1.8, label="PREM-stable")
        a.set(xlim=xlim, ylim=(bottom, 0), xlabel="Density (kg/m³)")
        a.set_title(title, loc="left", color=INK)
        style(a)
    ax[0].set_ylabel("Depth (km)")
    ax[1].legend(frameon=False, loc="lower left", labelcolor=INK)

    out = HERE / "prem_density.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out.name}")


if __name__ == "__main__":
    main()
