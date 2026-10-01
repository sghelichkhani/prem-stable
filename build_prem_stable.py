"""Build a gravitationally stable PREM from the published polynomials.

Starting point
--------------
PREM of Dziewonski and Anderson (1981), Table I, evaluated in `prem1981.py`:
isotropic wave speeds from the footnote of the table, the ocean replaced by the
upper crust, reference period 1 s. `check_prem1981.py` checks the transcription.

The node table
--------------
Each region of Table I gets evenly spaced nodes, at most 10 km apart from the
surface to the core-mantle boundary and at most 100 km apart in the core. The
spacing inside a region is the region thickness divided by the smallest number
of intervals that keeps it at or below that limit, so a node falls on every
region boundary. Every region boundary is written twice, once for each side,
with the lower side first. This is the convention of lovejx: a duplicated
radius is a declared interface. Some region boundaries (80 km, 600 km, 771 km,
2741 km) have no density jump in PREM, and the two rows there are equal.

What is changed
---------------
Only the density of the crust and the mantle, between the surface and the
core-mantle boundary at 2891 km depth. In that range PREM is super-adiabatic in
the crust (constant density in each layer, e = 1), between 24.4 and 220 km
(density increases upward, e = 1.13), between 220 and 400 km (e = 0.2), between
600 and 670 km (e = 0.63) and in parts of the lower mantle (e up to 0.029).
The stable transition zone between 400 and 600 km is inside the range, so that
the fit can move it if the constraints need that. It carries no constraint of
its own beyond staying stable. The core is PREM exactly. The outer core is a
fluid without shear strength, and GIA codes treat it separately.

The range includes the lithosphere and the crust. A table that is stable only
below a rigid lid is stable only for models whose lid is at least that thick,
and the GIAMIP models and the G-ADOPT models with a finite lid viscosity do not
all satisfy this. The reasons are in `DECISIONS.md`.

Inside the range the stable density is taken as linear between nodes. This is
the model that the file describes, and the mass and gravity below are computed
for it exactly. Outside the range, mass and gravity come from the PREM
polynomials, integrated exactly by Gauss-Legendre quadrature.

The rule
--------
The new density is the one closest to PREM at the nodes, in the mass-weighted
least-squares sense, subject to:

1. Stability. In every interval of the range the gradient is at most the
   adiabatic one, d rho/dr + rho_m^2 g_m / kappa_m <= 0, where the subscript m
   is the mean of the two end nodes, kappa is PREM's and g comes from the new
   density. This is the discrete measure that lovejx applies
   (`lovejx.lovenumbers.stratification.interval_measure`), so the check passes
   on the table exactly as it is written.
2. No inverted discontinuity: at every region boundary in the range the density
   below is at least the density above. At the core-mantle boundary the core,
   which is PREM, must be at least as dense as the base of the mantle. With
   PREM's core this condition is far from active.
3. The mass between the surface and the core-mantle boundary is PREM's. One
   constraint covers the whole range, so mass can move between the upper and
   the lower mantle. The total mass, the surface gravity and the gravity in the
   core are therefore PREM's.

Condition 1 is non-linear in the density because of rho^2 and g. It is solved by
fixed-point iteration: the coefficients rho_m g_m / kappa_m come from the
previous iterate, the fit is a quadratic programme with linear constraints, and
the loop stops when the density changes by less than 2e-5 kg/m^3. That is the
accuracy of the optimiser: below it, successive iterates move by 1e-6 to 1e-5
kg/m^3 without a trend.

Bulk and shear moduli stay PREM's. The wave speeds are recomputed from them and
the new density, so V_P and V_S change by about half the relative density
change, in the opposite direction.

Output
------
`PREM.csv` and `PREM-stable.csv`, on the same nodes, from the centre outward,
in SI units. Columns: radius (m), depth (km), region (index into
`prem1981.REGIONS`), density (kg/m^3), V_P and V_S (m/s), bulk and shear
modulus (Pa), gravity (m/s^2). lovejx's legacy decks carry radius, density,
V_P, V_S and gravity in the same units, so these columns convert directly.
"""

from __future__ import annotations

import pathlib

import numpy as np
from scipy.optimize import minimize, LinearConstraint

import prem1981 as P

HERE = pathlib.Path(__file__).resolve().parent
G = 6.6743e-11              # CODATA 2018, m^3 kg^-1 s^-2
R_BOT_KM = 3480.0          # core-mantle boundary radius, bottom of the adjusted range
STEP_MANTLE_KM = 10.0       # largest node spacing from the surface to the CMB
STEP_CORE_KM = 100.0        # largest node spacing in the core
R_CMB_KM = 3480.0

# Three-point Gauss-Legendre rule on [-1, 1]. It is exact for polynomials of
# degree 5, which covers r^2 times a cubic density (the PREM polynomials are at
# most cubic in x = r/a), so every shell mass below is exact up to round-off.
GL_X, GL_W = np.polynomial.legendre.leggauss(3)

CSV_HEADER = ("radius_m,depth_km,region,rho_kg_m3,vp_m_s,vs_m_s,"
              "kappa_Pa,mu_Pa,g_m_s2")


# --- PREM on the node table --------------------------------------------------

def prem_nodes():
    """PREM from the centre to the surface on the output node table, in SI.

    Returns a dict with radius (m), rho (kg/m^3), vp and vs (m/s), kappa and
    mu (Pa), the region index of each node, and a boolean mask `adj` of the
    nodes whose density the fit changes. Region boundaries appear as two nodes
    at one radius, the lower region's node first.
    """
    r_all, reg_all = [], []
    for k, (lo, hi, *_rest) in enumerate(P.REGIONS):
        # The core gets the coarse spacing, everything above the CMB the fine one.
        step = STEP_CORE_KM if hi <= R_CMB_KM + 1e-9 else STEP_MANTLE_KM
        # Smallest number of equal intervals with spacing <= step, so that the
        # first and last node sit exactly on the region boundaries.
        n_int = int(np.ceil((hi - lo) / step - 1e-9))
        r_all.append(np.linspace(lo, hi, n_int + 1))
        reg_all.append(np.full(n_int + 1, k))
    r_km = np.concatenate(r_all)
    reg = np.concatenate(reg_all)

    # Evaluate the polynomials region by region. Units: g/cm^3 and km/s.
    rho = np.empty_like(r_km)
    vp = np.empty_like(r_km)
    vs = np.empty_like(r_km)
    for k in np.unique(reg):
        m = reg == k
        rho[m], vp[m], vs[m] = P.evaluate(k, r_km[m])
    rho, vp, vs = rho * 1e3, vp * 1e3, vs * 1e3     # to kg/m^3 and m/s

    # Moduli from the wave speeds: mu = rho V_S^2 and kappa = rho V_P^2 - 4/3 mu.
    # These are the quantities the stable model keeps fixed.
    mu = rho * vs**2
    kappa = rho * vp**2 - 4.0 * mu / 3.0

    # Adjusted nodes: every node of the crust and the mantle. Of the two nodes
    # at the core-mantle boundary, the one that belongs to the outer core stays
    # PREM.
    adj = r_km >= R_BOT_KM - 1e-9
    adj &= ~((np.abs(r_km - R_BOT_KM) < 1e-9) & (reg == P.region_of(R_BOT_KM, "below")))
    return dict(r=r_km * 1e3, rho=rho, vp=vp, vs=vs, kappa=kappa, mu=mu, reg=reg, adj=adj)


def prem_shell_mass(r, reg):
    """Exact mass (kg) of PREM in every interval between consecutive nodes.

    The density inside an interval comes from the PREM polynomial of its
    region. A zero-width interval (a region boundary) has zero mass.
    """
    out = np.zeros(r.size - 1)
    for i in np.where(np.diff(r) > 0)[0]:
        r0, r1 = r[i], r[i + 1]
        # Map the Gauss points from [-1, 1] to [r0, r1].
        rq = 0.5 * (r1 - r0) * GL_X + 0.5 * (r1 + r0)
        rho_q = P.evaluate(reg[i + 1], rq / 1e3)[0] * 1e3
        out[i] = 0.5 * (r1 - r0) * np.sum(GL_W * 4.0 * np.pi * rq**2 * rho_q)
    return out


def linear_shell_weights(r):
    """Mass weights of a density that is linear between nodes.

    For the interval between nodes i and i+1, the mass of a density that is
    linear in r from rho_i to rho_{i+1} is a_i rho_i + b_i rho_{i+1}, with

        a_i = 4 pi int r^2 (r_{i+1} - r) / h dr,
        b_i = 4 pi int r^2 (r - r_i) / h dr,   h = r_{i+1} - r_i.

    Returns the arrays a and b, one entry per interval, in m^3. Zero-width
    intervals get zero weight.
    """
    a = np.zeros(r.size - 1)
    b = np.zeros(r.size - 1)
    for i in np.where(np.diff(r) > 0)[0]:
        r0, r1 = r[i], r[i + 1]
        h = r1 - r0
        rq = 0.5 * h * GL_X + 0.5 * (r1 + r0)
        jac = 0.5 * h * GL_W * 4.0 * np.pi * rq**2
        a[i] = np.sum(jac * (r1 - rq) / h)
        b[i] = np.sum(jac * (rq - r0) / h)
    return a, b


def shell_mass(r, rho, adj, prem_mass, a, b):
    """Mass of every interval: linear density where both nodes are adjusted,
    PREM's exact polynomial mass elsewhere."""
    lin = a * rho[:-1] + b * rho[1:]
    both = adj[:-1] & adj[1:]
    return np.where(both, lin, prem_mass)


def gravity(r, shells):
    """g(r) = G M(r) / r^2 at the nodes, from the interval masses.

    `shells` holds the mass of each interval, from the centre outward. Returns
    gravity (m/s^2) and the enclosed mass M(r) (kg) at every node.
    """
    m = np.concatenate([[0.0], np.cumsum(shells)])
    g = np.where(r > 0, G * m / np.where(r > 0, r, 1.0)**2, 0.0)
    return g, m


def stability_e(r, rho, kappa, g, mask):
    """e = 1 - Bullen parameter per interval whose two nodes are both in `mask`
    and that has finite width, on interval means as lovejx computes it.

    Returns the interval index and e. e > 0 is super-adiabatic (unstable where
    the material flows), e = 0 is the adiabatic gradient.
    """
    i = np.where(mask[:-1] & mask[1:] & (np.diff(r) > 0))[0]
    drdr = (rho[i + 1] - rho[i]) / (r[i + 1] - r[i])
    rm = 0.5 * (rho[i] + rho[i + 1])
    km = 0.5 * (kappa[i] + kappa[i + 1])
    gm = 0.5 * (g[i] + g[i + 1])
    return i, (drdr + rm**2 * gm / km) * km / (rm**2 * gm)


# --- the fit -----------------------------------------------------------------

def fit(model, tol=2e-5, max_iter=30):
    """Constrained least-squares fit of the adjusted density.

    See the module docstring for the rule. Returns the density at every node
    in kg/m^3, equal to PREM outside the adjusted range.
    """
    r, rho_p, kappa, reg, adj = model["r"], model["rho"], model["kappa"], model["reg"], model["adj"]
    idx = np.where(adj)[0]                        # model nodes that are unknowns
    n = idx.size
    assert np.all(np.diff(idx) == 1), "the adjusted nodes must be contiguous"

    # Interval masses: PREM's exact ones, and the weights of the linear density.
    prem_mass = prem_shell_mass(r, reg)
    a_w, b_w = linear_shell_weights(r)

    # Mass weight of each unknown node. Interval j of the unknowns is the
    # interval idx[j] of the model. Its mass is a x_j + b x_{j+1}.
    jj = idx[:-1]                                 # model interval index of each unknown interval
    w = np.zeros(n)
    w[:-1] += a_w[jj]
    w[1:] += b_w[jj]
    mass_prem = float(prem_mass[jj].sum())        # target: PREM's exact mass of the range
    wn = w / w.sum()                              # normalised objective weights

    # Interval and jump structure inside the adjusted range. Nodes are ordered
    # centre outward, so at a jump node j is the lower side and node j+1 the
    # upper side.
    dr = np.diff(r[idx])
    finite = np.where(dr > 0)[0]
    jumps = np.where(dr == 0)[0]

    # The outer-core node at the core-mantle boundary, just below the range,
    # is fixed at PREM.
    rho_core_top = rho_p[idx[0] - 1]

    rho = rho_p.copy()
    for it in range(max_iter):
        g, _ = gravity(r, shell_mass(r, rho, adj, prem_mass, a_w, b_w))
        gi, ki, xi = g[idx], kappa[idx], rho[idx]
        # Linearised stability: (x[j+1]-x[j])/dr + c_j (x[j]+x[j+1])/2 <= 0,
        # with c_j = rho_m g_m / kappa_m from the previous iterate. At a fixed
        # point this is exactly d rho/dr + rho_m^2 g_m / kappa_m <= 0.
        rm = 0.5 * (xi[finite] + xi[finite + 1])
        gm = 0.5 * (gi[finite] + gi[finite + 1])
        km = 0.5 * (ki[finite] + ki[finite + 1])
        c = rm * gm / km                          # 1/m
        A = np.zeros((finite.size, n))
        rows = np.arange(finite.size)
        A[rows, finite] = -1.0 / dr[finite] + 0.5 * c
        A[rows, finite + 1] = 1.0 / dr[finite] + 0.5 * c
        # Scale each row by dr so the constraint reads in kg/m^3, the units of
        # the unknowns. This keeps the rows well conditioned for the optimiser.
        A *= dr[finite][:, None]
        cons = [LinearConstraint(A, -np.inf, 0.0)]
        # Region boundaries inside the range: below (node j) >= above (node j+1).
        if jumps.size:
            B = np.zeros((jumps.size, n))
            B[np.arange(jumps.size), jumps] = 1.0
            B[np.arange(jumps.size), jumps + 1] = -1.0
            cons.append(LinearConstraint(B, 0.0, np.inf))
        # Core-mantle boundary: the lowest unknown (base of the mantle) must
        # not be denser than the core below it.
        lo_cmb = np.zeros((1, n)); lo_cmb[0, 0] = 1.0
        cons.append(LinearConstraint(lo_cmb, -np.inf, rho_core_top))
        # Mass of the range equal to PREM's.
        cons.append(LinearConstraint(w[None, :] / mass_prem, 1.0, 1.0))

        target = rho_p[idx]
        res = minimize(lambda x: float(wn @ (x - target)**2),
                       xi, jac=lambda x: 2.0 * wn * (x - target),
                       hess=lambda x: np.diag(2.0 * wn),
                       constraints=cons, method="trust-constr",
                       options=dict(gtol=1e-12, xtol=1e-12, maxiter=20000, verbose=0))
        new = rho.copy()
        new[idx] = res.x
        change = np.max(np.abs(new - rho))
        rho = new
        print(f"iteration {it}: max change {change:.2e} kg/m^3, optimiser status {res.status}")
        if change < tol:
            break
    return rho


# --- output --------------------------------------------------------------------

def write_csv(path, r, reg, rho, kappa, mu, g):
    """Write one model as CSV, centre outward, with V_P and V_S from the moduli."""
    vs = np.sqrt(mu / rho)
    vp = np.sqrt((kappa + 4.0 * mu / 3.0) / rho)
    depth = P.A_KM - r / 1e3
    # The depth of the centre node is written as exactly 6371 km, and the
    # surface as exactly 0, so the round-off of A_KM - r/1e3 does not show.
    depth = np.where(np.abs(depth) < 1e-9, 0.0, depth)
    data = np.column_stack([r, depth, reg, rho, vp, vs, kappa, mu, g])
    np.savetxt(path, data, delimiter=",", header=CSV_HEADER, comments="",
               fmt=["%.6f", "%.6f", "%d"] + ["%.15g"] * 6)


def main():
    model = prem_nodes()
    rho = fit(model)
    r, rho_p, kappa, mu, reg, adj = (model[k] for k in ("r", "rho", "kappa", "mu", "reg", "adj"))
    prem_mass = prem_shell_mass(r, reg)
    a_w, b_w = linear_shell_weights(r)
    g, m = gravity(r, shell_mass(r, rho, adj, prem_mass, a_w, b_w))
    g_p, m_p = gravity(r, prem_mass)

    _, e_new = stability_e(r, rho, kappa, g, adj)
    _, e_old = stability_e(r, rho_p, kappa, g_p, adj)
    d = rho - rho_p
    print(f"nodes: {r.size}, adjusted: {adj.sum()}")
    print(f"density change in the crust and mantle: {d[adj].min():+.2f} to {d[adj].max():+.2f} kg/m^3")
    lm = adj & (r / 1e3 < P.A_KM - 670.0 - 1e-9)
    print(f"density change below 670 km: {d[lm].min():+.2f} to {d[lm].max():+.2f} kg/m^3")
    print(f"e in the crust and mantle: PREM max {e_old.max():+.3f}, new max {e_new.max():+.2e}")
    print(f"total mass: PREM {m_p[-1]:.6e}, new {m[-1]:.6e}, relative change {(m[-1]-m_p[-1])/m_p[-1]:.1e}")
    print(f"surface gravity: PREM {g_p[-1]:.6f}, new {g[-1]:.6f} m/s^2")
    print(f"largest gravity change: {np.max(np.abs(g - g_p)):.2e} m/s^2")
    # Density jumps at the region boundaries in the range, and at the CMB.
    for k, (lo, hi, name, *_rest) in enumerate(P.REGIONS):
        if lo < R_BOT_KM - 1e-9:
            continue
        i = np.where(np.abs(r / 1e3 - lo) < 1e-6)[0]
        lo_i, hi_i = i[0], i[1]
        print(f"jump at {P.A_KM - lo:6.1f} km: below minus above "
              f"{rho[lo_i]-rho[hi_i]:+8.2f} kg/m^3 (PREM {rho_p[lo_i]-rho_p[hi_i]:+8.2f})")

    write_csv(HERE / "PREM.csv", r, reg, rho_p, kappa, mu, g_p)
    write_csv(HERE / "PREM-stable.csv", r, reg, rho, kappa, mu, g)
    print("wrote PREM.csv and PREM-stable.csv")


if __name__ == "__main__":
    main()
