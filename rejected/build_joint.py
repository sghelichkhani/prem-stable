"""Joint rule: jumps at 220 and 400 km keep PREM's size, density stays continuous
at 600 km, adiabatic gradient in the unstable layers (60-220, 220-400, 600-670 km),
PREM's gradient in the stable transition zone (400-600 km), and one free constant
chosen so that the mass between 60 and 670 km is PREM's. The lower mantle and core
are PREM, so gravity at 670 km is PREM's. Integrates upward from 670 km."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import prem1981 as P
from build_prem_stable import prem_si, shell_mass, mass_below, G, DR_KM

# (bottom radius km, top radius km, region index for PREM kappa/rho, mode)
SEG = [(5701.0, 5771.0, P.region_of(5750.0), "adiabat"),
       (5771.0, 5971.0, P.region_of(5900.0), "prem-gradient"),
       (5971.0, 6151.0, P.region_of(6000.0), "adiabat"),
       (6151.0, 6311.0, P.region_of(6200.0), "adiabat")]   # LID and LVZ share one polynomial
# PREM jumps (kg/m^3, value below minus value above) at 400 and 220 km.
def prem_rho(k, r): return prem_si(k, np.atleast_1d(r))[0][0]
def lvz(r):  # LID/LVZ polynomial, valid 6151-6346.6
    x = r / P.A_KM; return (2.6910 + 0.6924 * x) * 1e3
JUMP = {5971.0: prem_rho(P.region_of(5971.0, "below"), 5971.0) - prem_rho(P.region_of(5971.0, "above"), 5971.0),
        6151.0: prem_rho(P.region_of(6151.0, "below"), 6151.0) - lvz(6151.0)}

def profile(rho_670_top):
    """Integrate upward from just above 670 km. Returns list of (r_km, rho, M)."""
    m = mass_below(5701.0)
    rho = rho_670_top
    out = []
    for lo, hi, k, mode in SEG:
        n = int(round((hi - lo) / DR_KM)) + 1
        r_km = np.linspace(lo, hi, n); r_m = r_km * 1e3
        if lo in JUMP:                       # step across the discontinuity at the bottom of this segment
            rho = rho - JUMP[lo]
        if k == P.region_of(6200.0):
            x = r_km / P.A_KM; rp = (2.6910 + 0.6924 * x) * 1e3
            vp = (4.1875 + 3.9382 * x) * 1e3; vs = (2.1519 + 2.3481 * x) * 1e3
            kap = rp * vp**2 - 4 / 3 * rp * vs**2
        else:
            rp, kap, _ = prem_si(k, r_km)
        if mode == "prem-gradient":
            rh = rp + (rho - rp[0])
            f = 4 * np.pi * r_m**2 * rh
            M = m + np.concatenate([[0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(r_m))])
        else:
            kf = lambda r: np.interp(r, r_m, kap)
            sol = solve_ivp(lambda r, y: [-y[0]**2 * G * y[1] / r**2 / kf(r), 4 * np.pi * r**2 * y[0]],
                            (r_m[0], r_m[-1]), [rho, m], t_eval=r_m, rtol=1e-12, atol=[1e-9, 1.0], method="DOP853")
            rh, M = sol.y
        out.append((r_km, rh, rp, M))
        rho, m = rh[-1], M[-1]
    return out

target = mass_below(6311.0)
f = lambda x: profile(x)[-1][3][-1] - target
x0 = prem_rho(P.region_of(5701.0, "above"), 5701.0)
x = brentq(f, x0 - 300, x0 + 300, xtol=1e-9)
prof = profile(x)
names = ["600-670", "400-600", "220-400", "60-220"]
for nm, (r, rh, rp, M) in zip(names, prof):
    d = rh - rp
    print(f"{nm:8s} change {d.min():+7.2f} .. {d.max():+7.2f}  bottom {rh[0]:.2f} (PREM {rp[0]:.2f})  top {rh[-1]:.2f} (PREM {rp[-1]:.2f})")
print("jump at 670 km (lower mantle top minus new):", prem_rho(P.region_of(5701.0, 'below'), 5701.0) - prof[0][1][0])
print("jump at 60 km (lid minus new, negative = inverted):", lvz(6311.0) - prof[-1][1][-1])
print("mass residual", prof[-1][3][-1] - target)
np.save("joint_profile.npy", np.array([np.concatenate([p[0] for p in prof]), np.concatenate([p[1] for p in prof]), np.concatenate([p[2] for p in prof])]))
