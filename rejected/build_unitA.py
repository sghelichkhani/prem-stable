"""Unit 60-400 km: adiabat in 60-220 and 220-400, PREM's jump at 220 km kept, mass of
60-400 km equal to PREM's. Transition zone and below unchanged. Also reports what the
600-670 km layer does if made adiabatic from PREM's value at 600 km (continuous with
the stable transition zone above)."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import prem1981 as P
from build_prem_stable import prem_si, mass_below, G, DR_KM

def lvz_props(r_km):
    x = r_km / P.A_KM; rho = (2.6910 + 0.6924 * x) * 1e3
    vp = (4.1875 + 3.9382 * x) * 1e3; vs = (2.1519 + 2.3481 * x) * 1e3
    return rho, rho * vp**2 - 4 / 3 * rho * vs**2

def adiabat(r_km, kap, rho0, m0, up=True):
    r_m = r_km * 1e3
    kf = lambda r: np.interp(r, r_m, kap)
    span = (r_m[0], r_m[-1]) if up else (r_m[-1], r_m[0])
    te = r_m if up else r_m[::-1]
    sol = solve_ivp(lambda r, y: [-y[0]**2 * G * y[1] / r**2 / kf(r), 4 * np.pi * r**2 * y[0]],
                    span, [rho0, m0], t_eval=te, rtol=1e-12, atol=[1e-9, 1.0], method="DOP853")
    return (sol.y[0], sol.y[1]) if up else (sol.y[0][::-1], sol.y[1][::-1])

kU = P.region_of(6000.0)
r1 = np.linspace(5971.0, 6151.0, 361); rp1, kap1, _ = prem_si(kU, r1)
r2 = np.linspace(6151.0, 6311.0, 321); rp2, kap2 = lvz_props(r2)
jump220 = rp1[-1] - rp2[0]                     # PREM: below minus above, positive
m400 = mass_below(5971.0); target = mass_below(6311.0)

def run(rho400):
    a1, M1 = adiabat(r1, kap1, rho400, m400)
    a2, M2 = adiabat(r2, kap2, a1[-1] - jump220, M1[-1])
    return a1, a2, M2[-1]
x = brentq(lambda v: run(v)[2] - target, rp1[0] - 300, rp1[0] + 300, xtol=1e-9)
a1, a2, _ = run(x)
tz_top = prem_si(P.region_of(5971.0, "below"), np.array([5971.0]))[0][0]
print(f"220-400: change {np.min(a1-rp1):+.2f}..{np.max(a1-rp1):+.2f}; bottom {a1[0]:.2f} (PREM {rp1[0]:.2f})")
print(f"60-220 : change {np.min(a2-rp2):+.2f}..{np.max(a2-rp2):+.2f}; top {a2[-1]:.2f} (PREM {rp2[-1]:.2f})")
print(f"jump at 400 km, TZ top minus new: {tz_top - a1[0]:+.2f} (PREM {tz_top - rp1[0]:+.2f})")
print(f"jump at 60 km, lid minus new (positive = inverted under the lid): {rp2[-1] - a2[-1]:+.2f}")

# 600-670 anchored continuously at 600 km, integrated downward
kT = P.region_of(5750.0)
r3 = np.linspace(5701.0, 5771.0, 141); rp3, kap3, _ = prem_si(kT, r3)
m670 = mass_below(5701.0)
# downward integration needs M at 600; use PREM mass below 600 (layer mass changes slightly; iterate twice)
Mtop = mass_below(5771.0)
a3 = rp3.copy()
for _ in range(3):
    a3, M3 = adiabat(r3, kap3, rp3[-1], Mtop, up=False)
    f = 4*np.pi*(r3*1e3)**2*a3; dm = np.sum(0.5*(f[1:]+f[:-1])*np.diff(r3*1e3))
    Mtop = m670 + dm
f = 4*np.pi*(r3*1e3)**2*rp3; dm_prem = np.sum(0.5*(f[1:]+f[:-1])*np.diff(r3*1e3))
lm_top = prem_si(P.region_of(5701.0, "below"), np.array([5701.0]))[0][0]
print(f"600-670: change {np.min(a3-rp3):+.2f}..{np.max(a3-rp3):+.2f}; mass change {(Mtop-m670-dm_prem):.3e} kg = {(Mtop-m670-dm_prem)/mass_below(6371.0):.1e} of Earth; jump at 670 {lm_top-a3[0]:+.1f} (PREM {lm_top-rp3[0]:+.1f})")
