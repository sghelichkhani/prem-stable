# Rejected rules

These two scripts test rules for the stable density that we did not keep.
`DECISIONS.md` records their results and the reasons.

- `build_joint.py` keeps the mass of the whole 60 to 670 km range and PREM's
  density jumps at 220 and 400 km.
- `build_unitA.py` keeps the mass of 60 to 400 km as one unit and PREM's jump
  at 220 km.

Both scripts import helper functions (`prem_si`, `shell_mass`, `mass_below`,
`DR_KM`) from an earlier version of `build_prem_stable.py`. The current version
does not have these functions, so the scripts do not run as they are. They are
kept as a record of the method.
