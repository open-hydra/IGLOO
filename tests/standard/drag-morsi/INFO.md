# INFO — standard/drag-morsi and standard/drag-chang (drag laws along a slip relaxation, e2e)

## Reference
- **Chang:** the Clift-Gauvin correlation as Shimada et al. print it, `[Shim06]` eq. 18:
  `C_D = 24/Re (1 + 0.15 Re^0.687 + 0.0175 Re / (1 + 4.25e4 Re^-1.16))` — the drag law of every hydra JPL case.
- **Morsi-Alexander:** `C_D = K1/Re + K2/Re² + K3` on eight Re bands, the constants of Morsi and Alexander's
  appendix (`[MA72]`, printed p. 207); `check.py` types them as (a1, a2, a3) = (K3, K1, K2), and IGLOO's table
  (`Lib_Drag.f90`) equals the appendix in all eight bands.
- **Tags:** `[Shim06]`, `[MA72]` (see [../../REFERENCES.md](../../REFERENCES.md)).

## Fixture
drag-stokes's box and uniform gas (`INPUT/solfile.tec` is drag-stokes's: ρ_g = 1.2 kg/m³, u_g = 10 m/s along +x,
μ = 1.8×10⁻⁵ Pa·s, T_g = 600 K) with 100 µm particles: `INPUT/bc.txt` from
`tools/make_box_case.py INPUT --rp 5.0e-5` (κ_v = 0.1: injected along the gas at 1 m/s, T_p0 = T_g; ρ_p = 2950 kg/m³
from the common properties table). The slip Reynolds number falls from 60 to about 22 across the 0.15 m box
(t ≈ 0.033 s): the Morsi-Alexander path stays in the band 10 < Re ≤ 100. `drag-chang` symlinks `drag-morsi`'s
`INPUT/` and `check.py`; the two `input.ini` differ in `drag` only. `out-time = on`.

## What it verifies
`check.py` types the case's law (chosen by the case directory's name) and integrates the slip ODE
`du/dt = (3/4) (ρ_g / (ρ_p d)) C_D(Re) |u_g − u| (u_g − u)` (uniform gas, no transverse slip: one scalar ODE)
with RK4 at 2×10⁵ steps from u = 1 m/s, then a partial step to each printed time.

| id | check |
|---|---|
| D1 | every row: \|u − u_or(t)\| ≤ 5×10⁻⁷ (the F12.6 print) + \|du/dt\| 5×10⁻⁹ t (the ES16.8 time) + 10⁻⁹ |
| D2 | every row: v = w = 0 and T = 600 K to the print, d = 100 µm; each parcel's first row at x = 0, u = 1 m/s, t = 0 |
| D3 | 25 parcels with at least 50 rows each, 25 exits at x = 0.15 |

The check prints the Re range it verified.

## Discrimination
At Re = 60, Morsi-Alexander's C_D is 2.87 % below Schiller-Naumann's and Chang's 0.081 % above it; along the
box the u of a Schiller-Naumann run departs from Morsi-Alexander's by up to 5.3×10⁻² m/s and from Chang's by up
to 1.3×10⁻³ m/s, against a tolerance of about 5×10⁻⁷ m/s.
