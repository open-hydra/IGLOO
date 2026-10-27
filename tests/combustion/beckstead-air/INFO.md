# INFO — combustion/beckstead-air (burnout times of aluminium in air against Beckstead's Eq. 7, e2e)

## Reference paper
- **Title:** Correlating Aluminum Burning Times
- **Author(s):** M. W. Beckstead
- **Year:** 2005
- **Venue:** *Combustion, Explosions, and Shock Waves*, 41(5), pp. 533–546
- **DOI:** 10.1007/s10573-005-0067-2
- **Tag:** `[Beck05]` (see [../../REFERENCES.md](../../REFERENCES.md))
- **Reproduced:** Eq. 7 (p. 543), `t_b = a D^n/(X_eff p^0.1 T0^0.2)`, `a = 0.00735` for `n = 1.8` (D in µm, t in
  ms, p in atm, T0 in K, `X_eff = C_O2 + 0.6 C_H2O + 0.22 C_CO2`), at a stated condition; **Fig. 1**, the burning
  times measured by 12 sources, is digitized and drawn beside it (`reference/`, not gated).

## The condition
Fig. 1 pools many conditions; the case states one, the commonest of the laboratory sources: aluminium in **air at
1 atm and 300 K** (laser- or flame-ignited single particles in room air). `X_eff = C_O2 = 0.21`, a **mole fraction**
(IGLOO's `X-eff` lies in (0, 1]; with mole fractions Eq. 7 lies inside the Fig. 1 cloud, in percent 100× below it),
`p = 1 atm`, `T0 = 300 K` (Eq. 7's "initial temperature": the gas and the particles are both at 300 K, so either
reading holds). Then

    K_eff = X_eff p^0.1 T0^0.2 / a = 89.403847 µm^1.8/ms = 1.416955485e-6 m^1.8/s   (1 µm^1.8/ms = 10^-7.8 m^1.8/s)
    d^n(t) = D0^n − K_eff t,   t_b = D0^n/K_eff = 0.7057 ms (10 µm) … 2809.6 ms (1000 µm)

and IGLOO's `K-burn = K_eff/X_eff = 6.74740707e-6`. Eq. 7 at this condition lies 2.81× (10 µm), 1.82× (100 µm) and
1.17× (1000 µm) above Fig. 1's drawn fit, inside the scatter.

## Case construction
The uniform-gas box of `tools/make_box_case.py` (60×5×5 cells over 0.15 × 0.05 × 0.05 m), 25 diameters log-spaced
10 … 1000 µm, one per inlet cell (the reader's order: n outer, m inner):

    python3 tests/tools/make_box_case.py INPUT --kv 1.0 --kt 1.0 --tg 300 --ug 0.04 --rp-sweep 5.00000e-06,6.05764e-06,7.33900e-06,8.89140e-06,1.07722e-05,1.30508e-05,1.58114e-05,1.91559e-05,2.32079e-05,2.81171e-05,3.40646e-05,4.12702e-05,5.00000e-05,6.05764e-05,7.33900e-05,8.89140e-05,1.07722e-04,1.30508e-04,1.58114e-04,1.91559e-04,2.32079e-04,2.81171e-04,3.40646e-04,4.12702e-04,5.00000e-04

- `phase.txt`: `A 1 combustion=Beckstead K-burn=6.74740707e-06 n-burn=1.8 X-eff=0.21 T-ign=250`. No heat release
  (`beta-part`, `q-comb` at their 0 defaults), no oxide cap or residue (`xi-cap` 0, no `residue` key): ODE model 5.
  **`T-ign = 250 K` is a switch, not a physical claim**: Eq. 7's rate does not depend on the particle temperature
  above ignition, so the gate sits below the injection temperature and every particle burns from injection.
- `kv = kt = 1`: zero slip and zero thermal slip, the d^n law is the only dynamics; `u_g = 0.04 m/s` is a clock (the
  longest burn, 2.81 s at 1 mm, ends at x = 0.112 m inside the 0.15 m box; the 10 µm burn ends at x = 28 µm).
- **K-burn carries `K_eff/X_eff` to nine digits.** The late trajectory rows of the 121–1000 µm parcels resolve the
  burn rate to 1e-8: K-burn written to seven digits (6.747407e-6, 1.0e-8 below `K_eff/X_eff`) puts 6 of the 170 rows
  over their print budget (3.23×); nine digits (7.5e-12) leave every row inside it.
- `properties.dat` → `tests/common/properties.dat` (ρ_p = 2950 kg/m³, c_p = 1250 J/kg/K); the oracle reads ρ_p only
  for `d_floor` (it cancels in d^n(t)).
- `input.ini`: `out-time = on` (the time on every row and exit record), `print-dcell = 1`, `out-file = S`;
  `[IGLOO-ODE] H-sdirk4`, `relative-tol = 1e-10`, `absolute-tol = 1e-20` — the absolute tolerance sits below the
  smallest mass it controls (a 10 µm particle weighs 1.5e-12 kg); `time-end` absent (off).

**The burnout.** A parcel ends at the accepted step whose end mass is at or below 1e-15 kg, the solver's absolute
floor: `d_floor = (6 m_b/(π ρ_p))^{1/3} = 0.8651 µm`, `t_floor = (D0^n − d_floor^n)/K_eff` (0.98779 `t_b` at 10 µm,
0.99981 `t_b` at 100 µm). The exit record (bc 0) carries that step's end state, so its `t` and `d_p` are one state;
the overshoot `t_exit − t_floor` is the part of that last step past the floor: +1.8e-8 … +1.9e-6 s, +2.5e-7 …
+1.37e-3 of `t_floor` (the largest at 10 µm). The 13 sizes up to 100 µm (every D below 120.7 µm) burn out before
the first cell face (x = 2.5 mm, reached at t = 0.0625 s): their trajectory holds the injection row only, and the
burn time comes from the exit record. The longest solver call of the 10 µm parcel takes 62 computed steps (the
least `max-steps-ode` it runs with); the 1 mm parcel writes 45 rows and its longest call takes 331. A run takes
about 0.6 s; the log carries no revert, stuck, no-progress, outer-maxIter or solver-error line.

## What it verifies (`check.py`, Eq. 7 coded from the paper)

| id | gate | tolerance | measured |
|---|---|---|---|
| B0 | `phase.txt`: K-burn × X-eff = `K_eff` of Eq. 7; n-burn = 1.8, X-eff = 0.21 | 1e-6 relative | 7.5e-12 |
| B1 | 25 exit records, one per inlet cell (the exit (Y, Z) names the cell, its `bc.txt` record D0 = 2 rp), every `bc` = 0, `d_p ≤ d_floor`(1 + 1e-5), `t > 0` — the burnout witness no other bc-0 ending meets; each parcel's first row at d = D0, every row in its cell's (Y, Z) | half an E13.6 unit; F12.6 | `d_p` ≤ 0.8641 µm; 0.40 of the unit |
| B2 | `\|t_exit − (D0^n − d_p^n)/K_eff\| ≤ EPS t_exit`; the floor, one-sided: `t_exit ≥ t_floor (1 − EPS)`; the overshoot printed, not gated | EPS = 1e-5 (5 × measured, floored) | 3.9e-8; 25/25 |
| B3 | every row: `\|d^n − (D0^n − K_eff t)\|` within the row's print budget, `n d^{n−1}` × half an E13.6 unit of d + `K_eff` × half an ES16.8 unit of t | 1 | 0.9967 (170 rows) |
| B4 | `T_p = 300 K` and `u_p = (u_g, 0, 0)` on every row and exit record | half an F12.6 unit | exact |

## RED

| mutation | what fails |
|---|---|
| K-burn × 1.01 in `phase.txt` | B0 1.0e-2; B2 1.0e-2, the floor on 25 of 25 parcels; B3 (9.8e-3 of D0^n) |
| the X_eff exponent 1.1 in the solver (`aXeff`) | B2 0.14; B3 |
| `T-ign = 350` (above the injection temperature) | B1: nothing burns, 25 exits at the outlet (bc 400) |
| `beta-part=0.3 q-comb=1000000` (a heat release) | B4: T_p 882 K above 300 K |
| `kv 0.99` in `bc.txt` (a slip) | B4: 4.0e-4 m/s |
| the oracle without `d_floor` (t_floor = t_b) | B2's floor check: t_exit at 10 µm 1.09 % short of t_b |
| K-burn to seven digits (6.747407e-6) | B3: 6 of 170 rows over their budget, 3.23× |

## Comparison plot
`verify.py` writes `OUTPUT/beckstead-air.svg` (non-gating; no-op without matplotlib), on Fig. 1's axes (log-log,
`t_b` [ms] against D [µm], 10–1000 µm × 0.1–10⁴ ms): the 183 markers digitized from Fig. 1 as one grey glyph (the
source labels of the stacked markers are not recovered, `reference/PROVENANCE.md`), Fig. 1's fit **as drawn**
(`t = 0.00258 D^1.989`, dashed; its label prints the rounded `0.003 D^1.99`, which lies 16.5–16.8 % above the drawn
line), Eq. 7 at air, 1 atm, 300 K (solid), and IGLOO's 25 burnout times (the exit records' `t`, markers) at their
injection diameters. The parcels are taken by ascending ID; the figure does not depend on the order of the records
(`vv-order`). Fig. 9's reduced data are not drawn: the figure's own upper line is not Eq. 7.

## Reference data
- `reference/beckstead_fig1.csv` — the 183 markers of Fig. 1 (`D_um, t_ms, source, kind`), digitized from the page's
  300-ppi JPEG; 90 carry an NCC-confirmed source label. Method, calibration and the label rule in
  `reference/PROVENANCE.md`. Plotted, not gated.
