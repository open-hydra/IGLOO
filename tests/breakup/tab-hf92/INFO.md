# INFO — breakup/tab-hf92 (TAB products of Hsiang & Faeth's water series; e2e, the correlations drawn, not gated)

## Reference papers
- `[HF92]` L.-P. Hsiang, G. M. Faeth, "Near-limit drop deformation and secondary breakup", *Int. J. Multiphase
  Flow* 18(5) (1992) 635–652: the SMD correlation eq. [14] (p. 649), **SMD/d0 = 6.2 (ρ_L/ρ_G)^1/4 [μ_L/(ρ_L d0
  u0)]^1/2**, the standard deviation of its coefficient 20 %, for Oh < 0.1 and We < 10³ (bag, transition and shear
  breakup); Simmons' root-normal distribution with MMD/SMD = 1.2 after bag and multimode breakup (pp. 645–646); the
  breakup time t_b = 5 t*, t* = d0 (ρ_L/ρ_G)^1/2/u0 (eqs. [1]–[2], p. 642). We = ρ_G d0 u0²/σ, Oh = μ_L/(ρ_L d0
  σ)^1/2, Re = ρ_G d0 u0/μ_G, u0 the initial relative velocity.
- `[HF93]` L.-P. Hsiang, G. M. Faeth, "Drop properties after secondary breakup", *Int. J. Multiphase Flow* 19(5)
  (1993) 721–735: the water series this case reproduces (d0 = 1000 µm, We 15–600; Table 1: ρ_L 997, μ_L 8.94e-4,
  σ 70.8e-3, air ρ_G 1.18, μ_G 18.5e-6 at 297 K); the regimes from We = 13 (bag), 35 (multimode) and 80 (shear)
  (p. 724); the root-normal distribution with MMD/SMD = 1.2 after shear breakup with the core drop removed (§3.1,
  Fig. 1); the size/velocity correlation eq. [13] (p. 732), **u0/u_b − 1 = 2.7 [(ρ_G/ρ_L)^1/2 d0/d]^2/3**, u_b the
  relative velocity of a drop of diameter d at the end of breakup, Fig. 6's markers at (ρ_G/ρ_L)^1/2 d0/d ≈
  0.067–0.85 (its correlation line 0.058–1.07).
- `[ORA87]` the TAB model, as `tab-e2e` runs it (`breakup = TAB`, `method = 2`, `n = 3.5`).

## The case
One run, 25 classes: each inlet cell injects water drops of d0 = 1000 µm at its own slip, one diameter-based Weber
number per cell, so a class is an inlet cell and the series is the papers' (d0 fixed, u0 varied).

| quantity | value | against [HF92] / [HF93] |
|---|---|---|
| drop | water, d0 = 1000 µm in every cell; ρ_L 997, c_p 4182; `[IGLOO-Properties] sigma = 0.0708`, `mu = 8.94e-4` | [HF93] Table 1, water row |
| gas | air, ρ_G 1.18, μ_G 1.85e-5, T 297 K, uniform U = 200 m/s | [HF93] Table 1 |
| We | 25 classes, geometric 15.00 … 600.00, diameter-based | [HF93] We 15–600; bag 6 classes, multimode 5, shear 14 |
| slip | u0 = (We σ/(ρ_G d0))^1/2 = 30.0 … 189.7 m/s: the drop is injected at 170 … 10 m/s into the 200 m/s gas (kv 0.850 … 0.051) | the step change of relative velocity of a shock |
| Re, Oh, ρ_L/ρ_G | 1914 … 12 102; 0.00336; 844.9 | Re 1060–15 080, Oh 0.0025–0.039, ρ_L/ρ_G 579–985 |
| t*, t_b = 5 t* | 0.153 … 0.969 ms; 0.77 … 4.85 ms | [HF93] p. 726 |
| box | 1.2 × 0.05 × 0.05 m, 480 × 5 × 5 cells; faces 3–6 reflecting (code 300), outlet face 2 | every parcel is in the box at t_b |
| parcels | `[IGLOO-BC] ds = 0.2` (cm) in the 1-cm inlet cells: 714, 25 to 33 per class | |
| models | `breakup = TAB`, `method = 2`, `n = 3.5`; `drag = Schiller-Naumann`, `heat = Ranz-Marshall` (kt = 1: no thermal slip); `gas-order = 1`; H-sdirk4 at 1e-9 | `[ORA87]`, `[SN33]`, `[RM52]` |
| output | `print-dcell = 4` (a trajectory row per cm), `out-time = on` (the time column), `out-scatter = off`, `out-file = S`: `trajectories-A.dat` 11.7 MB | |

`INPUT/` (shared: `etab-hf92/INPUT` links it), from this directory:

    python3 ../../tools/make_pe_case.py INPUT --we-convention dia \
      --we-sweep 15.00,17.49,20.40,23.79,27.74,32.35,37.72,43.99,51.30,59.82,69.76,81.35,94.87,110.63,129.01,150.45,175.44,204.59,238.58,278.22,324.45,378.35,441.21,514.52,600.00 \
      --d0 1e-3 --u-gas 200 --sigma 0.0708 --t-gas 297 --kt 1.0 --rhog 1.18 --mug 1.85e-5 --lx 1.2 --nx 480 \
      --side-code 300 --title "IGLOO TAB/ETAB Hsiang-Faeth water series box"
    cp ../tab-e2e/INPUT/phase.txt INPUT/
    sed 's/ 1000\.000000 / 997.000000 /' ../tab-e2e/INPUT/properties.dat > INPUT/properties.dat

writes `solfile.tec` and `bc.txt` (`--d0` sets every cell's kv from its Weber number, `--rhog`/`--mug` the gas) and
the two property files (tab-e2e's, the density column at 997).

**The lateral faces reflect.** `--side-code 300` codes faces 3–6 as symmetry planes (specular reflection, the
grazing slide below |v_n| < 0.02 |v|, `obj_bc.f90`); a parcel leaves through the outlet only. The injection packing
puts the edge cells' parcels on the lateral faces, so with faces 3–6 coded 400 a parcel that drifts sideways leaves
the box: no TAB parcel does, and this case's output is the same under 400 and 300 (`.dat` files equal as sorted
multisets, `source.tec` equal at 1e-12). In a uniform gas parallel to the faces a reflection flips the sign of the
normal component of the slip and keeps its magnitude, which is what the drag and the breakup models read; the
grazing slide is the exception (`etab-hf92`, whose parcels reach the faces, measures its effect). P2's "no lateral
exit" therefore checks the input and the reflection, not the physics.

## What `check.py` gates (the preconditions of the figure)
A parcel's class is the inlet cell holding the (Y, Z) of its injection row (t = 0), cross-checked by that row's
|u_p| = kv·U (F12.6; a parcel on a cell edge belongs to the cell whose kv it carries); every class is injected at d0,
and the products move across cells. Inputs come from the case files (kv and d0 from `bc.txt`, the gas and the box
from `solfile.tec`, σ and μ_L from `input.ini`, ρ_L from `properties.dat`).

| gate | condition | measured |
|---|---|---|
| P0 population | every class ≥ 15 exit records, one per parcel, every parcel in one class | 714 parcels, 25 to 33 per class, 0 without a class |
| P1 every drop broke | d_exit < 0.999 d0 for every parcel (every class is above the onset We = 13) | max d_exit/d0 0.553539 |
| P2 the cascade ended inside the box | every exit at the outlet (X ≥ Lx − 1e-6 m), d constant over x ≥ Lx/2, the last row's radius Weber number ρ_G \|u_g − u_p\|² (d/2)/σ below 6 | 0 lateral exits, 0 parcels changing d over x ≥ 0.6 m, max We_r 3.828 |
| P3 the figure's size/velocity panel is defined | every parcel's rows bracket t_b = 5 t* of its class | 0 parcels leaving before t_b, min t_last/t_b 1.421 |

RED (scratch copies, each exit 1):

| mutation | what fails |
|---|---|
| every inlet record's kv 0.97 (`tools/set_kv.py`; slip 6 m/s, We 0.6) | P1 max d_exit/d0 1.000000; P3 658 of 658 leave before t_b; P0 56 parcels without a class (one kv for every class leaves the edge parcels ambiguous) |
| box `--lx 0.05 --nx 20` | P2 126 parcels change d over x ≥ 0.025 m, max We_r 34.481; P1, P3 |
| box `--lx 0.3 --nx 120` | P3 334 of 714 leave before t_b, min t_last/t_b 0.362; P2 52 parcels change d over x ≥ 0.15 m |
| `ds = 10` (one parcel per cell) | P0 1 parcel per class; P1–P3 pass |
| faces 3–6 coded 400 (`--side-code 400`) | nothing here (no TAB parcel reaches a lateral face); the same oracle on `etab-hf92`: P2 46 of 714 leave before the outlet |

## Measured (the reference run)
SMD = Σ ṁ / Σ (ṁ/d) and the cumulative volume fraction V(d) = Σ_{d_i ≤ d} ṁ_i / Σ ṁ_i over a class's exit records
(each parcel keeps its mass flow: the exit ṁ is the injected one); MMD the d at V = 0.5. "events" counts the d
changes between trajectory rows: a lower bound, the rows are 1 cm apart and none is written at an event.

| class | We | u0 [m/s] | Re | n | SMD/d0 | [HF92] eq. [14] | MMD/SMD | events | x(5 t*) [m] | We_r, last row |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 15.00 | 30.000 | 1914 | 27 | 0.3795 | 0.1827 | 1.036 | 1 | 0.829–0.837 | 3.220 |
| 2 | 17.49 | 32.394 | 2066 | 25 | 0.3379 | 0.1759 | 1.041 | 1 | 0.758–0.771 | 3.052 |
| 3 | 20.40 | 34.986 | 2232 | 25 | 0.3282 | 0.1692 | 1.081 | 1 | 0.691–0.704 | 3.828 |
| 4 | 23.79 | 37.781 | 2410 | 25 | 0.3149 | 0.1628 | 1.057 | 1 | 0.631–0.641 | 3.407 |
| 5 | 27.74 | 40.797 | 2602 | 25 | 0.2502 | 0.1567 | 1.117 | 1 | 0.576–0.605 | 2.626 |
| 6 | 32.35 | 44.057 | 2810 | 33 | 0.2467 | 0.1508 | 1.019 | 1 | 0.522–0.534 | 3.166 |
| 7 | 37.72 | 47.573 | 3034 | 30 | 0.2283 | 0.1451 | 1.058 | 1–2 | 0.474–0.487 | 2.963 |
| 8 | 43.99 | 51.375 | 3277 | 30 | 0.2017 | 0.1396 | 1.031 | 1–2 | 0.431–0.444 | 2.464 |
| 9 | 51.30 | 55.480 | 3539 | 30 | 0.1659 | 0.1344 | 1.103 | 1–2 | 0.391–0.411 | 1.831 |
| 10 | 59.82 | 59.910 | 3821 | 30 | 0.1491 | 0.1293 | 1.043 | 1–2 | 0.353–0.372 | 1.525 |
| 11 | 69.76 | 64.696 | 4127 | 33 | 0.1200 | 0.1244 | 1.084 | 1–2 | 0.320–0.339 | 1.125 |
| 12 | 81.35 | 69.864 | 4456 | 30 | 0.1018 | 0.1198 | 1.106 | 1–2 | 0.289–0.317 | 0.734 |
| 13 | 94.87 | 75.447 | 4812 | 30 | 0.0785 | 0.1152 | 1.206 | 1–2 | 0.262–0.312 | 0.447 |
| 14 | 110.63 | 81.473 | 5197 | 30 | 0.0738 | 0.1109 | 1.109 | 1–2 | 0.237–0.272 | 0.270 |
| 15 | 129.01 | 87.981 | 5612 | 30 | 0.0624 | 0.1067 | 1.078 | 1 | 0.213–0.247 | 0.186 |
| 16 | 150.45 | 95.011 | 6060 | 33 | 0.0499 | 0.1027 | 1.070 | 1–2 | 0.192–0.227 | 0.101 |
| 17 | 175.44 | 102.598 | 6544 | 30 | 0.0433 | 0.0988 | 1.046 | 1 | 0.174–0.219 | 0.054 |
| 18 | 204.59 | 110.794 | 7067 | 30 | 0.0389 | 0.0951 | 1.020 | 1 | 0.157–0.197 | 0.026 |
| 19 | 238.58 | 119.644 | 7631 | 30 | 0.0354 | 0.0915 | 1.039 | 1 | 0.142–0.179 | 0.012 |
| 20 | 278.22 | 129.202 | 8241 | 30 | 0.0295 | 0.0881 | 1.163 | 1 | 0.130–0.177 | 0.005 |
| 21 | 324.45 | 139.524 | 8899 | 28 | 0.0263 | 0.0847 | 1.222 | 1 | 0.118–0.154 | 0.002 |
| 22 | 378.35 | 150.669 | 9610 | 25 | 0.0188 | 0.0815 | 1.041 | 1 | 0.108–0.149 | 0.001 |
| 23 | 441.21 | 162.704 | 10378 | 25 | 0.0172 | 0.0785 | 1.021 | 1 | 0.099–0.135 | 0.000 |
| 24 | 514.52 | 175.702 | 11207 | 25 | 0.0152 | 0.0755 | 1.222 | 1 | 0.091–0.127 | 0.000 |
| 25 | 600.00 | 189.737 | 12102 | 25 | 0.0134 | 0.0727 | 1.189 | 1 | 0.084–0.114 | 0.000 |

- **SMD.** TAB's products are larger than eq. [14] in the bag and multimode regimes (2.08× at We 15), inside its
  ±20 % at We 59.8–81.4 (3 classes) and below it from We 95, 5.6× below at We 600: over the series TAB's SMD falls
  about as We^−0.96, eq. [14]'s as We^−1/4.
- **A class represents TAB's r_s.** A parcel draws its new radius from the *mass* distribution of the Rosin–Rammler
  law (spread n = 3.5, truncated to [0.01 r_s, r)) whose r32 is r_s, and keeps its mass flow, so the drops of a
  class have SMD = 1/E[1/d] over the draws = the law's r32 = r_s (`[KIVA2]` App. D and F draw the same way).
- **Distribution.** MMD/SMD is 1.019–1.222 over the classes, mean 1.088. TAB's own value after one event is the
  median of the draws over 1/E[1/d], 1.062 at n = 3.5, not the root-normal 1.2; with 25 to 33 one-event draws a class
  falls in 0.99–1.15 and the mean of 25 classes in 1.04–1.08 (95 %, sampled with this case's class sizes and MMD
  rule). The bag classes (We 15–32, mean 1.059) lie there; above them a product can break again (the events column,
  a lower bound) and the classes reach 1.22.
- **At t_b = 5 t*.** No parcel is still at d0, and no parcel has an event between the two rows bracketing t_b, so the
  interpolated d is a printed d; u_b = |u_g − u_p| is interpolated between the same rows. The ratio of IGLOO's
  u0/u_b − 1 to eq. [13] has the median 0.26 at We 15, 1.10 at We 239 and 2.8 at We 600: in the bag regime the
  products keep more of their relative velocity than H&F's drops (sphere drag against their deformed drops); from
  We 324 every product is smaller than Fig. 6's smallest drop (x > 0.85). 500 of the 714 parcels sit inside Fig. 6's
  x range; every parcel has a slip at t_b.
- **Run.** 13.9 s at 5 threads on a host at load 65 (e2e timeout 300 s).

## Comparison plot (`verify.py` → `OUTPUT/tab-hf92.svg`, not gated)
Parcels are taken by ascending ID; the figure does not depend on the order of the records (`vv-order`).
- (a) eq. [14] on [HF92] Fig. 10's axes, x = (ρ_L/ρ_G)^1/4 [μ_L/(ρ_L d0 u0)]^1/2 We, y = ρ_G SMD u0²/σ: the line
  y = 6.2 x and its ±20 % coefficient band, one marker per class, coloured by regime.
- (b) SMD/d0 against We with eq. [14] and the onsets We = 13, 35, 80.
- (c) MMD/SMD per class against 1.2.
- (d) V(d) against d/MMD per class, pooled by regime, on [HF93] Fig. 1's axes (x the normal quantile of V, 0.1–99 %;
  y linear in (d/MMD)^1/2, labels 0.0–4.0), with root-normal lines for MMD/SMD = 1.5, 1.2, 1.1. The lines come from
  the definition, (d/MMD)^1/2 normal by volume with mean 1 and spread s, cut at d/MMD = 0.01, s fixed by MMD/SMD =
  E_V[MMD/d]: s = 0.291, 0.219, 0.168. The papers name Simmons' function without its formula; the lines printed in
  [HF93] Fig. 1 are 5–15 % flatter than these. [HF93]'s text names its lines 1.0, 1.20 and 1.50 while its legend
  prints 1.5, 1.2, 1.1; the figure follows the legend.
- (e) eq. [13] per parcel at t_b: x = (ρ_G/ρ_L)^1/2 d0/d, y = u0/u_b − 1, the correlation solid over 0.058–1.07 and
  dotted beyond, Fig. 6's marker range shaded; the frame (x 0.01–10, y 0.01–10⁴) holds every parcel (the figure
  states the count of any parcel outside it or without a slip).

Not drawn: the papers' data points (not digitized; the correlations carry the comparison) and the core drop after
shear breakup ([HF93] §3.2–3.4: a TAB parcel stands for identical product drops and has no core).
