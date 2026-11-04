# INFO — erosion/wall-laws (the erosion laws in `OUTPUT/wall.tec`, e2e)

## Reference
- `[Oka05a]`, `[Oka05b]`: E(α) = g(α) E90 [mm³ of wall per kg of particles], g(α) = sin(α)^n1 (1 + Hv(1 − sin α))^n2
  (Part 1 eqs. 1–2), E90 = K (a Hv)^(k1 b) (v/v′)^k2 (D/D′)^k3 (Part 2 eq. 4); silica's constants (Part 2 Table 1:
  K 65, k1 −0.12, k2 = 2.3 Hv^0.038, k3 0.19, v′ 104 m/s, D′ 326 µm; Part 1 Table 3: n1 = 0.71 Hv^0.14,
  n2 = 2.4 Hv^−0.94).
- `[NG68]`: eqs. 3–5, Table I at 424 ft/s (φ 3571, ε 8375 ft·lbf/g, n 5).

The test is a plumbing gate: the laws themselves are `erosion/laws`'s; here the oracle (`check.py`, Python) applies
them, as the papers write them, to the impact speed, angle and diameter of each stream and checks what IGLOO sums
per wall cell and writes.

## Fixture
`erosion/wall-flux`'s `INPUT/` (symlinked): the 15×5×5 box of 0.01 m cells, face 2 (x = 0.15) coded 301, the
other faces 400; material A in two groups, B in one (`two-mat`'s properties). Eight assigned streams with
`NoDrag`/`NoHeat` start at x = 0.148 and move on straight lines at their injection velocities:

| stream | \|v\| [m/s] | angle from the wall | d [µm] | ṁ [kg/s] | wall cell (j, k) |
|---|---|---|---|---|---|
| 1 | 100 | 90° | 10 | 1.0e-4 | (1, 1) |
| 2 | 100 | 70° (+y) | 20 | 2.0e-4 | (2, 1) |
| 3 | 100 | 50° (+z) | 30 | 1.5e-4 | (3, 1) |
| 4 | 100 | 35° (+y +z) | 50 | 0.5e-4 | (4, 1) |
| 5 | 100 | 20° (−y) | 80 | 1.0e-4 | (5, 1) |
| 6 | 100 | 10° (+z) | 120 | 0.8e-4 | (1, 3) |
| 7 | 60 | 15° (+y) | 40 | 1.2e-4 | (2, 4) |
| 8 | 80 | 60° (−z) | 15 | 0.6e-4 | (1, 1) |

`[IGLOO-Erosion]`: `model = Oka Neilson-Gilchrist`, `wall-density = 1800`; the wall `Hv = 1.43` GPa, `a = 0.05`,
`b = 1.2` (a·Hv = 0.0715, inside Part 2's 0 ≤ aHv ≤ 1); A: silica's constants with `k2-s`/`k2-q`, `n1-s`/`n1-q`,
`n2-s`/`n2-q` (n1 = 0.7465, n2 = 1.7147, k2 = 2.3315 at Hv = 1.43); B: `phi = 4841625.893 J/kg`,
`epsilon = 11354975.317 J/kg` (Table I converted), `n = 5` (α0 = 18°), `k = 20 m/s` (the paper's K: the
Neilson-Gilchrist key is lower-case `k`, Oka's particle factor `K`, so `K = 65.0 0.0` and `k = 0.0 20.0`). Streams 6
and 7 hit below α0 and with V sin α ≤ k (cutting alone); streams 1 and 8 share a cell at 90° and 60°.

## What it verifies (`check.py`)

| id | gate | tolerance | measured |
|---|---|---|---|
| L0 | the run ends normally; the log names A: Oka and B: Neilson-Gilchrist with the constants it uses, n1, n2, k2 the correlations s·Hv^q | 6e-6 (ES12.5) | within |
| L1 | `er-A`, `er-B` in the 7 hit cells = Σ ṁ E / Σ ṁ of the oracle (the groups cancel); exactly 0 in the 18 others | 1e-12 relative | 3.9e-15, 1.8e-15 |
| L2 | `g-A`, `g-B` = groups × Σ ṁ / (dy dz); `rdot` = (g-A er-A + g-B er-B)/ρ_wall | 1e-12 relative | 1.5e-16; 1.9e-15 |
| L3 | the deck covers both branches of Neilson-Gilchrist and its threshold (streams 6, 7 below α0 and k; the others above), and cell (1,1)'s two impacts differ in angle | exact | — |
| L4 | `ref-direct/`: the deck with `n1`, `n2`, `k2` as values (the correlations, 17 digits): `er-A`, `er-B`, `rdot` the same | 1e-13 relative | 0 |

## RED (measured: each mutation built, this case and `wall-flux` run)

| mutation | what fails |
|---|---|
| the angle from the normal: acos for asin in `depositErosion` | L1 (er-A 0.61, er-B 1.55 relative), L2 (rdot 1.09) |
| every parcel's law of material 1 | L1 (er-B 0.88), L2 (rdot 0.71) |
| the correlations as q·Hv^s | L0 (n1 0.180, n2 −2.22, k2 0.087), L1 (er-A 0.86), L2, L4 |
| `rdot` from the incident flow instead of the eroded one | L2 (rdot 7.2e3 relative); `wall-flux` W5 |
| the tokens Oka and Neilson-Gilchrist swapped | the setup stops (A's `phi` 0 out of range): L0 |
| a fixed 10 µm for the parcel's diameter | L1 (er-A 0.38), L2 (rdot 0.097) |
| B's threshold read from Oka's `K` (B's value 0) instead of `k` | L0 (k 0 for 20), L1 (er-B 0.50), L2 (rdot 0.33) |

The impact speed is the parcel's at the crossing point (`vImp`); with `NoDrag` it equals the injection speed, so
this case cannot tell it from the end-of-step speed.
