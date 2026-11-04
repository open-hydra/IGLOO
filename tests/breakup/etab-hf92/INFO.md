# INFO — breakup/etab-hf92 (ETAB products of Hsiang & Faeth's water series; e2e, the correlations drawn, not gated)

## Reference papers
- `[HF92]`, `[HF93]`: the SMD correlation ([HF92] eq. [14]), the root-normal distribution with MMD/SMD = 1.2, the
  breakup time t_b = 5 t* and the size/velocity correlation ([HF93] eq. [13]), as `tab-hf92/INFO.md` states them.
- `[Tan97]`, `[Tan98]` the ETAB cascade, as `etab-e2e` runs it (`breakup = ETAB`, its defaults): an event resizes
  the parcel to r exp(−(K_br/ω) acos(1 − 1/We_Cr)), We_Cr = We_r/12, and gives the product a normal velocity of
  random azimuth.

## The case
`tab-hf92`'s case with `breakup = ETAB`: the same `INPUT/` (a link to `../tab-hf92/INPUT`, whose `INFO.md` gives
the generating command), 25 classes of d0 = 1000 µm water drops at We 15–600 (one inlet cell each, the slip set per
cell), air at 297 K, U = 200 m/s, the 1.2 m box with reflecting lateral faces, `ds = 0.2` (714 parcels, 25 to 33 per
class), `print-dcell = 4`, `out-time = on`, `out-scatter = off`. `check.py` runs `tab-hf92`'s oracle (P0–P3) on this
run and `verify.py` draws through `tab-hf92`'s, with this directory's output.

**The lateral faces reflect, and here parcels reach them.** The injection packing puts the edge cells' parcels on
the lateral faces and the ETAB kick drifts the products sideways (up to 12.4 mm whatever the width): with faces 3–6
coded 400, 46 of the 714 parcels leave laterally (26 of them before t_b, the last row's We_r up to 99.8), and 24
still do in a 0.2 m cross-section at `ds = 0.8` cm. Coded 300 (`--side-code 300`) the faces reflect specularly, or
let a grazing parcel slide along them (|v_n| < 0.02 |v|, `obj_bc.f90`). In a uniform gas parallel to the faces a
reflection flips the sign of the normal component of the slip and keeps its magnitude, which is what the drag and
the breakup models read, so it leaves the cascade as it is, up to the grazing slide. Measured: the 46 parcels that
leave laterally under 400 end, under 300, inside the final-d range of their class's other parcels or outside it by
at most 4.74e-3 of the range's lower end (class 21), against those classes' spreads (max − min)/min of 8.6e-4 to
5.1e-2, and each takes a number of d steps inside its class's range. P2's "no lateral exit" therefore checks the
input and the reflection, not the physics.

## What `check.py` gates
`tab-hf92`'s P0–P3 (the table and the definitions are in `tab-hf92/INFO.md`).

| gate | measured |
|---|---|
| P0 population | 714 parcels, 25 to 33 per class, 0 without a class |
| P1 every drop broke | max d_exit/d0 0.643858 |
| P2 the cascade ended inside the box | 0 lateral exits, 0 parcels changing d over x ≥ 0.6 m, max We_r on the last row 4.483 |
| P3 every parcel in the box at t_b | 0 parcels leaving before t_b, min t_last/t_b 1.444 |

RED (scratch copies, each exit 1):

| mutation | what fails |
|---|---|
| every inlet record's kv 0.97 (`tools/set_kv.py`; slip 6 m/s, We 0.6) | P1 max d_exit/d0 1.000000; P3 658 of 658 leave before t_b; P0 56 parcels without a class |
| box `--lx 0.05 --nx 20` | P2 274 parcels change d over x ≥ 0.025 m, max We_r 36.050; P1, P3 |
| box `--lx 0.3 --nx 120` | P3 313 of 714 leave before t_b, min t_last/t_b 0.364; P2 373 parcels change d over x ≥ 0.15 m |
| `ds = 10` (one parcel per cell) | P0 1 parcel per class; P1–P3 pass |
| faces 3–6 coded 400 (`--side-code 400`) | P2 46 of 714 leave before the outlet, 43 change d over x ≥ 0.6 m, max We_r 99.792; P3 26 leave before t_b, min t_last/t_b 0.174 |

## Measured (the reference run)
Columns as in `tab-hf92/INFO.md`; "events" (d changes between rows 1 cm apart) is a lower bound.

| class | We | n | SMD/d0 | [HF92] eq. [14] | MMD/SMD | events | x(5 t*) [m] | We_r, last row |
|---|---|---|---|---|---|---|---|---|
| 1 | 15.00 | 27 | 0.6058 | 0.1827 | 1.000 | 1 | 0.829 | 3.616 |
| 2 | 17.49 | 25 | 0.6439 | 0.1759 | 1.000 | 1 | 0.757 | 4.483 |
| 3 | 20.40 | 25 | 0.3850 | 0.1692 | 1.000 | 2 | 0.692 | 2.552 |
| 4 | 23.79 | 25 | 0.4405 | 0.1628 | 1.000 | 2 | 0.630 | 3.521 |
| 5 | 27.74 | 25 | 0.2659 | 0.1567 | 1.000 | 3 | 0.575 | 1.768 |
| 6 | 32.35 | 33 | 0.3258 | 0.1508 | 1.000 | 3 | 0.521 | 2.803 |
| 7 | 37.72 | 30 | 0.2084 | 0.1451 | 1.000 | 4 | 0.474 | 1.336 |
| 8 | 43.99 | 30 | 0.2606 | 0.1396 | 1.000 | 4 | 0.427 | 2.304 |
| 9 | 51.30 | 30 | 0.1781 | 0.1344 | 1.000 | 5 | 0.386 | 1.097 |
| 10 | 59.82 | 30 | 0.2172 | 0.1293 | 1.000 | 5 | 0.346 | 1.867 |
| 11 | 69.76 | 33 | 0.1506 | 0.1244 | 1.000 | 6 | 0.311 | 0.801 |
| 12 | 81.35 | 30 | 0.0902 | 0.1198 | 0.999 | 7 | 0.280 | 0.165 |
| 13 | 94.87 | 30 | 0.1066 | 0.1152 | 1.000 | 7 | 0.249 | 0.293 |
| 14 | 110.63 | 30 | 0.1037 | 0.1109 | 1.000 | 7 | 0.222 | 0.266 |
| 15 | 129.01 | 30 | 0.0852 | 0.1067 | 1.000 | 6–7 | 0.200 | 0.125 |
| 16 | 150.45 | 33 | 0.0904 | 0.1027 | 1.000 | 4–5 | 0.180 | 0.154 |
| 17 | 175.44 | 30 | 0.0686 | 0.0988 | 1.000 | 4 | 0.165 | 0.047 |
| 18 | 204.59 | 30 | 0.0468 | 0.0951 | 1.000 | 3–4 | 0.151–0.152 | 0.006 |
| 19 | 238.58 | 30 | 0.0552 | 0.0915 | 1.000 | 3–4 | 0.127–0.128 | 0.016 |
| 20 | 278.22 | 30 | 0.0342 | 0.0881 | 0.999 | 3 | 0.121 | 0.001 |
| 21 | 324.45 | 28 | 0.0346 | 0.0847 | 0.999 | 2–3 | 0.105 | 0.001 |
| 22 | 378.35 | 25 | 0.0295 | 0.0815 | 0.999 | 2–3 | 0.095–0.096 | 0.000 |
| 23 | 441.21 | 25 | 0.0235 | 0.0785 | 1.414 | 2 | 0.081–0.100 | 0.000 |
| 24 | 514.52 | 25 | 0.0218 | 0.0755 | 0.999 | 1 | 0.083 | 0.000 |
| 25 | 600.00 | 25 | 0.0143 | 0.0727 | 0.992 | 1 | 0.082–0.083 | 0.000 |

- **SMD.** ETAB's products are 3.3× eq. [14] at We 15 (3.7× at We 17.5), above it up to We 70, inside its ±20 % at
  We 95, 111 and 150 (3 classes) and 5.1× below it at We 600. Up to We 70 every second class adds an event: the SMD
  falls where the count grows and rises between two classes that share it (We 20.4: 0.385, We 23.8: 0.441, two
  events each), so it zig-zags down the series.
- **One size per class, except where the cascade ends at the onset.** The cascade is deterministic and the kick's
  random azimuth leaves the slip magnitude unchanged at the first event: the classes' (max − min)/SMD are 0 (classes
  1–4) to 5.1e-2 (class 25), MMD/SMD 0.992–1.000. Class 23 (We 441) is bimodal: its first event takes 23 parcels to
  90.9 µm, one to 116 µm, one to 335 µm; the cascade then reaches 33.3 µm at We_r 6.8 to 3.5 (the slips of the rows at
  1 and 2 cm), across the onset We_r = 6: 13 parcels end at 33.27–33.32 µm, 12 at 17.66–18.17 µm (MMD/SMD 1.414), a
  ratio of 0.530–0.546, which the ETAB law gives for one more event at We_r 6.12–6.27.
- **At t_b = 5 t*.** No parcel is still at d0 and none has an event between the rows bracketing t_b. The ratio of
  u0/u_b − 1 to eq. [13] is 0.19–0.27 up to We 111 and rises to 2.7 at We 600; it crosses 1 between We 441 and 515.
  504 of the 714 parcels sit inside Fig. 6's x range.
- **Run.** 12.2 s at 5 threads on a host at load 70 (e2e timeout 300 s).

## Comparison plot (`verify.py` → `OUTPUT/etab-hf92.svg`, not gated)
`tab-hf92`'s figure (parcels by ascending ID, `vv-order`): (a) eq. [14] on [HF92] Fig. 10's axes with its ±20 %
band; (b) SMD/d0 against We with the regime onsets; (c) MMD/SMD per class against 1.2; (d) instead of a
distribution, each parcel's d/SMD of its class against We, with the classes' spread stated (no curve: one size per
class; `tab-hf92`'s panel carries the root-normal lines); (e) eq. [13] per parcel at t_b, Fig. 6's marker range
shaded, every parcel inside the frame. Not drawn: the papers' data points and the
core drop (an ETAB parcel stands for identical product drops).
