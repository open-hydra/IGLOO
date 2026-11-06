# INFO — breakup/khrt-hf93 (Reitz KH-RT's parent drop against Hsiang & Faeth's core drop; e2e, drawn)

## What it verifies
Reitz KH-RT on Hsiang & Faeth's water series, the run behind a figure that sets the model's parent drop beside the
core drop `[HF93]` measured at the end of shear breakup. In KH-RT the parent is stripped by the KH rate until its first
RT event, which replaces it by drops of `λ_RT`; HF93's core drop is what remains of the original drop when the
stripping ends, at `t_b = 5 t*`. The comparison is drawn, never gated: a correlation of shock-tube drops is no closed
form of the model. `check.py` gates what the figure needs.

## The case
[tab-hf92](../tab-hf92/INFO.md)'s input (`INPUT` symlinked): d0 = 1000 µm water (ρ_L 997, σ 0.0708, μ_L 8.94e-4) in
air (ρ_G 1.18, 1.85e-5 Pa s, U = 200 m/s), 25 inlet cells, one diameter-based Weber number each, We 15 to 600; a 1.2 m
box of 480 cells whose lateral faces reflect. `breakup = Reitz-KHRT` at its defaults (B0 0.61, B1 20, Ctau 1, CRT 0.1,
WeLimit 6) with `mShedLim = 1`: no KH shed (a shed keeps the parent's d, so the parent's d(t) does not depend on it).
`print-dtime = 1e-12`, no `print-dcell`: a row at every segment end, so the first RT, which ends a segment, has its own
row and the record before it is the last cell face or plane crossed. `plane-x = 0.0005 0.001 0.0015 0.002`: four planes
in the first 2.5 mm cell. `ds = 0.5`: 126 parcels, 4 to 8 per class (KH-RT is deterministic, a class's parcels are
copies).

## Gates (`check.py`)
| id | gate | tolerance | measured |
|---|---|---|---|
| H0 | every class has a parcel, by tab-hf92's P0 rule (the inlet cell holding the (Y, Z) of the injection row, whose `|u_p|` is that cell's `kv U`); every parcel has a class and one exit record | F12.6 | 25 classes, 126 parcels, 4..8 per class |
| H1 | every parcel has a first RT: among its records (trajectory rows and plane records, by time) the first whose d is below 0.7 of the record before it; it is a trajectory row, and a record other than the injection row precedes it | exact on the printed values | 126 of 126; the RT step 0.564 (We 15) to 0.122 (We 600) of the record before it |

The threshold: an RT step takes d to 0.56 of the record before it at most; KH stripping alone moves d by 5.1 % at most
between two records before the first RT here, and by 14 % in a KH-only run (`CRT = 100`, class 25 between two face
rows): 0.7 separates the two. H1 proves that an RT record exists; that it is the parcel's first RT rests on two facts
it does not gate: every event writes a row (an event on a step of zero length would not), and every class's first RT
takes d below 0.7 of the record before it.

## RED
| change | what fails |
|---|---|
| `CRT = 100` (`λ_RT` about 0.5 m, far above d0: no RT) | H1: 0 of 126 parcels have a first RT |

## [HF93]'s core drop (pp. 725-728)
With an average drag coefficient `C̄_D` over the breakup period and the core drop's mass taken constant: [1]
`du/dt = −3 C̄_D ρ_G u²/(4 ρ_L d0)` for the relative velocity u; [2] `u = u0/[1 + (3 C̄_D t/4t*)(ρ_G/ρ_L)^{1/2}]`,
`t* = d0 (ρ_L/ρ_G)^{1/2}/u0`; [3] `(u'_b/u0)(ρ_L/ρ_G)^{1/2}(1 + 3C) = 3 C̄_D (t_b/t*)/4` for the velocity `u'_b` the
core drop gains, the relative velocity being `u0 − u'_b`; [4] `C = (3 C̄_D t_b/4t*)(ρ_G/ρ_L)^{1/2}`; [5] `Eo = a ρ_L
d²/σ`; [6] the acceleration, from [2]; [7] `Eo_cr = (3 C̄_D We/4)(d_b/d0)²/(1 + C)²` at the end of stripping; [8] the
core drop's Weber number `We_cr = (4 Eo_cr We/3 C̄_D)^{1/2}/(1 + C)`. The constants: `t_b/t* = 5`, `C̄_D = 5` (p. 726),
`Eo_cr = 16` (p. 727); the correlations are fitted over "the test range of the shear breakup regime, 100 < We < 600"
(p. 726). For this series (ρ_G/ρ_L = 1.18/997) `C` = 0.64506, the relative velocity at `t_b` `u_b/u0 = 1/(1 + C)` =
0.6079, and [7] gives `d_b/d0 = (1 + C)(4 Eo_cr/(3 C̄_D We))^{1/2}`: 0.8774 at We 15, 0.1387 at We 600.

Three of the printed forms disagree with the paper's own derivation and figures; the figure draws the derived ones:
- [3] prints `(1 + 3C)`. [2] at `t_b` gives `u_b/u0 = 1/(1 + C)`, i.e. `u'_b/u0 = C/(1 + C)`, which is [3] with `(1 +
  C)`; Fig. 2's ordinate (p. 727) reads `(u'_b/u0)[1 + 3C̄_D(ρ_G/ρ_L)^{1/2}(t_b/t*)/4]/(ρ_G/ρ_L)^{1/2}`, the `(1 + C)`
  form, as does Fig. 4's abscissa (p. 728), `We/[1 + 3C̄_D(ρ_G/ρ_L)^{1/2}(t_b/t*)/4]²`; and the text's "the relative
  velocity of the core drop at the end of breakup only was 30-40% lower than the initial relative velocity" is what
  `(1 + C)` gives here (39 %), not `(1 + 3C)` (22 %).
- [8]'s coefficient: "the coefficient in [8], (4Eo_cr/3C̄_D)^{1/2} = 2.3", where `(4·16/(3·5))^{1/2}` = 2.0656. Panel (a)
  draws both.
- [6] prints the factor `(1 + Ct/t_b)²`; the derivative of [2] carries it with the exponent −2, as [7]'s `/(1 + C)²`
  requires.

## The figure (`verify.py` → `OUTPUT/khrt-hf93.svg`, drawn, not gated)
One parcel per class (its lowest ID). (a) d/d0 against We: the parent before its first RT (the d of the record before
the RT row) and the parcel at `t_b` (hollow, interpolated between the rows that bracket it), beside [7]'s `d_b/d0` with
2.0656 and with the printed 2.3; (b) the relative velocity over u0 at the first RT and at `t_b` beside `1/(1 + C)`; (c)
the time of the first RT over t* beside `t_b/t* = 5`. HF93's lines are solid over We 100-600 and dotted below; the
regime onsets We 13, 35, 80 are marked.

The record before the first RT is a cell-face row, except for class 25 (We 600), whose first RT lies in the first
cell (x = 2.111 mm) and follows the 2-mm plane. The drawn d carries the stripping up to that record only: a run with a
plane every 0.5 mm puts the d just before the RT up to 1.7 % lower (classes 1-24), and the scatter of the drawn
parent above We 80 (0.880 to 0.910 d0) is the record's placement.

Drawn result: the parent before its first RT is 0.967 d0 at We 15, 0.880 at We 81, 0.903 at We 600, against HF93's
0.877 to 0.139: 1.10 times the core drop at We 15, 6.5 times at We 600. Its relative velocity at the first RT is 0.980
(We 15) to 0.993 (We 600) of u0 against 0.608. Its first RT comes at 2.22 t* (We 15) to 1.27 t* (We 600), before HF93's
`t_b = 5 t*`; the RT cascade then runs on, and at `t_b` the parcel is 0.30 d0 at We 15 and 0.002 d0 at We 600, its
relative velocity 0.90 of u0 at We 15 and below 1e-4 of it from We 324 on. KH-RT has no core drop: its parent, 0.88
d0 or more at the last record before its first RT, is replaced whole by the RT's products.

## Not covered
The KH shed and its children (`mShedLim = 1`); the other liquids of HF93's series (ethyl alcohol, n-heptane, the
glycerol mixtures); the size distribution of the stripped drops.
