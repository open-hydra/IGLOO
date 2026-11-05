# scatter-placement — where the scatter cloud's markers fall, how many there are and what they weigh

`khrt-stress`'s case: the `khrt-e2e` box, gas and inlet (`INPUT/` symlinked; 25 injected parcels whose diameters
sweep 0.36 to 12 mm over the inlet cells) under KH-RT breakup with `mShedLim = 0.01`, which sheds 701 children, and
Rayleigh-Taylor shatters that multiply a parcel's number rate in place. `out-file = ALL` writes the Eulerian number
moment (`euler1.tec`) and `out-time = on` stamps the trajectory rows with the parcel time.

## What it verifies
A parcel writes markers at a constant rate `r = ṅ_ref/dNscat`, `ṅ_ref` the number rate it entered with raised to
the material's floor `F·n̄` (`F` the `[IGLOO-General] scatter-floor`, 0.1 by default, `n̄ = Ṅ/nStreams` the mean of
the 25 injected streams' rates; a rate below `(1 − 10⁻⁶)·F·n̄` takes the floor; a breakup child keeps its
parent's), at the instants its accumulator, seeded at `frac(ID·φ)`, crosses a whole number:
`t_k = t_birth + (k − frac(ID·φ))/r`, inside the solver's steps. The position lies on the cubic Hermite interpolant
of the step's end positions and velocities, the other columns are linear in the step, and each marker carries its
weight `w = ṅ_p(t_k)/r` (the column `weight`), the droplets it stands for.

- **P1** markers within 1e-6 m (the `F12.6` print resolution) of an x-face plane are at most 3 × 2·10⁻⁶/Δx of all
  markers: a marker uniform in time is near a face as often as the face's share of the cell (0.0008), a marker
  at a step end crowds the faces, where the steps are cut.
- **P2** no two markers alike in every printed column.
- **P3** the markers' fractional x in their cells against U(0,1): Kolmogorov-Smirnov p ≥ 1e-3.
- **C** per parcel, the marker count `N = ⌊frac(ID·φ) + r·τ⌋`, `τ` its exit `t` less its birth `t`, `r` with the
  floor.
- **F** the floor's premise: `n̄/dNscat = fsample-traj · V_ref/L_ref` (the quantum's definition, `n̄` over the 25
  roots), `V_ref` the roots' mean first-row speed, `L_ref` the mesh's bounding-box diagonal.
- **W1** the Euler moment Σ n_p V (the ODE's integral of ṅ_p dt, deposited per cell) lies in the bracket
  [Σ I_lo, Σ I_hi] the parcels' printed states give.
- **W2** per parcel, the markers' Σ w lies in [I_lo − B, I_hi + B], `B = (TV + M)/r`.
- **W3** over all parcels, Σ w differs from the Euler moment by at most the sum of the parcels' bounds on
  |S_p − I_p|: I_hi for a parcel without a marker, B + I_hi − I_lo for the others.

C: the accumulator starts at `frac(ID·φ)` and adds `r·h` on every step the parcel keeps, and every crossing of a
whole number writes a marker, the exit face's included, so the count is exact whatever the steps. `r` is the
number rate of the parcel's injected ancestor at its entry, `ṅ_p = ṁ/m` from its birth flow and its first row,
raised to the floor `F·n̄` when below `(1 − 10⁻⁶)` of it, over `dNscat` (the 701 children are all of the first
generation, and the case holds no residue); `n̄` is the mean of the 25 roots' rates, an interval from theirs, and a
rate whose interval meets the floor's accepts either branch. Every printed
input enters with its half-ULP, and a parcel whose interval of `frac(ID·φ) + r·τ` holds an integer accepts either
floor; `dNscat`'s `ES10.3` dominates the interval (±8.5·10⁻⁵ relative, ±0.07 at 820 markers). The weight clauses
gate the `ṅ_p` factor of `w`, C gates `r`: a wrong rate leaves Σ w a Riemann sum of the same integral, a wrong
`ṅ_p` leaves the count.

The bracket and the bound come from the model and the parcel's printed states (its trajectory rows and its exit
record), not from the markers. With a constant density and no evaporation `ṅ_p = ṁ/m`: the flow `ṁ` falls only at
a shed, by the child's birth flow (a parcel's birth flow is its exit flow plus its children's), and the droplet
mass `m` never grows (KH stripping shrinks it at a constant flow, a shatter divides it at a constant flow, a shed
keeps it); `check.py` checks that premise on the printed `m`. Between two printed states `ṅ_p` therefore lies
between `ṁ(end)/m(start)` and `ṁ(start)/m(end)`, which times the interval give I_lo and I_hi, and its variation
is at most `ṁ(start)(1/m(end) − 1/m(start))` plus the interval's sheds over `m(end)`; TV is the sum over the
parcel's life, M the largest upper value. The markers are the Riemann sum, step 1/r, of the linear interpolant of
`ṅ_p` over the solver's steps: its integral lies in [I_lo, I_hi], its variation and maximum are at most TV and M,
and the sum's error is at most `(TV + M)/r` — TV for the samples inside their intervals (a shatter's jump
included), M for the seed and the unwritten tail. The printed inputs' rounding enters as factors (`dNscat`'s
`ES10.3` dominates B) and each weight with its half-ULP.

W3 is loose where parcels hold few markers, whose B exceeds their integral (its bound is 0.41 times the moment
here); W2 holds each parcel to its own.

## Pass criterion
P1-P3, C, F, W1-W3 and the premise all hold, and C's window holds one value on at least 90 % of the parcels.
Measured at the default floor (0.1; 18 of the 25 roots below it): 37 992 markers, 31 parcels without one (all of
them children, against 278 at `scatter-floor = 0`, 12 of them roots); P1 46 (0.0012, bound 0.0024); P2 none alike;
P3 D = 0.0055, p = 0.21; C 726/726, 721 windows of one value (two on IDs 33, 37, 46, 50 and 77), one quantum in
(0.585928, 0.58595] fitting every count; F `n̄/dNscat` 6.030594·10⁴ against 6.030227·10⁴ (tolerance 5.4); W1
Σ n_p V = 8.6544·10⁷ in [8.4663·10⁷, 8.8425·10⁷]; W2 726/726 parcels inside, 324 of them (69.8 % of the moment)
inside windows of ±50 %; W3 Σ w 0.32 % below the moment, bound 40.6 %.

## Falsification
At the default floor:
- Markers without the floor (every parcel at its own entry rate): C fails on 389 of 726 parcels, the 18 roots
  below the floor among them.
- The floor on the breakup children only, the roots at their own rate: C fails on the 18 roots below it.
- The floor `F·Ṅ`, not over the stream count: 83 841 markers, C fails on 663 of 726 parcels.
- The markers at the floored rate, the weights at the entry rate (`w = dNscat·ṅ_p/ṅ_entry`): W2 fails on 343 of
  726 parcels, W3 on their sum, 18.9 times the moment above it; C passes.

At `scatter-floor = 0` (the floor off):
- Markers at the end of the step in which they fall: P1 110 (0.0030), P2 4 158 alike, P3 D = 0.272.
- Markers at the step ends under a number-weighted accumulator (one per `dNscat` of `ṅ_p·dt`, no weight column):
  1 595 327 markers, P1 0.397, P2 0.215, P3 D = 0.368.
- A breakup child given its own number rate at the shed instead of its lineage's: 314 322 markers, C fails on
  675 of 726 parcels; P1-P3 and W1-W3 pass.
- A breakup child given its parent's `ṅ_p` at the shed instead of its `ṅ_ref`: 70 048 markers, C fails on 81 of
  726 parcels; P1-P3 and W1-W3 pass.
- Weights frozen at the parcel's entry (`w = ṅ_entry/r`): W2 fails on 312 of 726 parcels (70.9 % of the moment),
  W3 on their sum, 99.8 % below the moment against a bound of 92.4 %. `scatter-weight`, `scatter-weight-evap` and
  `scatter-weight-solid` pass under it, their number rates being constant.
