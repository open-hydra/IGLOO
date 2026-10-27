# repeatability/drag-stokes-time-end — the exit code's reset between sweeps

## What it verifies
A parcel's exit code is the boundary code of the face it leaves through, and 0 when it ends inside the domain (its
time-end record). `reset_state` zeroes it for every pinned parcel before a sweep. This case drives one population
through sweeps that end inside, at the outlet and inside again, so a code that survived a sweep would be printed by
the next sweep's time-end records.

`standard/drag-stokes`'s case (`INPUT` symlinked): a uniform 10 m/s, 600 K gas on the 0.15 m box, Stokes drag,
25 parcels of 11.89 µm (ρ_p 2950, τ_p = ρ_p d²/(18 μ) = 1.2872 ms) injected at `kV` 0.1 of the gas speed, with
`[IGLOO-General] time-end = 12e-3`. On x(t) = u_g t + (v_0 − u_g) τ_p (1 − e^(−t/τ_p)):

| sweep of the gas cycle | gas | the parcel | record |
|---|---|---|---|
| 0 | original, u_g = 10 m/s, v_0 = 1 m/s | reaches the outlet at 16.158 ms: at 12 ms it is at x = 0.108416 m | time-end, code 0 |
| 1 | U doubled, v_0 = 2 m/s (a 401 stream takes its speed from the gas at injection) | leaves through the outlet at 8.657 ms | exit, code 400 |
| 2 | original again | as sweep 0 | time-end, code 0 |

## Gates
The two-sweep harness's (`tools/check_twosweep.py`): **steady**, sweep 1 ≡ sweep 0 (both end inside), and
**gas-cycle**, sweep 1 differs from sweep 0 and sweep 2 ≡ sweep 0 (`outloc`, `trajectories`, `scatter` as exact
multisets, `source.tec` within 1e-12 of the field scale), with the premise that sweep 0's parcels moved. The MPI twin
`mpi-repeat-drag-stokes-time-end` runs the same at 2 ranks × 2 threads.

## Falsification
The gate fails on a `reset_state` that keeps the exit code, serial and at 2 ranks: sweep 2's 25 time-end records
print code 400 at x = 0.108416 m, and

    outloc-A: sweep2 (field restored) vs sweep0: 25 vs 25 records; 25 only in first, 25 only in second

The steady half and the premise still pass there, and so does `repeat-drag-stokes`, whose parcels leave through the
outlet in every sweep (each exit sets the code again). The trajectories and the scatter cloud stay equal: the exit row
of a trajectory is written only for a parcel gone through a face.
