# INFO — infrastructure/solver-fail-refine (a solver failure on a refinement call ends the parcel at its last accepted step, e2e)

## Reference
Not a physics case: no literature reference, no oracle curve. The case gates one exit of `Lib_Integration`'s
`ODEsystem`: a `Run_ODESolver` call that returns `err < 0` on a **refinement call**, the call that integrates a cut
step again after a crossing interrupted the previous one. The parcel reverts to the last step the solver accepted
and is marked gone there, so its exit record, its Eulerian deposit and its coupling source describe one state:
position and velocity at that step, the residence up to its time, the path and the moment integrals up to it.

## Fixture
`INPUT/` symlinks `solver-fail-midway/INPUT` (`tools/make_box_case.py --kv 1.0 --u2 5 --rp 1.0e-6`: 60×5×5 cells over
0.15 m, gas at 10 m/s up to x = 0.075 m and 5 m/s beyond, 25 parcels of d = 2 µm, ρ_p = 2950 kg/m³, Stokes drag,
model 1, injected at the gas velocity). `input.ini` is `solver-fail-midway`'s with `[IGLOO-ODE] max-steps-ode = 10`:
each parcel's steps grow until one crosses the face of the first cell (x = 0.0025 m); the crossing interrupts the
call, and the call that integrates the cut step again exceeds the 10 steps (SDIRK4 `EXIT OF SDIRK4`, err −1) before
it reaches the face. `out-file = ALL`, `out-time = on`, `mollify = off`, gas order 1.

**The limit's window** (one thread, `max-steps-ode` from 1 to 14): from 1 to 3 the first call fails, at injection
(t = 0); from 4 to 11 a refinement call fails in the first cell (from 4 to 9 the cut step from t = 2.312012e-4 s ends at
x = 0.002557 m; at 10 and 11 the one from t = 2.478175e-4 s at x = 0.002502 m); from 12 the first failure falls at the gas
velocity jump, on the first call of a later segment (`solver-fail-midway`).

**The sibling `ref-g/`**, run by `check.py` with the same executable: this `input.ini` with `[IGLOO-General]
body-accel = 0.1 0.0 0.0`, `max-steps-ode = 7` and `relative-tol = absolute-tol = 1e-6`. The body force along x starts
a slip transient (τ_p = 3.64e-5 s, terminal slip g τ_p = 3.6e-6 m/s) below the tolerance scale (1.1e-5 m/s), and the
window is again 4 to 11 with the same cut steps. At the default tolerances (1e-10) the window is 8 and 9; at
g = 1 m/s² and the tolerances at 1e-6 it is 4 to 6.

## What it verifies (`check.py`)

| id | gate | tolerance | measured |
|---|---|---|---|
| R0 | 25 exit records; 25 `Run_ODESolver err=` lines, one per parcel ID, each with `time=` > 0 (the time before the revert: the failing call follows latched steps; with R1, a refinement call of the first segment) | exact | 25; `time=` 2.4782e-4 s (ref-g 2.3120e-4 s) |
| R1 | every exit record inside the first cell, X ≤ 0.0025 m + 0.5e-6 m (the F12.6 half-unit), at t > 0 | exact | X = 0.002497 m, t = 2.49702669e-4 s, the solver's last accepted step (ref-g 0.002430 m, 2.42966031e-4 s) |
| R2 | Σ_cells ρ_p V_cell (`euler1.tec`, V_cell from its nodes) = Σ ṁ t of the exit records (model 1 deposits the flow over the residence) | 1e-5 relative (the E13.6 `mdot` column bounds the exit side at 8e-7) | −2.933e-7 (ref-g −2.930e-7) |
| R3 | ref-g: Σ_cells Fx (`source.tec`) = Σ_p [ṁ_p (u_inj,p − \|u_p\|) + ṁ_p g t_p], the drag reaction plus the impulse over the residence (each parcel's one segment starts at injection; gas and g along x, so \|u\| = u); u_inj,p the U of the parcel's t = 0 trajectory row, \|u_p\|, ṁ_p and t_p its exit record's | 3 Σ_p δ_p + 1e-12 Σ\|Fx\|, δ_p the prints: ṁ_p (2 × 0.5e-6 m/s) + h(ṁ_p) \|u_inj,p − \|u_p\| + g t_p\| + ṁ_p g h(t_p) (F12.6 U and \|u\|, E13.6 ṁ, ES16.8 t): 0.1235 of the impulse Σ ṁ g t = 3.754931e-7 N; the impulse exceeds twice the tolerance | +1.542e-2 of the impulse |
| R4 | both runs end normally (`Stop condition`, empty `run_err.txt`) | exact | |

## RED

| mutation | what fails |
|---|---|
| the revert to the last accepted step kept for model 8 alone | R1 on both runs (X = 0.002502 m and 0.002557 m: the cut step's end past the face, with the time latched before it), R2 −1.000 on both (no deposit), R3 −9.840e-1 of the impulse (the residence 0: the source misses ṁ T_stay g, while the exit record carries t) |
| `part%deltaL = 0` after the revert's `updatePart` | R2 −1.000 on both runs; R1 and R3 hold: the Euler deposit is the reverted state's tail, which `updatePart` sets |
