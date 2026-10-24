# drag-stokes-dopri5 — `ode-solver = H-dopri5` on the Stokes closed form

`drag-stokes` with the other ODE solver IGLOO exposes (`H-dopri5`, explicit Dormand–Prince 5(4)
via OSlo `wrap_dopri5`); `INPUT/` and `check.py` are symlinked, only the integrator differs.
`breakup/etab-e2e-dopri5` and the repeatability twins below run it too; every other fixture runs `H-sdirk4`.

**The DOPRI5 the case needs.** Both `SOLOUT` calls of `lib/hairer/dopri5.f` take the six-argument
`solout_if`, per step `SOLOUT(NACCPT+1,XOLDS,X,Y,N,IRTRN)`. **RED** on a DOPRI5 whose first call
is the six-argument form and whose per-step call keeps Hairer's original
eleven-argument form: every parcel dies at injection (`MORE THAN NMAX = 100000 STEPS ARE NEEDED`,
`Run_ODESolver err=-2 … marking gone`, 25/25) at every tolerance tried (1e-11, 1e-8, 1e-5). IGLOO's
`solout(NR,XOLD,X,Y,N,IRTRN)` then receives `CONT` where it expects `IRTRN` — `CONT`'s base address, an
unused slot inside the wrapper's `WORK` (`NRDENS = 0`), so the writes are benign and DOPRI5's own
`IRTRN` stays at the 1 set before the first call: it never interrupts, and every call runs to `XEND`
(= 20·Δt₀). NMAX itself comes from IGLOO's crossing-refinement loop (`Lib_Integration.f90`):
`IamOut`/`newGas` are still set by `solout`, so IGLOO rewinds to the last interior latch and shrinks
Δt (= DOPRI5's HMAX) against a fixed `t2` until (t2 − t)/hmax exceeds 10⁵.

Label 79 — reached only by a `SOLOUT`-requested stop (`IRTRN < 0`, a success, `IDID = 2` per Hairer's
header) — prints nothing, as in the SDIRK4 fork (Hairer's original prints `EXIT OF DOPRI5 AT X=` at that
label: 7850 lines per 1500 trajectory rows of this case). Both OSlo checkouts (IGLOO's `lib/OSLO`, hydra's
`lib/OSLO`) need this DOPRI5. Measured: 0 NMAX exits, 1525 rows, 25/25 parcels match the Stokes
closed form (worst resid/tol 0.014).

**The thread-private `XOLD`.** `dopri5.f` keeps the previous grid point `XOLD` in `COMMON /CONDO5/`
(for its dense-output routine), and IGLOO's `solout` accumulates the scatter weight on `x − xold`;
SOLOUT receives a thread-private copy. With `/CONDO5/` shared by every OpenMP thread a parcel's
cloud is timed by another thread's step: this case's `scatter-A.dat` differs between two runs of one
binary at OMP 5 (372 / 366 / 348 records; identical at OMP 1), and the breakup oscillator advance in
`solout`, which reads the same interval, makes breakup on DOPRI5 physically wrong under OpenMP. Two
twins gate it at OMP 5 — `repeatability/drag-stokes-dopri5` (scatter multiset sweep 1 ≡ sweep 0; 3/3
RED on the shared copy, 377 vs 345 records) and `repeatability/tab-dopri5` (`tab-e2e` on H-dopri5; 3/3
RED on the shared copy with 22/25 exits and 2855 vs 3783 trajectory records differing between sweeps) —
both 3/3 green, and OMP 5 reproduces the OMP 1 cloud (376 rows).

**What it pins.** That the explicit solver's cell-crossing interrupt reaches the integrator (the
failure mode above: a run that "completes" with every parcel lost and exit 0). It does not compare
the two solvers against each other beyond the shared oracle; `H-sdirk4` is the default and
the stiff-safe choice.
