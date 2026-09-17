# drag-stokes-dopri5 — `ode-solver = H-dopri5` on the Stokes closed form

`drag-stokes` with the other ODE solver IGLOO exposes (`H-dopri5`, explicit Dormand–Prince 5(4)
via OSlo `wrap_dopri5`); `INPUT/` and `check.py` are symlinked, only the integrator differs.
Every other fixture runs `H-sdirk4`, so until 2026-09-17 the token was accepted, wired, and
exercised by nothing (coverage audit G11 — the 09-14 audit first called it a dead token; it was
not, it was worse: live and broken).

**RED at OSlo `fb8255d` (ledger O25).** Every parcel died at injection: `MORE THAN NMAX = 100000
STEPS ARE NEEDED`, `Run_ODESolver err=-2 … marking gone`, 25/25, at every tolerance tried (1e-11,
1e-8, 1e-5). Cause: the OSlo fork's dopri5 import (`26383fa`) rewrote DOPRI5's *first* `SOLOUT`
call to the six-argument `solout_if` (`lib/hairer/dopri5.f:401`) but left the *per-step* call at
`:512` in Hairer's original eleven-argument form, so IGLOO's `solout(NR,XOLD,X,Y,N,IRTRN)`
received `CONT` where it expects `IRTRN` — `CONT`'s base address, an unused slot inside the
wrapper's `WORK` (`NRDENS = 0`), so the writes were benign and DOPRI5's own `IRTRN` stayed at the 1
set before the first call: it never interrupted, and every call ran to `XEND` (= 20·Δt₀). NMAX
itself came from IGLOO's crossing-refinement loop (`Lib_Integration.f90`): `IamOut`/`newGas` were
still set by `solout`, so IGLOO rewound to the last interior latch and shrank Δt (= DOPRI5's HMAX)
against a fixed `t2` until (t2 − t)/hmax exceeded 10⁵.

**Fixed in the OSlo fork, branch `fix/dopri5-solout` (2026-09-17), on top of `fb8255d`:** `:512`
calls `SOLOUT(NACCPT+1,XOLD,X,Y,N,IRTRN)`, and label 79 — reached only by a `SOLOUT`-requested stop
(`IRTRN < 0`, a success, `IDID = 2` per Hairer's header) — no longer prints `EXIT OF DOPRI5 AT X=`
(7850 lines per 1500 trajectory rows before; the SDIRK4 fork had already made the same change).
Both OSlo checkouts (IGLOO's `lib/OSLO`, hydra's `lib/OSlo`) must sit on that commit; IGLOO's
gitlink moved with it. With the fix: 0 NMAX exits, 1500 rows, 25/25 parcels match the Stokes
closed form (worst resid/tol 0.014), and the interim parse-time refusal (`refuse-dopri5`) and the
`DISABLED` flag this case carried went with it — this case is the gate now.

**Second defect on the same path (ledger O27, fixed in the fork's `fd9d178`).** The E A/B found
this case's `scatter-A.dat` differing between two runs of one binary at OMP 5 (372 / 366 / 348
records; identical at OMP 1). `dopri5.f` keeps the previous grid point `XOLD` in `COMMON
/CONDO5/` (for its dense-output routine), shared by every OpenMP thread, and IGLOO's `solout`
accumulates the scatter weight on `x − xold` — so a parcel's cloud was timed by another
thread's step. The breakup oscillator advance in `solout` reads the same interval (no fixture
runs breakup on DOPRI5). SOLOUT now receives a thread-private copy; the twin
`repeatability/drag-stokes-dopri5` (two sweeps at OMP 5, scatter multiset sweep 1 ≡ sweep 0) was
3/3 RED before and is 3/3 green after, and OMP 5 reproduces the OMP 1 cloud (376 rows).

**What it pins.** That the explicit solver's cell-crossing interrupt reaches the integrator (the
O25 failure mode: a run that "completes" with every parcel lost and exit 0). It does not compare
the two solvers against each other beyond the shared oracle; `H-sdirk4` remains the default and
the stiff-safe choice.
