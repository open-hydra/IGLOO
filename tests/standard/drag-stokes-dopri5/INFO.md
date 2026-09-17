# drag-stokes-dopri5 — `ode-solver = H-dopri5` on the Stokes closed form — **registered DISABLED (ledger O25)**

`drag-stokes` with the other ODE solver IGLOO exposes (`H-dopri5`, explicit Dormand–Prince 5(4)
via OSlo `wrap_dopri5`); `INPUT/` and `check.py` are symlinked, only the integrator differs.
Every other fixture runs `H-sdirk4`, so until 2026-09-17 the token was accepted, wired, and
exercised by nothing (coverage audit G11 — the 09-14 audit first called it a dead token; it is
not, it is worse: it is live and broken).

**RED at OSlo `fb8255d`.** Every parcel dies at injection: `MORE THAN NMAX = 100000 STEPS ARE
NEEDED`, `Run_ODESolver err=-2 … marking gone`, 25/25, at every tolerance tried (1e-11, 1e-8,
1e-5). Cause: the OSlo fork patched DOPRI5's *first* `SOLOUT` call to IGLOO's six-argument
form (`lib/OSLO/lib/hairer/dopri5.f:401`) but left the *per-step* call at `:512` in Hairer's
original eleven-argument form, so IGLOO's `solout(NR,XOLD,X,Y,N,IRTRN)` receives `CONT` where it
expects `IRTRN` — `CONT`'s base address, an unused slot inside the wrapper's `WORK` (`NRDENS = 0`),
so the writes are benign and DOPRI5's own `IRTRN` stays at the 1 set before the first (six-argument)
call: it never interrupts, and every call runs to `XEND` (= 20·Δt₀; the `EXIT OF DOPRI5 AT X =
0.03752` line is label 78's NMAX report, not an interrupt). NMAX itself comes from IGLOO's
crossing-refinement loop (`Lib_Integration.f90:549-581`): `IamOut`/`newGas` are still set by
`solout`, so IGLOO rewinds to the last interior latch and shrinks Δt (= DOPRI5's HMAX) against a
fixed `t2` until (t2 − t)/hmax exceeds 10⁵ — the failure mode the comment at `:566-569` names.

**Fix verified (one line, in OSlo, not in this repo).** With `:512` changed to
`CALL SOLOUT(NACCPT+1,XOLD,X,Y,N,IRTRN)` in a throwaway build: 0 NMAX exits, 1500 rows, 25/25
parcels match the Stokes closed form (worst resid/tol 0.014). The patch and the before/after
logs are in the 2026-09-17 session scratchpad (`dopri5/`). Two things to carry into the OSlo fix:
label 79 prints `EXIT OF DOPRI5 AT X=` on EVERY interrupt (7850 lines in the patched run's log
for 1500 trajectory rows) — silence the `IRTRN < 0` branch, as the SDIRK4 fork already does; and
hydra's own `lib/OSlo` gitlink is `2127869`, an ancestor that predates `dopri5.f` altogether
(only its dirty checkout is at `fb8255d`), so both gitlinks must move. The case is registered
`DISABLED` (ctest reports it "Not Run" on every run, a visible reminder, never a failure): delete
the `set_tests_properties(... DISABLED TRUE)` line when the gitlinks and the vendored copy carry
the fix; it then runs drag-stokes's own oracle. `H-dopri5` is otherwise still advertised by the
registry and accepted by the parser — a production run selecting it silently loses every parcel
with exit 0 (the run "completes" and writes its files).
