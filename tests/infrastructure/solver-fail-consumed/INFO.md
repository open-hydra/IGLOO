# INFO — infrastructure/solver-fail-consumed (`err < 0` exit hands the remnant to the gas, e2e)

## Reference
Not a physics case: no literature reference, no oracle curve. The physics is
`evaporation/tc-hexadecane`'s (Tonini–Cossali n-hexadecane drop) and is gated there. This
case gates **one branch of `Lib_Integration`'s exit contract** — see
[BUGS.md](../../BUGS.md) §E row **O30**, step 4b.

## Why it exists
`integrate` has two exits that both mean "this parcel died in flight", and they disagreed:

| exit | condition | `consumed` for models 2/5 |
|---|---|---|
| non-finite state | `any(y/=y)` after the solve | **set** — "remnant transfers to the gas, as on a clean burnout" |
| solver failure | `Run_ODESolver` returns `err < 0` | **not set** — the remnant silently left the balance |

`consumed` is read by the source block right after `ODEsystem` returns: it zeroes `massOut`,
so the segment-start mass is deposited in the cell instead of flowing out. Without it a
deleted parcel's remaining mass vanishes from the simulation — conservation is broken, and
in hydra-MI2 the closed-domain audit (`hydra_MI2.f90`, `massGas0` vs the end-of-run gas
mass) would read that as a coupling bug it cannot distinguish from a real one.

**No other fixture reaches the `err < 0` branch.** Before O30 the boiling/NaN parcels of
`tc-hexadecane` took the *non-finite* exit; after the `rhsEvaporation` scrub they burn out
cleanly through the normal branch. So this case is the branch's only coverage.

## Fixture
`INPUT/` symlinks `evaporation/tc-hexadecane/INPUT` whole. `input.ini` is tc-hexadecane's
own file with **one change**: the ODE tolerances go to `1e-30`.

    [IGLOO-ODE]
    ode-solver   = H-sdirk4
    relative-tol = 1e-30
    absolute-tol = 1e-30

`1e-30` is below anything double precision can resolve, so SDIRK4 cannot meet it on any
step and OSlo gives up immediately. The failure is **structural, not marginal** — no
compiler, optimisation level or OSlo revision can satisfy that tolerance — so the path
fires deterministically rather than by luck of a threshold.

Consequence, and it is the point: all 25 parcels die **at injection, carrying their full
mass**. The source total is therefore all-or-nothing, which makes the gate maximal.

⚠ This also shows what the contract costs. The deposited remnant is not necessarily small:
here it is **100 %** of the injected stream, dumped into the injection cells (spread over
225 of 1500 cells by the default 8-pass mollification). The alternative contract — drop the
remnant — is the single line `if (mod_model==2 .or. mod_model==5) consumed = .true.`
removed from the `err < 0` branch, and would make this gate read exactly 0.

## What it verifies

| gate | assertion |
|---|---|
| S1 | all 25 parcels took the `err < 0` exit — a **positive** witness that the branch ran |
| S2 | none took the non-finite exit (the other branch, gated in `tc-hexadecane` as H1) |
| S3 | 25 exit rows in `OUTPUT/outloc-A.dat` |
| S4 | `Σ wdot` over `source.tec` equals `Σ` outloc column 7 (`npdot·m_inj`) to `1e-5` |

S1 is deliberately positive. An oracle that only asserted the mass balance would pass on a
run where nothing failed at all, which is the whole point of the fixture.

## Red–green evidence (2026-09-22, release build, sprop2)

| | `bin/IGLOO.pre-C1` (no step 4b) | with step 4b |
|---|---|---|
| S1 | PASS — 25 `Run_ODESolver err=` lines (the branch fires either way) | PASS |
| S4 | **FAIL** — `Σ wdot = 0.000000e+00`, resid **1.000** | PASS — `1.545455e-02 kg/s`, resid **2.941e-07** |

`rc = 0` in both cases: the solver reports nothing wrong.

## Pass criterion
All four gates PASS. Tolerance `1e-5` sits 34× above the measured floor (`2.941e-7`, the
`E13.6` source print; identical at `OMP_NUM_THREADS` 1, 2 and 5) and five orders below the
RED signature. No tuning is involved — S1–S3 are exact integer counts.
