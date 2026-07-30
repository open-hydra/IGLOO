# INFO — breakup/pilch-erdman-e2e (B-VAL-1: paper-reproduction validation, e2e)

## Reference paper
- **Title:** Use of breakup time data and velocity history data to predict the maximum
  size of stable fragments for acceleration-induced breakup of a liquid drop
- **Author(s):** M. Pilch, C. A. Erdman
- **Year:** 1987 — *Int. J. Multiphase Flow*, 13(6), pp. 741–757
- **DOI:** 10.1016/0301-9322(87)90063-2 · **Tag:** `[PE87]` (`papers/PilchErdman1987.pdf`)
- **Reproduced result:** the dimensionless **total breakup time correlation `T*(We)`**
  (Fig. 7 / p. 748; five regimes, continuous at 18/45/351/2670) — a closed-form model
  output, coded directly (no digitization).

## Validation vs verification (why this case exists)
Tier **P** (paper-reproduction): `T*(We)` is an empirical correlation with no analytic
ODE solution across the regimes. The gate compares IGLOO's *integrated* output against
the **paper's published correlation** — not our re-coding of production (the A13/A16
lockstep trap) and not experiment. This case found **A16** (the `45<We≤351` branch coded
`+0.25` instead of PE87's `-0.25`, breakup ~10× too slow) and, once A16 was fixed,
immediately surfaced **A17** (all `[IGLOO-Models]` numeric breakup params zeroed at
runtime by the per-material `assign_breakup` re-allocation; **FIXED 2026-07-22** —
`bp` now allocated+filled once in `read_models()`).

## Case construction
- `tools/make_pe_case.py` writes `INPUT/solfile.tec` (uniform gas box 0.15×0.05×0.05 m,
  60×5×5 cells, `U=200 m/s`, `ρg=1.2`, `T=300 K`) and `INPUT/bc.txt`: the 25 face-1
  inlet cells (5×5) each inject drops of a DIFFERENT diameter — a Weber sweep
  `We ∈ {20…1000}` in one run. All drops share `κ_v=0.5` (`v0=100` ⇒ slip=100 m/s) and
  `σ=0.072`, so `We = d·ρg·slip²/σ` varies only through the per-cell `d`
  (0.12–6.0 mm). Drops are injected SLOWER than the carrier (shock-tube analog: the
  drop accelerates from near-rest, PE87's velocity-history premise).
- `INPUT/phase.txt` + `INPUT/properties.dat` are **ATLAS GPB artifacts** (regenerated
  byte-identical from `[GPB-Phase1]` with `material = A`, fixed ρ=1000/cp=4182); GPB
  does not touch `bc.txt`, so the per-cell diameter sweep is preserved.
- **Weber convention:** Pilch-Erdman is **diameter-based** (`Lib_Breakup.f90:118`);
  `check.py` uses the same convention. (Reitz-Diwakar is radius-based — do not reuse
  this kernel there.)

## What `check.py` gates (IGLOO vs the paper's correlation)
IGLOO's **initial breakup rate at the injection We** — the leading-window slope of the
recorded `d(t)` (`t` reconstructed as `Σ dx/ū`) — against the PE87 forward rate
`dd/dt = (d_stable−d)/τ`, `τ = T*·d/(slip·√(ρg/ρp))`, at the same `(d0, slip0)`.
This isolates `T*(We)` at a single known We per drop (no regime-crossing). The
reference applies the SAME leading-window LSQ to the kernel's RK4 samples, so the
finite-window bias of the decelerating rate cancels on both sides.
- **Gate 1 (rate):** all 22 drops with `We ≥ 45`, tol **1 %** (observed agreement
  ~1e-4 post-A16+A17 fixes — 100× margin), covering both the `-0.25` (45–351, the
  A16 branch) and `+0.25` (351–1000) regimes.
- **Gate 2 (d_stable, B-VAL-2):** plateau `d_end` vs PE87's closed form at initial
  conditions, tol **12 %** — see the B-VAL-2 section below.
- `We < 45` drops are excluded (break up inside the measurement window; the 18–45
  branch is covered by the `test_breakup_pe` unit test).

## Comparison plot (`verify.py` → `OUTPUT/pilch-erdman-e2e.svg`, non-gating)
PE87 `T*(We)` curve (log-x, Fig. 7 layout) vs IGLOO per-drop **implied `T*`** (the
measured initial rate inverted through the paper model); all `We≥45` points gated.
No experimental shock-tube scatter (Fig. 7 not digitized; the gate is the closed
form).

## Tier-2 (spray-level) documented rejection
Spray penetration `S(t)` / dense-spray axial SMD validations are **out of scope by
design**: they are governed by gas entrainment and two-way momentum coupling, need an
atomizing nozzle inlet + turbulent dispersion, and SMD-along-axis is a dense-spray
collective observable. IGLOO's carrier gas is steady and one-way and cannot entrain.
Revisit only if a transient/two-way IGLOO mode is added.

## B-VAL-2 — maximum stable fragment `d_stable(We)` (GATED since the A17 fix)
The paper's stated purpose is the closed form
`d_stable(We0) = Wec·σ/(ρg·U0²·(1−VdV)²)` at initial conditions. Post-A17-fix, drops
with `We0 ≤ 60` plateau in-domain (correct `bp` makes low-We breakup stop earlier at
larger `d`); `check.py::bval2_gate` gates the measured plateau against the closed
form, tol 12 % — observed 0.91–1.07. The tolerance covers the asymptotic-approach
tail plus the convention difference (production re-evaluates `d_stable` with the
*current* slip along the path; the paper folds the velocity history into `(1−VdV)`
at the *initial* `U0`). Discrimination: the pre-fix zero-bp state read 1.48–1.50
(printed as a diagnostic column), ~4× outside tolerance; the A16 sign bug far worse.

## Result (post-A16 + post-A17, 2026-07-22)
Both gates GREEN: rate 22/22 within 1 % (worst 1.1e-4); d_stable 5/5 within 12 %.
Implied `T*` rides the Fig. 7 curve including the V-minimum at the We=351 boundary.
