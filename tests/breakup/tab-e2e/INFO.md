# INFO — breakup/tab-e2e (TAB analytic-oscillator validation, e2e, **GREEN**)

## Reference paper
- **Title:** The TAB method for numerical calculation of spray droplet breakup
- **Author(s):** P. J. O'Rourke, A. A. Amsden
- **Year:** 1987 — SAE Technical Paper 872089
- **Tag:** `[ORA87]` (`papers/TAB-OrourkeAmsden1987.pdf`)
- **Reproduced result (Tier V):** the damped-oscillator deformation closed form
  (eq. 5): `y(t) = We_Cr·(1 − e^{−t/t_d}(cos ωt + sin ωt/(ω t_d)))` from rest,
  `We_Cr = We_r/12`, `1/t_d = C_d μ_l/(2ρ_l r²)` (C_d=5), `ω² = C_k σ/(ρ_l r³) − 1/t_d²`
  (C_k=8) — and its two testable consequences: the **onset** `We_r = 6` (max y =
  2·We_r/12 < 1 below it; the classic diameter-based `We_crit = 12`) and the
  **first-breakup time** `t_bu` (first `y = 1` crossing).

## Case construction
- `tools/make_pe_case.py --we-convention rad --u-gas 50 --lx 0.6 --nx 240` +
  25-point radius-based sweep `We_r ∈ {3 … 90}` dense around the onset. Slip = 25 m/s
  (`κ_v = 0.5`), σ=0.072, water drops (ATLAS-GPB `phase.txt`/`properties.dat`,
  ρ=1000/cp=4182). Box lengthened to 0.6 m (240 cells) so `t_bu` (1.7–17.5 ms) spans
  17–150 recording cells inside the 24 ms residence; slip drift over `t_bu` ≤ 2 %,
  `ω t_d ≥ 118` (near-undamped first crossing).
- **Weber convention:** TAB is **radius-based** in production (`Lib_Breakup.f90` `TABmodel`)
  and in ORA87 — this case and its kernel use the same convention (PE87 is
  diameter-based; do not mix).
- TAB constants at the ORA87 defaults: `Comega(C_k)=8`, `Cmu(C_d)=5`, `WeCrit=6`
  (doubled internally to 12, radius form); Rosin-Rammler product sampling
  (`method=2`, `n=3.5`, seeded RNG → deterministic).

## What `check.py` gates
- **A (onset, no-break):** `We_r ≤ 5.75` → `d` exactly constant.
- **B (must break):** `We_r ≥ 7` → a breakup event in-domain (analytic `t_bu` well
  inside the residence for the whole sweep).
- **C (Tier-V timing):** analytic `t_bu` inside the measured first-breakup bracket
  (recorded rows around the first `d` drop, `t = Σ dx/ū`), slack 5 % (slip drift +
  production's undamped-amplitude crossing detection + apply granularity).
- `5.75 < We_r < 7`: onset band, informational (drift-sensitive).

## Result (GREEN)
Gate passes: no-break onset clean (`We_r ≤ 5.75`), 16/16 supercritical drops break,
and the analytic `t_bu` lands inside every measured bracket to **sub-percent**
(e.g. `We_r=90`: analytic 17.481 ms vs bracket [17.471, 17.571]).

## What the gate guards: the event-breakup path
The event path has five coupled pieces. Without piece 1 no breakup event is ever called:
every TAB/ETAB run is physically breakup-free and KHRT spawns no children. Without piece 2
TAB/ETAB segfault; without 3 the oscillator runs on the wrong clock; without 4 it rolls back to
a stale value; without 5 the next segment reverts the breakup. Pieces 3 and 5 do not reach
KHRT's breakup (it takes the time, not the interval, and derives `d` from its ODE state). The
unit tests (`test_breakup_tab`, `test_tab_moments`) stay green without any of the five — they
call `breakupEvent` directly:
1. `assign_group2particle` copies `brkupEvent` to the particles: the only consumer
   (`eventType = part%brkupEvent`) otherwise reads the default `.false.`, and
   `breakupEvent` is never called.
2. `breakupEvent`'s non-optional `intent(out)` `addChild`/`childState` are integrate's own
   locals (`addChildLocal`, `childState`): handed an ABSENT optional dummy (a childless call
   site), TAB/ETAB **segfault**.
3. The oscillator advances on the accepted-step interval `x-xold` (the Hairer completed-step
   interval), not on `deltat` (the outer step cap).
4. `oldEvLocal = eventLocal` right after the pack at segment start (as for `oldLocal`);
   without it a first-`solout` crossing rolls the oscillator back to a STALE prior-segment
   value.
5. The per-droplet diameter/mass are frozen in `auxLocal` at injection (TAB = model 1,
   constant `d`), so the event finalize resyncs `auxLocal(ind_d)`/`(ind_m)` and `part%m`;
   without it the next segment recomputes `d` from the stale aux and REVERTS the breakup
   (the resync also keeps post-breakup drag consistent with the smaller, more numerous drops).
All five sit behind `eventType`/`eventLocal`: cases without a breakup event never reach them.

## Comparison plot (`verify.py` → `OUTPUT/tab-e2e.svg`, non-gating)
`t_bu(We_r)` analytic curve + shaded no-break region + IGLOO's measured brackets
(error bars), which ride the analytic curve to sub-percent.

## Spray-level rejection
As for the other breakup cases: spray penetration / dense-spray SMD validations are
out of scope by design (steady one-way carrier gas cannot entrain; no atomizing
nozzle inlet; SMD-on-axis is a dense-spray collective). Revisit only if a
transient/two-way IGLOO mode is added.
