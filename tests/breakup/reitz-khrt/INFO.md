# Reitz-KHRT (brkupSelect=3) — KH-path unit family

T-KHRT, unblocked 2026-07-07 by the A2 production fix (`omegaKH` was missing the
`sqrt` with `r³` inverted — dimensionally m⁶/s² (bug A2), and
`breakup/LITERATURE_TESTS.md` §5).

Oracle: the Reitz 1987 KH linear-stability chain (dimensionless growth-rate and
wavelength forms), assembled into the pure-KH `npdotDot` through the public
`breakupOde` entry with `acc=0` (RT branch inert: `force=0 → λ_RT=∞`).

| ID | check |
|---|---|
| KH1 | sub-critical `WeGas ≤ WeLimit` → rate exactly 0 |
| KH2 | moderate We (u=100 m/s, WeGas≈8.3, Re=1e4): code ≡ paper chain to 1e-12 |
| KH3 | high We (u=300 m/s, WeGas≈75, Re=1e3): code ≡ paper chain to 1e-12 |

Not covered here (deferred to a KHRT-specific pass, see LITERATURE_TESTS §5):
Ω_RT/Λ_RT/RT-trigger, child-drop cubic, the commented-out KH mass-shed branch,
and the `breakupEvent` (`ReitzKHRTevent`) path.

## test_khrt_interaction — the KH shed under evaporation (KS legs)

Under ODE model 4 the event variables carry `n0`, the number rate right after the last event; the
stripped mass flow is `(ṅ − n0)·m` (evaporation lowers `m`, never `ṅ`). One `breakupEvent` call per
leg on a 100 µm drop at 200 m/s slip (`We_g` = 33.3, `(d_stable/d)³` = 8.4e-3), `acc = 0` (RT inert),
`n0` = 1e6 /s:

| ID | state | oracle |
|---|---|---|
| KS1 | `ṅ = n0`, 10 % of the drop mass evaporated | no event, the event variables unchanged (exact) |
| KS2 | `ṅ = 1.1 n0`, 5 % of the stream mass evaporated | child flow `(ṅ − n0)·m` = 0.0864 `n0 m0`, `ṅ ← n0`, `m0 ← m`, `n0` kept (1e-14) |
| KS3 | the KS2 state in the model-3 layout (no `n0` slot) | the per-drop-deficit trigger and partition recomputed in the test (1e-14): child flow 0.1295 `n0 m0`, `ṅ ← 0.95 n0`; the result printed in hex |
| KS4 | 10 % of the stream mass evaporated, `ṅ = n0 (1 + 0.95 (d_stable/d)³)`, `mShedLim` 1e-3 | the trigger passes and `n_child` = 0.95 `n0` lies in `[ṅ m/m0, n0)`: no product parcel (exact) |

Falsification (measured, the per-drop-deficit shed): KS1 sheds 0.09 `n0 m0` of evaporated mass and
cuts `ṅ` to 0.9 `n0`; KS2's child carries 0.1295 `n0 m0` (1.5× the stripped mass); KS4 sheds.
