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
