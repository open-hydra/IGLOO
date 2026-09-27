# INFO — solidification (ODE model 6)

Verifies the solidification track: a molten droplet (alumina) cooling through its melting point,
**supercooling** to a nucleation temperature `T-nuc`, freezing abruptly with an adiabatic
**recalescence** back to `T-melt`, holding a **plateau** at `T-melt` while the rest of the latent heat
leaves by convection, then cooling as a **solid**; a solid heated to `T-melt` in a hotter gas **melts** on
the same plateau and continues as a liquid. Selected per material by the `solidification=on`
token on the material line of `INPUT/phase.txt` (with `h-fus`, `cp-solid`, optional `T-melt`,
`T-nuc`); closures in `IGLOO_Lib_Solidification`.

Physics, with `m` constant, `Q̇ = π d k_g Nu (T_g − T_p)` the convective heat received and `f` the
frozen mass fraction:

| phase | ODE | event that ends it |
|---|---|---|
| liquid / undercooled | `m c_l dT/dt = Q̇`, `df/dt = 0` | `T = T-nuc`: `f0 = c_l (T_m − T_nuc)/h_fus`; `f0 < 1` ⇒ plateau at `T_m` with `f = f0`; `f0 ≥ 1` ⇒ solid, `f = 1`, at `T_nuc + L(T_nuc)/c_s` with `L(T) = h_fus − (c_l − c_s)(T_m − T)` |
| plateau | `dT/dt = 0`, `m h_fus df/dt = −Q̇` | `f = 1` ⇒ solid at `T_m`; `f = 0` ⇒ liquid at `T_m` (molten in hot gas) |
| solid | `m c_s dT/dt = Q̇`, `df/dt = 0` | `T = T_m` ⇒ plateau at `T_m` with `f = 1 − c_s (T − T_m)/h_fus` (melting in hot gas) |

The specific enthalpy `h = c_l T + h_off` (liquid), `c_l T_m + h_off − f h_fus` (plateau),
`c_l T_m + h_off − h_fus − c_s (T_m − T)` (solid) is continuous across every event, so the
recalescence is adiabatic and the latent heat passes between droplet and gas only through the plateau's
`Q̇`. A step crosses a threshold when it ends strictly past it from a start at or before it; a threshold
found strictly passed at a segment start is applied there; every transition lands inside the new phase (a
freezing plateau at `T ≤ T_m`, a melting solid at `f ≤ 1`).

## Tests

Unit family `test_solidification.f90`: an alumina-like particle, `ρ = 3970`, `d = 50 µm`, `c_l = 1888`,
`c_s = 1420`, `h_fus = 1.09e6`, `T_m = 2327`, `T_nuc = 0.8 T_m`, `T_0 = 2500`, `T_g = 1500`, `Nu = 2`,
`k_g = 0.1` (`f0 = 0.806`, `t_n = 1.59e-2 s`, `t_plat = 2.11e-3 s`). SG7 heats the same particle from a solid
at 1500 K in a 3000 K gas (`T_m` at `t_m = 9.41e-3 s`, molten at `2.28e-2 s`).

| id | mode | locus | compared | tolerance | RED (mutation of the closure, measured) |
|---|---|---|---|---|---|
| SG0 | pin | metal slots `imTm/imHf/imTn/imCps = 7/8/9/10`, `nmp = nmetal = 11` | constants | exact | measured: `imTm`/`imHf` swapped fails SG0 alone |
| SG1 | analytic | `nucleationJump`: (a) `f0 = c_l (T_m − T_nuc)/h_fus = 0.806`, `T ← T_m`, plateau; (b) `h_fus = 4e5` (`f0 = 2.20`): solid, `f = 1`, `T = T_nuc + L(T_nuc)/c_s = 1989.9 K` | independent transcription | (a) 4 ULP of `c_l T_m/h_fus`; (b) 8 ULP of `c_l T_m/c_s` | measured: (a) the jump taken with `c_s` gives `f0 = 0.6063` (0.248 relative, as predicted), SG3 gap 2.18e5 J/kg, SG4 off by 6.0e-2 on `T` and 0.20 on `f`; (b) the single-capacity form `T_nuc + h_fus/c_l` gives 2073.46 K (83.6 K off), SG3 gap 1.19e5 J/kg (2.7e-2 of `c_l T_m`), SG5 whole freeze fails too |
| SG2 | analytic + RK4 | `plateauRate = −Q̇/(m h_fus)`: sign and magnitude; RK4 (4000 steps) from `f0` reaches `f = 1` at `t_plat`; heat received `= −m h_fus (1 − f0)` | closed form | 1e-14 / 1e-10 / 1e-12 relative | measured, rate sign flipped: SG2a 2.0, `f(t_plat)` off by 0.39, SG4 off by 0.44 on `T` and 1.0 on `f` |
| SG3 | analytic | `hSolid` continuous at nucleation (both branches), plateau end, re-melt, and through `solidTransition` from states past each threshold (`f = 1.001`, `f = −0.001`, `T = T_nuc − 5`, solid at `T = T_m + 5`); liquid branch `c_l T + h_off` to 1 ULP with `h_off = −1.7e7`; increasing in `T` within a phase; the melt from `T_m + 5` lands on the plateau at `T_m` with `f = 1 − 5 c_s/h_fus = 0.993486` (SG3d) | algebra | 4 ULP of `c_l T_m + |h_off|` | measured: the plateau branch without `− f h_fus` opens gaps of 8.79e5 J/kg at nucleation (0.250 of `c_l T_nuc`, as predicted) and 1.09e6 J/kg at the plateau end; a clamp on `f` opens 1.09e3 J/kg past each end; the single-capacity whole-freeze temperature 1.19e5 J/kg; without the melt transition SG3d fails (the state stays solid at 2332 K) |
| SG4 | independent integration | RK4 (`h = 1 µs`) of the three regimes with each phase end located by bisection on its event function and the transition applied at the crossing, vs the piecewise closed form at 14 liquid, 12 plateau and 14 solid epochs; the 12 plateau epochs at `T_m` exactly | closed form | 1e-10 relative on `T`, 1e-12 on `f` | measured, no nucleation event: `T` off by 0.218 (predicted 0.21 at the plateau midpoint), no epoch at `T_m`, `f` off by 1.0. The test applies its own event rule; the solver's is gated end to end |
| SG5 | corners | finite outputs on `f ∈ {−1e-3, 0, 1, 1+1e-3}` × `T ∈ {T_g, T_nuc, T_m, T_0}` × `Q̇ ∈ {−1, 0, 1} W`; event functions exactly 0 at `T_nuc`, `f = 1`, `f = 0` and `T_m` (the solid's) and of the right sign on either side; `hSolid` linear in `f` past both ends (no clamp on `f`); injection phases (solid at `T ≤ T_nuc`, undercooled below `T_m`, liquid from `T_m`); whole freeze with `L(T_nuc) < 0` (`h_fus = 1e5`) lands below `T_nuc` with the enthalpy kept | definition | bitwise / 4 ULP | measured, clamp on `f` in the plateau enthalpy: SG5c fails (with SG3); the solid without an event: SG5b fails |
| SG6 | production RHS | one `rhsSolidification` call per phase after `setupRHS(6, …)`: the mass, diameter, metal block and phase slots exist; `F(7) = Q̇/(m c)` (liquid, solid), `F(8) = −Q̇/(m h_fus)` (plateau), the unweighted euler tail; poisoning the nine metal slots the RHS must not read (all but `h-fus`, `cp-solid`) changes no bit of `F` | closed form at zero slip (`Nu = 2`) | 1e-13 relative; bitwise | measured: no metal block for model 6 fails SG6a; the plateau rate taking `cp-solid` fails SG6b by 766× |
| SG7 | independent integration | SG4's RK4 and event rule from a solid at 1500 K in a 3000 K gas: solid heating to `T_m` at `t_m = 9.412991e-3 s`, the melting plateau (`f` from 1 to 0 at 74.65 /s) to `2.280855e-2 s`, then liquid heating; vs the piecewise closed form at 14 solid, 12 plateau and 14 liquid epochs; the 12 plateau epochs at `T_m` exactly | closed form | 1e-10 relative on `T`, 1e-12 on `f` | measured, no melt event: a solid at 2493.9669 / 2619.5104 / 2713.9074 K where `f` should be 3/4, 1/2, 1/4, 2784.8851 K at the end of the melt, `T` off by 0.188, no epoch at `T_m`, `f` off by 1.0 |
| SG8 | landing | `solidTransition` from a state rounded past its threshold: freezing from `f = 1 − 1e-15` lands a solid at `T_m` exactly, melting from `T = T_m − 1e-12` the plateau at `f = 1` exactly; the new phase's event function `≥ 0`; the enthalpy kept | definition | bitwise; 4 ULP of `c_l T_m + |h_off|` | measured, landing clamps removed: `T − T_m = 9.1e-13 K` (a solid past its melt threshold) and `f − 1 = 1.1e-15` (a plateau past `f = 1`) |
| SG9 | rule | `eventCrossed(gOld, gEnd)`: (1, −1) and (0, −1) cross, (0, 1), (1, 0) and (−1, −2) do not; on a 5 × 5 grid of event values the 6 crossings interpolate at `s = gOld/(gOld − gEnd)` in [0, 1], at 0 only from a start on the threshold; `eventPassed`: 0 no, −1e-300 yes, 1e-300 no | definition | exact | measured, the non-strict rule (a crossing at `gEnd ≤ 0` from `gOld > 0`, the catch-up at `g ≤ 0`): SG9a and SG9b fail |

`dump_curve('solidification-recalescence', 'SG4', …)` writes the unit figure `unit-solidification-recalescence.svg`.

## End to end

| case | what it gates |
|---|---|
| [`solid-box`](solid-box/INFO.md) | the three regimes against their closed forms on the shared box, the event position, the energy deposit (global telescoping and the per-cell plateau heat) |
| [`solid-melt`](solid-melt/INFO.md) | a solid heated to `T-melt` in a hotter gas: the solid, melting-plateau and liquid regimes against their closed forms, the event positions, the heat taken from the gas per plateau cell, no parcel lost at a phase boundary |
| [`solid-box-euler`](solid-box-euler/INFO.md) | the equivalent-Eulerian field of model 6: coverage, `ρ_p/n_p = m`, `u_p`, `T_p` per cell |
| [`solid-box-2mat`](solid-box-2mat/INFO.md) | a model-6 material followed by a model-1 one in one run |
| [`scatter-weight-solid`](scatter-weight-solid/INFO.md) | the scatter cloud's marker count per parcel across the phase events: an event keeps the scatter weight of the part of the step it keeps |
| `repeatability/solid-box` | two sweeps reproduce each other: the phase and the frozen fraction are re-derived at injection |
| `infrastructure/refusals/solid-*` | the eight setup refusals of the input contract (`h-fus`, `cp-solid`, `T-nuc < T-melt`, no evaporation, combustion or breakup, constant `cp` and density) |
