# INFO — solidification (ODE models 6 and 7)

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

Breakup acts on the liquid only (phases 0 and 1). TAB and ETAB are event models: the material stays ODE model 6
and the breakup event is gated by the phase in force over the step. Pilch-Erdman, Reitz-Diwakar and KH-RT carry
a number-rate ODE and make **ODE model 7**: `Z(9) = ṅ`, the drop mass `m = ṁ/ṅ` from the stream flow,
`d = (6m/(πρ))^{1/3}`, the energy equation of the table above, and `dṅ/dt` the breakup model's rate in the liquid
phases, 0 on the plateau and in the solid. A KH-RT child is born in its parent's phase.

The drop's volume is its liquid part at the liquid density and its frozen fraction at the solid density, at
constant mass: `d = (6m ((1 − f)/ρ_l(T) + f/ρ_s)/π)^{1/3}`. `ρ_l` is the material's density, constant or the
`Density` column of the property table at the drop's temperature; `ρ_s` is the `rho-solid` token of the material
line, the liquid density when absent. The mass changes only at a breakup event.

## Tests

Unit family `test_solidification.f90`: an alumina-like particle, `ρ = 3970`, `d = 50 µm`, `c_l = 1888`,
`c_s = 1420`, `h_fus = 1.09e6`, `T_m = 2327`, `T_nuc = 0.8 T_m`, `T_0 = 2500`, `T_g = 1500`, `Nu = 2`,
`k_g = 0.1` (`f0 = 0.806`, `t_n = 1.59e-2 s`, `t_plat = 2.11e-3 s`). SG7 heats the same particle from a solid
at 1500 K in a 3000 K gas (`T_m` at `t_m = 9.41e-3 s`, molten at `2.28e-2 s`).

| id | mode | locus | compared | tolerance | RED (mutation of the closure, measured) |
|---|---|---|---|---|---|
| SG0 | pin | metal slots `imHf/imCps = 7/8`, `nmp = nmetal = 12` | constants | exact | measured: `imHf`/`imCps` swapped fails SG0, and SG6b by 767× (the plateau rate takes `cp-solid`) |
| SG1 | analytic | `nucleationJump`: (a) `f0 = c_l (T_m − T_nuc)/h_fus = 0.806`, `T ← T_m`, plateau; (b) `h_fus = 4e5` (`f0 = 2.20`): solid, `f = 1`, `T = T_nuc + L(T_nuc)/c_s = 1989.9 K` | independent transcription | (a) 4 ULP of `c_l T_m/h_fus`; (b) 8 ULP of `c_l T_m/c_s` | measured: (a) the jump taken with `c_s` gives `f0 = 0.6063` (0.248 relative, as predicted), SG3 gap 2.18e5 J/kg, SG4 off by 6.0e-2 on `T` and 0.20 on `f`; (b) the single-capacity form `T_nuc + h_fus/c_l` gives 2073.46 K (83.6 K off), SG3 gap 1.19e5 J/kg (2.7e-2 of `c_l T_m`), SG5 whole freeze fails too |
| SG2 | analytic + RK4 | `plateauRate = −Q̇/(m h_fus)`: sign and magnitude; RK4 (4000 steps) from `f0` reaches `f = 1` at `t_plat`; heat received `= −m h_fus (1 − f0)` | closed form | 1e-14 / 1e-10 / 1e-12 relative | measured, rate sign flipped: SG2a 2.0, `f(t_plat)` off by 0.39, SG4 off by 0.44 on `T` and 1.0 on `f` |
| SG3 | analytic | `hSolid` continuous at nucleation (both branches), plateau end, re-melt, and through `solidTransition` from states past each threshold (`f = 1.001`, `f = −0.001`, `T = T_nuc − 5`, solid at `T = T_m + 5`); liquid branch `c_l T + h_off` to 1 ULP with `h_off = −1.7e7`; increasing in `T` within a phase; the melt from `T_m + 5` lands on the plateau at `T_m` with `f = 1 − 5 c_s/h_fus = 0.993486` (SG3d) | algebra | 4 ULP of `c_l T_m + |h_off|` | measured: the plateau branch without `− f h_fus` opens gaps of 8.79e5 J/kg at nucleation (0.250 of `c_l T_nuc`, as predicted) and 1.09e6 J/kg at the plateau end; a clamp on `f` opens 1.09e3 J/kg past each end; the single-capacity whole-freeze temperature 1.19e5 J/kg; without the melt transition SG3d fails (the state stays solid at 2332 K) |
| SG4 | independent integration | RK4 (`h = 1 µs`) of the three regimes with each phase end located by bisection on its event function and the transition applied at the crossing, vs the piecewise closed form at 14 liquid, 12 plateau and 14 solid epochs; the 12 plateau epochs at `T_m` exactly | closed form | 1e-10 relative on `T`, 5e-12 on `f` (the harness's own rounding, the `f` sum over the plateau and the running time sum, is at most 1.3e-12 here and 2.1e-12 in SG7) | measured, no nucleation event: `T` off by 0.218 (predicted 0.21 at the plateau midpoint), no epoch at `T_m`, `f` off by 1.0. The test applies its own event rule; the solver's is gated end to end |
| SG5 | corners | finite outputs on `f ∈ {−1e-3, 0, 1, 1+1e-3}` × `T ∈ {T_g, T_nuc, T_m, T_0}` × `Q̇ ∈ {−1, 0, 1} W`; event functions exactly 0 at `T_nuc`, `f = 1`, `f = 0` and `T_m` (the solid's) and of the right sign on either side; `hSolid` linear in `f` past both ends (no clamp on `f`); injection phases (solid at `T ≤ T_nuc`, undercooled below `T_m`, liquid from `T_m`); whole freeze with `L(T_nuc) < 0` (`h_fus = 1e5`) lands below `T_nuc` with the enthalpy kept | definition | bitwise / 4 ULP | measured, clamp on `f` in the plateau enthalpy: SG5c fails (with SG3); the solid without an event: SG5b fails |
| SG6 | production RHS | one `rhsSolidification` call per phase after `setupRHS(6, …)`: the mass, diameter, metal block and phase slots exist; `F(7) = Q̇/(m c)` (liquid, solid), `F(8) = −Q̇/(m h_fus)` (plateau), the unweighted euler tail; poisoning the seven metal slots the RHS must not read (all but `h-fus`, `cp-solid`) changes no bit of `F` | closed form at zero slip (`Nu = 2`) | 1e-13 relative; bitwise | measured: no metal block for model 6 fails SG6a; the plateau rate taking `cp-solid` fails SG6b by 766× |
| SG7 | independent integration | SG4's RK4 and event rule from a solid at 1500 K in a 3000 K gas: solid heating to `T_m` at `t_m = 9.412991e-3 s`, the melting plateau (`f` from 1 to 0 at 74.65 /s) to `2.280855e-2 s`, then liquid heating; vs the piecewise closed form at 14 solid, 12 plateau and 14 liquid epochs; the 12 plateau epochs at `T_m` exactly | closed form | SG4's: 1e-10 relative on `T`, 5e-12 on `f` | measured, no melt event: a solid at 2493.9669 / 2619.5104 / 2713.9074 K where `f` should be 3/4, 1/2, 1/4, 2784.8851 K at the end of the melt, `T` off by 0.188, no epoch at `T_m`, `f` off by 1.0 |
| SG8 | landing | `solidTransition` from a state rounded past its threshold: freezing from `f = 1 − 1e-15` lands a solid at `T_m` exactly, melting from `T = T_m − 1e-12` the plateau at `f = 1` exactly; the new phase's event function `≥ 0`; the enthalpy kept | definition | bitwise; 4 ULP of `c_l T_m + |h_off|` | measured, landing clamps removed: `T − T_m = 9.1e-13 K` (a solid past its melt threshold) and `f − 1 = 1.1e-15` (a plateau past `f = 1`) |
| SG9 | rule | `eventCrossed(gOld, gEnd)`: (1, −1) and (0, −1) cross, (0, 1), (1, 0) and (−1, −2) do not; on a 5 × 5 grid of event values the 6 crossings interpolate at `s = gOld/(gOld − gEnd)` in [0, 1], at 0 only from a start on the threshold; `eventPassed`: 0 no, −1e-300 yes, 1e-300 no | definition | exact | measured, the non-strict rule (a crossing at `gEnd ≤ 0` from `gOld > 0`, the catch-up at `g ≤ 0`): SG9a and SG9b fail |
| SG10 | layout | `setupRHS(7, …)` for Pilch-Erdman, Reitz-Diwakar and KH-RT: the stream flow, `d`, `ρ`, `σ`, `μ`, the metal block and the phase slot (15 auxiliary slots), the KH-RT timer before the phase slot (3 state slots, 4 event variables; 1 and 0 for the other two); 9 equations, 15 with the Eulerian moments; a solidifying material is model 7 with a breakup ODE, model 6 without | definition | exact | — |
| SG11 | production RHS | one model-7 `rhsSolidification` call per phase and breakup model at a 100 m/s slip (Stokes, Ranz-Marshall), `σ = 5e-3`, `μ = 1e-3`, `Z(9) = 1.5 ṁ/m_0` (rebuilt `d = 43.68 µm` against the injection 50 µm); `F(1:8)` against a model-6 call on the rebuilt `m` and `d`; `F(9)` against a direct `breakupOde` call (the RHS's acceleration, the KH-RT timer started); the Eulerian tail; the timer in the auxiliary state | model 6; the breakup routine | 1e-13 relative on `F(1:8)` (measured 3.6e-16), 1e-12 on `F(9)` (2.4e-16): two compilation units; exact on the zeros, on `F(9) = 0` on the plateau and in the solid (the direct rate non-zero there), on the tail (model 6's shifted by one slot, `F(15) = Z(9)`) and on the timer (written back in the liquid phases, kept in the frozen ones) | measured: the rate not gated by the phase fails SG11d and SG11f; model 7's `d` taken from the injection slot fails SG11a (0.145 = 50/43.68 − 1), SG11b and SG11c |
| SG12 | analytic + production RHS | `rho-solid = 4600` on SG6's drop: `solidDiameter` at `f` ∈ {0, 0.3, 1} against `d0 (1 − f + f ρ_l/ρ_s)^{1/3}`; `setupRHS(6, …)` puts the solid-density slot right after the metal block (the emissivity after it with radiation); one `rhsSolidification` call per phase, `F(7)` or `F(8)` at that diameter, the other zero; a material set up without `rho-solid` after one with it has no such slot (13 auxiliary slots); with `emissivity = 0.5`, the solid's convective plus radiative heat at its diameter | closed form | 1e-14 relative (measured 5.7e-16 on `d`, 4.5e-16 on `F`, 5.3e-16 with radiation); exact on the layout | measured: the slot ignored in the RHS (`ρ_s = ρ_l`) fails SG12c by 5.0e-2; the reset removed fails SG12d; the radiation at the stored diameter fails SG12f by 1.1e-2 |
| SG12g–h | production RHS | model 7 (KH-RT) with `rho-solid = 4600`: `setupRHS(7, …)` puts the solid-density slot right after the metal block; one `rhsSolidification` call per phase at SG11's 100 m/s slip with `Z(9) = 1.5 ṁ/m_0`: `F(7)` or `F(8)` at the constant-mass diameter `(6 m ((1 − f)/ρ_l + f/ρ_s)/π)^{1/3}` of `m = ṁ/Z(9)`, Ranz-Marshall at its `Re`, the other zero (SG12g); `F(9)` against a direct `breakupOde` at that `d` in the liquid phases, exactly 0 on the plateau and in the solid (SG12h) | closed form; the breakup routine | 1e-14 relative on `F(7)`, `F(8)` (measured 1.8e-16), 1e-12 on `F(9)` (0); exact on the layout and the zeros | measured: model 7's `rho-solid` branch skipped (`d` from `ρ_l`) fails SG12g by 7.0e-2 (3.3e-2 on the plateau, the liquid phases unchanged) |
| SG13 | definition | `breakupStateAtSwitch(which, brkSelect, tNow, s)`, the breakup state the segment finalize and the segment-start catch-up set at a phase change, for `which` 0–4 and `brkSelect` 0, 3, 4, 5: `[0, 0]` (the TAB/ETAB oscillator at rest) at nucleation under TAB or ETAB, `[tNow, 0]` (the KH-RT timer restarted) at a re-melt under KH-RT, `s` otherwise | the table | exact (20 pairs) | measured: the nucleation branch dropped fails 2 pairs, the re-melt branch 1 |
| SB1 | birth state | `birthResidue` (`obj_condensed.f90`): a burnout residue of a freezing combustion product (`T_m = 2327`, `T_nuc = 1861.6 K`, `ρ = 3950`) from one parent record at 2500, 2000, 1500 K | phase liquid / undercooled / undercooled with `f = 0` (never solid at its own temperature); `m`, `d`, `ṅ`, `ṁ`, `m₀` for one residue per droplet (`m` = the cap, `d = (6m/(πρ))^{1/3}`) and at `d-residue` 5 µm (`m = ρπd³/6`, `ṅ` the parent's × cap/m); position, velocity, time, locators, parent ID, the birth row due | exact | the residue born by `solidPhaseAtInjection` (solid at 1500 K, `f = 1`), and sized with the burning material's density: both fail |

`dump_curve('solidification-recalescence', 'SG4', …)` writes the unit figure `unit-solidification-recalescence.svg`.

## End to end

| case | what it gates |
|---|---|
| [`solid-box`](solid-box/INFO.md) | the three regimes against their closed forms on the shared box, the event position, the energy deposit (global telescoping and the per-cell plateau heat), the frozen fraction of every row, the plane records before the nucleation point and on the plateau |
| [`solid-melt`](solid-melt/INFO.md) | a solid heated to `T-melt` in a hotter gas: the solid, melting-plateau and liquid regimes against their closed forms, the event positions, the heat taken from the gas per plateau cell, the frozen fraction of every row, no parcel lost at a phase boundary |
| [`solid-box-euler`](solid-box-euler/INFO.md) | the equivalent-Eulerian field of model 6: coverage, `ρ_p/n_p = m`, `u_p`, `T_p` per cell |
| [`solid-box-2mat`](solid-box-2mat/INFO.md) | a model-6 material followed by a model-1 one in one run |
| [`solid-box-rad`](solid-box-rad/INFO.md) | the three regimes under radiation alone (`NoHeat`, `emissivity = 0.8`) against their closed forms, the plateau rows, the energy deposit |
| [`scatter-weight-solid`](scatter-weight-solid/INFO.md) | the scatter cloud's marker count per parcel across the phase events: an event keeps the scatter weight of the part of the step it keeps |
| [`solid-brk-tab`](solid-brk-tab/INFO.md) | TAB on a solidifying material (model 6): a family that nucleates before its first breakup never breaks up, one that breaks up while liquid keeps its size once frozen; flows and energy |
| [`solid-brk-khrt`](solid-brk-khrt/INFO.md) | KH-RT on a solidifying material (model 7): stripping and RT on the liquid only, KH children born in their parent's phase and reaching the outlet; flows, energy, the Eulerian mass, the scatter rows |
| [`solid-box-rhos`](solid-box-rhos/INFO.md) | `rho-solid` on solid-box: the mass on every row, the diameter of every row from its `f`, the plateau `f(x)` and the solid cooling against their closed forms, the energy deposit, the scatter rows |
| [`solid-melt-rhos`](solid-melt-rhos/INFO.md) | `rho-solid` on solid-melt: the frozen-injected particle's mass at the solid density, the diameter of every row from its `f`, the solid heating, the melting plateau with its growing diameter and the liquid against their closed forms, the energy taken from the gas |
| [`solid-box-varrho`](solid-box-varrho/INFO.md) | a density table on solid-box: the mass on every row, the diameter of every row from the table at its `T`, the liquid and solid cooling against an independent integration, the plateau, the energy deposit |
| [`solid-brk-tab-rhos`](solid-brk-tab-rhos/INFO.md) | TAB with `rho-solid`: a frozen drop keeps the mass of its last breakup, the diameter of every row from its `m` and `f`; flows and energy |
| `repeatability/solid-box` | two sweeps reproduce each other: the phase and the frozen fraction are re-derived at injection |
| `infrastructure/refusals/solid-*`, `rho-solid` | the seven setup refusals of the input contract (`h-fus`, `cp-solid`, `T-nuc < T-melt`, no evaporation or combustion, a constant `cp`, `rho-solid` not negative) |
