# INFO — evaporation/tc-hexadecane (TC2012 Fig. 11 reproduction, e2e, **GREEN**)

## Reference
- **[TC2012]** Tonini, S.; Cossali, G. E. *"An analytical model of liquid drop evaporation
  in gaseous environment."* Int. J. Thermal Sciences **57** (2012) 45–53. **Figure 11** —
  non-dimensional drop size and temperature of an **n-hexadecane** drop, R_{d,0}=10 µm,
  T_{d,0}=300 K, vapour-free air at T_∞=600 K.
- Digitized present-model curves in `reference/` (WebPlotDigitizer; provenance in
  `reference/PROVENANCE.md`). **NON-gating overlay only.**

## What this case adds over `tc-box`
`tc-box` validates the TC Stefan-Fuchs rate for a synthetic **constant-density** heavy
fuel. This case is a **real fuel with temperature-dependent liquid density**
(`INPUT/properties.dat`, ρ_l 767→570 over 300→560 K, the paper's Table-1 anchors), so the
drop **swells** as it heats before the D²-law evaporation — the hallmark of Fig. 11. It runs
the T-dependent-property path end to end: with the density table (`rhoTab`) left
unallocated, the first variable-density read segfaults. The density and `Psat` columns are
read through `tableValue` (the end value held outside the table) at SDIRK4's Newton trial
temperatures; `Cp` is constant, so `lookupTab` is not called.

## Case construction
Same box + injection as `tc-box` (kv=1 ⇒ Re=0, Sh=2, coast at u_g ⇒ t=x/u_g; kt=0.5 ⇒
T_{d,0}=300 K; d0=20 µm). `evaporation=TC`, `interface=VLE`. The gas is air at 600 K and one
atmosphere, TC2012's state, with the gas-film properties at the film temperature
T_r = T_s + (T_∞ − T_s)/3 = 529.1 K (the 1/3 rule at the plateau, self-consistent): `INPUT/solfile.tec` is
`tests/tools/make_box_case.py --kv 1.0 --kt 0.5 --rp 1.0e-5 --rho 0.5884146341463415 --kl 0.04256 --gam 1.38328`
(ρ_∞ = p_atm/(R T_∞); k_g = 0.04256 W/m/K and c_pg = 1035.8 J/kg/K, `GAM` = c_pg/(c_pg − R)). Fuel
properties, all traceable to sources (`reference/PROVENANCE.md`):
- `Mv=226.45`, `Tboil=560` — TC2012 Table 1.
- `p_sat(T)` — TC2012's **own Table-1 curve**, the `Psat` column of `INPUT/properties.dat`:
  piecewise Clausius-Clapeyron through its six `(T, p/p_atm)` anchors
  (`tests/tools/make_psat_table.py`), read by the solver and by the oracle as the same column.
- `Lv=2.2695e5` — Table 1's latent heat at the boiling point, the **energy sink only**
  (`m c_p dT/dt = Q + ṁ Lv`); with the `Psat` column it does not set the saturation curve.
- `cp_l=2800` (`properties.dat`, constant) — NIST/Chemeo n-hexadecane liquid `Cp`, linearly
  extrapolated to the ~490 K operating point; no source of the liquid's `Cp` above 400 K is at hand.
- `Le=5.625` — the vapour diffusivity by Fuller-Schettler-Giddings at T_r and 1 atm,
  D_v = 1.2415e-5 m²/s, taken with ρ_∞ as in TC2012's eq. 1 (`ṁ = 4π R ρ_∞ D_v m̂`):
  `Le = k_g/(ρ_∞ c_pg D_v)`. TC2012 publishes no D_v (its gas properties are constant and unnamed).
- `cpv=2582` — n-hexadecane vapour at T_r, Burcat's C16H34 NASA-7 fit.
- `ρ_l(T)` variable from `properties.dat` (Table-1 boiling anchor 569.9 @ 560 K).

**Gotcha:** FiNeR mis-parses non-ASCII bytes in comments — keep `input.ini` ASCII-only.

## What `check.py` gates (Tier V, tight)
Along the MEASURED Tp(x): (1) the **variable-density mass rate** — the oracle integrates the
droplet MASS with the TC2012 present-model rate (the eq. 16 transcendental, bisection,
independent of production's Newton), its `X_s` from the `Psat` column of `INPUT/properties.dat`
(`tests/tools/proptab.py`, linear between the nodes as the solver reads it), then reconstructs
`d² = (6m/(π ρ_l(Tp)))^{2/3}`, coupling evaporative mass loss AND thermal swelling exactly as
production does, and matches the measured d²(x) (25/25 within 0.31 %, tol 2 %); (2)
**swelling** — measured max d²/d0² > 1.03 (a constant-density case can only shrink, so this
proves ρ_l(Tp) is live); (3) **PC0** — the run reported `p_sat tabulated from properties.dat`;
(4) **PC1** — the plateau (max T_p over the gated drops) within 1 K of the 0-D replica of the
kernel, `zero_d.py` on this fixture: 493.69 K. Exercises the density and `Psat` tables.

## IGLOO's kernel IS TC2012's present model (eq. 16), not Stefan-Fuchs
Proven from the paper: IGLOO's TC transcendental `m̂ + (T̃s−1)·Le_v·(f(m̂/Le_v)−1) = rhs0`
(`Lib_Evaporation.f90:TC_model`) has, in the small-rate limit, the reduction
`m̂ = (P̂vs − χ_v∞)/[½(T̃s+1)]` — **exactly** TC2012 eq. 16's stated small-rate form (the
`½(T̃s+1)` film-average denominator), NOT the Stefan-Fuchs eq. 2b (which lacks the `(T̃s−1)`
film term). So both production **and** `check.py`'s oracle are the *present model*, and the
overlay compares IGLOO to the correct one of Fig. 11's three curves.

## Saturation curve and latent heat, decoupled

`Lv` has two roles in a Clausius-Clapeyron closure: the slope of `ln p_sat(1/T)` and the energy
sink. For n-hexadecane one scalar cannot serve both: TC2012's Table-1 curve has an effective
latent heat of 258 kJ/kg over the operating range (219–269 kJ/kg segment by segment, falling
with T), while its latent heat at the boiling point is 226.95 kJ/kg. With the `Psat` column the
curve comes from the table and `Lv` is the sink alone.

`zero_d.py` integrates the kernel's equations for one drop (RK4, statement for statement) and
prints the plateau, the lifetime and `heat-frac = t_heat/t_life` (`t_heat` the first time T_p
reaches the plateau − 1 K, `t_life` the time d²/d0² reaches 0.02); `--figure` applies the same
definitions to the digitized Fig. 11 curves:

| `p_sat` | sink `Lv` | `cp_l` | plateau [K] | heat-frac |
|---|---|---|---|---|
| **`Psat` column** | **2.2695e5** | **2800** | **493.69** | **0.600** |
| `Psat` column | 2.2695e5 | 2900 | 493.69 | 0.612 |
| Clausius-Clapeyron, `Lv` 2.2695e5 | 2.2695e5 | 2800 | 487.13 | 0.625 |
| Clausius-Clapeyron, `Lv` 2.58e5 | 2.58e5 | 2800 | 491.44 | 0.575 |
| `Psat` column | 2.58e5 | 2800 | 491.41 | 0.568 |
| Clausius-Clapeyron, `Lv` 2.58e5 | 2.2695e5 | 2800 | 493.81 | 0.608 |
| TC2012 Fig. 11 (digitized) | | | 493.72 | 0.637 |

The run lands on its replica: plateau 493.69 K for all 25 drops (0.03 K below the figure),
swelling peak 1.089 (the replica's lifetime, to d²/d0² = 0.02, 4.26 ms). What moves the plateau is the sink, not the curve's shape: the
table and a Clausius-Clapeyron line with `Lv = 2.58e5` put it within 0.12 K of each other (rows 1, 6;
rows 4, 5). `cp_l` moves the heat-frac, not the plateau; 2800 is the sourced value (see
`reference/PROVENANCE.md`) and stays.

What each property does, by the replica: with ρ = 1.2 kg/m³ at 600 K (p = 2.04 atm) and air's
300 K properties (k_g 0.026, c_pg 1004.5, Le 2.5, cpv 2300) the plateau is 492.63 K, heat-frac 0.618;
the same gas at 1 atm, 472.84 K, 0.591; with D_v from Fuller there (Le 3.543), 481.88 K, 0.601; with
the film properties of this case, 493.69 K, 0.600. In Fig. 11's time τ = t D_v/R²_d,0 the drop lives τ = 528.6 against
the figure's 212.8: a factor 2.5 that no property of the case explains, since TC2012's eq. 1 and IGLOO
both scale the rate with ρ_∞ D_v and the figure's three curves share that scale; the overlay compares
on the normalized lifetime.

**RED** (each perturbation on the committed case, then restored):

| | perturbation | PC0 | PC1 | rate (1) |
|---|---|---|---|---|
| R1 | the `Psat` column removed (sink stays 2.2695e5): psat by Clausius-Clapeyron with the sink's `Lv` | RED, 0 reports | RED, 487.13 K | RED, no column for the oracle |
| R2 | column kept, oracle's psat back to Clausius-Clapeyron with `Lv` 2.2695e5 | GREEN | GREEN | RED, 2.42e-1 |
| R3 | column kept, oracle's psat Clausius-Clapeyron with `Lv` 2.58e5 | GREEN | GREEN | **GREEN, 4.4e-3** |
| R4 | the column zeroed (Clausius-Clapeyron kept, reported as such) | RED, 0 reports | RED, 487.13 K | RED, 1.10 |
| R5 | the reader of the table without the kernel connection (psat read, not used) | GREEN | RED, 487.13 K | RED, 2.75e-1 |

R3 is the blind spot of the rate gate: the table and the Clausius-Clapeyron fit to the same
anchors differ by at most a few percent along the heating path (4.4e-3 against the floor
of 3.1e-3, under the 2 % tolerance), so the rate gate cannot tell which one the run used. PC0
and PC1 can: PC0 sees the column read, PC1 sees it used (R5). The oracle's own constants, one at a time: `RHO_G`
1.2 (a 2.04-atm gas) reads 0.606 on the rate gate, `LE` 2.5 reads 0.640, `KG` 0.026 reads 0.441; a
vapour `CPV` of 2300 reads 3.7e-3 and passes, and in the solver it moves the replica's plateau by
0.92 K, inside PC1's 1 K: the vapour cp is the one input no gate of the case resolves.

## Comparison plot (`verify.py` → `OUTPUT/tc-hexadecane.svg`, NON-gating)
IGLOO vs the digitized TC2012 Fig. 11 **present-model** curves, on a **normalized lifetime**
axis (the absolute `τ = t D_v/R²` of the run is 2.5 times the figure's, above; the normalized
**shape** is the comparison): plateau 493.7 K against the digitized 493.7 K, heat-frac 0.600
against 0.637, swell peak 1.09. TC2012 takes its gas properties constant (p. 47: "Fluid properties
(Dv, μ, k, c) are assumed to be constant") and among the liquid's names only the density as
temperature dependent (p. 52); a constant `cp_l` at the operating point is the closer reading of its model.

## Exit-path and mass-balance gates

`check.py` gates three more things after the physics rows. They are about **how the drops
leave**, not about the rate law.

Without a finite scrub in `rhsEvaporation` all 25 parcels end on

    [WARNING] Particle   N non-finite state ==> reverting to last good step, marking gone
              last good m/m0 = 1.997E-03  (burnout fires at m = 1.000E-20 kg)

— the burnout threshold is never reached. The NaN is manufactured **inside `rhsEvaporation`**, not by
the boiling clamp: on a Newton trial state with `m < 0` the RHS evaluates
`d = (6m/(pi rho))^(1/3)`, which is NaN, and that poisons `Re`, `mdot` and `F(7:8)`; without
the scrub the Inf/NaN slope rides until `solout` catches the state. Measured, this case is the
suite's only instance.

With the scrub the bad trial becomes a `1e30` slope, the step is rejected, `deltat` shrinks,
and the drop reaches the burnout threshold — a **clean burnout** through the normal branch, not through
the solver-failure exit. That branch is gated separately in
[`infrastructure/solver-fail-consumed`](../../infrastructure/solver-fail-consumed/INFO.md),
which forces it with an unreachable ODE tolerance; its siblings `solver-fail-consumed-m4`
(model 4), `solver-fail-midway` (a step limit, mid-trajectory) and `solver-fail-refine` (a step
limit on a refinement call) are the branch's other cases.

| gate | assertion | RED without the scrub | GREEN |
|---|---|---|---|
| H1 | no parcel ends on a non-finite state | 25 | 0 |
| H2a | no parcel ends on a solver failure | 0 | 0 |
| H2b | all 25 parcels wrote an exit row | 25 | 25 |
| H3 | `Σ wdot` = `Σ npdot·m_inj` within `1e-5` | 1.022e-6 (passes) | 1.022e-6 |

With the scrub the drops burn out at `x = 0.043020`, downstream of the non-finite exit (`x = 0.042624`) —
the drops evaporate further before they are removed — and with it move the last trajectory
row, the `outloc` row and the burnout cells of `source.tec`. The step that crosses the burnout mass
(1e-20 kg) ends below zero on every parcel: the burnout sets their mass to 0, and their exit records
print `d = 0` (observed; without that clamp `d` prints NaN and `check_finite.py` fails). The physics
rows are read well before burnout and are the same either way (25 drops gated, 0 rate violations,
0 non-swelling).

**H3 is the mass books, not the scrub's witness.** Both the non-finite exit and the
burnout set `consumed`, so H3 reads the same with or without the scrub — H1 is the scrub's
witness. The burnout hands the gas at most 1e-20 kg per drop (3.1e-9 of m₀): with `consumed`
removed from the burnout branch H3 reads **1.022e-6**, the floor itself (identical
at `OMP_NUM_THREADS` 1, 2 and 5), and passes, so this case does not witness the burnout's hand-off;
`evaporation/evap-breakup-box`'s mass closure (1.7e-9 against 1e-12) and `evaporation/evap-breakup-box-varcp`'s
energy and mass books do. The `1e-5` tolerance sits 9.8x above the floor.

## Spray-level rejection
As for the other evaporation/breakup cases: dense-spray SMD / penetration is out of scope
(steady one-way carrier gas, no entrainment, no atomizing nozzle).
