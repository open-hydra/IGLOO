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
the T-dependent-property path end to end, and two failures of that path stop it:
- the variable-property tabs (`rhoTab`/`cpTab`/`hTab`) left unallocated: the first
  variable-density read segfaults;
- a `lookupTab` that indexes the table with `idint(T)` and no bounds guard: SDIRK4's
  Newton probes trial temperatures outside the tabulated range on stiff evaporation
  transients and reads out of bounds. The index is clamped.

## Case construction
Same box + injection as `tc-box` (kv=1 ⇒ Re=0, Sh=2, coast at u_g ⇒ t=x/u_g; kt=0.5 ⇒
T_{d,0}=300 K; d0=20 µm), gas air at 600 K. `evaporation=TC`, `interface=VLE`. Fuel
properties, all traceable to sources (`reference/PROVENANCE.md`):
- `Mv=226.45`, `Tboil=560` — TC2012 Table 1.
- `p_sat(T)` — TC2012's **own Table-1 curve**, the `Psat` column of `INPUT/properties.dat`:
  piecewise Clausius-Clapeyron through its six `(T, p/p_atm)` anchors
  (`tests/tools/make_psat_table.py`), read by the solver and by the oracle as the same column.
- `Lv=2.2695e5` — Table 1's latent heat at the boiling point, the **energy sink only**
  (`m c_p dT/dt = Q + ṁ Lv`); with the `Psat` column it does not set the saturation curve.
- `cp_l=2800` (`properties.dat`, constant) — NIST/Chemeo n-hexadecane liquid `Cp` at the
  ~490 K operating point.
- `Le=2.5`, `cpv=2300` — n-hexadecane vapour (not the air-default `Le=1`); kept physical.
- `ρ_l(T)` variable from `properties.dat` (Table-1 boiling anchor 569.9 @ 560 K).

**Gotcha:** FiNeR mis-parses non-ASCII bytes in comments — keep `input.ini` ASCII-only.

## What `check.py` gates (Tier V, tight)
Along the MEASURED Tp(x): (1) the **variable-density mass rate** — the oracle integrates the
droplet MASS with the TC2012 present-model rate (the eq. 16 transcendental, bisection,
independent of production's Newton), its `X_s` from the `Psat` column of `INPUT/properties.dat`
(`tests/tools/proptab.py`, linear between the nodes as the solver reads it), then reconstructs
`d² = (6m/(π ρ_l(Tp)))^{2/3}`, coupling evaporative mass loss AND thermal swelling exactly as
production does, and matches the measured d²(x) (25/25 within 0.19 %, tol 2 %); (2)
**swelling** — measured max d²/d0² > 1.03 (a constant-density case can only shrink, so this
proves ρ_l(Tp) is live); (3) **PC0** — the run reported `p_sat tabulated from properties.dat`;
(4) **PC1** — the plateau (max T_p over the gated drops) within 1 K of the 0-D replica of the
kernel, `zero_d.py` on this fixture: 492.63 K. Exercises the variable-property tabs and their
bounds guard.

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
| **`Psat` column** | **2.2695e5** | **2800** | **492.63** | **0.618** |
| `Psat` column | 2.2695e5 | 2900 | 492.63 | 0.629 |
| Clausius-Clapeyron, `Lv` 2.2695e5 | 2.2695e5 | 2800 | 486.00 | 0.644 |
| Clausius-Clapeyron, `Lv` 2.58e5 | 2.58e5 | 2800 | 490.17 | 0.591 |
| `Psat` column | 2.58e5 | 2800 | 490.19 | 0.583 |
| Clausius-Clapeyron, `Lv` 2.58e5 | 2.2695e5 | 2800 | 492.70 | 0.626 |
| TC2012 Fig. 11 (digitized) | | | 493.72 | 0.637 |

The run lands on its replica: plateau 492.63 K for all 25 drops (0.2 % below the figure),
swelling peak 1.088. What moves the plateau is the sink, not the curve's shape: the table and a
Clausius-Clapeyron line with `Lv = 2.58e5` agree to 0.07 % near 490 K (rows 4–6). `cp_l` moves the
heat-frac, not the plateau; 2800 is the sourced value (see `reference/PROVENANCE.md`) and stays.

**RED** (each perturbation on the committed case, then restored):

| | perturbation | PC0 | PC1 | rate (1) |
|---|---|---|---|---|
| R1 | the `Psat` column removed (sink stays 2.2695e5): psat by Clausius-Clapeyron with the sink's `Lv` | RED, 0 reports | RED, 486.00 K | RED, no column for the oracle |
| R2 | column kept, oracle's psat back to Clausius-Clapeyron with `Lv` 2.2695e5 | GREEN | GREEN | RED, 2.51e-1 |
| R3 | column kept, oracle's psat Clausius-Clapeyron with `Lv` 2.58e5 | GREEN | GREEN | **GREEN, 7.3e-3** |
| R4 | the column zeroed (Clausius-Clapeyron kept, reported as such) | RED, 0 reports | RED, 486.00 K | RED, 1.08 |
| R5 | the reader of the table without the kernel connection (psat read, not used) | GREEN | RED, 486.00 K | RED, 2.64e-1 |

R3 is the blind spot of the rate gate: the table and the Clausius-Clapeyron fit to the same
anchors differ by at most a few percent along the heating path (7.3e-3 against the floor
of 1.9e-3, under the 2 % tolerance), so the rate gate cannot tell which one the run used. PC0
and PC1 can: PC0 sees the column read, PC1 sees it used (R5).

## Comparison plot (`verify.py` → `OUTPUT/tc-hexadecane.svg`, NON-gating)
IGLOO vs the digitized TC2012 Fig. 11 **present-model** curves, on a **normalized lifetime**
axis (IGLOO's `D_v` is `Le`-set and differs from TC2012's n-hexadecane `D_v`, so absolute
`τ=tD_v/R²` is incomparable; the normalized **shape** is the valid comparison): plateau
492.6 K against the digitized ~493.7 K, swell peak 1.09. TC2012 uses *constant* gas-film
properties (only `ρ_l` T-dependent), so a **constant** `cp_l` is the faithful choice.

## Exit-path and mass-balance gates

`check.py` gates three more things after the physics rows. They are about **how the drops
leave**, not about the rate law.

Without a finite scrub in `rhsEvaporation` all 25 parcels end on

    [WARNING] Particle   N non-finite state ==> reverting to last good step, marking gone
              last good m/m0 = 1.357E-03  (burnout fires at m = 1.000E-15 kg)

— `mBurnTol` is never reached. The NaN is manufactured **inside `rhsEvaporation`**, not by
the boiling clamp: on a Newton trial state with `m < 0` the RHS evaluates
`d = (6m/(pi rho))^(1/3)`, which is NaN, and that poisons `Re`, `mdot` and `F(7:8)`; without
the scrub the Inf/NaN slope rides until `solout` catches the state. Measured, this case is the
suite's only instance.

With the scrub the bad trial becomes a `1e30` slope, the step is rejected, `deltat` shrinks,
and the drop reaches `mBurnTol` — a **clean burnout** through the normal branch, not through
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
| H3 | `Σ wdot` = `Σ npdot·m_inj` within `1e-5` | 2.941e-7 (passes) | 2.941e-7 |

The scrub moves the exit from `x = 0.070520` (the non-finite exit) to `x = 0.070854` — the
drops evaporate a little further before they are removed — and with it the last trajectory
row, the `outloc` row and the burnout cells of `source.tec`. The physics rows are read well
before burnout and are the same either way (25 drops gated, 0 rate violations, 0 non-swelling).

**H3 guards the consumption contract, not the scrub.** Both the non-finite exit and the
burnout set `consumed`, so H3 reads the same with or without the scrub — H1 is the scrub's
witness. H3's own RED, on a throwaway build with `consumed` removed from the
burnout branch: resid **2.048e-4** against a floor of **2.941e-7** that is identical at
`OMP_NUM_THREADS` 1, 2 and 5. The `1e-5` tolerance therefore sits 34x above the floor and 20x
below the signature; **do not loosen it past ~5e-5** or it stops seeing the dropped remnant.

## Spray-level rejection
As for the other evaporation/breakup cases: dense-spray SMD / penetration is out of scope
(steady one-way carrier gas, no entrainment, no atomizing nozzle).
