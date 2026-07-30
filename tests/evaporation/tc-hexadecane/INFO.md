# INFO — evaporation/tc-hexadecane (E-VAL-3: TC2012 Fig. 11 reproduction, e2e, **GREEN**)

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
drop **swells** as it heats before the D²-law evaporation — the hallmark of Fig. 11. It is
the suite's **first variable-property case**, and building it surfaced two production bugs
on the never-exercised T-dependent-property path:
- **A20** — the variable-property tabs (`rhoTab`/`cpTab`/`hTab`) were declared allocatable
  but never allocated; the first variable-density read segfaulted.
- **A21** — `lookupTab` indexed the table with `idint(T)` and no bounds guard; SDIRK4's
  Newton probes trial temperatures outside the tabulated range on stiff evaporation
  transients, reading out of bounds. Fixed by clamping the index.

## Case construction
Same box + injection as `tc-box` (kv=1 ⇒ Re=0, Sh=2, coast at u_g ⇒ t=x/u_g; kt=0.5 ⇒
T_{d,0}=300 K; d0=20 µm), gas air at 600 K. `evaporation=TC`, `interface=VLE`. Fuel
properties (Tier-1 reconstruction, 2026-07-28, all traceable to sources):
- `Mv=226.45`, `Tboil=560` — TC2012 Table 1.
- `Lv=2.58e5` — the effective latent heat of TC2012's **own Table-1 psat curve** (a CC slope
  through its `(T,Pvs/PN)` anchors ≈ Watson `L_v` at the ~490 K wet-bulb). The earlier `2.9e5`
  made IGLOO's CC-psat run **20–24 % low** vs Table 1, over-heating the drop.
- `cp_l=2800` (`properties.dat`, constant) — NIST/Chemeo n-hexadecane liquid `Cp` at the
  ~490 K operating point. The earlier `2200` was the **298 K** value (right property, wrong
  reference T) — confirmed: NIST `Cp(298 K)=499–500 J/mol·K = 2205 J/kg·K`.
- `Le=2.5`, `cpv=2300` — n-hexadecane vapour (not the air-default `Le=1`); kept physical.
- `ρ_l(T)` variable from `properties.dat` (Table-1 boiling anchor 569.9 @ 560 K).

**Gotcha:** FiNeR mis-parses non-ASCII bytes in comments — keep `input.ini` ASCII-only.

## What `check.py` gates (Tier V, tight)
Along the MEASURED Tp(x): (1) the **variable-density mass rate** — the oracle integrates the
droplet MASS with the TC2012 eq. 9 Stefan-Fuchs rate (bisection, independent of production's
Newton), then reconstructs `d² = (6m/(π ρ_l(Tp)))^{2/3}`, coupling evaporative mass loss AND
thermal swelling exactly as production does, and matches the measured d²(x) (25/25 within
0.2 %, tol 2 %); (2) **swelling** — measured max d²/d0² > 1.03 (a constant-density case can
only shrink, so this proves ρ_l(Tp) is live). Exercises A20 + A21.

## IGLOO's kernel IS TC2012's present model (eq. 16), not Stefan-Fuchs
Proven from the paper: IGLOO's TC transcendental `m̂ + (T̃s−1)·Le_v·(f(m̂/Le_v)−1) = rhs0`
(`Lib_Evaporation.f90:TC_model`) has, in the small-rate limit, the reduction
`m̂ = (P̂vs − χ_v∞)/[½(T̃s+1)]` — **exactly** TC2012 eq. 16's stated small-rate form (the
`½(T̃s+1)` film-average denominator), NOT the Stefan-Fuchs eq. 2b (which lacks the `(T̃s−1)`
film term). So both production **and** `check.py`'s oracle are the *present model*, and the
overlay compares IGLOO to the correct one of Fig. 11's three curves.

## Comparison plot (`verify.py` → `OUTPUT/tc-hexadecane.svg`, NON-gating)
IGLOO vs the digitized TC2012 Fig. 11 **present-model** curves, on a **normalized lifetime**
axis (IGLOO's `D_v` is `Le`-set and differs from TC2012's n-hexadecane `D_v`, so absolute
`τ=tD_v/R²` is incomparable; the normalized **shape** is the valid comparison). With the
Tier-1 property reconstruction IGLOO reproduces the **heating shape** (heat-frac ≈ 0.54,
matching Fig. 11) and the swell peak (1.09), with the plateau at 490.2 K (**0.7 % below**
the digitized ~493.6 K).

**Root cause of the shape mismatch — mis-set properties, corrected from TC2012's own data.**
The earlier attribution ("only the Lewis number") was incomplete: `Le` (2f881d0) tuned the
*plateau*, but the *heating shape / d²-decline* was wrong because two properties were mis-set
against TC2012's own Table 1: (1) **psat** — IGLOO's CC anchored with `Lv=2.9e5` ran 20–24 %
below Table 1's `(T,Pvs/PN)` curve, so the drop under-evaporated and over-heated (peak too
early, heat-frac 0.39 vs 0.54); (2) **cp_l** — `2200` is the 298 K value, but the drop
operates at ~490 K where `cp_l≈2800`. Correcting both (`Lv=2.58e5` → Table-1-consistent psat;
`cp_l=2800` → NIST operating-T value) moves heat-frac 0.39 → **0.536** — the shape now tracks
Fig. 11. TC2012 uses *constant* gas-film properties (only `ρ_l` T-dependent), so a **constant**
`cp_l` is the faithful choice — no variable-cp law (which would also trip an untested code
path). `Le=2.5` is kept (physical, `D_v`-justified).

**Deliberate trade — plateau vs shape (read alongside commit 2f881d0).** The prior state
(2f881d0) had a near-exact plateau (493.8 vs ~493.6 K) but the wrong shape (heat-frac 0.39).
That plateau was a **compensating-error coincidence**: a psat 20–24 % too low **and** `Le=2.5`
tuned against it happened to land the wet-bulb. With the psat corrected, the same physical
`Le=2.5` gives 490.2 K. Tier-1 therefore trades a hair of plateau accuracy (now **3.4 K / 0.7 %
low**) for the *shape* — which was the actual complaint. Tier-2 (decoupled psat) recovers the
plateau too (0-D: 492.7 K) without re-introducing the compensation.

**Residual (documented) → Tier 2, deferred.** The `d²` mid-decline still sits slightly below
TC2012's, and the plateau is 3.4 K low, because a **single scalar `Lv`** cannot be BOTH the
psat-curve slope (~258 kJ/kg) AND the optimal energy sink (~227 kJ/kg from Table-1 ΔHv). The
full match (0-D verified: plateau 492.7, heat-frac 0.538) needs the psat curve **decoupled**
from the sink `Lv` — i.e. the tabulated/variable-psat capability. This makes tc-hexadecane the
concrete demonstrator that justifies that feature; the plan is in
`plan-bucket/tc-hexadecane-tier2-decoupled-psat.md`. The tight gate stays IGLOO-vs-its-own-TC
kernel; the paper curves are a non-gating overlay.

## Tier-2 (spray-level) rejection
As for the other evaporation/breakup cases: dense-spray SMD / penetration is out of scope
(steady one-way carrier gas, no entrainment, no atomizing nozzle).
