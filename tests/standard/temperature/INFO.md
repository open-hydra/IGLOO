# INFO — temperature (family B)

Verifies the per-particle thermal ODE under various heat closures.
Same L1/L2 dual mechanism as family A.

Schema: see verification-phase1-plan.md §2.5.
Citation tags resolve in [../../REFERENCES.md](../../REFERENCES.md).

## Tests

| id | mode | source | doi | locus | compared | settings_provenance | status |
|---|---|---|---|---|---|---|---|
| B1 | analytic | `[RM52]` (Nu=2 limit) | n/a | `Nu=2` zero-slip limit ⇒ `dT/dt=(T_g−T_p)/τ_T`, `τ_T=ρ_p c_p d²/(12 k_g)`; exact exponential `T_p(t)=T_g+(T_p₀−T_g)e^{−t/τ_T}` | `T,vx: relL∞` (adaptive) | zero slip via `u_g=v₀=0` (no Nu=2 selector exists); Tg=400,T0=300 | **PASS** (3.6e-12) |
| BO | analytic | `[RM52]` (Nu=2) | n/a | order-of-accuracy on the LINEAR (Nu=2) heat RHS, fixed-step SDIRK4 | `p_obs ≈ 4 ± 0.3` over Δt/τ_T={0.4,0.2,0.1} | tend=2τ_T, fixed-step, tight Newton | **PASS** (err 3.5e-4→1.2e-6, p=4.12/4.08) |
| B2 | analytic + numeric | `[RM52]` full + `[CGW78]` | n/a | full `Nu=2+0.6 Re^½ Pr^⅓` COUPLED to Stokes drag; reference = independent RK4 of the same coupled system (`test_temperature.f90::rhs_heat`, 2e5 substeps) | `vx,T vs RK4 oracle: relL∞` | gas u=10 ⇒ relaxing slip; Tg=400,T0=300 | **PASS** (5.2e-11) |
| B3a | analytic | derived | n/a | equilibrium: T_p₀ = T_g ⇒ F(7) = 0 exactly, so every SDIRK4 stage equals the state | T = T_g **bit for bit** after every fixed step (one step per driver call, 20 steps at each Δt/τ_T = 0.1, 1, 10) | zero slip (Nu = 2), T_g = 400 K | **PASS** (60 of 60 steps exact). RED with the right-hand side's temperature taken 1e-12 high: 60 of 60 steps move |
| B3b | analytic | derived (SDIRK4's stability function lies in (0, 1) for every negative argument) | n/a | no overshoot, no growth: from T_p₀ = 300 K and 500 K (T_g = 400 K), at every fixed step T − T_g keeps its sign (or reaches 0) and \|T − T_g\| does not increase | inequalities, compared until \|T − T_g\| ≤ 1e3·ε·T_g = 8.9e-11 K (the roundoff floor) | zero slip (Nu = 2), Δt/τ_T = 0.1, 1, 10, one step per driver call | **PASS** (638 steps, 0 violations). RED with the heat sign flipped (`Lib_Equations.f90` `(Tg-temp)` → `(temp-Tg)`): \|T − T_g\| grows from the first step (5518 violations) |
| B3c | analytic | derived | n/a | B3b at a constant slip: `NoDrag` keeps v = 0 against u_g = 10 m/s, Ranz-Marshall Nu = 5.069 constant | as B3b, steps relative to τ_T = c_p ρ_p d²/(6 Nu k_g) | Δt/τ_T = 0.1, 1, 10 | **PASS** (638 steps, 0 violations); RED as B3b (5538 violations) |

## Notes

- B1 uses zero slip because IGLOO's heat module exposes no standalone "Nu=2"
  selector — `Nu=2` is realized via `Re=0` (gas velocity = particle velocity).
- B2 keeps the full momentum→heat coupling (no frozen velocity) and verifies it
  against an independent RK4 reimplementation of the same coupled system — a
  Tier-2 oracle that shares no code with production (uses production's `Pr`
  definition `μγR/((γ−1)k)` to match the closure, not its implementation).
- B3a–B3c step SDIRK4 one fixed step per driver call (`MODE_SDIRK_FIXEDSTEP`, nstep = 1, from the previous
  state), so the invariants are checked at every step, not only at the end.
- Status (T3): B1/BO/B2 all PASS; clean `build_verif`, ctest green. No findings.

## Unit tables (T5, 2026-07-02)

- `test_heat_nu.f90` (ctest `test_heat_nu`, GREEN): HN1 RM Nu(Re,Pr) table vs
  [RM52] formula (18 pts, machine-zero); HN2 conduction floor Nu(Re=0)=2 exact
  for RM/JAXA2/JAXA3 (deliberately NOT JAXA4/KD — their sub-2 low-Re limit is
  physical rarefaction); HN3 Kavanau-Drake Ma=0 collapse (transcription pin;
  D-HEAT-2/3 stay open flags).
- `test_heat_probes.f90` (ctest `test_heat_probes`, ordinary gate — bug-transcription pin): XH1
  JAXA1 missing +2 floor — Nu(Re=1e-6)=0.315 today. Flips the suite RED when
  production is fixed, signalling promotion to the gate (bug A7).
