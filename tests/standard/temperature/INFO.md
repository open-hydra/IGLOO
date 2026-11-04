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
  D-HEAT-2/3 stay open flags). HN4–HN8 pin the four other laws to their printed
  sources (production `heat()` against the forms typed in the test, with the
  printed exponents):

  | id | law | source (as printed) | grid | compared | tolerance | measured | RED |
  |---|---|---|---|---|---|---|---|
  | HN4 | JAXA2 | `[SP8039]` Table II ref. 82: 2 + 0.370 Re^0.6 Pr^⅓ | Re {0.5, 1, 10, 100, 1000} × Pr {0.7, 1, 7} × Ma {0, 0.01, 0.1, 0.5, 1} | Nu, relative | 1e-14 | 1.8e-16 | exponent 0.6 → 0.55: 0.28 |
  | HN5 | JAXA3 | `[SP8039]` Table II ref. 42: 2 + 0.459 Re^0.55 Pr^⅓ | as HN4 | Nu, relative | 1e-14 | 0 | 0.459 → 0.495: 7.5e-2 |
  | HN6 | JAXA4 | `[Shim06]` eq. 50: [(2 + 0.654 Re^½ Pr^⅓)^−1 + 3.42 Ma/(Re Pr)]^−1 (SP-8039's last row prints Pr^½; the code follows Shimada) | Re × Pr × Ma {0.01, 0.1, 0.5, 1} | Nu, relative | 1e-14 | 1.9e-16 | 0.654 → 0.645: 1.3e-2 |
  | HN7 | Kavanau-Drake | `[SP8039]` Table II ref. 72 = `[Shim06]` eq. 49: N₀/(1 + 3.42 (Ma/(Re Pr)) N₀), N₀ = 2 + 0.459 Re^0.55 Pr^0.33 | as HN6 | Nu, relative | 1e-14 | 0 | Pr^0.33 → Pr^⅓: 6.2e-3; 3.42 → 3.24: 5.3e-2 |
  | HN8 | JAXA4, Kavanau-Drake | the rarefaction identity 1/Nu(Re,Pr,Ma) − 1/Nu(Re,Pr,0) = 3.42 Ma/(Re Pr) (no transcription) | as HN6 | residual over the bound 8ε (1/Nu(Ma) + 1/Nu(0)) | 1 | 9.4e-2 | 3.42 → 3.24 in Kavanau-Drake: 2.9e13 |

  HN8's bound follows the conditioning: where 3.42 Ma/(Re Pr) ≪ 1/Nu the difference of the two reciprocals loses
  up to ≈ 5e3 in relative terms, so a relative tolerance would fail on correct code. The e2e twins of these
  legs are `conv-nu-jaxa2`, `-jaxa3`, `-jaxa4` and `-kd` (the laws through the solver's Re, Pr and Mach number).
- `test_heat_probes.f90` (ctest `test_heat_probes`, ordinary gate — bug-transcription pin): XH1
  JAXA1 missing +2 floor — Nu(Re=1e-6)=0.315 today. Flips the suite RED when
  production is fixed, signalling promotion to the gate (bug A7).
- `test_radiation.f90` (ctest `test_radiation`): the grey-body exchange `Q̇_r = ε σ π d² (T_g⁴ − T_p⁴)`
  against the form typed in the test (`σ = 5.670374419e-8`, CODATA) in quadruple precision; RU2 and RU4 at zero
  slip under `NoDrag` and `NoHeat`, so the radiative term is the only heat exchange:

  | id | compared | tolerance | measured | RED |
  |---|---|---|---|---|
  | RU1 | `radiativeHeat` at six points (both signs of `T_g − T_p`) against the closed form rounded once; the sign; `ε = 0` and `T_g = T_p` give exactly 0 | 4 ULP; exact | 2 ULP | `σ = 5.67e-8`: 5.6e11 ULP; area `π d²/2`: 4.3e15 ULP |
  | RU3 | the auxiliary layout, models 1–8 (nine layouts): with radiation one slot after every other block (`nAux + 1`, the slot at `nAux`); a material without radiation set up after one with it gets no slot and the plain `nAux` | exact | 9 of 9 | the slot index not reset per setup: it stays at the previous layout's index |
  | RU2 | per ODE model, `F(ε = 0.8) − F(ε = 0)` = `Q̇_r cpFactor/m` in the energy equation and 0 in every other component: models 1–5, model 1 on a varying-cp table (`cpFactor = 1`), models 6 and 7 in their four phases (`Q̇_r/(m c_l)` liquid and undercooled, `−Q̇_r/(m h_fus)` in `F(8)` on the plateau, `Q̇_r/(m c_s)` solid; model 7 with Reitz-Diwakar at `Z(9) = 1.5 ṁ/m_0`, its `m` and `d` rebuilt from the stream flow), model 2 with `blowing = LK` (`f₂ = 0.9116`: the term is not multiplied by `f₂`); model 8 (the oxide cap, `m_ox = m_Al`): `Q̇_r` on the equal-volume sphere over `m_Al c_Al + m_ox c_res` | 8 ULP of the larger `F` | worst 0.34 of tol (model 7 undercooled); model 8 0.05 | model 4 without the term: its leg (0 against 513.87 K/s); the term added before the blowing factor: the LK leg (468.45 = 513.87 `f₂`); model 7 without the term, or with the injection diameter: its four legs (`dF/ref` 0 and 1.31); model 8 without the term (`dF` 0) or on the metal's diameter (`dF/ref` 0.689, the area ratio 1.4504) |
  | RU4 | `getSourceTerms` (the Eulerian source branch): the energy component against `n_p ε σ π d² (T_g⁴ − T_p⁴)`, `d = (6 ρ̂_p/(π ρ n_p))^⅓`; `ε = 0` gives exactly 0 | 1e-14 relative | 9.5e-16 | the term absent from the branch: 1.0 |

  The e2e twins are `rad-cool`, `emissivity-cap`, `rad-foreign-key` and `solidification/solid-box-rad`.
