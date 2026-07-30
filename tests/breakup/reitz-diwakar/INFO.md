# INFO — breakup_reitz_diwakar (family D, test D2)

Verifies the Reitz–Diwakar breakup model (selector 2). This is an
ODE-rate model (`breakupOde` integrates `npdot`); since parcel mass·npdot
is constant, the reduced equation is `ḋ = (d_stab − d)/τ`, but
`τ = τ(d)` (bag `τ ∝ d^1.5`, strip `τ ∝ d`) — **not a clean exponential**.
⇒ needs a Tier-2 self-contained reference integrator.

Schema: see verification-phase1-plan.md §2.5.
Citation tags resolve in [../../REFERENCES.md](../../REFERENCES.md).

## Tests

| id | mode | source | doi | locus | compared | settings_provenance | status |
|---|---|---|---|---|---|---|---|
| D2 | hybrid | `[RD86]` Reitz & Diwakar 1986 | DOI:10.4271/860469 | bag regime: `τ ∝ d^1.5`, threshold `We_bag = WeBag`; strip regime: `τ ∝ d`, threshold from `Cs/√Re`; `npdot` evolution gives reduced ODE `ḋ = (d_stab − d)/τ(d)`; `bp = [WeBag, Cb, Cstrip, Cs]` | (a) regime selection vs (We, Re) sweep matches [RD86] regime map; (b) `d(t)` vs Tier-2 oracle of production `npdot` rate (re-implemented independently in `verif_oracle.f90`): `relL2 < tol`; (c) `τ` and `d_stab` values match closed-form expressions from [RD86] eq. set | constant gas state at (We, Re) points spanning bag-only / strip-only / threshold; single parcel; INI `bp` read in test driver | not started |

## Notes

- The "hybrid" mode reflects: (a) analytic regime map + closed-form `τ`, `d_stab`
  values from [RD86] eqs.; (b) literature-numeric character of `d(t)` (no clean
  closed form — Tier-2 oracle re-derives the production rate from [RD86]
  equations and integrates with `rk4_ref`).
- Test must **not** call any production breakup routine when constructing the
  oracle — independent re-implementation only.
- `npdot` carries dimension [1/s]; reduced `ḋ = (d_stab−d)/τ` derivation
  assumes parcel mass·npdot conservation (plan §1.4). Test asserts this
  invariant holds to ≤ 10·ε_mach.

## Status (2026-07-02)

`test_breakup_rd.f90` (ctest `test_breakup_rd`, GREEN): RD1 sub-critical zero;
RD2 bag + RD3 stripping vs OpenFOAM-form oracle (Cb=π/4, Cs=10=bp(4)/2), ~2e-16;
RD4 regime handoff at We_r=0.5√Re, each side matches its own branch (3.9e-13).
Re convention resolved: gas diameter Reynolds (Lib_Integration.f90:474), as the
stripping dStable requires. Algebra note: no dp<dStable guard needed — bag gives
dStable/d=6/We_r<1 and stripping dStable/d=Re/(4We_r²)<1 whenever entered.
