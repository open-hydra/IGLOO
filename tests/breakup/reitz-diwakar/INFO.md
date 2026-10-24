# INFO — breakup_reitz_diwakar (family D, test D2)

Verifies the Reitz–Diwakar breakup model (selector 2). This is an
ODE-rate model (`breakupOde` integrates `npdot`); since parcel mass·npdot
is constant, the reduced equation is `ḋ = (d_stab − d)/τ`, but
`τ = τ(d)` (bag `τ ∝ d^1.5`, strip `τ ∝ d`) — **not a clean exponential**.
⇒ needs a self-contained reference integrator.

Citation tags resolve in [../../REFERENCES.md](../../REFERENCES.md).

## Tests

| id | mode | source | doi | locus | compared | settings_provenance | status |
|---|---|---|---|---|---|---|---|
| D2 | hybrid | `[RD86]` Reitz & Diwakar 1986 | DOI:10.4271/860469 | bag regime: `τ ∝ d^1.5`, threshold `We_bag = WeBag`; strip regime: `τ ∝ d`, threshold from `Cs/√Re`; `npdot` evolution gives reduced ODE `ḋ = (d_stab − d)/τ(d)`; `bp = [WeBag, Cb, Cstrip, Cs]` | (a) regime selection: RD1 (sub-critical `We_r ≤ 6` → rate 0) and RD4 (the handoff at `We_r = 0.5√Re`, each side on its own branch); (b) `d(t)`: not here; [reitz-diwakar-e2e](../reitz-diwakar-e2e/INFO.md) gates its initial rate against an RK4 reference of the RD rate integrated along the measured slip (2 %); (c) `τ` and `d_stab` values: RD2 (bag) and RD3 (stripping) against the OpenFOAM-form oracle | constant gas state at (We, Re) points spanning bag-only / strip-only / threshold (`Re = 1e4`, so `0.5√Re = 50`); single parcel; `bp = [6, π, 0.5, 20]` (the INI defaults) set in the test driver | GREEN: `test_breakup_rd` (registered in `tests/CMakeLists.txt`) |

## Notes

- The "hybrid" mode reflects: (a) analytic regime map + closed-form `τ`, `d_stab`
  values from [RD86] eqs.; (b) literature-numeric character of `d(t)` (no clean
  closed form — the e2e case's reference integrates the RD rate with RK4).
- The oracle calls no production breakup routine — an independent re-implementation.
- `npdot` carries dimension [1/s]; the reduced `ḋ = (d_stab−d)/τ` assumes parcel
  mass·npdot conservation, which the ODE model 3 (breakup without evaporation) holds by
  construction: the per-drop mass is the stream flow over `npdot` (`m = aux(ind_m)/Z(8)` in
  `Lib_RHS.f90`), so `d ∝ npdot^{-1/3}`; the test does not assert it.

## Result

`test_breakup_rd.f90` (ctest `test_breakup_rd`, GREEN): RD1 sub-critical zero;
RD2 bag + RD3 stripping vs OpenFOAM-form oracle (Cb=π/4, Cs=10=bp(4)/2), ~2e-16;
RD4 regime handoff at We_r=0.5√Re, each side matches its own branch (3.9e-13).
Re convention: gas diameter Reynolds (`Re = gas(1)*slip*d/gas(6)` in `Lib_RHS.f90`, passed to `breakupOde`), as the
stripping dStable requires. Algebra note: no dp<dStable guard needed — bag gives
dStable/d=6/We_r<1 and stripping dStable/d=Re/(4We_r²)<1 whenever entered.
