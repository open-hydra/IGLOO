# INFO — evaporation/tc-analytic (F2: Tonini-Cossali analytical evaporation)

Verifies the `evaporation = TC` selector (`Lib_Evaporation.f90::TC_model`,
phase F2 of plan-bucket/evaporation-models-implementation-plan.md): the TC2012
variable-density Stefan-Fuchs analytical single-droplet rate, referenced as
"Eq (9)" in `[ATC24]`. The production Newton residual and the coded rate are
checked against an independently transcribed form of the pinned equation
(plan-bucket/f2-tc-derivation.md):

    rhs0 = (Mv/Minf)·ln[(1−Xinf)/(1−Xs)]                 (molar Stefan driver)
    G(m) = m + (Ts/Tg − 1)·Lev·(f(m/Lev) − 1) − rhs0,   f(x)=x/(1−e^−x)
    mdot = −π·d·ρg·Dv·Sh·m,   Sh = 2 + 0.6·Re^½·Sc^⅓     (Ranz-Marshall bolt-on)
    Qdot = −mdot·cpv·(Tg−Tp)/(e^χ − 1),  χ = −mdot·cpv/(π·d·kg·Nu)

`Lev = kg/(cpv·Dv·ρg)` is the vapor-cp Lewis number (= Le·cpg/cpv); the input
`Le` is the gas-cp film Lewis, converted internally. The e2e companion is
`tc-box/` (a uniform-gas box with `interface = VLE`).

Citation tags resolve in [../../REFERENCES.md](../../REFERENCES.md).

## Tests

| id | mode | source | locus | compared | status |
|---|---|---|---|---|---|
| TC1 | analytic | pinned Eq (9*) transcription | residual self-oracle: code m re-inserted into the transcribed G over a 27-pt (Tp,ρ,Le) grid; mdot < 0 everywhere | \|G\|/max(1,rhs0) ≤ 2e-13 | green |
| TC2 | analytic (corners) | design | boiling clamp (psat>p) and cpv=0 fallback stay finite and in tol | 2e-13 | green |
| TC3 | analytic (limit) | `[TC2012]` classical Stefan-Fuchs | isothermal (Tp=Tg) closed form m=rhs0 exact, Re=0 and Re=100 | rel 1e-14 | green |
| TC4 | analytic (limit) | derived | Mv=Mg (molar≡mass) + isothermal ⇒ TC ≡ CEM (Sh=2), Re=0 and Re=100 | rel 1e-13 | green |
| TC5 | analytic | own-BT transcription | heat: Qdot matches χ-form; override on; Stefan-blocked Qdot < π·d·kg·Nu·ΔT | rel 1e-14 | green |
| TC6 | analytic (theorem) | monotonicity bracket | cold drop rhs0 < m ≤ rhs0/Tts; hot drop rhs0/Tts ≤ m < rhs0 (variable-density direction) | exact bounds | green |
| TC7 | INFORMATIONAL | `[TC2012]` | TC vs ASM / CEM deviation rows | — | info |

## Notes

- The unit test and the tc-box oracle both solve the same implicit G, so the
  suite proves plumbing + limits, not the closure itself. The non-circular
  anchor for F2 is the independent Stefan-Fuchs re-derivation
  (plan-bucket/f2-tc-derivation.md) + the four analytic limits above, verified
  against the `[TC2012]`/`[ATC24]` primaries (source PDFs read 2026-07-13).
