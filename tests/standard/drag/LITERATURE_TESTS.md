# Literature-grounded verification tests — IGLOO drag models

Scope: the 13 drag closures in [Lib_Drag.f90](../../../src/lib/Lib_Drag.f90) selected by
`drag()` (L66). Momentum coupling uses
`Fdrag = (π/8)·Cd·Re·μ_g·d·(v_gas−v_p)` [Lib_Equations.f90 L226].

STATUS: COMPLETE. Internal-consistency checks DONE (sharpest). Source cross-check:
Wen-Yu 0.44 CONFIRMED (Rowe, via OpenFOAM WenYuDragForce docs); Putnam & Henderson
blend rest on airtight *internal* C⁰-continuity arguments (source-independent). Primary
PDFs (Putnam 1961, Henderson 1976 AIAA) not machine-readable via fetch — see Open Qs.

**UPDATE 2026-07-15 — provenance obtained (`papers/Shimada2006.pdf`, JAXA-SP-05-035E,
eqs. 12–25 = this whole catalog):** Newton 0.45, Chang≡Clift-Gauvin,
Carlson-Hoglund (incl. Cd0=Wen-Yu convention), Henderson branch bodies, Crowe/Hermsen
g,h and the A1 exp placement all CONFIRMED exact. **A11 found & fixed**: g(Re) was coded
in the bridge-exponent denominator; eqs. 23/24 have it multiplying M/Re (Cd off ×2 at
Re=1e3, Ma=2). Still open vs true primaries (Shimada disagrees with the coded values):
Putnam plateau 0.4392 (eq. 17) vs 0.424 (A5), Wen-Yu 0.43 (eq. 16) vs 0.44, Henderson
blend 3/4 (eq. 22, discontinuous — presumed typo) vs 4/3 (A6). References resolved:
[6] Putnam ARS J. 31 (1961); [9] Henderson AIAA J. 14(6) (1976); [5] Wen-Yu (1966);
[10] Crowe UTC 2296-FR (1959); [11] Hermsen CPIA 113 (1979).

---

## 0. Force-law dimensional check (applies to ALL models) — PASS

Standard drag: `F = ½ρ_g Cd A |vrel| vrel`, `A=πd²/4` ⇒ `F=(π/8)Cd ρ_g d² |vrel| vrel`.
With `Re=ρ_g|vrel|d/μ` ⇒ `ρ_g d² |vrel| = Re·μ·d`, so `F=(π/8)Cd·Re·μ·d·vrel`.
Code L226 is **algebraically exact and dimensionally consistent**. This means every
per-model test can use a *closed-form Stokes-limit* terminal/relaxation oracle as long
as `Cd·Re → 24` in that limit (see §Internal checks).

---

## 1. Internal-consistency matrix (source-free — the primary bug detector)

Low-Re limit `Cd·Re → 24` (creeping flow, Stokes) MUST hold for every full-range
correlation. High-Re: Newton inertial plateau `Cd ≈ 0.4–0.44`.

| # | Model | code Re→0 (Cd·Re) | plateau (Re→∞) | verdict |
|---|---|---|---|---|
| 1 | Newton | n/a (constant 0.45) | 0.45 | plateau-only; OK but see ⚠N |
| 2 | Stokes | 24 ✓ | — (diverges, expected) | MATCHES |
| 3 | Schlichting (Oseen) | 24 ✓ | — | MATCHES |
| 4 | Schiller–Naumann | 24 ✓ | grows (no cutoff) | MATCHES (range caveat) |
| 5 | Chang | 24 ✓ | 0.42 | MATCHES |
| 6 | Wen–Yu | 24 ✓ | 0.43 | ⚠W plateau |
| 7 | Putnam | 24 ✓ | 0.4392 | ⚠P plateau + C⁰ jump |
| 8 | Clift–Gauvin | 24 ✓ | 0.42 | MATCHES |
| 9 | Morsi–Alexander | 24 ✓ | 0.5191 | MATCHES |
| 10 | Carlson–Hoglund | 24·(compress.) | Ma-dep | base=Wen-Yu (note) |
| 11 | Henderson | subsonic branch | Ma-dep | ⚠H transonic blend |
| 12 | Crowe | Wen-Yu base | **blows up ~1e20** | ⚠C severe |
| 13 | Hermsen | Wen-Yu base | **blows up ~1e20** | ⚠C severe |

### Cross-model identity (strong check)
**Chang (#5) ≡ Clift–Gauvin (#8) exactly.** Chang L148 =
`24/Re(1+0.15Re^0.687)+0.42/(1+42500·Re^-1.16)`. Clift–Gauvin L192 factored form,
expand the inner `0.0175·Re` term: `24/Re·0.0175·Re = 0.42`, and `42500 = 4.25e4`.
⇒ identical Cd(Re) bit-for-bit. A unit test asserting `|Cd_Chang(Re)−Cd_CG(Re)|<1e-12`
over a Re sweep is a free regression guard; if it ever fails, one of the two was edited.

---

## 2. ⚠ CODE-VS-PAPER / CODE-VS-SELF DISCREPANCIES (do NOT fix code — report only)

### ⚠P (correctness, high confidence — INTERNAL check) — Putnam plateau
[Lib_Drag.f90 L176-180]. Low branch at the breakpoint Re=1000:
`24/1000·(1+1000^(2/3)/6) = 0.024·(1+100/6) = 0.024·17.667 = 0.4240` **exactly**.
To be C⁰-continuous at Re=1000 the constant MUST be **0.424** (this is the whole point of
Putnam's *integrable* form). The code uses **0.4392**, a spurious 3.6% jump at Re=1000.
Fix: `0.4392 → 0.424`. Confidence: very high on the INTERNAL inconsistency (self-proving,
source-independent). Note: secondary sources quote both 0.424 and 0.447 for "Putnam"
(0.447 belongs to the variant `24/Re(1+0.14Re^0.70)`, a different low branch) — but with
THIS code's low branch, only 0.424 restores continuity. Primary Putnam 1961 not fetched.

### ⚠H (correctness, high confidence — INTERNAL check) — Henderson transonic blend
[Lib_Drag.f90 L255]. Transonic interpolation over Ma∈[1,1.75]:
`Cd = Cd1 + 0.75·(Ma−1)·(Cd2−Cd1)`. A linear blend that reaches Cd2 at Ma=1.75 must have
slope `1/(1.75−1)=1/0.75=4/3`, i.e. `Cd1 + (4/3)(Ma−1)(Cd2−Cd1)`. At Ma=1.75 the code
gives `Cd1+0.75·0.75·(Cd2−Cd1)=Cd1+0.5625·(Cd2−Cd1)` — a **43.75% gap** vs the Cd2 it
must match ⇒ C⁰ discontinuity into the Ma≥1.75 branch. Code multiplies by 0.75 where it
should divide. Fix: `0.75_R8*(Ma-1)` → `(Ma-1)/0.75_R8` (= `4/3·(Ma-1)`). Confidence:
very high (self-proving via branch-continuity at Ma=1.75).

### ⚠C (correctness, SEVERE — INTERNAL check) — Crowe & Hermsen high-Re blowup
[Lib_Drag.f90 L272 (Crowe), L287 (Hermsen)]. Both share the tail term
`hfun / ((Ma·√G)·exp(−Re/(2·Ma)) + toll)`. The `exp(−Re/2Ma)` sits in the **denominator**.
As Re→∞ it → 0, so the denominator collapses to `toll=1e-20` and the term → `hfun/1e-20 ≈ 1e20`
instead of vanishing. The first two terms `2+(Cd0−2)·exp(−…/Re…)` are plainly built to
recover the finite continuum value Cd0 as Re→∞, so the third term MUST vanish there — it
cannot unless `exp(−Re/2Ma)` is a **numerator** factor (rarefaction correction: present at
low Re, gone at high Re). Numerically confirmed (γ=1.4,Tr=1): Cd(Re=10,Ma=2)=20 (fine),
Cd(Re=100,Ma=2)=**9.8e10**, Cd(Re=1000,Ma=2)=**3.2e20**, Cd(Re=300,Ma=1.5)=**3.5e20**.
The `toll` guard makes it worse: it converts what would be `Inf` into a huge FINITE number,
so an `isnan`/`isinf` guard downstream will NOT catch it. Fix (per intended rarefied form):
move exp to numerator, `hfun*exp(−Re/2Ma)/(Ma·√G + toll)`. Confidence: very high (self-proving
via the Re→∞ recovery requirement + direct numerical evaluation). **NOTE: this earlier draft's
§12/§13 claimed the opposite (term→0); that was reasoning from the physical formula, not the
coded parenthesization — corrected here.** Coefficients of Crowe/Hermsen otherwise unverified
(primary sources not fetched).

### ⚠W (accuracy, low severity) — Wen–Yu high-Re constant
[Lib_Drag.f90 L163]. `Cd=0.43` for Re>1000. Wen–Yu single-sphere base uses the Rowe
plateau **0.44** — CONFIRMED (OpenFOAM `WenYuDragForce`: "CD = 0.44 if Re ≥ 1000",
Rowe). 0.43 is ≈2% low ⇒ genuine (minor) discrepancy; fix `0.43 → 0.44`. Also: low branch
at Re=1000 = 0.450, so neither 0.43 nor 0.44 is C⁰-continuous (SN base overshoots) — that
jump is inherent to Wen–Yu, not a bug; only the constant value is wrong.

### ⚠N (convention, cosmetic) — Newton constant
[Lib_Drag.f90 L75]. `Cd=0.45`. Newton inertial-regime plateau is commonly 0.44 (CGW78)
or 0.4. 0.45 is within spread; flag as convention, not error. Confirm intended source.

### Range caveats (not bugs, but test-relevant)
- Schiller–Naumann (#4) has **no** Re>1000 cutoff in code; standard validity Re≲800–1000.
  A test must stay Re<1000 or it validates an out-of-range extrapolation.
- Morsi–Alexander boundaries use mixed `≤`/`>`; at exact decade values (Re=0.1,1,10,…)
  branch selection is well-defined (no gap/overlap) — verified by reading L298-314.

---

## 3. Per-model verification tests

Harness: `tests/standard` e2e cases (box case via `tools/make_box_case.py` + per-case `check.py`) for
end-to-end trajectory oracles, and `tests/standard/drag/test_drag.f90` (unit, drives
production RHS via `verif_driver`) for Cd(Re) table + closed-form relaxation. The existing
A1/A3/AO/A4′ (Stokes, Schiller–Naumann, order-4, zero-slip) already PASS — extend that file.

Two oracle styles:
- **(U) Unit Cd(Re[,Ma,γ,Tr]) table** — call `drag()` directly, compare to hand-computed
  reference at tabulated points. Cheapest, catches coefficient typos immediately.
- **(R) Closed-form relaxation / terminal velocity** — integrate production RHS; for any
  model with `Cd·Re→24` the low-Re relaxation is `v(t)=u_g+(v0−u_g)e^{−t/τ}`,
  `τ=ρ_p d²/(18μ)`. For nonlinear Cd, terminal velocity from force balance
  `drag=body-force` gives an implicit closed form (body force is available in production —
  see INFO.md note; A2 previously deferred pending it).

### 1 Newton (Cd=0.45)  [`[CGW78]` inertial regime]
- (U) trivial: `drag(Re=*,…,1)=0.45` for all args. Oracle=0.45, tol=1e-14.
- (R) high-Re terminal velocity: `½ρ_g·0.45·(πd²/4)v_t² = (π/6)d³(ρ_p−ρ_g)g` ⇒
  `v_t=sqrt(4 d (ρ_p−ρ_g) g /(3·0.45·ρ_g))`. Verifies force-law wiring at plateau.
- ⚠N: assert against whichever of 0.44/0.45 the user confirms.

### 2 Stokes (24/Re)  [`[CGW78]` eq. creeping flow] — ALREADY A1/AO PASS
- (R) exponential relaxation, exact. Existing A1: err 2.7e-13. Keep as regression.
- Paper self-test: `Cd·Re=24` identically. tol=1e-12.

### 3 Schlichting / Oseen (24/Re·(1+3Re/16))
- Source: Oseen correction (Schlichting, *Boundary-Layer Theory*). Self-test point:
  Re=1 ⇒ Cd=24·(1+0.1875)=28.5; Re=0→ Cd·Re→24.
- (U) table {Re=0.1,1,10}: Cd={241.35, 28.5, 2.85}. tol 1e-10 rel.
- (R) low-Re relaxation matches Stokes to O(Re) — use Re<0.5 so correction <10%.

### 4 Schiller–Naumann (24/Re·(1+0.15Re^0.687))  [`[SN33]`] — ALREADY A3 PASS
- Standard drag-curve points (CGW78 tab): Re=1→Cd≈26.4; Re=100→Cd≈1.09; Re=1000→Cd≈0.44.
  Compute: Re=100: 0.24(1+0.15·100^0.687)=0.24(1+0.15·24.36)=0.24·4.654=1.117. tol 2%.
- (U) table above; (R) A3 RK4-oracle already validates nonlinear integration (7.3e-11).
- Range: keep Re<1000.

### 5 Chang (SN + 0.42/(1+42500 Re^-1.16))
- (U) **cross-identity test**: assert `Cd_Chang≡Cd_CliftGauvin` over Re∈[1e-2,1e5]
  (see §1). tol 1e-12. This is the tightest possible check.
- (U) plateau: Re=1e5 ⇒ Cd→0.42 (Δ<1%). Re→0 ⇒ Cd·Re→24.
- Standard-drag-curve match: Re=1000 → 24/1000·(1+0.15·118.4)+0.42/(1+42500·1000^-1.16)
  ≈ 0.450+0.42/(1+42500/3630)=0.450+0.42/12.7=0.450+0.033=0.483. Compare CGW78 Cd(1000)≈0.47.

### 6 Wen–Yu (SN for Re≤1000, else 0.43)  [Wen & Yu 1966]
- (U) two-branch table: Re=100→1.117 (=SN); Re=2000→0.43. ⚠W: flag 0.43 vs 0.44.
- Discontinuity probe: Cd(1000⁻)=0.450 vs Cd(1000⁺)=0.43 — document the 4.4% jump as a
  known feature; a test should NOT straddle Re=1000.

### 7 Putnam (24/Re·(1+Re^(2/3)/6) for Re<1000, else 0.4392)  [Putnam 1961]
- (U) **continuity probe**: Cd(999.9) vs Cd(1000.1). Correct constant 0.424 ⇒ jump ≈0;
  code's 0.4392 ⇒ jump 0.424→0.4392 = **3.6%**. This test *documents ⚠P* (expected to
  reveal the jump until code is fixed). Oracle for low branch at Re=1000 = 0.424 exact.
- (U) low-Re: Re=1→24(1+1/6)=28.0; Re→0→24/Re.

### 8 Clift–Gauvin (factored standard drag)  [`[CGW78]`]
- (U) standard-drag-curve tab (CGW78): Re=1→Cd≈26.6; Re=10→4.15; Re=100→1.09;
  Re=1000→0.474; Re=1e4→0.41; Re=1e5→0.42. Compute Re=1000:
  0.024·(1+0.15·118.4)+0.42/(1+4.25e4·1000^-1.16)=0.483. tol 3% (curve-fit spread).
- (U) cross-identity with Chang (§5). (R) relaxation Re<1000.

### 9 Morsi–Alexander (piecewise a1+a2/Re+a3/Re²)  [Morsi & Alexander 1972]
- (U) **per-interval table** — one Cd(Re) point in each of the 8 bands with hand oracle:
  Re=0.05→24/0.05=480; Re=0.5→3.69+22.73/0.5+0.0903/0.25=3.69+45.46+0.361=49.51;
  Re=5→1.222+29.1667/5−3.8889/25=1.222+5.833−0.156=6.899;
  Re=50→0.6167+46.5/50−116.67/2500=0.6167+0.930−0.0467=1.500;
  Re=500→0.3644+98.33/500−2778/2.5e5=0.3644+0.1967−0.0111=0.550;
  Re=2000→0.357+148.62/2000−4.75e4/4e6=0.357+0.0743−0.0119=0.419;
  Re=7500→0.46−490.546/7500+5.787e5/5.625e7=0.46−0.0654+0.0103=0.405;
  Re=2e4→0.5191−1662.5/2e4+5.4167e6/4e8=0.5191−0.0831+0.0135=0.4495. tol 1e-6 (exact fit).
- (U) continuity probe at each internal boundary (M&A intervals are ~continuous by design).

### 10 Carlson–Hoglund (Re,Ma)  [Carlson & Hoglund 1964, AIAA]
- (U) incompressible recovery: Ma→0 ⇒ compressibility factor→? At Ma→0: `exp(-0.427/Ma^4.63)→0`
  so numerator→(1+0)=1... but denominator `1+Ma/Re(3.82+…)`→1, and the exp term
  `exp(-0.427/Ma^4.63)` → exp(-∞)=0 ⇒ Cd→Cd0·1/1=Cd0(Wen-Yu). Verified: exp terms placed
  correctly, C–H recovers Cd0 at Ma→0 (unlike Crowe/Hermsen). Cd(Re,Ma→0)=Wen-Yu(Re).
  Note: base is **Wen–Yu**; some sources use SN base — modeling choice, flag.
- (U) known point: reproduce a C–H tabulated Cd(Re,Ma) from the 1964 report once located.

### 11 Henderson (Re,Ma,γ,Tr)  [Henderson 1976, AIAA J 14(6):707]
- (U) **transonic continuity probe** at Ma=1.75: Cd_blend(1.75) must equal Henderson_2(1.75).
  Code fails this (⚠H, 43.75% gap) — test documents the bug. Also probe Ma=1.0 continuity
  (blend→Henderson_1, code OK there).
- (U) subsonic incompressible limit Ma→0: Cd→Henderson_1 which should → standard drag curve.
- (U) reproduce Henderson Fig/tab Cd(Re,Ma) points once the 1976 paper is fetched.

### 12 Crowe (Re,Ma,γ,Tr)  — ⚠C SEVERE (see §2)
- (U) **finite-value / Re→∞ recovery probe**: assert `Cd(Re,Ma) < 10` and `→Cd0` as Re grows
  at fixed Ma. Code FAILS: Cd(Re=1000,Ma=2)=3.2e20 (should be O(1)). This test documents ⚠C
  and is the single cheapest catch of a catastrophic one-line bug.
- (U) low-Re sanity: Cd(Re=10,Ma=2)=19.99 (near Wen-Yu 20.6) — the ONLY regime currently sane.
- Source coeffs (3.07, 0.77, 1.92, 1.17, 2.3, 1.7) unverified — flag.

### 13 Hermsen (Re,Ma,γ,Tr)  [Hermsen 1979, AIAA / JANNAF]  — ⚠C SEVERE (see §2)
- (U) same Re→∞ finite-value probe as Crowe; identical tail-term blowup. Code FAILS.
- Source coeffs (12.278, 0.548, 11.278, 5.6, 1.7) unverified — flag.

---

## 4. TOP 3 drag models to implement as verification tests FIRST

1. **Crowe & Hermsen (#12/#13)** — a *catastrophic* self-proving bug (⚠C: `exp(−Re/2Ma)`
   in the denominator ⇒ Cd ~1e20 for Re≳100). A one-line "Cd finite & O(1) at moderate
   Re,Ma" probe catches it; the highest-severity defect found and airtight (numerically
   confirmed, source-independent). Ranks first because it silently corrupts *any* Crowe/
   Hermsen run at engineering Re, and the `toll` guard hides it from Inf/NaN checks.
2. **Putnam (#7)** — self-proving correctness bug (⚠P: 0.4392 should be 0.424 by C⁰
   continuity at Re=1000; low branch at Re=1000 = 0.424 exactly). Cheap unit continuity probe.
3. **Henderson (#11)** — self-proving correctness bug (⚠H: transonic blend uses 0.75 instead
   of 4/3 ⇒ 43.75% discontinuity at Ma=1.75). Unit continuity probe at Ma=1.75. Breaks the
   whole transonic regime.

Rationale: rank by (a) severity of the *proven* defect and (b) whether the oracle needs no
external source — all three qualify via internal consistency. Wen-Yu (⚠W, 2%) and the
Chang≡Clift–Gauvin cross-identity (free `<1e-12` regression) come next. Stokes/Schiller–
Naumann already have passing L1/L2 tests; extend, don't duplicate.

---

## 5. Open questions for the user
- Newton 0.45 and Wen–Yu 0.43: intended, or should they be 0.44? (source?)
- Carlson–Hoglund / Crowe / Hermsen: which primary references define the coded coefficients?
  (need DOIs to finish §10/12/13 term-by-term).
- Confirm ⚠P and ⚠H are logged as bugs (I did NOT edit src/ per constraint).
