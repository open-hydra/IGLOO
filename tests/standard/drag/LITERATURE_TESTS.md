# Literature-grounded verification tests — IGLOO drag models

Scope: the 13 drag closures in [Lib_Drag.f90](../../../src/lib/Lib_Drag.f90) selected by
`drag()`. Momentum coupling uses
`Fdrag = (π/8)·Cd·Re·μ_g·d·(v_gas−v_p)` [Lib_Equations.f90, `interphase`].

Sources. `[Shim06]` (JAXA-SP-05-035E) eqs. 12–25 hold this whole catalog: Newton 0.45, Chang≡Clift-Gauvin,
Carlson-Hoglund (incl. Cd0=Wen-Yu convention), Henderson branch bodies, Crowe/Hermsen g,h and their
`exp(−Re/2Ma)` factor match it exactly, and so does the place of g(Re), multiplying M/Re (eqs. 23/24; in the
bridge-exponent denominator Cd is off ×2 at Re=1e3, Ma=2). The primaries of Putnam, Henderson and Wen-Yu are
unavailable, so `[Shim06]` is authoritative: Putnam plateau 0.4392 (eq. 17; C⁰ continuity would give 0.424, §2),
Wen-Yu 0.43 (eq. 16; Rowe's plateau, via OpenFOAM's WenYuDragForce docs, is 0.44). Henderson's blend is the
exception: the code's 4/3 against eq. 22's 3/4 (discontinuous — a typo, §2). The C⁰-continuity arguments of §1-§2
are *internal* (source-independent). References: [6] Putnam ARS J. 31 (1961); [9] Henderson AIAA J. 14(6) (1976);
[5] Wen-Yu (1966); [10] Crowe UTC 2296-FR (1959); [11] Hermsen CPIA 113 (1979).

---

## 0. Force-law dimensional check (applies to ALL models) — PASS

Standard drag: `F = ½ρ_g Cd A |vrel| vrel`, `A=πd²/4` ⇒ `F=(π/8)Cd ρ_g d² |vrel| vrel`.
With `Re=ρ_g|vrel|d/μ` ⇒ `ρ_g d² |vrel| = Re·μ·d`, so `F=(π/8)Cd·Re·μ·d·vrel`.
The code's form is **algebraically exact and dimensionally consistent**. This means every
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
| 9 | Morsi–Alexander | 24 ✓ | 0.5191 | MATCHES (`[MA72]` appendix, all eight bands) |
| 10 | Carlson–Hoglund | 24·(compress.) | Ma-dep | base=Wen-Yu (note) |
| 11 | Henderson | subsonic branch | Ma-dep | MATCHES (blend 4/3, ⚠H) |
| 12 | Crowe | Wen-Yu base | Ma-dep, finite | MATCHES `[Shim06]` eq. 23 (⚠C) |
| 13 | Hermsen | Wen-Yu base | Ma-dep, finite | MATCHES `[Shim06]` eq. 24 (⚠C) |

### Cross-model identity (strong check)
**Chang (#5) ≡ Clift–Gauvin (#8).** Chang (`drag_Chang`) =
`24/Re(1+0.15Re^0.687)+0.42/(1+42500·Re^-1.16)`. Clift–Gauvin's factored form,
expand the inner `0.0175·Re` term: `24/Re·0.0175·Re = 0.42`, and `42500 = 4.25e4`.
⇒ the same Cd(Re), the two forms placing `toll` differently. DL2 (`test_drag_lit`) asserts
`|Cd_Chang(Re)−Cd_CG(Re)|<1e-12` (relative) over a Re sweep: if it fails, one of the two has been edited.

---

## 2. Code vs paper and code vs self

### ⚠P — Putnam plateau (source-faithful, C⁰ jump at Re=1000)
[`drag_Putnam`]. Low branch at the breakpoint Re=1000:
`24/1000·(1+1000^(2/3)/6) = 0.024·(1+100/6) = 0.024·17.667 = 0.4240` **exactly**.
C⁰ continuity at Re=1000 needs the constant **0.424** (this is the whole point of Putnam's *integrable* form).
The code takes `[Shim06]` eq. 17's **0.4392**: a 3.6% jump at Re=1000, a property of the published form. Note:
secondary sources quote both 0.424 and 0.447 for "Putnam" (0.447 belongs to the variant `24/Re(1+0.14Re^0.70)`,
a different low branch) — but with THIS low branch, only 0.424 is continuous. XD3 pins 0.4392.

### ⚠H — Henderson transonic blend (4/3, continuous)
[`drag_Henderson`]. Transonic interpolation over Ma∈[1,1.75]: a linear blend that reaches Cd2 at Ma=1.75 has
slope `1/(1.75−1)=1/0.75=4/3`, and the code's is `Cd = Cd1 + (4/3)(Ma−1)(Cd2−Cd1)`. `[Shim06]` eq. 22 prints
`Cd = Cd1 + 0.75·(Ma−1)·(Cd2−Cd1)`: at Ma=1.75 that gives `Cd1+0.75·0.75·(Cd2−Cd1)=Cd1+0.5625·(Cd2−Cd1)` — a
**43.75% gap** vs the Cd2 it must match ⇒ C⁰ discontinuity into the Ma≥1.75 branch, so eq. 22's 0.75 is read
as a typo. XD4 pins the continuity at Ma=1.75.

### ⚠C — Crowe & Hermsen rarefaction term (a factor, vanishing at high Re)
[`drag_Crowe`, `drag_Hermsen`]. Both share the tail term `hfun/(Ma·√G)·exp(−Re/(2·Ma))` (`[Shim06]` eqs.
23/24). The first two terms `2+(Cd0−2)·exp(−…/Re…)` are plainly built to recover the finite continuum value
Cd0 as Re→∞, so the third term MUST vanish there, as it does with `exp(−Re/2Ma)` a **numerator** factor
(rarefaction correction: present at low Re, gone at high Re). With the exp in the **denominator**,
`hfun / ((Ma·√G)·exp(−Re/(2·Ma)) + toll)`, the denominator collapses to `toll=1e-20` as Re→∞ and the term →
`hfun/1e-20 ≈ 1e20` instead of vanishing: such a build (γ=1.4,Tr=1) gives Cd(Re=10,Ma=2)=20,
Cd(Re=100,Ma=2)=**9.8e10**, Cd(Re=1000,Ma=2)=**3.2e20**, Cd(Re=300,Ma=1.5)=**3.5e20**, and the `toll` guard
converts what would be `Inf` into a huge FINITE number, so an `isnan`/`isinf` guard downstream does NOT catch
it. XD1/XD2 pin a finite Cd at Re=1000, Ma=2. The coefficients match `[Shim06]`; the primaries are not checked.

### ⚠W — Wen–Yu high-Re constant (source-faithful)
[`drag_Wen_Yu`]. `Cd=0.43` for Re>1000, `[Shim06]` eq. 16. Wen–Yu single-sphere base uses the Rowe
plateau **0.44** (OpenFOAM `WenYuDragForce`: "CD = 0.44 if Re ≥ 1000", Rowe); the primary is unavailable, so
the code follows `[Shim06]`. The low branch at Re=1000 = 0.438, so neither 0.43 nor 0.44 is C⁰-continuous —
that jump is inherent to Wen–Yu. XD5 pins 0.43 and the step.

### ⚠N (convention, cosmetic) — Newton constant
[`drag()`, model 1]. `Cd=0.45` (`[Shim06]` eq. 13; DL4). Newton inertial-regime plateau is commonly 0.44
(CGW78) or 0.4. 0.45 is within spread: a convention, not an error.

### Range caveats (not bugs, but test-relevant)
- Schiller–Naumann (#4) has **no** Re>1000 cutoff in code; standard validity Re≲800–1000.
  A test must stay Re<1000 or it validates an out-of-range extrapolation.
- Morsi–Alexander boundaries: the paper's inequalities are strict on both sides of every edge; the code takes each
  edge (Re=0.1,1,10,…,10⁴) in the lower band, so branch selection is well-defined (no gap/overlap).

---

## 3. Per-model verification tests

Harness: `tests/standard` e2e cases (box case via `tools/make_box_case.py` + per-case `check.py`) for
end-to-end trajectory oracles; `tests/standard/drag/test_drag.f90` (unit, drives the production RHS via
`verif_driver`: A1/A3/AO/A4′ for Stokes, Schiller–Naumann, order 4, zero slip) for closed-form relaxation;
`test_drag_lit.f90` (DL1–DL8) for the Cd(Re) tables; `test_drag_probes.f90` (XD1–XD5) for the pins of §2.

Two oracle styles:
- **(U) Unit Cd(Re[,Ma,γ,Tr]) table** — call `drag()` directly, compare to hand-computed
  reference at tabulated points. Cheapest, catches coefficient typos immediately.
- **(R) Closed-form relaxation / terminal velocity** — integrate production RHS; for any
  model with `Cd·Re→24` the low-Re relaxation is `v(t)=u_g+(v0−u_g)e^{−t/τ}`,
  `τ=ρ_p d²/(18μ)`. For nonlinear Cd, terminal velocity from force balance
  `drag=body-force` gives an implicit closed form (the body force is in production:
  `[IGLOO-General] body-accel`).

### 1 Newton (Cd=0.45)  [`[CGW78]` inertial regime]
- (U) trivial: `drag(Re=*,…,1)=0.45` for all args. Oracle=0.45, tol=1e-14: DL4.
- (R) high-Re terminal velocity: `½ρ_g·0.45·(πd²/4)v_t² = (π/6)d³(ρ_p−ρ_g)g` ⇒
  `v_t=sqrt(4 d (ρ_p−ρ_g) g /(3·0.45·ρ_g))`. Verifies force-law wiring at plateau.
- ⚠N: 0.45 is `[Shim06]` eq. 13's (§2).

### 2 Stokes (24/Re)  [`[CGW78]` eq. creeping flow] — A1, AO, DL1
- (R) exponential relaxation, exact: A1 (`test_drag`), measured err 2.7e-13.
- Paper self-test: `Cd·Re=24` identically. tol=1e-12; DL1 holds it to 1e-15.

### 3 Schlichting / Oseen (24/Re·(1+3Re/16))
- Source: Oseen correction (Schlichting, *Boundary-Layer Theory*). Self-test point:
  Re=1 ⇒ Cd=24·(1+0.1875)=28.5; Re=0→ Cd·Re→24.
- (U) table {Re=0.1,1,10}: Cd={244.5, 28.5, 6.9}. tol 1e-10 rel. DL5 holds the form at Re 0.1, 1, 4 to 1e-14.
- (R) low-Re relaxation matches Stokes to O(Re) — use Re<0.5 so correction <10%.

### 4 Schiller–Naumann (24/Re·(1+0.15Re^0.687))  [`[SN33]`] — A3, DL3
- Standard drag-curve points (CGW78 tab): Re=1→Cd≈26.4; Re=100→Cd≈1.09; Re=1000→Cd≈0.44.
  Compute: Re=100: 0.24(1+0.15·100^0.687)=0.24(1+0.15·23.66)=0.24·4.549=1.092. tol 2%.
- (U) table above; (R) A3's RK4 oracle validates the nonlinear integration (7.3e-11).
- Range: keep Re<1000.

### 5 Chang (SN + 0.42/(1+42500 Re^-1.16))
- (U) **cross-identity test**: assert `Cd_Chang≡Cd_CliftGauvin` over Re∈[0.1,1e5]
  (see §1). tol 1e-12. This is the tightest possible check: DL2, on Re 0.1..1e5; DL8 against `[Shim06]` eq. 18.
- (U) plateau: Cd→0.42 as Re→∞. Re→0 ⇒ Cd·Re→24.
- Standard-drag-curve match: Re=1000 → 24/1000·(1+0.15·115.1)+0.42/(1+42500·1000^-1.16)
  ≈ 0.438+0.42/(1+42500/3020)=0.438+0.42/15.07=0.438+0.028=0.466. Compare CGW78 Cd(1000)≈0.47.
- (R) e2e `drag-chang`: the slip relaxation against a quadrature of the form typed in its check.

### 6 Wen–Yu (SN for Re≤1000, else 0.43)  [Wen & Yu 1966; `[Shim06]` eq. 16]
- (U) two-branch table: Re=100→1.092 (=SN, DL3); Re=2000→0.43. ⚠W: 0.43 vs Rowe's 0.44 (§2).
- Discontinuity probe: Cd(1000⁻)=0.438 vs Cd(1000⁺)=0.43 — a 1.9% jump, a feature of the published form,
  which XD5 pins; a test of the low branch should NOT straddle Re=1000.

### 7 Putnam (24/Re·(1+Re^(2/3)/6) for Re<1000, else 0.4392)  [Putnam 1961; `[Shim06]` eq. 17]
- (U) **continuity probe**: Cd(999.9) vs Cd(1000.1). A constant of 0.424 gives a jump ≈0;
  the code's 0.4392 gives 0.424→0.4392 = **3.6%**, the published form's (⚠P); XD3 pins the plateau.
  Oracle for low branch at Re=1000 = 0.424 exact.
- (U) low-Re: Re=1→24(1+1/6)=28.0; Re→0→24/Re.

### 8 Clift–Gauvin (factored standard drag)  [`[CGW78]`]
- (U) standard-drag-curve tab (CGW78): Re=1→Cd≈26.6; Re=10→4.15; Re=100→1.09;
  Re=1000→0.474; Re=1e4→0.41; Re=1e5→0.42. Compute Re=1000:
  0.024·(1+0.15·115.1)+0.42/(1+4.25e4·1000^-1.16)=0.466. tol 3% (curve-fit spread).
- (U) cross-identity with Chang (§5): DL2; DL8 against `[Shim06]` eq. 18. (R) relaxation Re<1000.

### 9 Morsi–Alexander (piecewise a1+a2/Re+a3/Re²)  [`[MA72]`]
- **Source.** `[MA72]`'s appendix (printed p. 207) states the eight bands as C_d = K1/R_N + K2/R_N² + K3 (eq. 4,
  p. 195); the code's (a1, a2, a3) are (K3, K1, K2) and equal the print in all eight bands. (K1, K2, K3):
  (24.0, 0, 0) for R_N < 0.1; (22.73, 0.0903, 3.69) 0.1–1; (29.1667, −3.8889, 1.222) 1–10; (46.5, −116.67, 0.6167)
  10–100; (98.33, −2778, 0.3644) 100–1000; (148.62, −4.75×10⁴, 0.357) 1000–5000; (−490.546, 57.87×10⁴, 0.46)
  5000–10⁴; (−1662.5, 5.4167×10⁶, 0.5191) 10⁴–5×10⁴.
- (U) **per-interval table** — one Cd(Re) point in each of the 8 bands with hand oracle:
  Re=0.05→24/0.05=480; Re=0.5→3.69+22.73/0.5+0.0903/0.25=3.69+45.46+0.361=49.51;
  Re=5→1.222+29.1667/5−3.8889/25=1.222+5.833−0.156=6.899;
  Re=50→0.6167+46.5/50−116.67/2500=0.6167+0.930−0.0467=1.500;
  Re=500→0.3644+98.33/500−2778/2.5e5=0.3644+0.1967−0.0111=0.550;
  Re=2000→0.357+148.62/2000−4.75e4/4e6=0.357+0.0743−0.0119=0.419;
  Re=7500→0.46−490.546/7500+5.787e5/5.625e7=0.46−0.0654+0.0103=0.405;
  Re=2e4→0.5191−1662.5/2e4+5.4167e6/4e8=0.5191−0.0831+0.0135=0.4495. tol 1e-6 (exact fit).
- **Table 1** (p. 195: C_d experimental and calculated at 29 R_N from 0.1 to 5×10⁴) does not reproduce the
  appendix at 13 of its 29 points beyond half a unit of its last printed digit (worst R_N 7000: 0.391 against
  0.401732, 2.7 %; R_N 200: 0.771 against 0.786600, 2.0 %). The appendix is the paper's statement of the constants
  (p. 195), so it is the pin; Table 1 is not an oracle.
- **Range and edges.** The paper fits 0.1 < R_N < 5×10⁴ with strict inequalities on both sides of every edge; the
  code takes each edge in the lower band (the fit jumps there by 5.3e-5 to 2.3e-2 relative, the largest at 10⁴; not
  gated) and continues the last band past 5×10⁴ (C_d 0.488017 at 5×10⁴, 0.503017 at 10⁵, 0.517443 at 10⁶).

### 10 Carlson–Hoglund (Re,Ma)  [Carlson & Hoglund 1964, AIAA; `[Shim06]` eq. 19]
- (U) incompressible recovery: at Ma→0 `exp(-0.427/Ma^4.63)`→0, so the numerator `1+exp(…)`→(1+0)=1 and the
  denominator `1+Ma/Re(3.82+…)`→1 ⇒ Cd→Cd0·1/1=Cd0(Wen-Yu): the exp terms are placed correctly,
  Cd(Re,Ma→0)=Wen-Yu(Re). Note: base is **Wen–Yu** (`[Shim06]`'s convention); some sources use SN base.
- (U) DL6: against `[Shim06]` eq. 19 on Re 1..1000 × Ma 0.1..2, and Cd = Wen-Yu at Ma = 0.

### 11 Henderson (Re,Ma,γ,Tr)  [Henderson 1976, AIAA J 14(6):707]
- (U) **transonic continuity probe** at Ma=1.75: Cd_blend(1.75) must equal Henderson_2(1.75): XD4 (§2, ⚠H).
  At Ma=1.0 the blend meets Henderson_1 by construction.

### 12 Crowe (Re,Ma,γ,Tr)  [`[Shim06]` eq. 23]
- (U) **finite-value probe**: `Cd(Re=1000,Ma=2) < 10` (O(1), recovering Cd0 as Re grows at fixed Ma): XD1
  (§2, ⚠C).
- Source coeffs (3.07, 0.77, 1.92, 1.17, 2.3, 1.7) as `[Shim06]` prints them; the primary is not checked.

### 13 Hermsen (Re,Ma,γ,Tr)  [Hermsen 1979, AIAA / JANNAF; `[Shim06]` eq. 24]
- (U) the same finite-value probe as Crowe (the same tail term): XD2.
- Source coeffs (12.278, 0.548, 11.278, 5.6, 1.7) as `[Shim06]` prints them; the primary is not checked.
