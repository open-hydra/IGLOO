# LITERATURE_TESTS — heat-transfer (Nusselt) models (paper-grounded verification)

Companion to [INFO.md](INFO.md) (B1/BO/B2/B3 schema). LITERATURE side of the convective
heat closure: per production Nu-correlation — original paper, the paper's own verification
test, a reproducible-in-IGLOO design, and a term-by-term CODE-vs-PAPER cross-check of
[../../../src/lib/Lib_Heat.f90](../../../src/lib/Lib_Heat.f90) (selector `heat()` L45).

**Confirmed vs inferred** stated per claim. Coefficients re-derived / re-verified here, not
recalled. **DO NOT modify src/** — coded≠paper is flagged, never patched.

**UPDATE 2026-07-15 — provenance obtained (`papers/Shimada2006.pdf`, JAXA-SP-05-035E,
eqs. 40–50 = this whole catalog; the "JAXA1–4" names resolve to NASA-SP-8039 (1971), its
ref [12]):** JAXA2 (eq. 46), JAXA3 (eq. 47), Ranz-Marshall (eq. 48) and
Kavanau-Drake (eq. 49) CONFIRMED exact — eq. 49 is EXPLICIT and prints `Pr^0.33`, refuting
the D-HEAT-2/3 flags. **A7 closed as source-faithful**: eq. 45 has no `+2` floor — a
documented limitation of the published correlation, not a transcription bug (xfail probe
promoted to a transcription pin). **A12 found & fixed**: JAXA4 base 0.645 → 0.654
(eq. 50, transposed digits). Kavanau-Drake primary: Univ. of California Rept. HE-150-108
(1953), ref [14].

---

## 0. How Nu enters the ODE (shared oracle skeleton) — VERIFIED

[Lib_Equations.f90#L228](../../../src/lib/Lib_Equations.f90#L228):

    Qdot = Nu * k_g * pi * d * (T_g - T_p) * cpFactor          (cpFactor = 1/c_p,p, L399)
    F(7) = Qdot / m ,   m = rho_p (pi/6) d^3

⇒ **dT_p/dt = (6 Nu k_g)/(c_p,p rho_p d^2)·(T_g − T_p) = (T_g − T_p)/tau_T**, with

    tau_T = c_p,p rho_p d^2 / (6 Nu k_g)          [GENERAL — re-derived]

Re-derivation: `Qdot/m = Nu k_g pi d (T_g−T_p)/c_p ÷ rho_p(pi/6)d^3 = 6 Nu k_g(T_g−T_p)/(c_p rho_p d^2)`.
This is standard lumped-capacitance: `m c_p dT/dt = h A (T_g−T_p)`, `h=Nu k_g/d`, `A=pi d^2`.
At Nu=2 ⇒ `tau_T = c_p rho_p d^2/(12 k_g)`, matching e2e temp-relax
([check.py#L17](../temp-relax/check.py) `TAU=CP*RHO_P*D**2/(12 K_G)`) and INFO.md B1. **CONFIRMED.**

> ⚠ The assignment text stated `tau_T = c_p rho_p d^2/(2 Nu k_g)`. **Off by 3×** — correct
> denominator is **6 Nu k_g** (the only form consistent with the /12 k_g at Nu=2 the message
> itself cites, and with the passing e2e case). NOT a code bug; the code is correct — the
> shorthand formula in the brief was. Every design below uses `6 Nu k_g`.

`Pr` inside the closure ([#L230](../../../src/lib/Lib_Equations.f90#L230)):
`Pr = mu_g gamma R_g/((gamma−1)k_g) = c_p,g mu_g/k_g` (since `c_p,g=gamma R_g/(gamma−1)`). Standard.
`Ma = slip/sqrt(gamma R_g T_g)`.

---

## Internal (source-free) checks — with the compressible-model caveat

- **Stagnant-sphere conduction floor:** steady conduction from an isothermal sphere radius `a`
  into an infinite quiescent medium gives `h=k/a` ⇒ **Nu=2**. So `Nu(Re→0)→2` is the physical
  anchor. **BUT this applies only to the incompressible/continuum forms** (JAXA2, JAXA3,
  Ranz-Marshall). See next bullet.
- **⚠ The blanket "Nu→2 as Re→0" is a FALSE-POSITIVE generator for the two rarefied models
  (JAXA4, Kavanau-Drake).** In a coupled trajectory the rarefaction group is
  `Ma/Re = mu_g/(rho_g d sqrt(gamma R_g T_g)) ≈ Kn` — **fixed by particle size + gas state,
  independent of slip.** As the particle decelerates, `Re→0` but `Ma/(Re Pr)` stays constant, so
  `Nu → 2/(1 + 6.84·Ma/(Re Pr)) < 2`, never 2. That is physical (continuum→slip corner), NOT a
  bug. `Nu→2` for these models requires BOTH `Re→0` AND `Kn→0` (large `d`, high `rho_g`).
  *Edge note:* at EXACTLY slip=0 the code's `Ma` numerator is 0, so `Ma/(Re Pr+toll)=0/toll=0` ⇒
  Nu snaps to the continuum value 2 — a genuine discontinuity vs the slip→0⁺ limit, from the
  `toll=1e-20` guard. The e2e B1 zero-slip case sits exactly on this point, so it (correctly)
  sees Nu=2 for all six models.
- **JAXA1 →0 at Re→0 IS a genuine flag** (it is Ma-free, so no rarefaction alibi) — see D-HEAT-1.
- **Monotonic in Re** at fixed Pr; **Nu>0** for all admissible inputs.

---

## Provenance map (what the labels REALLY are)

| selector | code Nu | true origin (confidence) |
|---|---|---|
| JAXA1 (1) | `2.5 Re^0.15 + 0.04 Re` | **UNIDENTIFIED** JAXA report; high-speed/rarefied fit, no conduction floor. UNVERIFIABLE — do not cite |
| JAXA2 (2) | `2 + 0.37 Re^0.6 Pr^(1/3)` | NOT Whitaker (his is `2+(0.4Re^0.5+0.06Re^(2/3))Pr^0.4`), NOT Ranz-Marshall. Origin UNIDENTIFIED — UNVERIFIABLE |
| JAXA3 (3) | `2 + 0.459 Re^0.55 Pr^(1/3)` | **Kavanau (1955) continuum limit** — same base as Kavanau-Drake, rarefaction off (inferred, high confidence via internal algebra) |
| JAXA4 (4) | `[(2+0.645 Re^0.5 Pr^(1/3))^-1 + 3.42 Ma/(RePr)]^-1` | **Kavanau-Drake bridging, algebraically identical** to model 6 (proven below); base 0.645 Re^0.5 UNVERIFIED |
| Ranz-Marshall (5) | `2 + 0.6 Re^0.5 Pr^(1/3)` | **Ranz & Marshall 1952 `[RM52]`** — 0.6 CONFIRMED (textbook + evap review); valid 0<Re<200 |
| Kavanau-Drake (6) | `(2+0.459 Re^0.55 Pr^0.33)/(1+3.42 Ma/(RePr)·Nu_0)` | **Kavanau 1955 continuum + Drake rarefaction**, the plasma-spray form; 0.459/0.55/3.42 INFERRED (literature-consensus, not primary-verified this session) |

Suggested REFERENCES.md additions (verify DOIs before citing):
`[Kav55]` Kavanau, L.L., "Heat transfer from spheres to a rarefied gas in subsonic flow,"
*Trans. ASME* 77 (1955) 617–623 (no DOI, historical). `[RM52]` already listed.

---

## Source-free spine (strongest findings — need NO external citation)

**S1 — JAXA4 ≡ Kavanau-Drake bridging (algebraically identical).** Re-derived:
JAXA4 `Nu = [Nu_0^-1 + k]^-1 = Nu_0/(1+k·Nu_0)` with `k=3.42 Ma/(RePr)`; Kavanau-Drake
[#L134](../../../src/lib/Lib_Heat.f90#L134) `Nu = Nu_0/(1+k·Nu_0)` — **same functional form**,
differing ONLY in the continuum base (`0.645 Re^0.5 Pr^(1/3)` vs `0.459 Re^0.55 Pr^0.33`).
Both are **explicit** (single eval, continuum `Nu_0` in the denominator). *The brief called KD
"implicit"; as coded it is explicit.* ⚠ open question: does the source bridging put the **total**
Nu or the **continuum** Nu in the denominator? Code uses continuum — flag for owner.

**S2 — JAXA3 ≡ Kavanau-Drake continuum base** (`2+0.459 Re^0.55·Pr^{1/3 or 0.33}`), rarefaction
off, Pr exponent `1/3` vs `0.33` the only difference. See D-HEAT-2.

**S3 — Cross-model continuum-corner consistency (Re→0, large d ⇒ Kn→0):** JAXA2, JAXA3,
Ranz-Marshall, and (in the Kn→0 corner) JAXA4 & KD must ALL give Nu=2. JAXA1 gives 0 ⇒ isolates D-HEAT-1.
No citation needed — the oracle is inter-model agreement + the analytic Nu=2 floor.

---

## Per-model detail

### Model 5 — Ranz-Marshall (heatSelect=5) — the reference case
**(a) Paper.** Ranz, W.E. & Marshall, W.R., "Evaporation from drops I & II," *Chem. Eng. Prog.*
48 (1952) 141–146, 173–180 `[RM52]`. Form `Nu = 2 + 0.6 Re^0.5 Pr^(1/3)`, `Sh` analogue with Sc.
**(b) Paper's own test.** Water + aniline/naphthalene drops in forced convection; they plot
`(Nu−2)/Pr^(1/3)` (equivalently `(Nu−2)` for gases) vs `Re^0.5` and fit slope **0.6**; anchor
`Nu(Re=0)=2`. Validity 0<Re<200. The **non-tautological oracle** is RM52's tabulated
experimental drop data (naphthalene sublimation, water evaporation), NOT a re-evaluation of the
formula. To weaponize: transcribe ≥3 RM52 data points into the unit test.
**(c) Reproducible in IGLOO.**
- *Unit* — `heat(Re,Pr,·,5)`: Re={0,10,100}; assert Nu=2 at Re=0 (analytic, source-free), and vs
  RM52 data at Re>0. Machine-precision on the Nu=2 anchor.
- *ODE* — generalize B2 ([INFO.md](INFO.md)): full `Nu=2+0.6 Re^0.5 Pr^(1/3)` coupled to Stokes
  drag, oracle = independent RK4 of the SAME coupled system (already PASS 5.2e-11). **Zero
  approximation error** — dominates the frozen-slip route.
**(d) ⚠ CODE-VS-PAPER.** 0.6 **CONFIRMED**; Pr^(1/3) exact **CONFIRMED**. Nu→2 floor present. Correct.

### Model 6 — Kavanau-Drake (heatSelect=6)
**(a) Paper.** Kavanau 1955 continuum `Nu_0=2+0.459 Re^0.55 Pr^(1/3)` + Drake temperature-jump
rarefaction `Nu=Nu_0/(1+3.42·Ma/(RePr)·Nu_0)` (plasma-spray / rarefied-particle form).
**(b) Paper's own test.** Kavanau: sphere heat-transfer data in a subsonic rarefied wind tunnel,
Nu vs Re at several Kn; the rarefaction collapses Nu below the continuum curve as Kn↑. Anchor:
Kn→0 recovers Nu_0; Nu_0(Re→0)=2.
**(c) Reproducible in IGLOO.** TWO regimes needed to cover both terms:
- *Continuum base* — large `d` / high `rho_g` (Kn→0) ⇒ rarefaction term ≈0 ⇒ test `Nu_0=2+0.459 Re^0.55 Pr^0.33`.
- *Rarefaction* — small `d` / low `rho_g` (finite Ma/Re) ⇒ exercise the 3.42 factor; oracle from Drake's bridging value at chosen (Re,Ma,Pr).
- A single test leaves 3.42 unverified.
**(d) ⚠ CODE-VS-PAPER.** Base `0.459/0.55` INFERRED (Kavanau) — not primary-verified this
session, mark `no` in REFERENCES. `3.42` INFERRED (Drake). **Pr^0.33 (not 1/3)** — see D-HEAT-2.
Explicit not implicit — see S1. `toll=1e-20` guards Re=0 (good) but creates the slip=0 discontinuity above.

### Model 3 — JAXA3 (heatSelect=3)
= Kavanau continuum limit (S2). `Nu=2+0.459 Re^0.55 Pr^(1/3)`. Nu→2 at Re=0 ✓.
**Test:** unit-check equal to KD's continuum base at Kn→0 (source-free, S3) + Nu=2 anchor.
**(d)** 0.459/0.55 same provenance/confidence as KD; Pr^(1/3) here vs Pr^0.33 in KD — internal inconsistency, D-HEAT-2.

### Model 2 — JAXA2 (heatSelect=2)
`Nu=2+0.37 Re^0.6 Pr^(1/3)`. Nu→2 at Re=0 ✓ (floor OK). Provenance UNIDENTIFIED (not Whitaker,
not RM). **Test degrades to REGRESSION only** (no independent oracle) + the Nu=2 anchor + RK4-coupled
ODE consistency. State explicitly it is not literature-verified.

### Model 4 — JAXA4 (heatSelect=4)
Bridging ≡ KD (S1), base `2+0.645 Re^0.5 Pr^(1/3)`. Continuum-corner Nu→2 only at Kn→0 (caveat above).
Base 0.645 UNVERIFIED. **Test:** same two-regime plan as KD; the bridging algebra (S1) is the
source-free part; the 0.645 base is regression-only.

### Model 1 — JAXA1 (heatSelect=1)
`Nu=2.5 Re^0.15 + 0.04 Re`. **No +2 conduction floor ⇒ Nu(Re→0)=0.** Ma-free, so no rarefaction
excuse. Genuine physical anomaly for a general closure — D-HEAT-1. Provenance UNIDENTIFIED
(likely a high-Re/rarefied reentry-demise fit where the conduction floor is irrelevant).
**Test:** regression only + explicit low-Re WARNING that this model is invalid as Re→0.

---

## Discrepancies (severity-ranked)

### D-HEAT-1 — JAXA1 lacks the Nu=2 conduction floor (CONFIRMED, correctness-breaking at low Re)
[Lib_Heat.f90#L76](../../../src/lib/Lib_Heat.f90#L76) `Nu = 2.5 Re^0.15 + 0.04 Re` → **Nu→0 as Re→0**.
Physical floor is Nu=2 (isothermal sphere conduction). Any low-slip particle under JAXA1 gets
**zero convective heat** — unphysical. Unlike JAXA4/KD this is NOT rarefaction (no Ma term).
*Do not fix* (may be an intentional high-Re-only fit) — **flag + guard the test domain**. If this
model is ever used at low Re it under-heats badly. Owner must confirm the intended Re range /
whether a `max(·,2)` floor was dropped.

### D-HEAT-2 — Pr exponent inconsistency 0.33 vs 1/3 (CONFIRMED, accuracy-cosmetic)
KD [#L133](../../../src/lib/Lib_Heat.f90#L133) uses `Pr^0.33`; JAXA3 [#L98](../../../src/lib/Lib_Heat.f90#L98)
uses `Pr^(1/3)` for the *same* Kavanau base. `0.33≠0.33333…`. For Pr≈0.71: `0.71^0.33=0.8893`
vs `0.71^0.3333=0.8890` — ~0.03% on the convective term, **cosmetic**. Fix (owner): make KD use
`Pr**(1._R8/3._R8)` for consistency with JAXA3 and the literature `1/3`.

### D-HEAT-3 — "Kavanau-Drake" is explicit, brief/paper may intend implicit (INFERRED, needs owner)
S1: code puts the **continuum** Nu_0 in the bridging denominator (explicit). If the source form
is truly implicit (`Nu=Nu_0/(1+k·Nu)`, total Nu), values differ at finite Kn. Flag, do not fix
blind — resolve against Kavanau/Drake original once available.

### D-HEAT-4 — JAXA4 base coefficient 0.645 unverified (INFERRED)
`0.645 Re^0.5` base has no confirmed provenance (not RM 0.6, not Frössling 0.552, not Kavanau
0.459/0.55). Regression-only until the JAXA report is located. Do not cite.

---

## ⚠ KEY METHOD FINDING — the brief's "very heavy particle freezes slip" route DOES NOT WORK

The proposed route (heavy/velocity-frozen particle ⇒ quasi-constant Re ⇒ exponential T) rests on
`tau_T << tau_v`. But both relaxation times scale identically:

    tau_v = rho_p d^2/(18 mu_g)   [Stokes];   tau_T = c_p,p rho_p d^2/(6 Nu k_g)
    tau_T/tau_v = 3 (c_p,p/c_p,g) Pr_g / Nu        ← INDEPENDENT of rho_p and d

Making the particle heavy or large scales BOTH times by the same `rho_p d^2` — the ratio is fixed
by material/gas properties (`~O(1)`: e.g. water/air c_p ratio ~4, Pr~0.7, Nu~2 ⇒ tau_T/tau_v~4).
**You cannot freeze slip by adding mass.** Slip decays on the same (or faster) timescale as
temperature ⇒ Re is NOT quasi-constant ⇒ the exponential-T oracle carries an *uncontrolled*
approximation error, not a small foldable one.

**Fixes (present, owner picks):**
1. **RK4-coupled oracle (generalize B2) — RECOMMENDED.** Integrate the TRUE coupled (v,T) system
   with an independent RK4 (production shares no code with the oracle), swapping in each Nu
   function. **Zero approximation error** — the frozen-Re error term vanishes entirely. Already
   proven for RM (B2 PASS 5.2e-11). This is strictly better than any frozen-slip design.
2. **Body-force terminal-velocity constant-Re — analytic alternative.** Apply a body force
   (existing body-force e2e pattern); inject AT terminal velocity `v_t=tau_v g` so slip = `v_t` is
   **exactly constant from t=0** ⇒ Re constant ⇒ Nu constant ⇒ **exact** `T_p(t)=T_g+(T_p0−T_g)e^{−t/tau_T}`,
   `tau_T=c_p,p rho_p d^2/(6 Nu(Re_t) k_g)`. Closed-form, no frozen-slip approximation. `Nu(Re_t)`
   for the oracle taken from the PAPER (unit-table value), not the code.

The frozen-slip closed form is only a *degenerate* special case and should not be the primary route.

---

## Shortlist — which Nu model & test to implement first

1. **Unit Nu(Re,Pr[,Ma]) table on `heat()` — RANK 1.** Cheapest, deterministic, sharpest
   coefficient-bug detector (mirrors evap Rank-1). Per model: Nu=2 continuum anchor (source-free,
   analytic) for models 2/3/5 and 4/6-at-Kn→0; Ranz-Marshall Re>0 point vs **RM52 tabulated data**
   (the one truly independent oracle); JAXA1/2/4 = regression + explicit "not literature-verified"
   tag. Immediately exposes D-HEAT-1 (JAXA1→0) and D-HEAT-2 (Pr 0.33 vs 1/3), and confirms S1/S2/S3
   inter-model identities with no citation risk.
2. **Generalize B2 RK4-coupled-oracle to all six Nu — RANK 2.** Zero-approximation ODE test of the
   full momentum→heat coupling; swap the Nu function, reuse the proven RK4 harness. Covers the
   convective term the brief targets, correctly (no frozen-slip error). Start with Ranz-Marshall
   (the trusted 0.6 form), then JAXA2/3.
3. **Body-force terminal-velocity exact-exponential — RANK 3.** The analytic convective-Nu test;
   use only where a closed form is wanted for a report. Requires the constant-Re (terminal)
   injection, NOT a heavy particle.

**Implement Ranz-Marshall first** (only fully-CONFIRMED coefficient, existing B2 oracle). The five
JAXA/KD models should be gated as regression + internal-consistency (S1–S3) until their JAXA/Kavanau
primary sources are located and their DOIs verified — do not present them as literature-verified.

---
*Residual uncertainties:* Kavanau 0.459/0.55, Drake 3.42, JAXA1/2/4 bases, and all JAXA report
identities are INFERRED/UNVERIFIABLE from open sources this session — flagged, not asserted. Nu=2
floor, the tau_T=6Nu derivation, S1/S2 algebra, and the tau_T/tau_v mass-independence are
independently re-derived and CONFIRMED.
