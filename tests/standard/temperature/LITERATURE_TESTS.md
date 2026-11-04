# LITERATURE_TESTS — heat-transfer (Nusselt) models (paper-grounded verification)

Companion to [INFO.md](INFO.md) (B1/BO/B2/B3 schema). LITERATURE side of the convective
heat closure: per production Nu-correlation — original paper, the paper's own verification
test, a reproducible-in-IGLOO design, and a term-by-term CODE-vs-PAPER cross-check of
[../../../src/lib/Lib_Heat.f90](../../../src/lib/Lib_Heat.f90) (selector `heat()`).

**Confirmed vs inferred** stated per claim. Coefficients re-derived / re-verified here, not
recalled.

Sources. `[Shim06]` (JAXA-SP-05-035E) eqs. 40–50 hold this whole catalog; the "JAXA1–4" names resolve to
NASA-SP-8039 (1971), its ref [12]. JAXA2 (eq. 46), JAXA3 (eq. 47), Ranz-Marshall (eq. 48) and
Kavanau-Drake (eq. 49) match it exactly — eq. 49 is EXPLICIT and prints `Pr^0.33`. Eq. 45 (JAXA1) has no
`+2` floor — a documented limitation of the published correlation, not a transcription error
(`test_heat_probes` pins it). JAXA4's base is 0.654 (eq. 50). Kavanau-Drake primary: Univ. of California
Rept. HE-150-108 (1953), ref [14].

---

## 0. How Nu enters the ODE (shared oracle skeleton) — VERIFIED

[Lib_Equations.f90](../../../src/lib/Lib_Equations.f90), `interphase`:

    Qdot = Nu * k_g * pi * d * (T_g - T_p) * cpFactor          (cpFactor = 1/c_p,p)
    F(7) = Qdot / m ,   m = rho_p (pi/6) d^3

⇒ **dT_p/dt = (6 Nu k_g)/(c_p,p rho_p d^2)·(T_g − T_p) = (T_g − T_p)/tau_T**, with

    tau_T = c_p,p rho_p d^2 / (6 Nu k_g)          [GENERAL — re-derived]

Re-derivation: `Qdot/m = Nu k_g pi d (T_g−T_p)/c_p ÷ rho_p(pi/6)d^3 = 6 Nu k_g(T_g−T_p)/(c_p rho_p d^2)`.
This is standard lumped-capacitance: `m c_p dT/dt = h A (T_g−T_p)`, `h=Nu k_g/d`, `A=pi d^2`.
At Nu=2 ⇒ `tau_T = c_p rho_p d^2/(12 k_g)`, matching e2e temp-relax
([check.py](../temp-relax/check.py) `TAU=CP*RHO_P*D**2/(12 K_G)`) and INFO.md B1. **CONFIRMED.**
`Pr` inside the closure (`interphase`):
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
  `toll=1e-20` guard. The unit case B1 (zero slip, `test_temperature`) sits exactly on this point,
  where every model but JAXA1 gives Nu=2.
- **JAXA1 →0 at Re→0** (it is Ma-free, so no rarefaction alibi): the published form's, see Discrepancies.
- **Monotonic in Re** at fixed Pr; **Nu>0** for all admissible inputs.

---

## Provenance map (what the labels REALLY are)

| selector | code Nu | true origin (confidence) |
|---|---|---|
| JAXA1 (1) | `2.5 Re^0.15 + 0.04 Re` | **`[Shim06]` eq. 45**, quoting `[SP8039]`; high-speed/rarefied fit, no conduction floor as published; pinned by `test_heat_probes` XH1 |
| JAXA2 (2) | `2 + 0.37 Re^0.6 Pr^(1/3)` | **`[SP8039]` Table II, ref. 82**, printed `2 + 0.370 Re^0.6 Pr^⅓` = the code (`[Shim06]` eq. 46); checked by `test_heat_nu` HN4 and end to end by `conv-nu-jaxa2` |
| JAXA3 (3) | `2 + 0.459 Re^0.55 Pr^(1/3)` | **`[SP8039]` Table II, ref. 42**, printed `2 + 0.459 Re^0.55 Pr^⅓` = the code (`[Shim06]` eq. 47) — the Kavanau continuum base of model 6, rarefaction off; checked by HN5 and `conv-nu-jaxa3` |
| JAXA4 (4) | `[(2+0.654 Re^0.5 Pr^0.5)^-1 + 3.42 Ma/(RePr)]^-1` | **`[SP8039]` Table II, the last row** (p. 26), printed `[(2 + 0.654 Re^½ Pr^½)^−1 + 3.42 M/(Re Pr)]^−1` = the code; the Kavanau-Drake bridging of model 6 on a `0.654 Re^½ Pr^½` base. **Two prints:** `[Shim06]` eq. 50, which cites SP-8039, prints the same law with **Pr^⅓**, the exponent of refs. 82/42 (whose ⅓ is a glyph distinct from the ½ of the last row on the same page); the code takes SP-8039's Pr^½ (with Pr^⅓: +6.7e-3 of Nu at `conv-nu`'s state). Checked by HN6, the rarefaction identity HN8 and `conv-nu-jaxa4` |
| Ranz-Marshall (5) | `2 + 0.6 Re^0.5 Pr^(1/3)` | **Ranz & Marshall 1952 `[RM52]`** — 0.6 CONFIRMED (textbook + evap review); valid 0<Re<200 |
| Kavanau-Drake (6) | `(2+0.459 Re^0.55 Pr^0.33)/(1+3.42 Ma/(RePr)·Nu_0)` | **`[SP8039]` Table II, ref. 72** = `[Shim06]` eq. 49, both printing `Pr^0.33` = the code: Kavanau's continuum base, Drake's rarefaction term, explicit; checked by HN7 (the 0.33 pinned), HN8 and `conv-nu-kd` |

Kavanau's paper: Kavanau, L.L., "Heat transfer from spheres to a rarefied gas in subsonic flow,"
*Trans. ASME* 77 (1955) 617–623 (no DOI, historical). `[RM52]` is listed in `tests/REFERENCES.md`.

---

## Source-free spine (strongest findings — need NO external citation)

**S1 — JAXA4 ≡ Kavanau-Drake bridging (algebraically identical).** Re-derived:
JAXA4 `Nu = [Nu_0^-1 + k]^-1 = Nu_0/(1+k·Nu_0)` with `k=3.42 Ma/(RePr)`; Kavanau-Drake
[#L146](../../../src/lib/Lib_Heat.f90#L146) `Nu = Nu_0/(1+k·Nu_0)` — **same functional form**,
differing ONLY in the continuum base (`0.654 Re^0.5 Pr^0.5` vs `0.459 Re^0.55 Pr^0.33`).
Both are **explicit** (single eval, continuum `Nu_0` in the denominator), as `[Shim06]` eq. 49 prints
Kavanau-Drake.

**S2 — JAXA3 ≡ Kavanau-Drake continuum base** (`2+0.459 Re^0.55·Pr^{1/3 or 0.33}`), rarefaction
off, Pr exponent `1/3` vs `0.33` the only difference, each as its source prints it (see Discrepancies).

**S3 — Cross-model continuum-corner consistency (Re→0, large d ⇒ Kn→0):** JAXA2, JAXA3,
Ranz-Marshall, and (in the Kn→0 corner) JAXA4 & KD must ALL give Nu=2. JAXA1 gives 0 (its published form).
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
formula; HN1 and HN2 (`test_heat_nu`) re-evaluate the formula and its floor, no RM52 data point.
**(c) Reproducible in IGLOO.**
- *Unit* — `heat(Re,Pr,·,5)`: Re={0,10,100}; assert Nu=2 at Re=0 (analytic, source-free), and vs
  RM52 data at Re>0. Machine-precision on the Nu=2 anchor.
- *ODE* — B2 ([INFO.md](INFO.md)): full `Nu=2+0.6 Re^0.5 Pr^(1/3)` coupled to Stokes
  drag, oracle = independent RK4 of the SAME coupled system (5.2e-11). **Zero
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
**(d) ⚠ CODE-VS-PAPER.** Base `0.459/0.55`, `3.42` and **Pr^0.33 (not 1/3)** as `[SP8039]` Table II ref. 72 and
`[Shim06]` eq. 49 print them (HN7); Kavanau's and Drake's primaries are not checked.
Explicit, as printed — see S1. `toll=1e-20` guards Re=0 (good) but creates the slip=0 discontinuity above.

### Model 3 — JAXA3 (heatSelect=3)
= Kavanau continuum limit (S2). `Nu=2+0.459 Re^0.55 Pr^(1/3)`. Nu→2 at Re=0 ✓.
**Test:** unit-check equal to KD's continuum base at Kn→0 (source-free, S3) + Nu=2 anchor.
**(d)** 0.459/0.55 and Pr^(1/3) as `[SP8039]` Table II ref. 42 prints them (HN5); KD's ref. 72 prints Pr^0.33.

### Model 2 — JAXA2 (heatSelect=2)
`Nu=2+0.37 Re^0.6 Pr^(1/3)`. Nu→2 at Re=0 ✓ (floor OK). `[SP8039]` Table II ref. 82 (`[Shim06]` eq. 46; not
Whitaker, not RM). **Test:** HN4 against the print, HN2 the Nu=2 anchor, `conv-nu-jaxa2` end to end.

### Model 4 — JAXA4 (heatSelect=4)
Bridging ≡ KD (S1), base `2+0.654 Re^0.5 Pr^0.5`. Continuum-corner Nu→2 only at Kn→0 (caveat above).
The base 0.654 is the one `[SP8039]` Table II's last row and `[Shim06]` eq. 50 both print; the exponent Pr^½ is
SP-8039's (Shimada prints Pr^⅓). **Test:** HN6 pins the law to SP-8039's print on Re × Pr × Ma, HN8 the
bridging algebra (S1, source-free), `conv-nu-jaxa4` the law end to end.

### Model 1 — JAXA1 (heatSelect=1)
`Nu=2.5 Re^0.15 + 0.04 Re`. **No +2 conduction floor ⇒ Nu(Re→0)=0.** Ma-free, so no rarefaction
excuse. A physical anomaly for a general closure, and the published form's (`[Shim06]` eq. 45, quoting
`[SP8039]`; likely a high-Re/rarefied reentry-demise fit where the conduction floor is irrelevant).
**Test:** XH1 (`test_heat_probes`) pins the form, Nu(Re→0)→0 included; the model is invalid as Re→0.

---

## Discrepancies (severity-ranked)

### JAXA1 lacks the Nu=2 conduction floor (as published; under-heats at low Re)
`heat_JAXA_1` ([Lib_Heat.f90](../../../src/lib/Lib_Heat.f90)) `Nu = 2.5 Re^0.15 + 0.04 Re` → **Nu→0 as Re→0**.
Physical floor is Nu=2 (isothermal sphere conduction). Any low-slip particle under JAXA1 gets
**zero convective heat** — unphysical. Unlike JAXA4/KD this is NOT rarefaction (no Ma term).
`[Shim06]` eq. 45 prints the form without the floor and the code transcribes it (XH1): used at low Re,
the model under-heats badly.

### The Kavanau base's Pr exponent: 0.33 (KD) vs 1/3 (JAXA3), each as printed (accuracy-cosmetic)
KD (`heat_Kavanau_Drake`) uses `Pr^0.33`; JAXA3 (`heat_JAXA_3`) uses `Pr^(1/3)` for the *same* Kavanau base.
`0.33≠0.33333…`. For Pr≈0.71: `0.71^0.33=0.8931` vs `0.71^0.3333=0.8921` — ~0.11% on the convective term,
**cosmetic**. Each is its source's print (`[SP8039]` Table II refs. 72 and 42; `[Shim06]` eqs. 49 and 47),
pinned by HN7 and HN5.

### "Kavanau-Drake" is explicit, as printed
S1: the code puts the **continuum** Nu_0 in the bridging denominator (explicit), as `[SP8039]` ref. 72 and
`[Shim06]` eq. 49 print it. An implicit form (`Nu=Nu_0/(1+k·Nu)`, total Nu) would differ at finite Kn;
Kavanau's and Drake's originals are not checked.

### JAXA4's Prandtl exponent: the two sources print different ones
`[SP8039]` Table II's last row prints `[(2 + 0.654 Re^½ Pr^½)^−1 + 3.42 M/(Re Pr)]^−1`; `[Shim06]` eq. 50, which
cites SP-8039, prints the same base 0.654 with Pr^⅓. The code takes SP-8039's Pr^½, which gives a Nu up to 12 %
below Pr^⅓'s for gases (Pr 0.42–0.8, Re ≤ 1000, Ma = 0); HN6 and `conv-nu-jaxa4` pin it to that print.

---

## ⚠ Method: the "very heavy particle freezes slip" route DOES NOT WORK

A heavy/velocity-frozen particle route (⇒ quasi-constant Re ⇒ exponential T) rests on
`tau_T << tau_v`. But both relaxation times scale identically:

    tau_v = rho_p d^2/(18 mu_g)   [Stokes];   tau_T = c_p,p rho_p d^2/(6 Nu k_g)
    tau_T/tau_v = 3 (c_p,p/c_p,g) Pr_g / Nu        ← INDEPENDENT of rho_p and d

Making the particle heavy or large scales BOTH times by the same `rho_p d^2` — the ratio is fixed
by material/gas properties (`~O(1)`: e.g. water/air c_p ratio ~4, Pr~0.7, Nu~2 ⇒ tau_T/tau_v~4).
**You cannot freeze slip by adding mass.** Slip decays on the same (or faster) timescale as
temperature ⇒ Re is NOT quasi-constant ⇒ the exponential-T oracle carries an *uncontrolled*
approximation error, not a small foldable one.

**Two routes avoid it:**
1. **RK4-coupled oracle** (B2, `test_temperature`). Integrate the TRUE coupled (v,T) system
   with an independent RK4 (production shares no code with the oracle), with the Nu
   function swapped in. **Zero approximation error** — the frozen-Re error term vanishes entirely:
   B2 measures 5.2e-11 for RM. This is strictly better than any frozen-slip design.
2. **Body-force terminal-velocity constant-Re — analytic** (`conv-nu`, `conv-nu-jaxa2/-jaxa3/-jaxa4/-kd`).
   Apply a body force; inject AT terminal velocity `v_t=tau_v g` so slip = `v_t` is
   **exactly constant from t=0** ⇒ Re constant ⇒ Nu constant ⇒ **exact** `T_p(t)=T_g+(T_p0−T_g)e^{−t/tau_T}`,
   `tau_T=c_p,p rho_p d^2/(6 Nu(Re_t) k_g)`. Closed-form, no frozen-slip approximation. `Nu(Re_t)`
   for the oracle taken from the PAPER (typed in the check), not the code.

The frozen-slip closed form is only a *degenerate* special case and should not be the primary route.

---

## Checks

- **Unit Nu(Re,Pr[,Ma]) tables on `heat()`** (`test_heat_nu`): HN1 Ranz-Marshall, HN2 the Nu=2 continuum
  anchor (source-free) for models 2/3/5, HN3 KD at Ma=0, HN4–HN7 JAXA2/JAXA3/JAXA4/KD against `[SP8039]`'s
  print, HN8 the rarefaction term of the two bridged forms (S1). JAXA1: XH1 (`test_heat_probes`).
- **RK4-coupled ODE** (`test_temperature` B2, Ranz-Marshall); B1/BO/B3 the Nu=2 relaxation, its order and
  the temperature invariants.
- **Body-force terminal-velocity exact exponential** (e2e): `conv-nu` (Ranz-Marshall) and
  `conv-nu-jaxa2/-jaxa3/-jaxa4/-kd`.

---
*Residual uncertainties:* the primaries behind `[SP8039]`'s rows (Kavanau's 0.459/0.55, Drake's 3.42, the JAXA
reports) are not checked. Nu=2 floor, the tau_T=6Nu derivation, S1/S2 algebra, and the tau_T/tau_v
mass-independence are independently re-derived and CONFIRMED.
