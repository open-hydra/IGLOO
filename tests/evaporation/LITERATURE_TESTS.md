# LITERATURE_TESTS — evaporation models (paper-grounded verification)

Companion to [INFO.md](INFO.md) (which holds the C0–C4 test schema). This file is
the LITERATURE side: for each production closure it gives the original paper, the
verification test the *paper itself* presents, a reproducible-in-IGLOO design, and a
term-by-term CODE-vs-PAPER cross-check of [../../src/lib/Lib_Evaporation.f90](../../src/lib/Lib_Evaporation.f90).

**Confirmed vs inferred** is stated per claim. Coefficients re-derived here, not
recalled. Assumes the known call-site scramble (bug A4) is fixed;
subroutine-level formulas checked independently of that bug.

---

## Canonical spherical result (the yardstick for all four)

Quasi-steady, spherically-symmetric droplet vaporization, constant properties. The
heat-conduction-controlled mass rate (Turns, *An Introduction to Combustion*, Eq. 3.57;
Law, *Combustion Physics*, ch. 13; Lefebvre, *Atomization and Sprays* `[Lef89]`):

    mdot = 4π (k_g/c_pg) r_s ln(1+B_T)  =  2π (k_g/c_pg) d ln(1+B_T)          (radius r_s = d/2)

Mass-diffusion form (Le=1 identical): `mdot = 4π ρ_g D_v r_s ln(1+B_M) = 2π ρ_g D_v d ln(1+B_M)`.
Convective generalization: `mdot = π d (k_g/c_pg) Nu* ln(1+B_T)`, with the stagnant
film `Nu* = Sh* = 2` recovering the line above. Evaporation constant (d² linear law):

    d²(t) = d0² − K t,   K = −d(d²)/dt = 8 (k_g/(ρ_l c_pg)) ln(1+B_T)          (factor 8)

Derivation check (from mdot, m=ρ_l(π/6)d³): `dm/dt = ρ_l(π/2)d² dd/dt = −mdot` ⇒
`d(d²)/dt = 2d·dd/dt = −8(k_g/(ρ_l c_pg))ln(1+B_T)`. Factor **8** confirmed. This is the
number every d²-law textbook prints (Lefebvre λ; Turns Eq. 3.58). **Hold this: the
code, INFO.md C1, and check_evap.py all use 4, i.e. half — see Discrepancy D-EVAP-1.**

---

## Model 1 — d²-law (evapSelect=1)

**(a) Paper.** Godsave, G.A.E., *Proc. Combust. Inst.* 4 (1953) 818–830 (suspended-droplet
combustion, first statement of the d²-law); Spalding, D.B., *Proc. Combust. Inst.* 4 (1953)
847–864 (transfer-number formulation). No DOI (1953 symposia). Compact modern
statement: `[Lef89]`.

**(b) Verification test the papers present.** THE canonical experiment: a suspended
single droplet (n-heptane, ethanol, water) in a quiescent hot/burning environment;
measured D² falls linearly in time; slope = evaporation/burning constant K.
Reproducible numbers (Godsave/Spalding-era, order of magnitude): fuel droplets K ≈
0.5–1.0 mm²/s in combustion; pure evaporation in hot air K ≈ 10⁻⁷–10⁻⁶ m²/s. The exact
oracle is the **closed form** `d²(t)=d0²−Kt` with `K=8(k_g/(ρ_l c_pg))ln(1+B_T)`.

**(c) Reproducible in IGLOO.**
- **Cleanest = quiescent, constant-B_T.** Re=0 (set inlet v = u_g so slip→0), hold gas
  state uniform. B_T is constant *only if T_p is constant*; the coupled T_p ODE
  (F(7), [Lib_RHS.f90](../../src/lib/Lib_RHS.f90#L516)) drives T_p toward the
  **wet-bulb** T_wb where net heating→0 and B_T freezes. Two ways to get constant B_T:
  (i) inject at T_p0 = T_wb (pre-solve T_wb from `psat`/energy balance) so T_p sits still
  from t=0 ⇒ exact linear d²(t); or (ii) accept the short heat-up transient and fit
  the linear tail. Map t = x/u_g (coast at u_g) ⇒ closed-form d²(x)=d0²−K·x/u_g.
  Oracle: slope match to K and R²>0.999. **NOTE:** this exposes D-EVAP-1 — with the
  current code the fitted slope is HALF the literature K.
- **Measured-T_p integral oracle** — already prototyped at
  [d2law/check_evap.py](d2law/check_evap.py): integrate the rate
  kernel along the measured T_p(x). Independent of the T_p equation. **But its
  constant is wrong (see D-EVAP-2).**

**(d) ⚠ CODE-VS-PAPER.** See **Discrepancy D-EVAP-1 (factor-2, CONFIRMED, NEW)** below.

---

## Model 2 — CEM, Classical Evaporation Model (evapSelect=2)

**(a) Paper.** Spalding (1953) transfer-number theory + Ranz & Marshall convective
correlation, `[RM52]` (*Chem. Eng. Prog.* 48, Parts I&II). Code line
[Lib_Evaporation.f90#L144](../../src/lib/Lib_Evaporation.f90#L144):
`Sh = 2 + 0.6 Re^0.5 Sc^(1/3)`, `mdot = −π d ρ_g D_v Sh ln(1+B_M)`, `D_v=k_g/(ρ_g c_pg Le)`.

**(b) Verification test the paper presents.** Ranz–Marshall Parts I & II: water and
aniline drops, forced convection; correlation `Nu = 2 + 0.6 Re^0.5 Pr^(1/3)`,
`Sh = 2 + 0.6 Re^0.5 Sc^(1/3)`, validated 0 < Re < 200. Verification = reproduce Sh(Re)
and the resulting evaporation rate at tabulated (Re, B_M).

**(c) Reproducible in IGLOO.** Re-sweep unit test on `CEM_model`: feed Re=0,10,100,
check Sh=2, 2+0.6·√10·Sc^⅓, … to machine precision. **Cross-model consistency:** at
Re=0, Le=1, CEM must equal d²-law (both → stagnant film). It does NOT in the current
code (D-EVAP-1): CEM(Re=0) gives `mdot=−2π d (k_g/c_pg) ln(1+B_M)`, exactly **twice**
d²-law. This inconsistency is itself the proof of the d²-law bug — no external source
needed.

**(d) ⚠ CODE-VS-PAPER.** Ranz–Marshall coefficient **0.6 CONFIRMED correct**. Prefactor
`π d ρ_g D_v Sh` with Sh→2 gives `2π d ρ_g D_v = 4π r_s ρ_g D_v` = canonical. **CEM is
correct.**

---

## Model 3 — CEM-B, film-corrected CEM (evapSelect=3)

**(a) Paper.** Hubbard, Denny & Mills, *Int. J. Heat Mass Transfer* 18 (1975) 1003–1008
(the "1/3 rule": reference properties at `T_f = T_s + (T∞−T_s)/3`,
`Y_f = Y_s + (Y∞−Y_s)/3`). DOI:10.1016/0017-9310(75)90215-9 *(inferred from title/journal —
verify)*.

**(b) Verification test.** Hubbard-Denny-Mills compare the 1/3-rule reference-property
scheme against full variable-property numerical solutions; agreement in Nu/Sh and rate
across a T∞ range. Verification = property-freeze vs 1/3-rule rate at a chosen (T_p,T_g).

**(c) Reproducible in IGLOO.** Unit test on `CEMB_model`: check `T_f=T_p+(T_g−T_p)/3`,
`ρ_f=ρ_g T_g/T_f`, `μ_f=μ_g(T_f/T_g)^0.7`. At Re=0, Le=1 CEM-B → stagnant film ⇒ must
match CEM(Re=0) and (once fixed) d²-law.

**(d) ⚠ CODE-VS-PAPER.** 1/3-rule `T_f` CONFIRMED. Film density `ρ_f=ρ_g T_g/T_f`
(ideal gas, p const) CONFIRMED. Power-law `μ,k ∝ T^0.7` is a *reasonable but not unique*
1/3-rule companion (INFERRED convention — 0.7 is a common air fit, not from HDM itself;
acceptable). `Re_f = Re·(ρ_f/ρ_g)·(μ_g/μ_f)`
([#L170](../../src/lib/Lib_Evaporation.f90#L170)) is the correct film-Reynolds rescale.
Prefactor `π d ρ_f D_f Sh` → **canonical (2π d at Sh=2). CEM-B correct.**

---

## Model 4 — ASM, Abramzon–Sirignano (evapSelect=4)

**(a) Paper.** Abramzon, B. & Sirignano, W.A., "Droplet vaporization model for spray
combustion calculations," *Int. J. Heat Mass Transfer* **32**(9) (1989) 1605–1618.
DOI:10.1016/0017-9310(89)90064-7 (suggest tag `[AS89]` for REFERENCES.md).

**(b) Verification tests the paper presents.** A–S give: the modified `Sh*`, `Nu*` with
the film-thickness correction `F(B)`; and figure comparisons of `Nu*/Nu0`, `Sh*/Sh0`,
droplet lifetime vs the classical model and vs the "extended film" reference. Reference
limit checks: `B→0 ⇒ F(B)→1 ⇒ Sh*→Sh0, Nu*→Nu0`; `Re→0 ⇒ Sh0,Nu0→2`.

**(c) Reproducible in IGLOO.**
- Unit test on `ASM_model`: `F(0)=1` (limit), `F(B)=(1+B)^0.7 ln(1+B)/B`; verify the
  `φ–B_T–F_T–Nu*` fixed-point converges; `Sh0=2+0.552 Re^0.5 Sc^⅓`.
- **Cross-model limit:** Re=0 ⇒ Sh0=Nu0=2, and with B_M small `F→1` ⇒ ASM → stagnant
  film `mdot=−2π d ρ_g D_v ln(1+B_M)` = CEM(Re=0) = (fixed) d²-law. Good consistency
  gate spanning all four.

**(d) ⚠ CODE-VS-PAPER (term by term).**
| A–S eq | quantity | code | verdict |
|---|---|---|---|
| 17 | `F(B)=(1+B)^0.7 ln(1+B)/B` | [#L251](../../src/lib/Lib_Evaporation.f90#L251) | **CONFIRMED**, exponent 0.7 correct |
| 19/Frössling | `Sh0=2+0.552 Re^0.5 Sc^⅓` | [#L200](../../src/lib/Lib_Evaporation.f90#L200) | **CONFIRMED** 0.552 |
| 19/Frössling | `Nu0=2+0.552 Re^0.5 Pr^⅓` | [#L214](../../src/lib/Lib_Evaporation.f90#L214) | **CONFIRMED** |
| 18 | `Sh*=2+(Sh0−2)/F_M` | [#L201](../../src/lib/Lib_Evaporation.f90#L201) | **CONFIRMED** |
| 8 | `mdot=−π d ρ_g D_v Sh* ln(1+B_M)` | [#L204](../../src/lib/Lib_Evaporation.f90#L204) | **CONFIRMED** (=canonical 2π d at Sh*=2) |
| 22 | `φ=(c_pv/c_pg)(Sh*/Nu*)(1/Le)` | [#L220](../../src/lib/Lib_Evaporation.f90#L220) | **CONFIRMED** |
| 23 | `B_T=(1+B_M)^φ−1` | [#L221](../../src/lib/Lib_Evaporation.f90#L221) | **CONFIRMED** |
| 24 | heat into liquid | [#L233](../../src/lib/Lib_Evaporation.f90#L233) | ⚠ **see D-EVAP-3 (Lv sign, INFERRED)** |

ASM mass path is faithful to A–S. **The ASM mdot is correct (unlike d²-law).**

---

## Discrepancies (severity-ranked, NEW — beyond the known call-site bug)

### D-EVAP-1 — d²-law mdot is HALF the Godsave/Spalding rate (CONFIRMED, correctness-breaking)
**Code** [Lib_Evaporation.f90#L127](../../src/lib/Lib_Evaporation.f90#L127):
`mdot = −π d (k_g/c_pg) ln(1+B_T)`.
**Canonical** (Turns Eq. 3.57; Lefebvre λ): `mdot = 2π d (k_g/c_pg) ln(1+B_T)`
(= `4π r_s (k_g/c_pg) ln(1+B_T)`). The code uses `π d`, i.e. effectively **Nu=1**; the
stagnant film is **Nu=2**. Missing factor 2.
**Two independent proofs:** (i) external — canonical spherical solution, factor 8 in K
(above); (ii) INTERNAL, source-free — CEM/CEM-B/ASM all carry the explicit `Sh`/`Nu`
multiplier and give `2π d …` at Sh=2, so `d²-law = ½·CEM(Re=0,Le=1)` in the very limit
where they must be identical. A cross-model Re=0 consistency test fails by exactly 2.
**Fix (owner):** `mdot = −2π d (k_g/c_pg) ln(1+B_T)` — insert the Nu=2 factor (equivalently
`π d k_g/c_pg · Nu` with Nu=2). Then K→8 and d²-law = CEM(Re=0).
**Impact:** every d²-law run under-evaporates 2×; droplet lifetime 2× too long.

### D-EVAP-2 — the verification oracle encodes the buggy half-constant (CONFIRMED, blocks detection)
Both the planned test and the prototype bake in K=4, not the literature K=8, so they
would **rubber-stamp** D-EVAP-1 (and *fail* a correct fix):
- [INFO.md](INFO.md) row C1: `K = 4 k_g ln(1+B_T)/(ρ_l c_pg)` — half; should be `8`.
- [check_evap.py](d2law/check_evap.py#L73): `APRE = 4.0*KG/(CP_G*RHO_P)` and the
  docstring "d(d²)/dt = −4 kg/(cp_g rho_p) ln(1+BT) … the classic d2-law" — the classic
  constant is **8**. The derivation is arithmetically consistent *with the coded (half)
  mdot*, i.e. it is a code tautology, not a literature oracle.
**Fix:** set the oracle constant to `8 k_g/(ρ_l c_pg)` (APRE=8·KG/…) so it tests the
LITERATURE. With correct-8 oracle vs current-half code, the e2e case fails at ~2× —
which is the desired signal. Do this **together with** D-EVAP-1 (fixing one without the
other flips the pass/fail).

### D-EVAP-3 — ASM heat-to-liquid Lv sign (INFERRED, needs owner convention check — NOT a confirmed bug)
**Code** [#L233](../../src/lib/Lib_Evaporation.f90#L233):
`Qdot_evap = mdot (L_v + c_pv (T_g−T_p)/B_T)`.
**A–S eq 24** (heat conducted into the liquid interior, their mdot>0):
`Q_L = mdot ( c_pv (T∞−T_s)/B_T − L_v )` — the L_v term is **subtracted**, and A–S mdot
is positive whereas IGLOO's `mdot<0`. The relative sign of L_v differs. This *may* be
absorbed correctly by IGLOO's F(7) assembly
([Lib_RHS.f90#L516](../../src/lib/Lib_RHS.f90#L516),
`F(7)=(Qdot − mdot(Z(7) − L_v·cpFactor))/m`) and the `override_Qdot` path — the sign
convention there is non-obvious. **Flag, do not fix blind:** trace the ASM Qdot end-to-end
against A–S eq 24 once ASM is exercised (currently moot — call-site bug + ASM untested).
Wet-bulb sanity: as T_p→T_wb, `Q_L→0` in A–S; verify the code's Qdot→0 there.

---

## Best evaporation verification test to stand up first (call-site bug fixed)

1. **Cross-model quiescent consistency (Re=0) unit test — RANK 1.** Call
   `evaporation` with Re=0, Le=1, identical gas/droplet state, for evapSelect 1/2/3/4;
   assert `mdot` equal across models to machine precision. Cheapest (no integration, no
   mesh), deterministic, and it **immediately exposes D-EVAP-1** (d²-law = ½ of the other
   three). Highest signal-to-effort. Source-free — the oracle is the code's own
   consistency, so no citation dispute.
2. **Analytic d²(t) linear-decay, quiescent constant-B_T — RANK 2.** Inject at T_wb (or
   fit linear tail), map t=x/u_g, fit slope to `K=8(k_g/(ρ_l c_pg))ln(1+B_T)`, R²>0.999,
   lifetime `d0²/K`. This is THE canonical d²-law test. **Requires fixing the oracle
   constant to 8 (D-EVAP-2) first**, else it certifies the bug.
3. **Measured-T_p integral oracle (existing check_evap.py) — RANK 3.** Already built;
   only needs `APRE 4→8` (D-EVAP-2) and the call-site fix, then rename to `check.py` for
   the gate. Verifies the rate law along the real temperature history without modelling
   T_p.

Consistency-first (Rank 1) is recommended: it needs no chosen constant, no wet-bulb
setup, and it is the sharpest detector of the factor-2 error that all three other
artifacts currently hide.
