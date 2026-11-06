# LITERATURE_TESTS — evaporation models (paper-grounded verification)

Companion to [INFO.md](INFO.md) (which holds the C0–C4 test schema). This file is
the LITERATURE side: for each production closure it gives the original paper, the
verification test the *paper itself* presents, a reproducible-in-IGLOO design, and a
term-by-term CODE-vs-PAPER cross-check of [../../src/lib/Lib_Evaporation.f90](../../src/lib/Lib_Evaporation.f90).

**Confirmed vs inferred** is stated per claim. Coefficients re-derived here, not
recalled. The formulas are checked at the subroutines; the production call site that
feeds them is pinned by `test_evap_probes` XE2.

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
number every d²-law textbook prints (Lefebvre λ; Turns Eq. 3.58), and the one the code
(`d2law`'s stagnant film), [INFO.md](INFO.md) C1 and [d2law/check.py](d2law/check.py) use.

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
  (F(7), [Lib_RHS.f90](../../src/lib/Lib_RHS.f90#L733)) heats the drop: with the
  Ranz–Marshall heat at Re=0 (`π d k_g Nu (T_g−T_p)`, Nu=2) and this rate the net heating
  is `2π d k_g (L_v/c_pg)(B_T − ln(1+B_T)) > 0`, so T_p climbs toward T_g and B_T falls.
  [d2law-line](d2law-line/INFO.md) holds B_T at its injection value with a particle c_p
  ×100 and density 8000 (ρ_p c_p ≈ 271× d2law's; T_p drifts ~0.8 K over the path) and
  maps t = x/u_g (coast at u_g) ⇒ the
  closed-form d²(x)=d0²−K·x/u_g from the case inputs alone; the drift enters the
  tolerance, not the reference.
- **Measured-T_p integral oracle** — [d2law/check.py](d2law/check.py): integrates the
  rate kernel, with the literature 8, along the measured T_p(x). Independent of the T_p
  equation.

**(d) ⚠ CODE-VS-PAPER.** `d2law` codes `mdot = −2π d (k_g/c_pg) ln(1+B_T)`
([#L340](../../src/lib/Lib_Evaporation.f90#L340)): the Nu=2 stagnant film, K = 8,
**CONFIRMED** (the cross-checks below).

---

## Model 2 — CEM, Classical Evaporation Model (evapSelect=2)

**(a) Paper.** Spalding (1953) transfer-number theory + Ranz & Marshall convective
correlation, `[RM52]` (*Chem. Eng. Prog.* 48, Parts I&II). Code line
[Lib_Evaporation.f90#L353](../../src/lib/Lib_Evaporation.f90#L353):
`Sh = 2 + 0.6 Re^0.5 Sc^(1/3)`, `mdot = −π d ρ_g D_v Sh ln(1+B_M)`, `D_v=k_g/(ρ_g c_pg Le)`.

**(b) Verification test the paper presents.** Ranz–Marshall Parts I & II: water and
aniline drops, forced convection; correlation `Nu = 2 + 0.6 Re^0.5 Pr^(1/3)`,
`Sh = 2 + 0.6 Re^0.5 Sc^(1/3)`, validated 0 < Re < 200. Verification = reproduce Sh(Re)
and the resulting evaporation rate at tabulated (Re, B_M).

**(c) Reproducible in IGLOO.** Re-sweep unit test on `CEM_model`: feed Re=0,10,100,
check Sh=2, 2+0.6·√10·Sc^⅓, … to machine precision. **Cross-model consistency:** at
Re=0, Le=1 and equal transfer numbers CEM equals d²-law (both → stagnant film,
`mdot=−2π d (k_g/c_pg) ln(1+B)`): `test_evaporation` EV1 pins CEM's rate,
`test_evap_probes` XE1 d²-law's.

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
`ρ_f=ρ_g T_g/T_f`, `μ_f=μ_g(T_f/T_g)^0.7`. At Re=0 CEM-B is the stagnant film at the
1/3-rule film properties: ρ_f cancels in ρ_f·D_f, so its rate is CEM(Re=0)'s times
(T_f/T_g)^0.7 (EV3 pins the film chain, EV5 the factor).

**(d) ⚠ CODE-VS-PAPER.** 1/3-rule `T_f` CONFIRMED. Film density `ρ_f=ρ_g T_g/T_f`
(ideal gas, p const) CONFIRMED. Power-law `μ,k ∝ T^0.7` is a *reasonable but not unique*
1/3-rule companion (INFERRED convention — 0.7 is a common air fit, not from HDM itself;
acceptable). `Re_f = Re·(ρ_f/ρ_g)·(μ_g/μ_f)`
([#L376](../../src/lib/Lib_Evaporation.f90#L376)) is the correct film-Reynolds rescale.
Prefactor `π d ρ_f D_f Sh` → **canonical (2π d at Sh=2). CEM-B correct.**

---

## Model 4 — ASM, Abramzon–Sirignano (evapSelect=4)

**(a) Paper.** Abramzon, B. & Sirignano, W.A., "Droplet vaporization model for spray
combustion calculations," *Int. J. Heat Mass Transfer* **32**(9) (1989) 1605–1618.
DOI:10.1016/0017-9310(89)90043-4, tag `[AS89]` ([../REFERENCES.md](../REFERENCES.md)).

**(b) Verification tests the paper presents.** A–S give: the modified `Sh*`, `Nu*` with
the film-thickness correction `F(B)`; and figure comparisons of `Nu*/Nu0`, `Sh*/Sh0`,
droplet lifetime vs the classical model and vs the "extended film" reference. Reference
limit checks: `B→0 ⇒ F(B)→1 ⇒ Sh*→Sh0, Nu*→Nu0`; `Re→0 ⇒ Sh0,Nu0→2`.

**(c) Reproducible in IGLOO.**
- Unit test on `ASM_model`: `F(0)=1` (limit), `F(B)=(1+B)^0.7 ln(1+B)/B`; verify the
  `φ–B_T–F_T–Nu*` fixed-point converges; `Sh0=2+0.552 Re^0.5 Sc^⅓`.
- **Cross-model limit:** Re=0 ⇒ Sh0=Nu0=2, and with B_M small `F→1` ⇒ ASM → stagnant
  film `mdot=−2π d ρ_g D_v ln(1+B_M)` = CEM(Re=0) = d²-law: a consistency check
  spanning d²-law, CEM and ASM (EV6 and EV6b pin ASM = CEM at Re=0 within 9 ε, in the boiling clamp and below it, for B_M ≥ 1e-2); CEM-B
  sits the factor (T_f/T_g)^0.7 below CEM (EV5).

**(d) ⚠ CODE-VS-PAPER (term by term).**
| A–S eq | quantity | code | verdict |
|---|---|---|---|
| 17 | `F(B)=(1+B)^0.7 ln(1+B)/B` | [#L522](../../src/lib/Lib_Evaporation.f90#L522) | **CONFIRMED**, exponent 0.7 correct |
| 20/Frössling | `Sh0=2+0.552 Re^0.5 Sc^⅓` | [#L402](../../src/lib/Lib_Evaporation.f90#L402) | **CONFIRMED** 0.552 |
| 20/Frössling | `Nu0=2+0.552 Re^0.5 Pr^⅓` | [#L415](../../src/lib/Lib_Evaporation.f90#L415) | **CONFIRMED** |
| 10 | `Sh*=2+(Sh0−2)/F_M` | [#L403](../../src/lib/Lib_Evaporation.f90#L403) | **CONFIRMED** |
| 8 | `mdot=−π d ρ_g D_v Sh* ln(1+B_M)` | [#L406](../../src/lib/Lib_Evaporation.f90#L406) | **CONFIRMED** (=canonical 2π d at Sh*=2) |
| 22 | `φ=(c_pv/c_pg)(Sh*/Nu*)(1/Le)` | [#L420](../../src/lib/Lib_Evaporation.f90#L420) | **CONFIRMED** |
| 23 | `B_T=(1+B_M)^φ−1` | [#L422](../../src/lib/Lib_Evaporation.f90#L422) | **CONFIRMED** |
| 24 | heat into liquid | [#L432](../../src/lib/Lib_Evaporation.f90#L432) | **CONFIRMED**: `Q_G=−mdot c_pv (T_g−T_p)/B_T`, the latent term added by the F(7) assembly (below) |

ASM mass and heat paths are faithful to A–S. **The ASM mdot is correct.**

---

## Cross-checks across the models

### The d²-law rate — the Godsave/Spalding stagnant film (CONFIRMED)
**Code** [Lib_Evaporation.f90#L340](../../src/lib/Lib_Evaporation.f90#L340):
`mdot = −2π d (k_g/c_pg) ln(1+B_T)`.
**Canonical** (Turns Eq. 3.57; Lefebvre λ): `mdot = 2π d (k_g/c_pg) ln(1+B_T)`
(= `4π r_s (k_g/c_pg) ln(1+B_T)`): the stagnant film, **Nu=2**; `π d` would be
effectively **Nu=1**, half the rate.
**Two independent proofs:** (i) external — canonical spherical solution, factor 8 in K
(above); (ii) INTERNAL, source-free — CEM/CEM-B/ASM all carry the explicit `Sh`/`Nu`
multiplier and give `2π d …` at Sh=2, so `d²-law = CEM(Re=0,Le=1)` in the very limit
where they must be identical (equal transfer numbers). A half-rate d²-law fails both by
exactly 2: lifetime 2× too long.

### The oracles carry the literature constant 8
- [INFO.md](INFO.md) row C1: `K = 8 k_g ln(1+B_T)/(ρ_l c_pg)`.
- [d2law/check.py](d2law/check.py): `APRE = 8.0 * KG / (CP_G * RHO_P)`, the Nu=2 film: it
  tests the LITERATURE, not the coded rate, and a half-rate `d2law` fails the e2e case at
  ~2×.

### The ASM heat to the liquid — A–S eq. 24 (CONFIRMED)
**Code** [#L432](../../src/lib/Lib_Evaporation.f90#L432):
`Qdot_evap = −mdot c_pv (T_g−T_p)/B_T`, the gas-side heat, positive into the droplet,
latent sink excluded; `override_Qdot` makes it the heat in place of the convective one
([Lib_RHS.f90#L718](../../src/lib/Lib_RHS.f90#L718)), and the F(7) assembly
([Lib_RHS.f90#L733](../../src/lib/Lib_RHS.f90#L733), `F(7)=(Qdot + mdot·L_v·cpFactor)/m`)
adds the latent term.
**A–S eq 24** (heat conducted into the liquid interior, their mdot>0):
`Q_L = mdot ( c_pv (T∞−T_s)/B_T − L_v )`. IGLOO's `mdot<0` is −mdot_AS, so the code's
`Qdot_evap + mdot·L_v` is A–S's `Q_L`. Wet-bulb sanity: as T_p→T_wb, `Q_L→0`. Pinned by
`test_evap_probes` XE3 (the sign and the override), `test_evaporation` EV10/EV11 (the
value near equilibrium; the Stefan-free limit) and `asm-box` K2 (the temperature end to
end).

---

## The tests, strongest oracle first

1. **Cross-model quiescent consistency (Re=0)**, unit, no integration, no mesh:
   `test_evap_probes` XE1 (d²-law against the literature rate), `test_evaporation` EV1
   (CEM against the Spalding chain), EV3 (CEM-B's film chain), EV6 and EV6b (ASM equals CEM at
   Re=0 within 9 ε, the rounding of the two sequences, in the boiling clamp and below it, for
   B_M ≥ 1e-2, where both take log(1 + B_M)). Deterministic; the oracle is the literature or the
   code's own consistency, so no citation dispute.
2. **Analytic d²(t) linear decay, constant B_T** — [d2law-line](d2law-line/INFO.md): the
   slope `K=8(k_g/(ρ_l c_pg))ln(1+B_T)` from the case inputs alone, t=x/u_g. THE canonical
   d²-law test.
3. **Measured-T_p integral oracle** — [d2law/check.py](d2law/check.py): the rate law along
   the real temperature history without modelling T_p.

---

## Evaporation with ODE breakup (ODE model 4)

**Model.** Per droplet of a stream with number rate `n`, `dm/dt = ṁ_evap − m·(dn/dt)_brk/n` and
`dn/dt = (dn/dt)_brk`, so `d(n m)/dt = n ṁ_evap`: breakup conserves the stream mass (`[RD87]`, "with each
change in drop size the drop number was changed correspondingly to conserve liquid mass"; `[RD86]` App. A1,
`n_f r_f³ = n_i r_i³`; `[Reitz87]` eq. 11 with `N a³ = N0 a0³`) and evaporation shrinks each drop at
constant `n`. In diameter, `dd/dt = (d_s − d)/τ + 2 ṁ_evap/(ρ_l π d²)`.

**No published frozen-gas benchmark.** IGLOO integrates non-interacting parcels through a frozen, one-way gas,
so a spray reduces to independent single drops and every spray-level observable is out of reach:

| candidate | what it is | why rejected |
|---|---|---|
| Hsiang & Faeth, IJMF 18(5):635–652, 1992 | shock-tube secondary breakup, room temperature | no evaporation (the TAB/ETAB product sizes against `[HF92]`/`[HF93]`: [tab-hf92](../breakup/tab-hf92/INFO.md), [etab-hf92](../breakup/etab-hf92/INFO.md)) |
| Hiroyasu & Kadota, SAE 740715, 1974 | diesel spray drop-size distribution | measured distribution of a transient spray; no evaporation oracle |
| Hiroyasu & Arai, SAE 900475, 1990 | diesel spray tip penetration | set by two-way entrainment |
| Wu, Santavicca, Bracco & Coghe, AIAA J 22:1263–1270, 1984 | LDV in non-vaporizing dense diesel sprays | non-vaporizing, dense, two-way |
| Sommerfeld & Qiu, IJHFF 19:10–22, 1998 | isopropyl-alcohol hollow-cone spray in a heated pipe expansion | two-way turbulent; no breakup model in play |
| Chen & Pereira, IJHMT 39:441–454, 1996 | RANS computation of the Sommerfeld case | as above |
| O'Rourke & Amsden 1987, Tanner 1997 validation sprays | TAB/ETAB spray penetration and SMD | non-evaporating, transient, two-way (the kernels are gated by `tab-e2e`/`etab-e2e`) |
| Sirignano / Sazhin monodisperse spray | analytic two-way coupled dilute spray | with a frozen gas it reduces to a single drop |

**Chosen:** the composed single-drop oracle of [evap-breakup-box](evap-breakup-box/INFO.md) — the two kernels,
each source-verified in the suite (Reitz-Diwakar: `reitz-diwakar-e2e`, `test_breakup_rd`; d²-law: `d2law`,
`test_evap_probes` XE1), composed by the stream-mass rule above and integrated along the run's measured `u(x)`
and `T_p(x)`. The oracle is independent of the right-hand-side assembly, which is what model 4 adds. Twins
isolate one half each: [d2law-brk-dormant](d2law-brk-dormant/INFO.md) (breakup rate zero),
[../breakup/rd-evap-frozen](../breakup/rd-evap-frozen/INFO.md) (evaporation rate zero); the right-hand side at
one state is pinned by [evap-breakup](evap-breakup/INFO.md).
