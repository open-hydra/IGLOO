# Breakup — Literature-Grounded Verification Tests

Scope: secondary breakup models in [Lib_Breakup.f90](../../src/lib/Lib_Breakup.f90).
Selector `assign_breakup` ([L58](../../src/lib/Lib_Breakup.f90#L58)):
Pilch-Erdman(1), Reitz-Diwakar(2), Reitz-KHRT(3), TAB(4), ETAB(5).
Defaults read in [Lib_INI.f90 L308-402](../../src/lib/Lib_INI.f90#L308) (`read_breakup_constants`).

A theory review: each model's coded form against its paper; each section names its unit test.

**Verified against the primary PDFs** (`papers/`): Pilch-Erdman and TAB fully verified (every
constant confirmed, incl. PE87's 0.75/3 velocity poly,
4.5/1.2/1.64 Oh-override, and TAB's K=10/3 child energy balance). KHRT: KH curve fits and RT
structure confirmed; the Ohnesorge Z and the Taylor T on the paper's fluid, length and Weber
(Reitz-87 Eqs 4-5) and the RT child count with its cube (Beale-Reitz Eq. 11). ETAB:
oscillator/Kbr/size law confirmed; the stripping branch uses k2, with the Tanner-1998 defaults
0.2222/We_t=80; the child v⊥=A·ẋ spray-cone component is in the code (Tanner Eqs 8–10,
gated ET4a/b/c: closed-form magnitude, orthogonality, A²<0 no-kick). Reitz-Diwakar: criteria, forms,
stable sizes confirmed; time-scale defaults follow the Reitz-1987 lineage, whose primary
("Structure of High-Pressure Fuel Sprays", 1987, `[RD87]`) gives the bag D=π (Eq. 7) and Cs=20.

---

## 0. DOI verification

| Tag | Paper | DOI | Resolves? |
|---|---|---|---|
| `[ORA87]` | O'Rourke & Amsden, TAB Method, SAE 872089 (1987) | 10.4271/872089 | **verified** — SAE page + OSTI biblio 6118786 confirm title/authors |
| `[RD86]` | Reitz & Diwakar, Effect of Drop Breakup on Fuel Sprays, SAE 860469 (1986) | 10.4271/860469 | **verified** — SAE Mobilus lists 860469 w/ matching title; the source PDF is read (`../REFERENCES.md`) |
| `[PE87]` | Pilch & Erdman, IJMF 13(6) 741-757 (1987) | 10.1016/0301-9322(87)90063-2 | **verified** — ScienceDirect PII 0301932287900632 = this DOI, title/pages match |

Cross-check implementations used as secondary authorities (not citable primaries):
ANSYS/CFX theory (TAB Eq. 24, reproduced in Bhandarkar-Manna-Chakraborty, *Int. J. Fluid Mech. Res.* 27(1) 2017); OpenFOAM v7 `ReitzDiwakar.C` (RD reference form). Both confirmed independently below.

---

## 1. TAB — model (4) — STRONGEST ORACLE

### Paper / equation
[ORA87]. Distortion `y` (equator displacement `x` normalized: `y = x/(Cb·r)`, breakup at `y>1`) obeys a linear damped driven oscillator (ANSYS/CFX Eq. 24, = ORA87):

    ÿ = (Cf/Cb)·(ρg/ρl)·(u²/r²)  −  (Ck·σ/(ρl·r³))·y  −  (Cd·μl/(ρl·r²))·ẏ

ORA87 constants: **Cb=1/2, Ck=8, Cd=5, Cf=1/3, K=10/3**.

Grouped forms:
- damping   `1/td = (Cd/2)·μl/(ρl·r²)`
- frequency `ω² = Ck·σ/(ρl·r³) − 1/td²`
- equilibrium (steady, ÿ=ẏ=0) `y_eq = (Cf/(Cb·Ck))·We_r` with radius Weber `We_r = ρg·u²·r/σ`.
  With the constants `Cf/(Cb·Ck) = (1/3)/((1/2)·8) = 1/12` ⇒ **y_eq = We_r/12**.

> ⚠ Grouping: `We_eff = Cf/Ck·(…)` omits the `1/Cb`. Correct grouping is `Cf/(Cb·Ck)=1/12`, NOT `Cf/Ck=1/24`; the code reproduces 1/12 (below).

Closed form at constant forcing (heavy/fixed-velocity drop ⇒ `u`,`r` const ⇒ constant coeffs), with `We_eff ≡ y_eq`:

    y(t) = We_eff + e^(−t/td)·[ (y0−We_eff)·cos ωt + ( (y0−We_eff)/td + ẏ0 )·sin(ωt)/ω ]

### Code cross-check — MATCHES
[TABmodel L440-441](../../src/lib/Lib_Breakup.f90#L440):
`rdt = 0.5·bp(2)·mup/(rhop·radius²)` = `(Cd/2)μl/(ρl r²) = 1/td`, with `bp(2)=Cmu=5=Cd` ✓
`omega = bp(1)·sigma/(rhop·radius³) − rdt²` = `Ck σ/(ρl r³) − 1/td²`, with `bp(1)=Comega=8=Ck` ✓ (stored SQUARED).
[L445-446](../../src/lib/Lib_Breakup.f90#L445): `We = radius·vel²/sigma·rhog` = We_r ✓; `WeCr = We/bp(3)`.
[YupdateTAB L644-645](../../src/lib/Lib_Breakup.f90#L644):
`y = WeCr + exp(−dt·rdt)·(y1·cost + (yDot+y1·rdt)·sint/omega)` with `y1=y−WeCr`, `cost=cos(ω dt)`, `sint=sin(ω dt)`.
⇒ term-for-term identical to the closed form with `We_eff=WeCr`. **EXACT MATCH.** `yDot` line is the analytic derivative — also exact.

Product-drop size [L473](../../src/lib/Lib_Breakup.f90#L473): `k = 1 + 4/3 + ρl r³/(8σ)·ẏ²`, `rs=r/k`.
ORA87 energy balance: `r/r_new = 1 + (8K/20) + (ρl r³ ẏ²/σ)·(6K−5)/120`. With **K=10/3**: `8K/20=4/3` ✓ and `(6K−5)/120=1/8` ✓. **K=10/3 confirmed, MATCHES.**

`bp(3)` doubling — CORRECT DESIGN, not a bug. [Lib_INI.f90 L367-368](../../src/lib/Lib_INI.f90#L367): `WeCrit` default 6, then `bp(3)=2·bp(3)` (the `if(err) x=6; x=2*x` semicolon-guard makes the ×2 unconditional — INTENTIONAL). So `bp(3)=12` = the equilibrium denominator `Cb·Ck/Cf`. `WeCr = We_r/12 = y_eq` ✓. A drop from rest (y0=ẏ0=0) has amplitude `a=WeCr`, first peak `y_max = WeCr+a = 2·WeCr`; breakup gate `a+WeCr>1` ⇒ `We_r > bp(3)/2 = 6` = canonical TAB rest-drop critical Weber. User-facing `WeCrit` IS that rest threshold; ×2 maps it to the equilibrium factor. Self-consistent — no discrepancy.

### Verification test T-TAB (`tab/test_breakup_tab.f90`, TAB1-TAB4)
Non-tautological triangulation (YupdateTAB is literally the closed form, so code-vs-formula alone only checks transcription), through the public `breakupEvent` entry at a fixed slip `u` and radius `r`:
- **Oracle A**: the analytic damped oscillator, with `ω, td, We_eff=We_r/12` derived *independently* from {Ck=8,Cd=5,Cf=1/3,Cb=1/2}.
- **Oracle B**: RK4 of the raw ODE `ÿ = (Cf/Cb)(ρg/ρl)(u²/r²) − (Ck σ/ρl r³)y − (Cd μl/ρl r²)ẏ`. Confirms A *solves* the ODE (guards against a mistranscribed constant that would agree between code and A but not B).
- **TAB1**: the code's stepping vs A over 500 steps at `We_r = 4` (sub-critical, no event), `y` and `ẏ/ω` within 1e-9 (round-off bound).
- **TAB2**: B vs A on the same drop, within 1e-8 (RK4 truncation plus round-off).
- **TAB3**: rest drop (`y0=ẏ0=0`), `We_r ∈ {5, 5.9, 6.1, 12}` ⇒ the event fires iff `We_r>6` (checks the `bp(3)` design above).
- **TAB4**: breakup time at `We_r = 9`, weak damping: the code's `tb` ([L466](../../src/lib/Lib_Breakup.f90#L466)), bisected on the event flag, vs B's first `y=1` crossing, within the damping-neglect bound plus the bisection width.

---

## 2. ETAB — model (5)

[ORA87 §ETAB extension / Tanner 1997]. Same oscillator (`bp(5)=Comega=8`, `bp(6)=Cmu=5`, `bp(3)=WeCrit×2=12`), but product size from exponential mass-rate law instead of energy sampling.
[L543](../../src/lib/Lib_Breakup.f90#L543): `AWe = (k2/k1·√WeTrans − 1)/WeTrans⁴`, defaults `k1=k2=0.2222`, `WeTrans=80` ([Lib_INI.f90 L389-396](../../src/lib/Lib_INI.f90#L389)).
[L580-584](../../src/lib/Lib_Breakup.f90#L580): `Kbr = k1·ω·(AWe·We⁴+1)` for `We≤WeTrans`, `Kbr = k2·ω·√We` for `We>WeTrans`; `r_new = r·exp(−(Kbr/ω)·acos(1−1/WeCr))`.

### Verification test T-ETAB (`etab/test_breakup_etab.f90`, ET1-ET4)
**Continuity of `Kbr` at `We=WeTrans`.** Below: `Kbr = k1·ω·(AWe·WeTrans⁴+1) = k1·ω·(k2/k1·√WeTrans) = k2·ω·√WeTrans`. Above: `k2·ω·√WeTrans`. Continuous for any k1, k2: ET2 asserts it at the default k1=k2, ET3 for k1≠k2 (the ln(r_new/r) ratio across WeTrans is 1 to 1e-6); ET1 is the sub-threshold drop, ET4 the product kick. Oscillator half re-uses T-TAB.

---

## 3. Pilch-Erdman — model (1)

### Paper
[PE87]. Diameter Weber `We = ρg·u²·d0/σ`. Dimensionless total breakup time `T` (nondim by `t* = (d0/u)·√(ρl/ρg)`), piecewise:

| We range | T |
|---|---|
| 12 < We ≤ 18 | 6·(We−12)^(−0.25) |
| 18 < We ≤ 45 | 2.45·(We−12)^(0.25) |
| 45 < We ≤ 351 | 14.1·(We−12)^(−0.25) |
| 351 < We ≤ 2670 | 0.766·(We−12)^(0.25) |
| We > 2670 | 5.5 |

Critical Weber (Oh-corrected): `Wc = 12·(1 + 1.077·Oh^1.6)`, `Oh = μl/√(ρl d σ)`.
Max stable fragment: `d_stable = Wc·σ / (ρg·((1−V)·u)²)` where `V` = accumulated dimensionless drop velocity at breakup.

### Code cross-check — T(We), Wc(Oh) and d_stable MATCH
[PilchErdman L142-163](../../src/lib/Lib_Breakup.f90#L142): `We = dp·rhog·vel²/sigma` (diameter ✓, unlike the other models which use radius). `Wec = 12(1+1.077·Oh^1.6)` ✓ [L144]. Branch table [L149-161] = paper table **exactly** (5.5 / 0.766 / 14.1 / 2.45 / 6, with (We−12)^±0.25, breakpoints 2670/351/45/18/12) ✓.
`tau = TBT·dp/(vel·√(ρg/ρl))` [L169] = `T·t*` ✓.
Oh-override [L163]: `Oh>0.1 & We<228 ⇒ T = 4.5(1+1.2·Oh^1.64)` — the PE87 viscous-drop breakup-time correction; coefficients (4.5, 1.2, 1.64) verified against the PE87 PDF.
`VdV = √(ρg/ρl)·(0.75·Cd·T + 3·B·T²)` [L164], `dStable = Wec·σ/(ρg v²(1−VdV)²+toll)` [L165]. The `d_stable` structure matches PE87, and so does the truncated 2-term velocity-history fit (0.75, 3), with the defaults `Cd=bp(1)=1`, `B=bp(2)=0.116` ([Lib_INI.f90 L322-326](../../src/lib/Lib_INI.f90#L322)).

### Verification tests (`pilch-erdman/test_breakup_pe.f90`, PE1-PE4)
- **T-PE-time** (PE2): the rate at 12 Weber numbers, inside the branches and just inside each break (17.9/18.1, 44.9/45.1, 350/352, 2669/2671), each matching the table to 1e-12 (with `B` = 0.01, so that the low branches break); continuity is not asserted (the table's rounded constants meet to 0.2 % at the breaks). Confirms transcription + branch bounds (`>` vs `≥`). Oracle = the table above.
- **T-PE-Wec**: `Wc(Oh)` enters the PE2-PE3 rates through the oracle's own Wec (Oh ≈ 0.01 and 0.2); no separate sweep at Oh ∈ {0, 0.1, 1} (expect {12, 12·(1+1.077·0.1^1.6), 12·2.077}).
- **T-PE-dStable**: PE4 pins the `dp < dStable` guard (We = 20, default B); the d_stable plateau against PE87's closed form is gate 2 of [pilch-erdman-e2e](pilch-erdman-e2e/INFO.md) (tol 12 %).

---

## 4. Reitz-Diwakar — model (2)

### Paper / reference form
[RD86], reference implementation = OpenFOAM v7 `ReitzDiwakar.C`. OpenFOAM `We = 0.5·ρg·u²·d/σ` = **radius Weber** (= IGLOO convention). Two regimes:
- **Bag**: `We > Cbag` (Cbag=6). `d_stable = 2·Cbag·σ/(ρg u²)`. `τ_bag = Cb·d·√(ρl d/σ)`, OpenFOAM Cb=0.785.
- **Stripping**: `We > Cstrip·√Re` (Cstrip=0.5). `d_stable = (2 Cstrip σ)²/(ρg u³ μg)`. `τ_strip = Cs·d·√(ρl/ρg)/u`, OpenFOAM Cs=10.

### Code cross-check — MATCHES (constant repackaging)
[ReitzDiwakar L182-192](../../src/lib/Lib_Breakup.f90#L182): `We = radius·ρg·vel²/σ` = OpenFOAM's `0.5·ρg u² d/σ` ✓.
Bag: `We>bp(1)=6` ✓; `dStable = 2σ·bp(1)/(ρg vel²) = 12σ/(ρg u²)` ✓; `τ = bp(2)·dp·√(0.0625·ρl dp/σ) = (π/4)·dp·√(ρl dp/σ) = 0.785·dp√(ρl dp/σ)` (bp(2)=π, √0.0625=1/4 ⇒ π/4=0.785=OpenFOAM Cb) ✓.
Stripping: `We>bp(3)√Re`, bp(3)=0.5 ✓; `τ = bp(4)·0.5·dp·√(ρl/ρg)/u = 10·dp√(ρl/ρg)/u` (bp(4)=Cs=20, ×0.5 ⇒ 10=OpenFOAM Cs) ✓; `dStable=(2·bp(3)σ)²·Re/(dp ρg² u⁴)` = OpenFOAM `(2 Cstrip σ)²/(ρg u³ μg)` **iff `Re = ρg u dp/μg`** (gas Reynolds), which is the `Re` `breakupOde` receives ([Lib_RHS.f90 L804](../../src/lib/Lib_RHS.f90#L804): `Re = gas(1)·slip·d/gas(6)`). **All RD constants match once repackaged.**

### Verification test T-RD (`reitz-diwakar/test_breakup_rd.f90`, RD1-RD4)
Unit checks: (i) bag regime `We_r∈(6, 0.5√Re)` ⇒ `dStable=12σ/(ρg u²)`, `τ` per formula; (ii) stripping `We_r>0.5√Re` ⇒ stripping branch; (iii) sub-critical `We_r<6` ⇒ `npdotDot=0`; (iv) both sides of the handoff at `We_r=0.5√Re`. Oracle = OpenFOAM formulas above.

---

## 5. Reitz-KHRT — model (3)

### Paper
`[Reitz87]` Reitz 1987 KH linear stability: max growth rate and wavelength
`Ω_KH·√(ρl a³/σ) = (0.34 + 0.38 We_g^1.5)/((1+Oh)(1+1.4 T^0.6))` ⇒ `Ω_KH = […]·√(σ/(ρl a³))` [units 1/s];
`Λ_KH/a = 9.02(1+0.45 Oh^0.5)(1+0.4 T^0.7)/(1+0.87 We_g^1.67)^0.6`.
RT: `Ω_RT = √[2·(g_t(ρl−ρg))^1.5 / (3√(3σ)·(ρl+ρg))]`, `Λ_RT = 2π C_RT/√(g_t(ρl−ρg)/(3σ))`.

### Code cross-check
`lambdaKH` [L225](../../src/lib/Lib_Breakup.f90#L225): matches Reitz KH wavelength **exactly** ✓.
`omegaRT` [L232](../../src/lib/Lib_Breakup.f90#L232), `lambdaRT` [L233](../../src/lib/Lib_Breakup.f90#L233), `tauRT=bp(3)/omegaRT` (Ctau=1): match ✓.
`omegaKH` [L223-224](../../src/lib/Lib_Breakup.f90#L223) (and its copy on the event path, [L331-332](../../src/lib/Lib_Breakup.f90#L331)):
`(0.34+0.38·WeGas^1.5)/((1+Oh)(1+1.4·Tay^0.6)) * sqrt(sigma/(rhop*radius**3))` = the Reitz growth rate, units 1/s ✓.
Consequently `tauKH = 3.726·B1·r/(Λ_KH·omegaKH)` [L227] and `dStable=2·B0·Λ_KH` [L228] are dimensionally sound ✓.

### Verification test T-KHRT (`reitz-khrt/`)
`test_breakup_khrt.f90` (KH1-KH3: the pure-KH rate against the Reitz-87 paper chain, Z and T as the paper defines them, 1e-12), `test_kh_rayleigh_limit.f90` (the fit's `We_g→0, Z→0` corner against Rayleigh's 1878 dispersion relation), `test_khrt_interaction.f90` (the ODE and event copies of the KH and RT blocks agree on the RT predicate `Λ_RT<d`, `tc>τ_RT`; the RT branch lands the parcel on `Λ_RT`; the KH shed's stripped-mass partition). The cubic child count has no pointwise test. The KH-mass-shed branch is live code ([L357-401](../../src/lib/Lib_Breakup.f90#L357)).

---

## The tests, strongest oracle first

1. **T-TAB (TAB closed-form oscillator)** — exact analytic oracle, all constants (Ck,Cd,Cf,Cb,K) verified vs [ORA87]/ANSYS Eq.24, triangulated code↔analytic↔RK4. Cleanest, strongest.
2. **T-PE-time** — PE87 table verified exactly; pure branch/transcription unit check, trivial oracle.
3. **T-PE-Wec** — `Wc(Oh)` closed form, inside the PE rate checks.
4. **T-RD** — matches OpenFOAM after repackaging, with the gas Reynolds number on the diameter.
5. **T-ETAB** — Kbr continuity at WeTrans (for any k1, k2) and the product kick.
6. **T-PE-dStable** — the `dp < dStable` guard; the d_stable plateau end to end.
7. **T-KHRT** — the paper chain, the Rayleigh limit and the KH↔RT path agreement.

---

## Breakup with evaporation (ODE model 4) and events under a mass-state model

- **ODE breakup with evaporation** (Pilch-Erdman, Reitz-Diwakar, Reitz-KHRT with `evaporation` set) runs as
  ODE model 4: the droplet mass and the number rate are both ODE states, and the breakup share
  `−m·(dn/dt)/n` of `dm/dt` keeps `n d³` constant under breakup. Literature judgment and the composed oracle:
  [../evaporation/LITERATURE_TESTS.md](../evaporation/LITERATURE_TESTS.md), "Evaporation with ODE breakup".
  Gates: [rd-evap-frozen](rd-evap-frozen/INFO.md) (the `[RD87]` rate oracle under model 4),
  [khrt-evap-frozen](khrt-evap-frozen/INFO.md) (KH rate, RT persistence, children, source budgets, the exit
  flows of parents and children adding up to the injected flow) and its Eulerian twin
  [khrt-evap-frozen-euler](khrt-evap-frozen-euler/INFO.md) (children keep their model-4 state; the Eulerian mass);
  model 3's Eulerian twin: [khrt-e2e-euler](khrt-e2e-euler/INFO.md).
- **Events under a mass-state model** (TAB/ETAB with evaporation run as model 2; KH-RT with evaporation as
  model 4): the event's diameter and number rate become the ODE state together, `N^(n+1) = N^n (r^n/r^(n+1))³`
  (`[ORA87]`). Gates: [tab-evap-frozen](tab-evap-frozen/INFO.md), [etab-evap-frozen](etab-evap-frozen/INFO.md)
  (ETAB's product velocity kick included).
- **KH shed and the gas source** (`[Reitz87]` product-parcel rule): the stripped mass leaves the parent as a
  child parcel born at the shed point, so the parent hands the child's birth flux on instead of depositing
  it. Gates: `khrt-e2e/check.py::check_source_budget`, [khrt-shed-noexchange](khrt-shed-noexchange/INFO.md).

---

## TAB and ETAB product sizes against Hsiang & Faeth

`[HF92]` eq. [14] (SMD/d0 = 6.2 (ρ_L/ρ_G)^1/4 [μ_L/(ρ_L d0 u0)]^1/2 after bag, multimode and shear breakup, the
standard deviation of its coefficient 20 %) and `[HF93]` (the root-normal size distribution with MMD/SMD = 1.2 in all
three regimes once the core drop is removed; the size/velocity correlation eq. [13] at the end of breakup, t_b = 5 t*)
correlate shock-tube drops; they are not closed forms of TAB or ETAB, so they are drawn, not gated.
[tab-hf92](tab-hf92/INFO.md) and [etab-hf92](etab-hf92/INFO.md) run `[HF93]`'s water series (d0 = 1000 µm, We 15–600,
the air and water of its Table 1) and gate the figure's preconditions only (every drop breaks, the cascade ends inside
the box, every parcel is in it at t_b). TAB's SMD falls about as We^−0.87 against eq. [14]'s We^−1/4: 1.54× eq. [14]
at We 15, inside its ±20 % at We 37.7–59.8, 6.5× below it at We 600; its MMD/SMD, 0.968–1.375 (mean 1.164), is
the scatter of 25 to 33 draws about TAB's own (ln 2)^1/n Γ(1 − 1/n) = 1.149 (n = 3.5), not 1.2.
ETAB's cascade gives one size per class (two where it ends at the onset We_r = 6), 3.3× eq. [14] at We 15 and 5.1×
below it at We 600. Neither model has `[HF93]`'s core drop.
