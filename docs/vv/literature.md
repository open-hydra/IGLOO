# Literature Tests

The literature-grounded unit tests compile isolated Fortran drivers against the
IGLOO library and exercise each physical model against results derived from the
source papers.  Each family is independent: no ODE solver, no box mesh, no
injection — just a targeted call to the production physics routine compared against
a hand-computed or RK4-oracle reference.  All families are registered in CTest
under `-DBUILD_VERIFICATION=ON` and run by `tests/test.sh`.

Every family on this page is a **routine-level implementation check** — kind 1 in
the [taxonomy](index.md#kinds-of-verification).  Because the call is direct, these
tests are the sharpest instrument available for the *formula*: a wrong coefficient
shows up pointwise, at machine precision, with nothing to hide behind.  For the
same reason they cannot see whether the routine is reached at all, whether its
inputs arrive intact, or whether the surrounding integration uses its result — the
A16/A18/A19 finds all passed their unit family and were caught only by the
end-to-end paper reproductions in [End-to-End Cases](e2e.md).  Where the test's
oracle is a re-coding of the same paper the production code was written from, the
two can also lockstep on the same error (A13 did); the fix for that is an external
reference, not a tighter tolerance.

One family here is not a physics closure: `test_gas_reconstruction`
(§ Infrastructure) checks the trilinear interpolation and hex-inverse machinery,
and is the routine-level member of the base-feature kind.

```bash
./tests/test.sh unit         # all unit families
./tests/test.sh standard     # drag + temperature families only
./tests/test.sh breakup      # breakup families only
./tests/test.sh infrastructure
```

The full bibliography (citation tags `[CGW78]`, `[RM52]`, etc.) is in
`tests/REFERENCES.md`.

---

## Drag — family A (`tests/standard/drag`)

Tests the 13 drag closures in `Lib_Drag.f90`, exercised through the production
`rhsStandard` path via `verif_driver`.

**What is tested.** The Stokes limit `Cd·Re → 24` for every full-range correlation
(internal self-check, no external citation needed); the relaxation integral for
Stokes (`[CGW78]`) and Schiller–Naumann (`[SN33]`) against independent oracles; the
SDIRK4 tableau order ($p \approx 4$) on the linearized Stokes RHS; and probes XD1–XD5
for the model-formula bugs A1/A5/A6 (fixed 2026-07-07) and the Wen-Yu value
correction.

**Oracle type.** Closed-form exponential relaxation for Stokes (A1); independent
Python RK4 re-implementation of the Schiller–Naumann RHS (A3); direct algebraic
evaluation at tabulated Re for the probe tests.

**Literature references.** `[CGW78]` (Clift, Grace & Weber 1978 — all correlations);
`[SN33]` (Schiller & Naumann 1933).

**Status.** `test_drag_lit` GREEN (DL1–DL3); `test_drag_probes` **GATED GREEN**
(XD1–XD5 — bugs A1/A5/A6 + Wen-Yu C-flag fixed 2026-07-07).

<figure>
  {% include "vv/images/unit-standard-drag.svg" ignore missing %}
</figure>

---

## Temperature — family B (`tests/standard/temperature`)

Tests the six Nusselt correlations in `Lib_Heat.f90`, exercised through the
production `rhsStandard` path.

**What is tested.** The Nu = 2 conduction floor at Re = 0 for Ranz–Marshall, JAXA2,
and JAXA3 (analytic, source-free); the full Ranz–Marshall correlation `[RM52]`
against a table of (Re, Pr) values; SDIRK4 tableau order on the linearized (Nu = 2)
thermal RHS; and the full coupled (velocity + temperature) ODE under Ranz–Marshall
against an independent RK4 re-implementation.  JAXA4 ≡ Kavanau–Drake bridging is
verified algebraically (S1 in `standard/temperature/LITERATURE_TESTS.md`).

**Oracle type.** Analytic Nu = 2 anchor (source-free); Ranz–Marshall table values
from `[RM52]`; independent RK4 of the coupled system for B2; machine-zero inter-model
identity checks.

**Literature references.** `[RM52]` (Ranz & Marshall 1952).

**Status.** `test_heat_nu` GREEN (HN1–HN3); `test_heat_probes` GREEN — A7 closed
2026-07-15/16 as SOURCE-FAITHFUL (the published correlation, Shimada2006 eq. 45
quoting NASA-SP-8039 Table II, has no +2 floor; the probe now pins the exact
transcription).  JAXA4 base coefficient 0.645 → **0.654** (A12, both sources).

<figure>
  {% include "vv/images/unit-standard-heat.svg" ignore missing %}
</figure>

---

## Evaporation — family C (`tests/evaporation`)

Tests the four evaporation closures in `Lib_Evaporation.f90` at function level.

**What is tested.** The Sherwood / Nusselt prefactor structure at Re = 0 for all
four models (d²-law, CEM, CEM-B, ASM); cross-model Re = 0 consistency (CEM/CEM-B/ASM
must agree to machine precision; d²-law agrees at $\mathrm{Nu}=2$); the
Abramzon–Sirignano $F(B)$ correction factor `[AS89]`; Ranz–Marshall Sh coefficient
in CEM `[RM52]`.  Probes: XE1 verifies the d²-law / CEM ratio ($= 1$ after A3 fix);
XE2 verifies $\dot{m} > 0$ from the corrected call sites (A4 fix); XE3 verifies
$\dot{Q}_G > 0$ for hot gas (A9 fix).

**Oracle type.** Algebraic (function-level direct calls); cross-model consistency
(source-free).

**Literature references.** `[RM52]`; `[Lef89]` (Lefebvre 1989 — d²-law canonical
constant $K = 8 k_g \ln(1+B_T)/(\rho_l c_{pg})$); Abramzon & Sirignano (1989) eq. 20/24.

**Status.** `test_evap_probes` **GATED GREEN** (XE1/XE2/XE3 — bugs A3/A4/A9 fixed
2026-07-07).  The end-to-end d²-law case (`evaporation/d2law`) is also GATED GREEN —
see [Evaporation](evaporation.md).

<figure>
  {% include "vv/images/unit-evaporation-cem.svg" ignore missing %}
</figure>

### Non-equilibrium interface — LK (`tests/evaporation/interface-neq`)

Langmuir–Knudsen interface depression (MHB98 model M2) around the CEM gas-side
rate: fixed-point self-consistency, deep-Picard oracle agreement, VLE limit,
the exact $1/(p\,d)$ scaling inversion, and a 27-point convergence envelope
(LK1–LK6).  **GATED GREEN.**

<figure>
  {% include "vv/images/unit-evaporation-lk.svg" ignore missing %}
</figure>

### Tonini–Cossali analytical (`tests/evaporation/tc-analytic`)

TC2012 eq. 9 Stefan–Fuchs rate: residual self-oracle on a 27-point grid, boiling
clamp, exact isothermal Stefan–Fuchs limit, CEM-collapse identity, heat coupling,
and the bracket theorem (TC1–TC7).  **GATED GREEN.**

<figure>
  {% include "vv/images/unit-evaporation-tc.svg" ignore missing %}
</figure>

### Stefan-blowing reduction $f_2$ — MHB98 Fig. 4 decane

`tests/evaporation/test_mhb98_decane.f90` — E-VAL-2b, the **model-discriminating**
MHB98 case: Fig. 4(b) spans ~200 K across
their eight models, against 0.47 K in the Fig. 2(b) water case gated by
[mhb98-water](e2e.md#mhb98-water).  What separates them here is *not*
non-equilibrium — at $D_0 = 2$ mm the Langmuir–Knudsen depression is
$L_K\beta/(D/2) = 4.9\cdot10^{-5}$ against $\chi_{s,eq}\approx0.63$, as vacuous as in
Fig. 2 — but the **evaporative heat-transfer reduction** $f_2$ (their eq. 19):

$$
f_2 = \frac{\beta}{e^{\beta}-1}, \qquad
\beta = -\tfrac{3}{2}\,\mathrm{Pr}_G\,\tau_d\,\frac{\dot m}{m_d} \quad\text{(eq. 17)}
$$

Here $\beta \approx 1.6 \Rightarrow f_2$ spans **0.40–0.88** over the droplet's life, i.e.
the convective heat is cut by up to 60 %.  IGLOO carries it as the *opt-in* global
selector `[IGLOO-Models] blowing = none|LK`, **off by default** (it would move every
evaporating case), applied to the convective $\dot Q$ only — never to $\dot m L_v$,
which eq. (3) keeps unreduced — and **skipped when the evaporation model already
returns a blowing-corrected $\dot Q$** (ASM's $1/B_T$, TC's $1/(e^{\chi}-1)$), so the
two reductions cannot compound.

A **unit** test, not an e2e box case: MHB98's Fig. 4 droplet is at **fixed slip**
(suspended drop, $\mathrm{Re}\propto d$), established by measurement — free flight
gives $\beta = 1.31$–$1.44$ for *every* $\mathrm{Sc}$ against their stated $\approx1.6$.
Three legs: **DC1** production `evaporation` (CEM+LK) $\dot m$ vs an independently
re-coded M7 chain ($3.5\cdot10^{-16}$); **DC2** `blowingFactor` vs eq. 19
($9.4\cdot10^{-15}$); **DC3** the re-code at MHB98's *own* Sh/Nu coefficient 0.552 vs
three pixel-measured Fig. 4 quantities — $\beta$ 1.3 %, $T_d(3.5\,\mathrm{s})$ 3.3 K,
$D^2(4.0\,\mathrm{s})$ 1.2 %.  **GATED GREEN** (5/5).

**Documented deviations, quantified not hidden.** IGLOO's Ranz–Marshall carries
$0.600$ where `[MHB98]` eq. (6) prints $0.552$ (8.7 %), so DC1/DC2 run the oracle at
0.600 ("does IGLOO integrate its own model correctly") and DC3 at 0.552 ("does the
model reproduce the figure"), with IGLOO-at-0.600 reported non-gating.  The case runs
$\mathrm{Le}=1$ **deliberately**: MHB98's own Appendix-A $\Gamma_V$ for decane gives
$\mathrm{Sc}=3.018$, $\mathrm{Le}=4.40$, which misses the measured plateau by +40 K —
their Appendix and their Fig. 4 are mutually inconsistent.

<figure>
  {% include "vv/images/unit-evaporation-mhb98-decane.svg" ignore missing %}
</figure>

---

## Breakup — family D (`tests/breakup`)

Tests the five secondary-breakup models in `Lib_Breakup.f90`.  Each sub-family has
its own directory and `INFO.md`.

### TAB / ETAB (`tests/breakup/tab`, `tests/breakup/etab`)

`[ORA87]` (O'Rourke & Amsden 1987).

**TAB.** The distortion ODE is a linear damped driven oscillator.  At constant forcing
the closed form is:

$$
y(t) = y_\mathrm{eq} + e^{-t/t_d}
  \!\left[(y_0 - y_\mathrm{eq})\cos\omega t
  + \frac{(y_0 - y_\mathrm{eq})/t_d + \dot{y}_0}{\omega}\sin\omega t\right]
$$

with constants $C_k=8$, $C_d=5$, $C_f=1/3$, $C_b=1/2$, $K=10/3$ verified
term-for-term against `[ORA87]` (and the ANSYS/CFX theory guide transcription).
`YupdateTAB` is an exact closed-form transcription — all constants confirmed.

**ETAB.** The `Kbr` mass-rate coefficient is continuous at `We = WeTrans` iff
$k_1 = k_2$ (default); the test is a continuity probe sweeping `We` through
`WeTrans`.

**Status.** **GREEN.** `test_breakup_tab` (TAB1 oscillator vs the closed form,
TAB2 independent RK4 of the raw `[ORA87]` ODE, TAB3 rest-drop critical Weber,
TAB4 breakup time via bisection); `test_tab_moments` (child-size $E[r]$, $E[r^2]$
over 20 000 sampled breakups vs a truncated Rosin–Rammler moment oracle,
CLT-band tolerance); `test_breakup_etab` (ET1–ET4 product-size ratio
$\ln(r_\mathrm{new}/r)$ vs the Tanner $K_\mathrm{br}$ branches, continuity across
`WeTrans`, ET4 child $v_\perp = A\dot{x}$).

<figure>
  {% include "vv/images/unit-breakup-tab.svg" ignore missing %}
</figure>

<figure>
  {% include "vv/images/unit-breakup-etab.svg" ignore missing %}
</figure>

### Pilch–Erdman (`tests/breakup/pilch-erdman`)

`[PE87]` (Pilch & Erdman 1987).

**What is tested.** The piecewise breakup-time table $T(\mathrm{We})$ at five
representative Weber numbers plus just inside each breakpoint — branch values and
discontinuities confirmed against the paper table; $W_c(\mathrm{Oh})$ critical Weber
at Oh = 0, 0.1, 1.

**Status.** `test_breakup_pe` **GATED GREEN** (PE1–PE4: the $T(\mathrm{We})$
breakup-time table with the Oh-override branch and the $d<d_\mathrm{stable}$
stability guard; all constants confirmed against `[PE87]`, 2026-07-15).

<figure>
  {% include "vv/images/unit-breakup-pe.svg" ignore missing %}
</figure>

### Reitz–Diwakar (`tests/breakup/reitz-diwakar`)

`[RD86]` (Reitz & Diwakar 1986); reference implementation cross-checked against
OpenFOAM v7 `ReitzDiwakar.C`.

**What is tested.** Bag and stripping branch stable-diameter and breakup-time
formulas, confirmed by repackaging the `bp` coefficients to match the OpenFOAM
reference form.

**Status.** `test_breakup_rd` **GATED GREEN** (RD1–RD4: bag/stripping stable-diameter
and breakup-time branches plus the regime handoff; the Reynolds convention was
resolved — `breakupOde` receives $Re = \rho_g\,|\mathrm{slip}|\,d/\mu_g$, which is
what RD needs).  The $C_s = 20$ stripping constant is **closed source-faithful**
(2026-07-23): its primary source `[RD87]` (SAE 870598) obtains $C=20$ explicitly
by curve-fitting, and the bag constant $D=\pi$ — so `[RD86]`'s earlier
"of order unity" wording does not apply to the calibrated form production uses.
The e2e companion is `reitz-diwakar-e2e` — see [End-to-End Cases](e2e.md).

<figure>
  {% include "vv/images/unit-breakup-rd.svg" ignore missing %}
</figure>

### ETAB — see TAB above.

### Reitz-KHRT (`tests/breakup/reitz-khrt`)

**What is tested.** The corrected $\omega_\mathrm{KH}$ chain against the Reitz (1987)
dimensional form: sub-critical $\omega_\mathrm{KH} = 0$ (KH1); two representative
$(We_g, Oh, Ta)$ points code ≡ Reitz-1987 formula to 1e-12 (KH2/KH3); pure-KH path
with acc = 0.  RT quantities, event path, and child-drop moments are deferred.

**Oracle type.** Direct algebraic evaluation against the published Reitz-1987 chain.

**Literature references.** `[RD86]` (Reitz & Diwakar / Reitz 1987 KH formulation).

**Status.** `test_breakup_khrt` **GATED GREEN** (bug A2 fixed 2026-07-07; KH1/KH2/KH3
pass to 1e-12; RT/event path deferred).

<figure>
  {% include "vv/images/unit-breakup-khrt.svg" ignore missing %}
</figure>

---

## Combustion — family M (`tests/combustion`)

Tests the Beckstead $d^n$ aluminium burn law (`Lib_Combustion.f90::becksteadRate`,
model 5): closed-form mass evolution $m(t)$ with complex-step derivative checks
(CB1/CB2), the bitwise ignition gate at $T_\mathrm{ign}$ (CB3), and the
$X_\mathrm{eff}$-exponent discriminator calibrated at the paper anchor (CB4).

**Literature references.** `[Beck05]` (Beckstead 2005); exponent $X_\mathrm{eff}^{1.0}$
primary-verified 2026-07-13.

**Status.** `test_combustion` **GATED GREEN** (CB1–CB4); e2e companion `burn-box`
GATED GREEN — see [End-to-End Cases](e2e.md).

<figure>
  {% include "vv/images/unit-combustion-beckstead.svg" ignore missing %}
</figure>

---

## Infrastructure — family E + T9

### Gas reconstruction (`tests/infrastructure/gas_reconstruction`)

Verifies `IGLOO_RayFaceIntersection3D`, `interp2ndOrder`, and the Newton inverse hex
map — the cell-tracking and interpolation routines used on every time step.

**Tests.** E1 constant field (machine eps); E2 multilinear field on parallelepiped
(machine eps); E3 smooth field convergence order ($p_\mathrm{obs} \approx 2 \pm 0.3$
over four refinement levels); E4 forward/inverse hex map round-trip (residual
$< 10^{-4}$); E5 partition-of-unity and node recovery.  E6 (Stokes-in-linear-field,
matrix-exponential oracle) deferred.

**Status.** E1–E5 GREEN; E6 scaffolded.

<figure>
  {% include "vv/images/unit-infrastructure-interp.svg" ignore missing %}
</figure>

### INI pipeline (`tests/infrastructure/ini_pipeline`)

Pins the `input.ini` → module-variable round-trip through `read_IGLOO_input` (FiNeR
parsing) for five representative parameter groups: TAB defaults, model selectors,
body force, ODE tolerances, and RNG seed.

**Status.** IP1–IP5 GREEN.  No overlay figure — a configuration contract has
nothing physical to plot.

### DB injection (`tests/infrastructure/db-injection`)

End-to-end case for the assigned-position (`AP`) injection path.  Places 5 particles
on a uniform box with `up=1`; checks placement fidelity, first-step $u = u_p$
hand-off via `obj_particle%vInj`, Stokes velocity relaxation, and domain exits.

**Status.** **GATED GREEN** (bug B1 fixed 2026-07-08).

---

## Production bugs found

The verification campaign catalogued **fifteen model-formula bugs (A1–A15)**, one
combustion-constant fix (M1), and six infrastructure bugs (B1–B6); as of 2026-07-16
**all are fixed, refuted, or closed as source-faithful** — every fix is gated by a
CTest probe or e2e oracle.  Highlights by physics area:

- **Drag** — A1 (Crowe/Hermsen exp placement) and A11 (their $g(Re)$ bridge factor,
  off by $g^2$) fixed and confirmed at the true primary (NASA-SP-8039 Table I);
  A6 Henderson blend $4/3$; Putnam/Wen-Yu plateaus source-faithful (0.4392 / 0.43,
  Shimada2006 authoritative per the 2026-07-16 decision).
- **Heat** — A12 (JAXA4 0.654) fixed; A7 (JAXA1 floor) closed *not-a-bug*: the
  published correlation has no +2 term (NASA-SP-8039 Table II).
- **Evaporation & droplet energy** — A3 (d²-law factor 2), A4 (scrambled gas args),
  A8 (LEB stub), A9 (ASM heat), A10 (spurious sensible-carry in F(7), ~65 % of the
  latent sink for water) all fixed.
- **Breakup** — A2 (Ω_KH sqrt), A13 (KHRT Z/T definitions, T was ~29× high),
  A14 (RT child count missing cube), A15 (ETAB stripping k₂ + Tanner defaults)
  fixed against the primary PDFs; Pilch-Erdman and TAB verified with zero
  discrepancies.

Each fix is gated by the test named in its row above; `tests/VERIFICATION_MATRIX.md`
maps every CTest entry to its reference source and oracle.
