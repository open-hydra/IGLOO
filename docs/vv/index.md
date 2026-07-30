# Verification & Validation

IGLOO's V&V suite lives in a single model-first tree under `tests/`, categorized
by physics: `standard/` (drag + heat), `evaporation/`, `breakup/`, and
`infrastructure/`.  Each category holds tests of two kinds, all registered in
CTest (opt-in via `-DBUILD_VERIFICATION=ON`).

**End-to-end cases** run the full production solver on a canonical box geometry,
then compare the trajectory output against an independent analytic oracle.  The
pass criterion is quantitative: the oracle's `check.py` must exit 0.

**Literature-grounded unit tests** compile isolated test drivers against the IGLOO
library and exercise each physical model against results derived from the source
papers.  Pass/fail is reported by CTest; a Python aggregation step writes
per-family CSV reports.

!!! note "Legacy test suite — retired"
    The legacy `test/` directory (zero-byte-stderr pass criterion, no physics
    verification) was retired on 2026-07-16 and parked out-of-git at
    `~/Desktop/Software/toBeRemoved/IGLOO-legacy-test/` as a manual byte-level
    A/B provider.  Its only physics-bearing case lives on in this suite as
    `infrastructure/db-2daxi`.  See [Testing](../development/testing.md).

---

## Kinds of verification

The tree is organized by *physics* (`standard/`, `evaporation/`, `breakup/`,
`combustion/`, `infrastructure/`), but the tests are not all the same **kind** of
verification, and a green gate licenses very different claims depending on which
kind it is.  The discriminator is a single question: **what is the reference the
gate compares against?**

| # | Kind | Reference is… | Independence |
|---|---|---|---|
| 1 | **Routine-level implementation check** | the literature expression re-coded in the test and evaluated at swept inputs, or an independent/deep-converged oracle, against a *direct call* to the production closure | full — but only of one function |
| 2 | **Tier V — closed-form solution** | the exact analytic solution of the model's governing equation, from case inputs only | full |
| 3 | **Run-conditioned kernel** | the model's rate law integrated along a quantity *read back from the run* (the measured $T_p$) | partial (semi-independent) |
| 4 | **Tier P — published-model reproduction** | something the paper itself published: its empirical correlation coded directly, or its own figure digitized | full, and external |
| 5 | **Base-feature / infrastructure exercise** | a feature contract — placement, transport, finiteness, accumulator identity — with either a Tier-V closed form or *no* reference curve at all (behavioral) | varies; stated per case |

Three further kinds — **6**, **7** and **8** — produce **no figure**, and are easy
to overlook when counting coverage:

- **6 — Conservation / consistency audits** — not standalone tests but secondary
  gates *inside* other cases: the mass-telescoping audits $\sum\dot w =
  \sum\dot n\,\Delta m$ in `lk-neq` and `burn-box`, the source-deposit identity in
  `coupled-body` (the deposit must be the drag reaction only), and the two-phase
  mass-closure gate in `mhb98-water` (injected vs deposited mass flow to 1 %, added
  with the A22 burnout-remnant fix).  These are the only checks that
  test IGLOO *against itself* rather than against a model, and they catch a class
  nothing else does — a trajectory can be right while the Eulerian feedback is
  not.
- **7 — Bug-transcription pins** — `test_drag_probes`, `test_heat_probes`,
  `test_evap_probes`.  They assert a *known* value (including deliberately
  source-faithful transcriptions with documented limitations), so they stay green
  while a fixed bug stays fixed and turn RED on regression.  They verify no model;
  they pin history.  Each was *written* as an expected-failure probe — exit 1 while
  the bug it documents is still live — but **every probed bug is now fixed** (A1/A5/A6
  + the Wen-Yu C-flag, A3/A4/A9, and A7 closed as *source-faithful*), so all three
  exit 0 and run as **ordinary gates**: no `WILL_FAIL` property is set on any of them,
  and a regression turns them red the normal way.  Their "promote probes to the gate"
  message is a leftover instruction to fold them into the parent families — they are
  already gating.
- **8 — Contract and harness meta-tests** — `test_ini_pipeline` (the
  `input.ini` → module-variable round-trip through the production reader) and
  `self_test` (which exercises the verification support library itself).  No
  physics, no curve, by design.

### Figure inventory by kind

Every SVG in `docs/vv/images/`, the test that emits it, and the reference its
gate is held to.  The file layout is deliberately **flat** — the sync loop in
`tests/test.sh` copies `OUTPUT/*.svg` and `curves-svg/*.svg` by basename, so
adding a case needs no harness edit and category membership stays a
documentation concern.

| Figure | Test (CTest name) | Kind | Gated against | Why this kind |
|---|---|---|---|---|
| `unit-standard-drag.svg` | `test_drag`, `test_drag_lit`, `test_drag_probes` | 1 routine-level | `[CGW78]`/`[SN33]` expressions at swept Re; `rk4_ref` of the Schiller–Naumann RHS | production `Lib_Drag::drag` called directly, no mesh/BC |
| `unit-standard-heat.svg` | `test_heat_nu`, `test_temperature` | 1 routine-level | `[RM52]` $\mathrm{Nu}=2+0.6\,\mathrm{Re}^{1/2}\mathrm{Pr}^{1/3}$ table; Kavanau–Drake collapse | production `Lib_Heat::heat` Nu vs coded correlation, pointwise |
| `unit-evaporation-cem.svg` | `test_evaporation` | 1 routine-level | independent Spalding chain $p_{sat}\to X_s\to Y_s\to B_M$; Abramzon–Sirignano $F(B)$ | algebraic function-level comparison |
| `unit-evaporation-lk.svg` | `test_interface_lk` | 1 routine-level | deep-Picard (1e-15) fixed point of the `[MHB98]` LK interface, 27-pt envelope | independent iteration, tighter than production's |
| `unit-evaporation-tc.svg` | `test_tc_analytic` | 1 routine-level | independent bisection root of the `[TC2012]`/`[ATC24]` $G(m)$ residual | different root-finder than production's Newton |
| `unit-evaporation-mhb98-decane.svg` | `test_mhb98_decane` | 1 routine-level **+ 4 Tier P (pixel-measured)** | DC1/DC2: an independently re-coded `[MHB98]` M7 chain and eq. 19 $f_2$, called against production at 25 sampled states (1e-12).  DC3: the re-code at MHB98's own 0.552 vs $\beta$, $T_d(3.5\,\mathrm{s})$ and $D^2(4.0\,\mathrm{s})$ **measured off Fig. 4 in pixels** | kernels called directly at prescribed fixed slip (kind 1), but the DC3 leg is the paper's own figure (kind 4); a unit test because the Fig. 4 drop is suspended, not free-flying |
| `unit-combustion-beckstead.svg` | `test_combustion` | 1 routine-level | closed-form $m(t)$ + complex-step $\mathrm{d}m/\mathrm{d}t$ vs `becksteadRate` `[Beck05]` | direct rate-function check |
| `unit-breakup-tab.svg` | `test_breakup_tab`, `test_tab_moments` | 1 routine-level | `[ORA87]` damped-oscillator closed form + `rk4_ref` of the raw ODE; Rosin–Rammler moment oracle | `breakupEvent` called directly |
| `unit-breakup-etab.svg` | `test_breakup_etab` | 1 routine-level | `[Tan97]`/`[Tan98]` $K_{br}$ branches, continuity across $We_\mathrm{trans}$ | product-size ratio, pointwise |
| `unit-breakup-pe.svg` | `test_breakup_pe` | 1 routine-level | `[PE87]` breakup-time table + Oh-override branch | `breakupOde` rate vs the paper table |
| `unit-breakup-rd.svg` | `test_breakup_rd` | 1 routine-level | `[RD86]`/`[RD87]` bag/stripping branches and the regime handoff | swept-$We_r$ pointwise formula check |
| `unit-breakup-khrt.svg` | `test_breakup_khrt` | 1 routine-level | `[Reitz87]` KH linear-stability chain to 1e-12 (pure-KH, acc = 0) | `breakupOde` vs coded $\Lambda/\Omega$ chain |
| `drag-stokes.svg` | `drag-stokes` | 2 Tier V | closed-form $x(v)$ Stokes relaxation, $\tau=\rho_p d^2/18\mu$ | exact solution of the momentum ODE, inputs only |
| `temp-relax.svg` | `temp-relax` | 2 Tier V | closed-form $T(x)$, $\mathrm{Nu}=2$ lumped capacitance | case built so $\mathrm{Re}=0$ ⇒ the ODE has a closed form |
| `body-force.svg` | `body-force` | 2 Tier V | closed-form terminal drift $v(x)=v_\infty(1-e^{-x/L})$ | `toll=1e-20` makes drag exactly linear ⇒ analytic |
| `conv-nu.svg` | `conv-nu` | 2 Tier V | closed-form $T(x)$ at frozen $\mathrm{Re}$, $\mathrm{Nu}=2.240$ `[RM52]` | injection at terminal velocity freezes Re ⇒ analytic |
| `d2law-line.svg` | `d2law-line` | 2 Tier V | input-only line $d^2=d_0^2-K_0(x-x_a)/u_g$, theory tolerance budget | $c_p\times100$ freezes $B_T$ ⇒ the textbook analytic case |
| `tab-e2e.svg` | `tab-e2e` | 2 Tier V | `[ORA87]` eq. 5 oscillator: onset $We_r=6$ + first-breakup time $t_{bu}$ | $t_{bu}$ is a root of the model's own analytic solution |
| `etab-e2e.svg` | `etab-e2e` | 2 Tier V | ORA87 $t_{bu}$ + the `[Tan97]` cascade ratio $\exp(-(K_{br}/\omega)\arccos(1-1/We_{Cr}))$ | $\omega$ cancels analytically ⇒ closed form in $We$ alone |
| `burn-box.svg` | `burn-box` | 2 Tier V | closed-form $d^n(x)=d_0^n-K_\mathrm{eff}x/u_g$ `[Beck05]` + independent RK4 energy balance | $d$ is $T_p$-independent above ignition ⇒ analytic |
| `d2law.svg` | `d2law` | 3 run-conditioned | Godsave–Spalding `[God53]`/`[Spa53]` kernel integrated along the **measured** $T_p(x)$ | $T$ is coupled ⇒ no closed form; oracle rides the run's $T_p$ |
| `lk-neq.svg` | `lk-neq` | 3 run-conditioned | `[MHB98]` LK-corrected CEM $d^2$-ODE, RK4 along measured $T_p$; VLE curve as the gated regression | same, plus a gated discrimination margin |
| `tc-box.svg` | `tc-box` | 3 run-conditioned | `[TC2012]` Stefan–Fuchs $d^2$-ODE along measured $T_p$; CEM as the gated regression | same |
| `tc-hexadecane.svg` | `tc-hexadecane` | 3 run-conditioned | variable-$\rho$ TC mass rate along measured $T_p$, $d^2$ reconstructed from mass + swelling assertion | the gate is the kernel; the digitized Fig. 11 overlay is **non-gating** |
| `pilch-erdman-e2e.svg` | `pilch-erdman-e2e` | 4 Tier P (correlation) | `[PE87]` $T^*(We)$ correlation (Fig. 7, five regimes) coded directly + the $d_\mathrm{stable}$ closed form | an empirical fit with no analytic solution — the paper's own numbers |
| `reitz-diwakar-e2e.svg` | `reitz-diwakar-e2e` | 4 Tier P (correlation) | `[RD87]` bag ($D=\pi$) and stripping ($C=20$, curve-fit) $\tau$ and $d_\mathrm{stable}$ | constants are the paper's calibration, not derived |
| `khrt-e2e.svg` | `khrt-e2e` **+** `khrt-e2e-rt` (two CTest entries, one figure) | 4 Tier P (correlation) **+** 5/audit | `khrt-e2e`: `[Reitz87]` KH $\lambda_{KH}$/$\Omega_{KH}$ correlations in $\dot d=(d_\mathrm{stable}-d)/\tau_{KH}$.  `khrt-e2e-rt`: **no external reference** — eight drops must shatter to the fragment scale, persist there, and stay mass-self-consistent | the KH rate is an analytic *reduction* of *empirical* correlations (the matrix tags it Tier V/P); the RT gate is behavioral + a conservation audit |
| `mhb98-water.svg` | `mhb98-water` | 4 Tier P (digitized), + 3 + audit | three gates, classified by the most external one: the **digitized `[MHB98]` Fig. 2 M7 line** ($\beta$ slope ±10 %, wet-bulb plateau ±0.5 K); also the LK kernel along measured $T_p$ (kind 3), the $K=(T_G-T_\mathrm{wb})8k_g/\rho_\ell L_v$ identity, and the two-phase mass closure (audit) | the reference is the paper's own plotted result for the paper's own case |
| `unit-infrastructure-interp.svg` | `test_gas_reconstruction` | 5 base feature (routine) | analytic constant/multilinear fields at machine eps; observed order $p\approx2$; hex Newton round-trip | interpolation, not a closure — the routine-level base-feature check |
| `vie-plait.svg` | `vie-plait` | 5 base feature (Tier-V oracle) | `[Vie15]` §5.1 closed-form damped oscillator per strand, scipy-validated to 1e-13 | drag already gated elsewhere; what is unique is the in-solver non-uniform gas rebuild (`gas-order = 2`) |
| `periodic-y.svg` | `periodic-y` | 5 base feature (Tier-V oracle) | the `body-force` closed form folded **modulo $L_y$**; velocity unchanged across each wrap | physics reused from `body-force`; what is unique is the bcdef-201 transport |
| `coupled-body.svg` | `coupled-body` | 5 base feature (Tier-V oracle + audit) | `body-force` $v(x)$ verbatim **plus** source totals vs closed forms (drag reaction only) | the only case exercising euler + source + body together; carries a conservation audit |
| `db-injection.svg` | `db-injection` | 5 base feature (**behavioral**) | no reference curve — placement fidelity, `vInj` hand-off, monotone relaxation, domain exit | assigned-position injection plumbing; bug-B1 regression gate |
| `db-2daxi.svg` | `db-2daxi` | 5 base feature (**behavioral**) | no reference curve — both parcels integrate, exit the nozzle outlet, all fields finite and in range | a real MOSE nozzle field has no analytic solution; deliberately md5-free |

Note that one entry is deliberately **hybrid**: `test_mhb98_decane` (E-VAL-2b,
MHB98 Fig. 4, $T_G=1000$ K decane) is the strong discriminator Fig. 2 is not — its
panel (b) spans ~200 K across the paper's eight models against 0.47 K in Fig. 2(b) —
and it carries both a kind-1 leg (production kernels vs an independently re-coded M7
chain, 1e-12) and a kind-4 leg (the same re-code at the paper's own coefficient
against three quantities *pixel-measured* off Fig. 4).  It is a **unit** test rather
than an e2e case because MHB98's droplet is suspended at fixed slip, which the box
solver could only emulate with a no-drag production flag.  See
[Literature-grounded families](literature.md).

---

## End-to-end case status

| Case | Physics exercised | Oracle | Status |
|------|-------------------|--------|--------|
| [Stokes drag](drag-stokes.md) | Stokes drag, velocity relaxation | Closed-form $x(v)$ | **GREEN** |
| [Temperature relaxation](temp-relax.md) | Lumped-capacitance heat, $\mathrm{Nu}=2$ (conduction limit) | Closed-form $T(x)$ | **GREEN** |
| [Body-force drift](body-force.md) | Transverse body acceleration + exact linear Stokes drag | Closed-form $v(x)$ | **GREEN** |
| [Convective Nu](../vv/e2e.md#conv-nu) | Ranz–Marshall heat at constant non-zero slip ($\mathrm{Re}=0.20$, $\mathrm{Nu}=2.24$) | Closed-form $T(x)$ | **GREEN** |
| [Vié plait](../vv/e2e.md#vie-plait) | Inertial particles in a compressive field $v_g=-\epsilon(y-1)$; trajectory crossing at $\mathrm{St}=5$ (first non-uniform-gas e2e) | Closed-form damped oscillator $y(x)$ per strand | **GREEN** (2026-07-20; Vié 2015 Eq. 5.4 sin sign-typo caught) |
| [Evaporation](evaporation.md) | d²-law mass transfer | Godsave kernel integrated along measured $T_p(x)$ | **GREEN** (bugs A3/A4 fixed 2026-07-07; 25/25 particles) |
| [d²-law analytic line](../vv/e2e.md#d2law-line) | d²-law in the frozen-$B_T$ textbook regime ($c_{p,p}\times100$, $T_p$ drift 0.76 K) | Input-only closed form $d^2 = d_0^2 - K_0 x/u_g$, theory tolerance budget | **GREEN** (E-VAL-1, 2026-07-22; worst residual $3\cdot10^{-4}$ vs $7\cdot10^{-4}$ budget) |
| LK non-equilibrium (`evaporation/lk-neq`) | Langmuir-Knudsen interface correction on CEM (MHB98 M2) | LK-corrected kernel RK4-integrated along measured $T_p(x)$; LK-vs-VLE discrimination; mass telescoping | **GREEN** (F1, 2026-07-11; telescoping tail-extrapolation 2026-07-15) |
| Tonini–Cossali (`evaporation/tc-box`) | TC analytical evaporation (Stefan–Fuchs) | Own re-derivation kernel along measured $T_p(x)$ | **GREEN** (F2, 2026-07-11) |
| [MHB98 water (validation)](../vv/e2e.md#mhb98-water) | Single water droplet vs a paper's own case ($D_0$=1.1 mm, quiescent); CEM+LK = MHB98 **M7** | IGLOO-vs-LK-kernel + wet-bulb identity + **IGLOO β vs the digitized M7 slope β_M7≈6.35e-3 (±10%)** + **wet-bulb vs M7 plateau (±0.5 K)**, all gated | **GREEN** (E-VAL-2; 2026-07-29 the case was set to MHB98's **own Appendix-A properties** at their T_R=293.97 K — IGLOO lands on M7, β **0.8%**, wet-bulb **0.10 K**; eyeball exp dropped — re-digitize follow-up) |
| [TC n-hexadecane (validation)](../vv/e2e.md#tc-hexadecane) | Variable-liquid-density TC evaporation (TC2012 Fig. 11: n-hexadecane swells then $D^2$-law); first T-varying-property case | IGLOO-vs-TC-kernel variable-$\rho$ mass rate + swelling (gated) + digitized TC2012 Fig. 11 (non-gating overlay) | **GREEN** (E-VAL-3, 2026-07-23; 25/25 <0.2%; found+fixed A20/A21) |
| [Pilch–Erdman breakup (validation)](../vv/e2e.md#pilch-erdman-e2e) | 25-drop Weber sweep ($We$ 20–1000) vs PE87's published $T^*(We)$ correlation (Fig. 7) + $d_\mathrm{stable}(We)$ closed form (B-VAL-2) | Windowed initial breakup rate vs the paper kernel (22 drops, tol 1 %) + plateau $d_\mathrm{stable}$ (5 drops, tol 12 %) | **GREEN** (B-VAL-1+2, 2026-07-22; gated the A16 sign fix, found+gated the A17 bp-wipe fix) |
| [TAB breakup (validation)](../vv/e2e.md#tab-e2e) | 25-drop radius-based We sweep across the ORA87 onset $We_r=6$; damped-oscillator $t_{bu}$ closed form | Onset (no-break) + first-breakup time vs ORA87 eq. 5 (Tier V) | **GREEN** (B-VAL-3, 2026-07-22; found + gated the A18 fix — event path was dead; 16/16 breaks, $t_{bu}$ sub-percent) |
| [ETAB breakup (validation)](../vv/e2e.md#etab-e2e) | 25-drop radius-based We sweep across onset $We_r=6$ and $We_\mathrm{trans}=80$; ETAB cascade product size + shared TAB $t_{bu}$ | Onset + $t_{bu}$ (ORA87) + child-size ratio vs Tanner eq. 6/8 (Tier V) | **GREEN** (B-VAL-5, 2026-07-22; 20/20 breaks, cascade size within 0.03 % over bag+strip branches) |
| [KHRT breakup (validation)](../vv/e2e.md#khrt-e2e) | 25-drop radius-based We sweep ($We_r$ 30–1000); continuous KH-stripping rate + RT/shed shatter vs Reitz-87 | Initial $\mathrm{d}d/\mathrm{d}t$ vs the Reitz-87 KH rate ($We_r\ge340$) + `khrt-e2e-rt` RT-shatter persistence gate | **GREEN** (B-VAL-6, 2026-07-22; KH stripping <0.03 %; found + gated the A19 fix 2026-07-23 — RT shatter now persists) |
| [Reitz–Diwakar breakup (validation)](../vv/e2e.md#reitz-diwakar-e2e) | 25-drop radius-based We sweep ($We_r$ 8–1000) across the bag$\to$stripping handoff vs RD 1987 (SAE 870598) | Initial $\mathrm{d}d/\mathrm{d}t$ per drop vs the RD closed form for its regime (bag/stripping) | **GREEN** (B-VAL-4, 2026-07-23; 25/25 <0.1 %; $C_s{=}20$ confirmed curve-fit in SAE 870598 — no bug) |
| Al combustion (`combustion/burn-box`) | Beckstead $d^n$ burn law (model 5) | Closed-form burn-time kernel | **GREEN** (M1, 2026-07-11) |
| DB injection | Assigned-position particle placement, `vInj` hand-off, Stokes velocity relaxation, domain exit | Placement fidelity + analytic Stokes relaxation | **GREEN** (bug B1 fixed 2026-07-08) |
| Coupled outputs (`infrastructure/coupled-body`) | euler + source + body-force accumulators together (model 1) | Body-force $v(x)$ + source totals vs closed forms (drag-reaction-only deposit, agreement ~4·10⁻⁷) | **GREEN** (2026-07-15) |
| 2Daxi + DB (`infrastructure/db-2daxi`) | Axisymmetric wedge, real MOSE gas field, DB injection, euler-only output | Behavioral (both particles integrate to the outlet, fields finite) | **GREEN** (B5/B6 fixed 2026-07-15) |
| Periodic BC (`infrastructure/periodic-y`) | Translational periodic pair (bcdef 201): transport, velocity-unchanged contract, relocation | Body-force $v(x)$ closed form across 2–3 wraps + $y(x)$ modulo $L_y$ (residual ~$5\cdot10^{-7}$ m) | **GREEN** (2026-07-16, first exercise of the path) |

---

## Verification family status

| Family (tag) | What it grounds | Pass criterion | Status |
|---|---|---|---|
| Drag (A) | Full 13-model catalog vs Shimada2006 (eqs. 12–24, the provenance) + NASA-SP-8039 Table I | CTest; `test_drag_lit` GREEN, `test_drag_probes` GATED GREEN (XD1–XD5; A1/A6/A11 fixed; Putnam/Wen-Yu plateaus source-faithful per the 2026-07-16 decision) | **GREEN** — no inferred constants remain |
| Temperature (B) | Full 6-model Nu catalog vs Shimada2006 (eqs. 45–50) + NASA-SP-8039 Table II (the true primary) | CTest; `test_heat_nu` GREEN, `test_heat_probes` GREEN (A7 closed source-faithful; A12 fixed) | **GREEN** — no inferred constants remain |
| Evaporation (C) | `evaporation()` function-level; d²-law, CEM, LK, TC; droplet energy F(7) | CTest; `test_evap_probes` GATED GREEN (XE1/XE2/XE3 — A3/A4/A9 fixed); A10 sensible-carry fixed 2026-07-15; d²-law/lk-neq/tc e2e GATED GREEN | **GREEN** |
| Gas reconstruction (E) | `ord2` interpolation, E1–E5 sub-tests | CTest GREEN (E6 deferred) | **GREEN** |
| INI pipeline (T9) | INI → module-variable round-trip, IP1–IP5 | CTest GREEN | **GREEN** |
| Breakup (TAB/ETAB/RD/PE/KHRT) | All five models vs the primary papers (PE87, TAB-872089, Tanner 97/98, Reitz-87, RD-860469, Beale-Reitz-99) | CTest; 7 gated tests GREEN (A2/A13/A14/A15 fixed; PE87+TAB fully closed; ETAB child $v_\perp = A\dot{x}$ implemented 2026-07-16, gated ET4) | **GREEN** (RD Cs=20 confirmed curve-fit in SAE 870598, 2026-07-23) |
| Combustion (M1) | Beckstead $d^n$ Al burn law | CTest; `test_combustion` + burn-box e2e GREEN | **GREEN** |

---

## Running the suite

```bash
# Full suite: builds build/verif/ (and bin/IGLOO for the e2e cases), runs ctest
./tests/test.sh all

# By kind or category
./tests/test.sh e2e          # end-to-end solver cases only
./tests/test.sh unit         # literature/unit tests only
./tests/test.sh standard     # a single model category
```

Step-by-step walkthrough (reading reports, diagnosing failures, xfail promotion):
[Using the Suite](guide.md). Layout and conventions: [Testing](../development/testing.md).

---

## Production bugs found

The verification effort has catalogued **twenty-one model-formula bugs (A1–A21)**, one
combustion-constant fix (M1), and six infrastructure bugs (B1–B6).  As of 2026-07-23
**every one is fixed, refuted, or closed as source-faithful**, each gated by a CTest
probe or e2e oracle: A7 (JAXA1 floor) and B2 closed *not-a-bug*; A16/A17 (`pilch-erdman-e2e`),
A18 (the event-breakup path had never executed; `tab-e2e`), and A19 (the KHRT
Rayleigh–Taylor child/shed event reverted because its droplet count was never written to
the ODE state; `khrt-e2e-rt`) all fixed.  The A16/A18/A19 finds — and the fixes they gate —
came out of the paper-reproduction validation campaign.
Three reproducers survive as transcription pins — `test_drag_probes`, `test_heat_probes`,
`test_evap_probes` — holding the fixed values (and the source-faithful JAXA1 limitation) so a
regression turns them RED.  They are **ordinary gates**, not expected failures: every bug they
probe is fixed, so each exits 0 and no `WILL_FAIL` property is set.
See [Literature tests](literature.md#production-bugs-found)
for the model-formula table, and `tests/VERIFICATION_MATRIX.md` (one row per CTest
entry, reference-source and oracle provenance) for the source of truth.
