# Verification matrix

One row per registered CTest entry (42 total: 22 e2e + 19 unit + `self_test`).
Companion to [`REFERENCES.md`](REFERENCES.md) (bibliography — all `[tags]`
below resolve there). Regenerate the reconciliation with
`ctest --test-dir build/verif -N`.

Columns:

- **Reference source** — the paper/standard the model is taken from.
- **Source-data generation** — how the *reference* the code is compared against
  is produced. Five kinds recur:
  - *closed form* — exact analytic solution of the model equation, computed from
    the test inputs only (fully independent of production output);
  - *analytic formula @ swept points* — the literature expression coded directly
    and evaluated at swept inputs, compared to the production function pointwise;
  - *independent integration* — a separate, easy-to-check numerical integrator
    (`support/verif_oracle.f90::rk4_ref`, or a local RK4) run on the model RHS;
  - *deep-converged fixed point* — an independent iterative solve to a tighter
    tolerance than production;
  - *run-conditioned* — the reference ODE is integrated **along a quantity read
    back from production output** (semi-independent; used where the rate law is
    driven by the measured particle temperature).
- **Gas solfile (e2e only)** — mesh topology + gas state of the field the full
  `IGLOO` executable reads. IGLOO reads geometry from the solfile, not `MESH/`.
- **Condensed-phase BC (e2e only)** — injection mode + face boundary conditions
  for the parcels.

---

## Shared box fixture (13 of 22 e2e cases)

Every e2e case except `db-2daxi`, `vie-plait`, `mhb98-water`,
`pilch-erdman-e2e`, `tab-e2e`, `etab-e2e`, `khrt-e2e`, and `reitz-diwakar-e2e` runs on the
**same** axis-aligned uniform-gas box emitted by
[`tools/make_box_case.py`](tools/make_box_case.py). The gas field is held *fixed*
across all of them; cases differ only by particle injection properties, enabled
models, body acceleration, and boundary conditions — never by the gas. Stated
once here, referenced as **"shared box"** below. The eight exceptions: `db-2daxi`
(spatially varying, real MOSE nozzle solution) and `vie-plait` (spatially varying,
analytic `−ε(y−1)`) are the two non-uniform-gas cases; `mhb98-water` keeps the box
*geometry* but swaps in a different **uniform** gas (water in air at T_G=298 K with
a near-static U=1.9×10⁻⁴ clock) to match the Miller-Harstad-Bellan Fig-2 conditions;
`pilch-erdman-e2e` likewise keeps the geometry but runs U=200 m/s water drops
(`tools/make_pe_case.py`) for the PE87 Weber sweep; `tab-e2e`/`etab-e2e` stretch the box to
0.6 m × 240 cells at U=50 m/s and `khrt-e2e`/`reitz-diwakar-e2e` run the 0.15 m box at U=200 m/s
(all via `make_pe_case.py`) for the TAB/ETAB onset and the KHRT / Reitz-Diwakar sweeps.

| Property | Value |
|---|---|
| Mesh topology | single Cartesian block, 60×5×5 cells (nodes 61×6×6), 3-D (no axisymmetric fold) |
| Domain | x∈[0,0.15], y∈[0,0.05], z∈[0,0.05] m; Δx=2.5 mm |
| Gas state | **spatially uniform**: ρ=1.2 kg/m³, **U**=(10,0,0) m/s, T=600 K, μ=1.8×10⁻⁵ Pa·s, k=0.026 W/m/K, γ=1.4, R=287 J/kg/K (P=101325 Pa, unused) |
| Injection | inlet-face (bcdef 401) on face 1 (x=0): one parcel stream per inlet cell, from `bc.txt` — number density κ_ρ·ρ (κ_ρ=0.34), speed v₀=κ_v·|U|, temperature T_p0=κ_t·T_g, monodisperse (Dirac) d=2r_p |
| Faces | face 2 (x=L_x) = outlet (parcel exits); faces 3–6 = wall. Pure +x motion only ever reaches the outlet |

Because the gas is uniform, the per-particle trajectory is an exact function of
the injection state and the enabled models — which is what makes a closed-form
oracle possible.

---

## End-to-end cases (full `IGLOO` executable)

| Case | Reference source | Source-data generation | Gas solfile | Condensed-phase BC |
|---|---|---|---|---|
| **drag-stokes** | `[CGW78]`,`[SN33]` (Stokes limit) | *closed form* — exact Stokes velocity relaxation x(u)=x_a−u_g·τ·ln[(u−u_g)/(u_a−u_g)]−τ(u−u_a), τ=ρ_p d²/18μ; precomputed to `reference/stokes.txt`, recomputed in `check.py` | shared box | inlet-face 401; κ_v=0.1 (90% axial slip), κ_t=1.0 (no thermal slip), d=11.89 µm, ρ_p=2950; outlet + walls |
| **temp-relax** | `[RM52]` | *closed form* — lumped-capacitance Nu=2 relaxation T(x)=T_g+(T_p0−T_g)e^{−x/L}, L=u_g·ρ_p·c_p·d²/(12k_g); `reference/nu2.txt` | shared box | inlet-face 401; κ_v=1.0 (Re=0, drag off), κ_t=0.5 (T_p0=300 K, ΔT=300), d=11.89 µm; outlet + walls |
| **body-force** | analytic (body-accel model) | *closed form* — transverse terminal drift v(x)=v_∞(1−e^{−x/L}), v_∞=g_y·τ, L=u_g·τ (linear Stokes, `toll=1e-20`); `reference/drift.txt` | shared box + `body-accel=(0,−200,0)` m/s² | inlet-face 401; κ_v=1.0, κ_t=1.0, d=11.89 µm; outlet + walls |
| **conv-nu** | `[RM52]` | *closed form* — Ranz–Marshall Nu=2+0.6 Re^{1/2}Pr^{1/3} at **constant** slip (parcel injected at body-force terminal velocity ⇒ Re, Nu fixed), T(x)=T_g+(T_a−T_g)e^{−(x−x_a)/L}, L=u_g·τ_T; `reference/nu_re.txt` | shared box + `body-accel=(0,−200,0)` | inlet-face 401 at terminal velocity: κ_v=1.00033, α=−0.02574 rad, κ_t=0.5, d=11.89 µm; outlet + walls |
| **vie-plait** | `[Vie15]` | *closed form* (derived from ODE + IC, not the paper's typeset Eq. 5.4 which has a sin sign typo; scipy-validated to 1e-13) — anchored damped oscillator y(x)=1+e^{−t/2τ}[Y_a cos ωt+((V_a+Y_a/2τ)/ω)sin ωt], t=(x−x_a)/u_g0, ω=√(ε/τ−1/4τ²); recomputed in `check.py` | **NON-UNIFORM** box (first such e2e): x∈[0,6], y∈[0,2], z thin; 120×50×5; U=u_g0=0.2 uniform, **V=−ε(y−1)** (ε=1, compressive, linear in y), W=0. **`gas-order=2` required** (2nd-order rebuild; order 1 = staircase) | **assigned-position (DB)** injection: 6 parcels at x=0.125, z=0.025, y∈{0.5,0.7,0.9,1.1,1.3,1.5} (cell centres), up=0.2/vp=0 (U_p0=u_g0, V_p0=0); τ=5 s via ρ_p=1620 (local `properties.dat`), d=1 mm ⇒ St=5; outlet + walls (never hit) |
| **d2law** | `[God53]`,`[Spa53]` (+`[RM52]`,`[Lef89]`) | *run-conditioned* — Godsave–Spalding kernel d(d²)/dt=−(8k_g/c_pg·ρ_p)·ln(1+B_T) integrated in `check.py` **along the measured T_p(x)** from production output (no static file) | shared box | inlet-face 401; `evaporation=d2-law`; κ_t=0.5, d=30 µm, ρ_p=8000; fuel L_v=2e5, M_v=100, T_boil=900, p_sat=1e4; outlet + walls |
| **d2law-line** | `[God53]`,`[Spa53]` | *closed form* (Tier V) — input-only line d²(x)=d0²−K0·(x−x_a)/u_g, K0=8k_g·ln(1+B_T0)/(cp_g·ρ_p), B_T0 at T_p0=κ_t·T_g=300 K; tolerance = theory budget (monotone K-drift bound from measured T_p span + E13.6 truncation + floor); measured T_p bounds the budget, never the reference | shared box | inlet-face 401; `evaporation=d2-law`; κ_v=1 (Re=0), κ_t=0.5, d=30 µm; ρ_p=8000 with **cp=125000 (100×)** via local ATLAS-GPB `properties.dat` (freezes T_p: drift 0.76 K, d² loss 39 %); fuel L_v=2e5; outlet + walls |
| **lk-neq** | `[MHB98]` | *run-conditioned* — LK-corrected CEM d²-ODE RK4-integrated along measured T_p; VLE-equilibrium curve as the regression the case must out-distance | shared box | inlet-face 401; `evaporation=CEM`+`interface=LK`, α_e=1.0; κ_t=0.5, d=20 µm, ρ_p=8000; outlet + walls |
| **tc-box** | `[TC2012]`,`[ATC24]` | *run-conditioned* — Tonini–Cossali Stefan–Fuchs d²-ODE RK4-integrated along measured T_p; CEM-Spalding as the silent regression | shared box | inlet-face 401; `evaporation=TC`+`interface=VLE`; κ_t=0.5, d=20 µm; fuel L_v=1e6, T_boil=600; outlet + walls |
| **burn-box** | `[Beck05]` | *mixed* — diameter d^n(x)=d₀^n−K_eff·(x−x_a)/u_g is a **closed form** (T_p-independent above ignition); temperature T_p(x) is an **independent RK4** energy balance (Nu=2 conduction + combustion release); `reference/beckstead.txt` | shared box | inlet-face 401; `combustion=Beckstead` (K-burn=4.5e-7, n=1.8, X-eff=0.5, T-ign=550, β=0.3, q-comb=1e6); κ_t=1.0, d=30 µm, ρ_p=8000; outlet + walls |
| **test_mhb98_decane** (unit) | `[MHB98]` Fig. 4 **M7** (decane, T_G=1000 K, T_d,0=315 K, D0=2.0 mm, Re_d,0=17; exp Wong & Lin 1992) | *validation* (E-VAL-2b) — THE model-discriminating case: Fig. 4b spans ~200 K across M1–M8 vs 0.47 K in Fig. 2b. GATES: **DC1** IGLOO `evaporation()` CEM+LK mdot vs an independently re-coded M7 chain (eqs 10, 15–18) at 25 sampled states, 1e-12 (measured **3.5e-16**) + **DC2** IGLOO `blowingFactor()` vs eq. 19 `f2=b/(e^b−1)`, 1e-12 (measured **9.4e-15**; f2 spans **0.40–0.88**, i.e. convective heat cut up to 60 %) + **DC3** the re-coded M7 at MHB98's OWN 0.552 vs three PIXEL-MEASURED Fig. 4 quantities: β 1.579 vs 1.60 (**1.3 %**), T_d(3.5 s) 393.7 vs 397 K (**3.3 K**), D²(4.0 s) 1.908 vs 1.885 mm² (**1.2 %**) | **FIXED SLIP** (suspended drop, Re ∝ d): free flight was REFUTED by measurement (β 1.31–1.44 for every Sc). All properties = MHB98 App. A at their T_R=420.19 K; **Le=1 deliberately, NOT their App. Γ_V** (which gives Le=4.40 and misses the plateau by +40 K — paper is self-inconsistent) | needs `blowing = LK`; trajectory RK4'd in-test over IGLOO's kernels — a unit test, not e2e, since fixed slip would need a no-drag flag and `t=x/U` is invalid with real slip. **Confound quantified, not hidden:** IGLOO's Ranz-Marshall carries 0.600, MHB98 eq. 6 prints 0.552 (8.7 %); DC1/DC2 use 0.600 (implementation), DC3 uses 0.552 (paper) |
| **mhb98-water** | `[MHB98]` Fig. 2 **M7** model line (WPD-digitized) | *validation* (paper-reproduction) — GATES: IGLOO-vs-LK-kernel RK4-integrated along measured T_p (tight) + the `K=(T_G−T_wb)·8k_g/(ρ_l L_v)` identity + **IGLOO β vs the digitized M7 slope β_M7≈6.35×10⁻³ (Eq.28), ±10%** (IGLOO=6.30e-3, rel 0.8%; guarded vs paper's 6×10⁻³) + **wet-bulb vs M7 plateau 282.4 K, ±0.5 K** (IGLOO 282.33, ΔT 0.10 K); overlay is absolute D²[mm²]. IGLOO runs CEM+LK = uniform-temp LK = **M7**. Eyeball `[RM52]` exp dropped (reads high; re-digitize follow-up) | box, but **water @ T_G=298 K**, U=1.9×10⁻⁴ (pure clock); **all properties = MHB98 Appendix A (Harpole 1981) at their T_R=293.968 K**: ρ=1.184727, μ=1.8735×10⁻⁵, k=0.0261731, c_p,g=989.45 (via GAM) | inlet-face 401; `evaporation=CEM`+`interface=LK`; single water droplet d=1.05 mm (the FIGURE's D0², not the caption's), κ_t=0.9463 (T_p0=282 K), κ_v=1 (Re=0), ρ_p=997, L_v=2.462478e6, dry air Y∞=0; outlet + walls |
| **tc-hexadecane** | `[TC2012]` Fig. 11 | *run-conditioned* (Tier V) + digitized overlay (non-gating) — variable-density TC mass rate: the oracle integrates droplet MASS with the TC2012 eq.9 Stefan-Fuchs rate (bisection) along measured Tp, reconstructs d^2=(6m/(pi rho_l(Tp)))^{2/3} (evaporation + thermal swelling), matches measured d^2 (25/25 <0.2%, tol 2%) + asserts swelling. First VARIABLE-property case: found+fixed A20 (tabs never allocated) + A21 (lookupTab OOB) | shared box but **water->n-hexadecane, T-dependent rho_l(T)** in `properties.dat` (Mv=226.45, Tboil=560, Lv=2.9e5); air at T_G=600 K | inlet-face 401; `evaporation=TC`+`interface=VLE`; d0=20 um, kt=0.5 (T_p0=300 K), kv=1 (Re=0); outlet + walls |
| **pilch-erdman-e2e** | `[PE87]` Fig. 7 / p. 748 | *validation* (paper-reproduction) — GATE is IGLOO's initial breakup rate per drop vs PE87's closed-form T\*(We) forward rate `dd/dt=(d_stable−d)/τ` at the injection We (t reconstructed as Σdx/ū); **diameter-based** We (PE convention); the reference applies the same leading-window LSQ to the kernel's RK4 samples (window-bias cancels). Gates: rate at We≥45 (22 drops, tol 1 %) + B-VAL-2 d_stable plateau vs the paper's initial-conditions closed form (5 drops, tol 12 %) | shared box geometry but **U=200 m/s**, ρ_g=1.2, T=300 K (`tools/make_pe_case.py`); σ=0.072, μ_p=1e-3 via `[IGLOO-Properties]` | inlet-face 401 on face 1: **25 inlet cells, one diameter each** (We sweep 20…1000 at slip=100), κ_ρ=0.34, κ_v=0.5 (drop SLOWER than gas — shock-tube analog), κ_t=1.0, Dirac; `breakup=Pilch-Erdman`, Cd=0.5, B=0.0758; ρ_p=1000 (ATLAS-GPB `phase.txt`/`properties.dat`); outlet + walls |
| **tab-e2e** | `[ORA87]` eq. 5 (SAE 872089) | *closed form* (Tier V) — onset (no-break `We_r<=5.75`) + first-breakup time vs the damped-oscillator crossing `y(t)=We_Cr(1−e^{−t/t_d}(cos ωt+sin ωt/(ω t_d)))`; **radius-based** We (TAB/ORA87 convention); measured bracket = rows around the first d drop, slack 5 % (t_bu inside bracket to sub-percent). Found + gates the A18 fix (event path was dead) | box geometry stretched: 0.6×0.05×0.05 m, 240×5×5, **U=50 m/s** (`make_pe_case.py --we-convention rad`) | inlet-face 401: 25 cells, one diameter each (`We_r ∈ {3…90}` at slip=25), κ_ρ=0.34, κ_v=0.5, κ_t=1.0, Dirac; `breakup=TAB` (Comega=8, Cmu=5, WeCrit=12 internal, method=2 n=3.5); water (ATLAS-GPB fixtures); outlet + walls |
| **etab-e2e** | `[Tan97]`,`[Tan98]` (SAE 970050/980808) | *closed form* (Tier V) — onset (no-break `We_r<=5.75`) + first-breakup time (ORA87 oscillator, shared with tab-e2e) + the ETAB deterministic cascade product size `d_child/d_parent = exp(-(Kbr/ω)acos(1-1/WeCr))` (bag `k1(AWe We^4+1)`, strip `k2 sqrt(We)`, WeTrans=80); **radius-based** We; reference at the We AT BREAKUP, tol 1% (worst 0.03%) | box geometry stretched: 0.6×0.05×0.05 m, 240×5×5, **U=50 m/s** (`make_pe_case.py --we-convention rad`) | inlet-face 401: 25 cells, one diameter each (`We_r ∈ {3…105}` at slip=25), κ_v=0.5; `breakup=ETAB` (k1=k2=0.2222, WeCrit=12 internal, WeTrans=80, Comega=8, Cmu=5); water (ATLAS-GPB fixtures); outlet + walls |
| **khrt-e2e** | `[Reitz87]` (Atomisation & Spray Tech. 1987) | *closed form* (Tier V/P) — the continuous KH-stripping rate `dd/dt=(dStable-d)/tauKH`, dStable=2·B0·λ_KH, τ_KH=3.726·B1·r/(λ_KH·Ω_KH) with the Reitz-87 Λ/Ω correlations; **radius-based** We; initial-rate (same-window LSQ) at We_r>=340, tol 2% (worst 2e-4). `khrt-e2e-rt` = RT-shatter persistence gate (found+gated the A19 fix: 8 drops shatter discontinuously to the fragment scale and persist, mass-consistent) | shared 0.15 m box but **U=200 m/s** (`make_pe_case.py --we-convention rad`) | inlet-face 401: 25 cells, one diameter each (`We_r ∈ {30…1000}` at slip=100), κ_v=0.5; `breakup=Reitz-KHRT` (B0=0.61, B1=20, Ctau=1, CRT=0.1, mShedLim=0.03, WeLimit=6); water (ATLAS-GPB fixtures); outlet + walls |
| **reitz-diwakar-e2e** | `[RD87]` (SAE 870598) | *closed form* (Tier P) — the RD breakup rate `dd/dt=(dStable-d)/τ`, branch per drop: **bag** (`We_r>6`, `We_r<=0.5√Re`) τ=π√(ρ_l·r³/2σ), dStable=12σ/(ρ_g u²); **stripping** (`We_r>0.5√Re`) τ=C·(r/u)√(ρ_l/ρ_g) C=20, dStable=σ²/(ρ_g u³ μ_g). **radius-based** We; initial-rate (same-window LSQ), tol 2% (worst 8e-4). All four constants + both stable sizes verified vs the [RD87] PDF (bag D=π Eq.7, stripping C=20 curve-fit Eqs 6/8) — no production change | shared 0.15 m box but **U=200 m/s** (`make_pe_case.py --we-convention rad`) | inlet-face 401: 25 cells, one diameter each (`We_r ∈ {8…1000}` at slip=100, handoff at We_r=20), κ_v=0.5; `breakup=Reitz-Diawakar` (WeBag=6, Cb=π, Cstrip=0.5, Cs=20); water (ATLAS-GPB fixtures); outlet + walls |
| **db-injection** | infrastructure (no paper; bug-B1 regression) | *behavioral* — no reference curve; asserts placement fidelity, `vInj` hand-off (u starts at u_p=1), Stokes relaxation, and domain exits | shared box | **assigned-position (DB)** injection: 5 parcels at explicit x=0.01, y=0.005…0.045, z=0.025, u_p=1.0 m/s, d=11.89 µm, ṁ=1e-4, T=300 K; outlet + walls |
| **coupled-body** | analytic (euler+source+body) | *mixed* — body-force v(x) closed form (as body-force) **plus** source-deposit totals vs closed forms: the deposit must be the drag reaction only (body-gained momentum/energy stripped) | shared box + `body-accel=(0,−200,0)`; euler+source accumulators ON (default `out-file`) | inlet-face 401; κ_v=1.0, κ_t=1.0, d=11.89 µm; outlet + walls |
| **db-2daxi** | infrastructure (no paper; promoted legacy `assigned-pos`) | *behavioral*, md5-free — both parcels integrate, exit the nozzle outlet, all fields finite; no analytic oracle (real flow field) | **NON-UNIFORM real MOSE nozzle**: `common/solfile_mose.tec`, single block I=201×J=181×K=2 (2-D axisymmetric wedge, one cell thick), ~6×10⁵ distinct field values; spatially varying ρ/U/T (~300–3600 K). Drag=Morsi-Alexander, heat=Kavanau-Drake | **assigned-position (DB)** injection: 2 parcels x=[−0.49,0.064], y=[0.55,0.66], d=[1e-4,2e-5], ṁ=[0.1,1.0]; face1 inlet, face2 outlet, **face3 symmetry** (axisymmetric fold), face4 wall; euler-only output |
| **periodic-y** | infrastructure (translational periodic path) | *closed form* — body-force closed form folded **modulo L_y**; velocity unchanged across each wrap (transport = exact ±L_y translation), residual ~5e-7 m | shared box + `body-accel=(0,−8000,0)` (drives 2–3 y-wraps) | inlet-face 401; κ_v=1.0, κ_t=1.0, d=11.89 µm; **faces 3/4 translational-periodic (bcdef 201)**; face2 outlet; faces 5/6 wall |

---

## Unit families (compiled Fortran, function-level)

These link the production static library `IGLOOL` and call the production
closure **directly** — no time loop, no mesh, no gas file, no boundary
conditions. **Gas-solfile and condensed-phase-BC columns are therefore N/A** by
construction. `test_drag`/`test_temperature` are the exception: they drive the
production `rhsStandard` through `verif_driver` + OSlo's SDIRK4, still with no
box/BC. Sub-test IDs in parentheses.

| Test | Reference source | Source-data generation |
|---|---|---|
| **test_drag** | `[CGW78]`,`[SN33]` | *closed form* exp relaxation (`exp_relax`, A1 Stokes); SDIRK4 tableau order (AO); *independent integration* `verif_oracle::rk4_ref` of Schiller–Naumann RHS (A3); zero-slip straight line (A4). Production path: `Lib_Drag::drag` via `rhsStandard`+SDIRK4 |
| **test_drag_lit** | `[CGW78]`,`[SN33]` | *analytic formula @ swept Re* — production `Lib_Drag::drag` vs coded literature expressions: Stokes 24/Re (DL1), Chang≡Clift-Gauvin identity (DL2), Schiller-Naumann≡Wen-Yu low branch (DL3) |
| **test_temperature** | `[RM52]` | *closed form* exp relaxation (`exp_relax`, Nu=2, B1); SDIRK4 order (BO); *independent integration* `rk4_ref` of the coupled velocity+temperature system under Ranz-Marshall (B2). Production `Lib_Heat::heat` via `rhsStandard` |
| **test_heat_nu** | `[RM52]` | *analytic formula @ swept Re×Pr* — production `Lib_Heat::heat` Nu vs Ranz-Marshall 2+0.6 Re^{1/2}Pr^{1/3} table (HN1), Nu=2 conduction floor (HN2), Kavanau-Drake collapse (HN3) |
| **test_evaporation** | `[RM52]`,`[Lef89]`, Abramzon-Sirignano 1989 | *analytic/algebraic* — production `Lib_Evaporation::evaporation` CEM/CEM-B vs independent Spalding chain psat→X_s→Y_s→B_M (EV1/EV3); convective ratio vs Ranz-Marshall Sh (EV2) |
| **test_evap_probes** | — | *bug-transcription pin* — probes XE1/XE2/XE3 for bugs A3/A4/A9 (green while a fixed bug stays fixed; RED signals a regression). **Ordinary gate, not xfail:** written as an expected-failure probe, but A3/A4/A9 are all fixed, so it exits 0 and no `WILL_FAIL` property is set |
| **test_interface_lk** | `[MHB98]` | *deep-converged fixed point* — production LK interface vs independent deep-Picard (1e-15) oracle (LK2/LK5); exact 1/(p·d) chain inversion (LK4); VLE limit (LK3); 27-pt T_p×d×ρ envelope |
| **test_tc_analytic** | `[TC2012]`,`[ATC24]` | *deep-converged fixed point + closed form* — production TC vs independent bisection root of the G(m) residual (TC1, 27-pt grid); exact isothermal Stefan-Fuchs limit m=rhs0 (TC3); CEM-collapse identity (TC4); bracket theorem (TC6) |
| **test_combustion** | `[Beck05]` | *closed form + independent integration* — closed-form m(t)=(ρπ/6)(d₀^n−K_eff·t)^{3/n} with complex-step dm/dt vs `becksteadRate` (CB1); local RK4 of the production rate vs the closed form (CB2); bitwise ignition gate (CB3); X_eff-exponent discriminator (CB4) |
| **test_breakup_khrt** | `[Reitz87]` (KH), `[BR99]` (RT) | *analytic formula @ swept We* — production `Lib_Breakup::breakupOde` (pure-KH, acc=0) vs the coded Reitz-1987 KH linear-stability chain to 1e-12 (KH1 sub-critical, KH2/KH3) |
| **test_breakup_tab** | `[ORA87]` | *closed form + independent integration* — production `breakupEvent` oscillator vs the analytic damped-driven-oscillator solution (TAB1); *independent* `rk4_ref` of the raw ORA87 ODE (TAB2); rest-drop critical Weber (TAB3); breakup time via bisection (TAB4) |
| **test_tab_moments** | `[ORA87]` | *stochastic moments* — child-size E[r], E[r²] over 20 000 sampled breakups vs a truncated Rosin-Rammler moment oracle (CLT-band tolerance) |
| **test_breakup_rd** | `[RD86]` (+`[RD87]` lineage) | *analytic formula @ swept We_r* — production `breakupOde` vs the OpenFOAM-form bag/stripping oracle branches and the regime handoff (RD1–RD4) |
| **test_breakup_etab** | `[Tan97]`,`[Tan98]` | *analytic formula @ swept We* — production `breakupEvent` product-size ratio ln(r_new/r) vs the Tanner closed-form K_br branches, checking continuity across We_trans (ET1–ET4) |
| **test_breakup_pe** | `[PE87]` | *analytic formula @ swept We* — production `breakupOde` rate vs the PE87 breakup-time (TBT) table oracle, incl. Oh-override branch and the d<d_stable stability guard (PE1–PE4) |
| **test_gas_reconstruction** | internal (trilinear hex interpolation) | *analytic field + convergence order* — production `Lib_Equations::interp2ndOrder` vs analytic constant/multilinear fields (E1/E2, machine eps); observed order p≈2 under uniform refinement (E3); Newton forward/inverse hex round-trip (E4); partition-of-unity / node recovery (E5) |
| **test_ini_pipeline** | internal (config contract) | *config round-trip* — `IGLOO_IO_INI::read_IGLOO_input` parsed globals vs known `input.ini` values (TAB defaults, model selectors, body force, tolerances, RNG seed; IP1–IP5). No physics/curve |
| **test_heat_probes** | `[Shim06]`,`[SP8039]` | *bug-transcription pin* — probe XH1, originally for A7 (JAXA1 Nu floor), now pinning the source-faithful transcription (documented limitation, no +2 term). **Ordinary gate, not xfail:** A7 closed 2026-07-15 as source-faithful, so it exits 0 and no `WILL_FAIL` property is set |
| **test_drag_probes** | `[CGW78]`,`[Shim06]` | *bug-transcription pin* — probes XD1–XD5 for bugs A1/A5/A6 + the Wen-Yu C-flag correction. **Ordinary gate, not xfail:** all five now report `[FIXED]`, so it exits 0 and no `WILL_FAIL` property is set |
| **self_test** | — (harness) | *meta-test* — exercises the support library itself (`verif_norms`, `verif_oracle` `rk4_ref`/`exp_relax`, `assert_lt`); no production code |

---

## Reconciliation

22 e2e + 19 unit + 1 `self_test` = **42 CTest entries** (incl. `khrt-e2e-rt`, the KHRT RT-shatter persistence gate). Unit families that emit
a production-vs-reference overlay (`verif_dump` → `tools/plot_curves.py` →
`docs/vv/images/unit-*.svg`): drag, heat, evaporation-CEM, LK, TC, combustion,
the five breakup models, and interp — 12 figures. `ini_pipeline` (config
contract) and the three `*_probes` (scalar bug pins) emit no curve by design.
