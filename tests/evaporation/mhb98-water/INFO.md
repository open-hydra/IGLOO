# INFO — evaporation/mhb98-water (E-VAL-2: paper-reproduction validation, e2e)

## Reference paper
- **Title:** Evaluation of equilibrium and non-equilibrium evaporation models for
  many-droplet gas–liquid flow simulations
- **Author(s):** R. S. Miller, K. Harstad, J. Bellan
- **Year:** 1998 — *Int. J. Multiphase Flow*, 24(6), pp. 1025–1055
- **DOI:** 10.1016/S0301-9322(98)00028-7 · **Tag:** `[MHB98]`
- **Experimental data:** Ranz & Marshall (1952b), `[RM52]` — the points MHB98 plot.
- **Reproduced case:** **Fig. 2** — single water droplet, `D₀=1.1 mm`, `T_d,0=282 K`,
  `T_G=298 K`, `Re_d=0` (quiescent, 1 atm). This fulfils the "digitized overlay
  pending" note in [../lk-neq/INFO.md](../lk-neq/INFO.md).

## Validation vs verification (why this case exists)
The rest of the suite is *verification* (IGLOO vs re-coded model equations). This is
*validation*: the full integrated IGLOO trajectory is compared, under MHB98's exact
Fig-2 conditions, against **MHB98's own model result** (their reported `β≈6×10⁻³`) and
the re-coded LK kernel. Per the validation-plan gating rule, both the **kernel** and the
**paper model β** are gated (tight, deterministic); the **experimental points are a
NON-gating overlay** drawn by `verify.py`.

## Case construction
`tools/make_box_case.py` variant: uniform air at `T_G=298 K`; **`U=1.9e-4 m/s` is a pure
clock** (`κ_v=1` ⇒ zero slip ⇒ `Re=0` ⇒ `Nu=Sh=2`, the quiescent limit), so `t = x/U`.
`evaporation=CEM`+`interface=LK`; dry air `Yinf=0`.

**All properties are MHB98's OWN Appendix-A correlations (2026-07-29).** The paper publishes
every correlation *and* its source — air and water both **Harpole (1981)** — and states that
properties are frozen **once** at `T_R = T_WB` from their eq. (27), which for water at
`T_G=298 K` gives **`T_R = 293.9676 K`**. Evaluated there:

| | MHB98 App. A @ `T_R` | (was) | note |
|---|---|---|---|
| `k_g`   | **0.0261731** | 0.026    | air `λ_C(T)` |
| `μ_g`   | **1.8735e-5** | 1.85e-5  | air `μ_C(T)` |
| `c_p,g` | **989.45**    | 1004.5   | `= Pr_C·λ_C/μ_C`; carried in `GAM=1.4085719` because IGLOO builds `cp_g=γR/(γ−1)` — γ has no other effect at zero slip |
| `ρ_g`   | **1.184727**  | 1.184    | chosen so `p=ρRT_G` is exactly 1 atm = their `P_G`; under `Le=1` the flux uses `ρ_g D_v = k_g/c_p,g`, so ρ cancels there |
| `L_v`   | **2.462478e6**| 2.45e6   | water `2.257e6+2595(373.15−T)` |
| `ρ_l`   | **997**       | 1000     | `properties.dat` |
| `c_l`   | **4184**      | 4182     | |
| `c_p,v` | **2366.95**   | 1900     | INERT on the CEM path (only ASM/TC read `cpv`) |

Two IGLOO choices turn out to be the **paper's own**, not stand-ins:
- **`Le=1`** — water has no tabulated `Γ_V` or `Sc` in Appendix A, and the Appendix's stated
  default for exactly that case is *"the Lewis number is equal to unity"*.
- **`p_sat` by single-anchor Clausius–Clapeyron at `T_B` with constant `L_V`** — that is their
  eq. (10) verbatim.

**`D0` follows the FIGURE, not the caption (2026-07-28).** MHB98 Fig. 2a is internally
inconsistent: the caption says `D0=1.1 mm` (⇒ `D0²=1.21 mm²`) but its own plotted curve
starts at **`D0²≈1.1027 mm²`** — and the panel's y-axis tops at 1.2, so 1.21 could not even
be drawn. The plot is self-consistent (1.1027 ÷ 1.347e-3 ⇒ dies at ~818 s, matching the
drawn ~800 s). Since the overlay compares against *that curve*, the case uses the plot's
value: `bc.txt rp=5.25048e-4` ⇒ **`D0=1.050095 mm`**. This is a pure bookkeeping choice —
`β` (Eq. 28) is `D0`-independent and was **numerically unchanged** by this choice — but it
collapses the IGLOO-vs-M7 **lifetime** gap from **+10.9 % to +0.8 %** (825.3 s vs 818.5 s),
which is exactly the residual `β` difference.

**Box lengthened to `Lx=0.17 m` (68 cells, `dx` unchanged).** At `Lx=0.15` the drop left the
domain still at 15 % of `d0²`, truncating the curves at `t/t_life≈0.86`. The longer box lets
it evaporate to **1.2 % of `d0²`** (`t/t_life=0.989`), so the tight kernel gate now validates
**98.9 % of the droplet lifetime** instead of ~85 %. Doing this exposed production bug
**A22** (full evaporation poisoned the Eulerian source with NaN) — found here, fixed, and
documented with the A22 fix.

## What `check.py` gates
Three gates, tightest first:
1. **IGLOO-vs-kernel (tight).** IGLOO `d²(x)` vs the LK-corrected CEM rate law (written
   from the literature, not the code) RK4-integrated along the MEASURED `Tp(x)`; plus the
   coupled-physics identity `K = (T_G − T_wb)·8k_g/(ρ_l L_v)` (Re=0) and the mass-source
   telescoping audit. At `D₀=1.1 mm` the LK correction is small (`L_K/d≈2e-4`) but the
   kernel carries it — so it is *not* a VLE-discrimination case (unlike `lk-neq`).
2. **Two-phase MASS CONSERVATION (A22 burnout deposit).** Every droplet here is consumed
   *inside* the domain, so all injected droplet mass must reach the gas:
   `sum(wdot) == sum(mdot_inj)`, tol 1 % (measured closure **3.1e-6**). This is the direct
   assertion that the burnout remnant is handed to the Eulerian source rather than vanishing —
   without that deposit a coupled run leaks one remnant per parcel. (Only meaningful because the
   drops fully evaporate in-domain; a case whose droplets exit alive legitimately carries mass out.)
3. **Tier-P paper reproduction (the E-VAL-2 point) — now built on the DIGITIZED M7 model
   line.** IGLOO runs CEM+LK = the uniform-temperature Langmuir-Knudsen model = MHB98's
   **M7** (M8 is the finite-conductivity variant; superimposed on M7 in Fig. 2, so only M7
   was digitized). Two sub-gates:
   - (a) IGLOO's integrated `β = (ρ_d Pr/8μ_G)·|dD²/dt|` (Eq. 28) vs **`β_M7`** = the LSQ
     slope of the WebPlotDigitizer M7 curve (`reference/mhb98_fig2_M7_d2.csv`), **±10%**.
     `β_M7 ≈ 6.35e-3` is guarded to lie within ±20% of the paper's stated `β ≈ 6×10⁻³`
     (p. 1038, all M1–M8 at `Re_d=0`; the fallback if the file is absent). IGLOO `β = 6.30e-3`
     → **0.8 %**. (Both numbers shifted from the pre-2026-07-29 6.44/6.51e-3 because Eq. 28's
     coefficient `ρ_d Pr_G/8μ_G` is now built from MHB98's own `ρ_l`, `Pr_G`, `μ_G`; the gate is
     a slope ratio, in which the coefficient cancels. Bonus: in the paper's own constants the
     digitized slope reads 6.35e-3, i.e. *closer* to their stated ≈6×10⁻³.)
   - (b) IGLOO wet-bulb vs the **M7 temperature plateau** 282.43 K (`..._Td.csv`), **±0.5 K**.
     IGLOO 282.33 K → **|ΔT| = 0.10 K**.
   This upgrades the earlier analytic-`β` proxy to a genuine IGLOO-vs-digitized-model-**curve**
   gate in both panels.

## Comparison plot (`verify.py` → `OUTPUT/mhb98-water.svg`, non-gating)
Two panels vs the **normalized** `t/t_life` (t_life = D0²/K from the LSQ slope, extrapolated to
D²→0 rather than taken as the last sample, which lands at d²/d0²≈0.012): (a) `d²/d0²` — IGLOO, its LK kernel (gated oracle), and the
digitized **M7 model line**; (b) `T_d` — IGLOO wet-bulb history and the digitized M7 line.
Normalized (not absolute) because the paper's caption `D0=1.1 mm`⇒`D0²=1.21` disagrees with its
own plot (M7 starts at D0²≈1.10; axis tops at 1.2); since β (Eq. 28) is D0-independent and
matches to <1 %, the two normalized d²-law lines **coincide** (an honest presentation — the
coincidence is partly tautological, so the β annotation + panel (b) carry the real evidence).
The old eyeball Ranz–Marshall overlay is **no longer plotted** (see below).

**Lifetime now matches too.** With the figure's `D0²` and the lengthened box, `t_life = D0²/K`
is **825.3 s (IGLOO) vs 818.5 s (M7) — +0.8 %**, i.e. exactly the residual `β` difference; both
curves span the full `0→1` normalized lifetime. (Historical note: with the caption's `D0²=1.21`
the gap was **+10.9 %** and IGLOO's curve stopped at `t/t_life=0.86` because the drop left the
short box still at 15 % of `d0²` — that was **truncation, not a shorter lifetime**; IGLOO's
lifetime was in fact *longer*. Both artefacts are gone.)

**Panel (b) — why IGLOO's plateau sits 0.10 K below the digitized M7.** IGLOO's `T_p=282.33 K`
is the **equilibrium wet-bulb of the MHB98 parameterisation as IGLOO implements it** — an
independent conduction=latent balance with the same inputs gives the same value, i.e. IGLOO
solves its own model correctly (the K↔wet-bulb identity gate confirms it, 25/25). MHB98 never
state a numerical wet-bulb for Fig. 2 (the text only says `T_d,0=282 K` is "approximately equal
to the predicted wet bulb"). At `D0=1.1 mm` the LK correction is vacuous (`L_K/d≈2×10⁻⁴`), so
IGLOO correctly behaves as the classical/equilibrium model (non-equilibrium only matters
`<50 µm` — MHB98's own point). The overlay is against **M7**, which sits at the **top** of
MHB98's 8-model cluster, and IGLOO's value lies *inside* that band.

*The gap is REAL, not a digitization artifact:* the Td digitization reproduces the exactly-known
initial condition `T_d,0=282 K` to **+0.003 K**, with 0.013 K tail scatter; and a direct pixel
measurement of Fig. 2b (2026-07-29, 400 dpi render, axis calibrated on the 280/282/284/286 major
ticks ⇒ **100.7 px/K**) puts the **M7/M8 marker symbols at 282.44 ± 0.03 K**. The digitized 282.43
is a faithful read.

***…but the figure cannot resolve 0.1 K between models.*** The same measurement shows the eight
plotted curves occupy **282.08 → 282.55 K (≈0.47 K)**, not the ~282.2–282.43 stated here before:
- M7/M8 (markers) at the top, **282.44**;
- the rest of the M3–M8 stack in **282.28–282.50**, individual strokes **0.04 K** thick and the
  M7/M8 symbols **0.09–0.19 K** tall;
- two low outliers, **282.16** (solid) and **282.10** (dashed) — consistent with **M1 and M2**,
  the only two models that use the *time-dependent 1/3-rule mixture* properties (Table 1 footnote)
  rather than the pure far-field gas.

IGLOO's 282.33 K therefore lands **on one of the plotted curves**, ~0.10 K under M7 — i.e. about
*one* intra-cluster line spacing. Meanwhile the paper's own text for this figure says the models
are indistinguishable here ("all of the mathematical expressions for `f₂` approach unity; while
all expressions for `H_M` approach `B_M,eq`"). A 0.47 K drawn spread under a claim of identity is
the honest statement of Fig. 2b's resolving power: **0.10 K carries no model information.**

*What the gap is NOT (settled 2026-07-29, once Appendix A was read).* Sensitivity: a **1 % change
in any of k_g, D_v, L_v or p_sat moves the wet-bulb ≈0.076 K**, so 0.10 K ≈ 1.3 % of one property.
Every property was then set to **MHB98's own Appendix-A values at their `T_R`** (table above),
which moved the wet-bulb 282.30 → 282.33 K and β from 1.1 % → 0.8 %. So the residual is **not**:
- **`Le`** — MHB98 also run water at `Le=1` (Appendix A gives no `Γ_V`/`Sc` for water and states
  `Le=1` as the default in that case). The old note that IGLOO's `Le=1` is "13 % below real
  H₂O–air" describes a deviation from *reality* that the paper **shares**; it cannot open an
  IGLOO-vs-M7 gap. Contrast tc-hexadecane, where `Le` genuinely *was* the lever.
- **the `p_sat` closure** — their eq. (10) is the same single-anchor CC at `T_B` with constant
  `L_V` (both ~10 % below Magnus at 282 K).
- **`k_g`, `μ_g`, `c_p,g`, `L_v`, `ρ_l`** — now the paper's own numbers, worth +0.025 K in total.
- **`f₂`** — M7's eq. (3) carries the evaporative heat-transfer correction `f₂ = β/(e^β−1)`
  (eq. 19) which IGLOO's CEM path lacks (`f₂≡1`), but here `β≈6.3e-3 ⇒ f₂=0.9969`: 0.31 %, and
  its sign *lowers* the wet-bulb, i.e. adding it moves IGLOO ~0.02 K **further** from M7. (It is
  a first-order term only at high evaporation rate — see the Fig. 4 decane scope note.)

*Honest limit (revised 2026-07-29 after the pixel measurement):* ~0.10 K remains unattributed, and
**it is not decidable against Fig. 2b** — it is smaller than the figure's own model-to-model
spacing (see above). In rate terms 0.10 K ≈ **1.3 %** of the group `ρ_g D_v L_v/k_g`, i.e. the
size of a single Table-1 *form* difference at this condition: `B_M,eq` vs `ln(1+B_M,eq)` is 0.3 %,
`(Y_s−Y_G)` vs `B_M,eq` is 0.6 %, `f₂=(1+B_T)⁻¹` (M4/M6) is 1.5 %. IGLOO's own known form
difference (`f₂≡1` instead of eq. 19) is 0.3 % and the wrong way. So the residual is a
model-**form**/bookkeeping difference of the same order as the spread MHB98 themselves call
negligible — **not a fixable IGLOO error**, and well inside the ±0.5 K gate.

*How it WOULD be decided:* only on a case where the models actually separate — **Fig. 4 decane**
(`T_G=1000 K`, `Re_d=17`, E-VAL-2b), where `B_T=O(1)` so `f₂` is genuinely `≠1`, the LK
non-equilibrium term is no longer vacuous, and M1/M2/M5 visibly over-predict while only M7/M8
track the experiment. At Fig. 2's condition every discriminating term is ~1 %.

## Result & finding
Wet-bulb **T_p = 282.33 K vs the M7 plateau 282.43 K (|ΔT| 0.10 K)**; D²-law reproduced.
**IGLOO's integrated β = 6.30×10⁻³ matches the digitized M7 line β_M7 = 6.35×10⁻³ to 0.8%**
(K_M7 = 1.347×10⁻³ mm²/s) — IGLOO lands *on* MHB98's M7 curve, not just our re-coded kernel.
This is a genuine IGLOO-vs-paper-**model-curve** reproduction (both panels), tighter than the
earlier analytic-β proxy (which used the paper's rounded 6×10⁻³, giving ~7%).

**Caveat on the earlier "~13% below EXPERIMENT" (E-VAL-2, 2026-07-23).** That gap rested on
the *eyeball* Ranz–Marshall CSV (β_exp≈7.4e-3, start 1.21). The WebPlotDigitizer M7 line
extrapolates to D0²≈1.10 and in the figure the exp circles sit essentially **on M7**
(β≈6.5e-3), so the "13% below exp" is likely **largely an eyeball-digitization artifact**, not
a real model-vs-exp offset. This does **not** overturn the committed finding on eyeball
evidence alone — confirming it needs a proper exp re-digitization (follow-up). The
`k_g`-film-rule hypothesis stays refuted either way. Not a bug.

**Scope note (out of this case):** Fig. 2 (water, Re=0) is a *weak* discriminator — all
eight MHB98 models coincide. The strong model discriminator is **Fig. 4 decane**
(`T_G=1000 K`, `Re_d=17`), where M1/M2/M5 visibly over-predict and only M7/M8 (non-eq)
track the experiment. That is a new *convective* case, beyond this water retrofit; see the
campaign follow-up.

## Reference data
- `reference/mhb98_fig2_M7_d2.csv`, `..._Td.csv` — **WebPlotDigitizer** exports of the M7
  model line (Fig. 2a/2b). The **gated** references: β_M7 (LSQ slope) and the T_d plateau.
- `reference/mhb98_fig2.csv` — Ranz–Marshall **experiment**, eyeball-digitized (±2–3%). Reads
  high (β_exp≈7.4e-3 vs the figure's ≈6.5e-3); **no longer plotted or gated** — a proper
  re-digitization is a follow-up. Full provenance + the M7-vs-eyeball discrepancy in
  `reference/PROVENANCE.md`.
