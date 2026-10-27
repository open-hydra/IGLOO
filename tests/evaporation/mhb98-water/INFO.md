# INFO — evaporation/mhb98-water (paper-reproduction validation, e2e)

## Reference paper
- **Title:** Evaluation of equilibrium and non-equilibrium evaporation models for
  many-droplet gas–liquid flow simulations
- **Author(s):** R. S. Miller, K. Harstad, J. Bellan
- **Year:** 1998 — *Int. J. Multiphase Flow*, 24(6), pp. 1025–1055
- **DOI:** 10.1016/S0301-9322(98)00028-7 · **Tag:** `[MHB98]`
- **Experimental data:** Ranz & Marshall (1952b), `[RM52]` — the points MHB98 plot.
- **Reproduced case:** **Fig. 2** — single water droplet, `D₀=1.1 mm`, `T_d,0=282 K`,
  `T_G=298 K`, `Re_d=0` (quiescent, 1 atm). It is the digitized MHB98 overlay that
  [../lk-neq/INFO.md](../lk-neq/INFO.md) points to.

## Validation vs verification (why this case exists)
The rest of the suite is *verification* (IGLOO vs re-coded model equations). This is
*validation*: the full integrated IGLOO trajectory is compared, under MHB98's exact
Fig-2 conditions, against **MHB98's own model result** (their reported `β≈6×10⁻³`) and
the re-coded LK kernel. Both the **kernel** and the
**paper model β** are gated (tight, deterministic); the **experimental points are not
gated**: `verify.py` draws them, re-digitized from the figure (see Reference data).

## Case construction
`tools/make_box_case.py` variant: uniform air at `T_G=298 K`; **`U=1.9e-4 m/s` is a pure
clock** (`κ_v=1` ⇒ zero slip ⇒ `Re=0` ⇒ `Nu=Sh=2`, the quiescent limit), so `t = x/U`.
`evaporation=CEM`+`interface=LK`; dry air `Yinf=0`.

**All properties are MHB98's OWN Appendix-A correlations.** The paper publishes
every correlation *and* its source — air and water both **Harpole (1981)** — and states that
properties are frozen **once** at `T_R = T_WB` from their eq. (27), which for water at
`T_G=298 K` gives **`T_R = 293.9676 K`**. Evaluated there:

| | MHB98 App. A @ `T_R` | generic value | note |
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

**`D0` follows the FIGURE, not the caption.** MHB98 Fig. 2a is internally
inconsistent: the caption says `D0=1.1 mm` (⇒ `D0²=1.21 mm²`) but its own plotted curve
starts at **`D0²≈1.1027 mm²`** — and the panel's y-axis tops at 1.2, so 1.21 could not even
be drawn. The plot is self-consistent (1.1027 ÷ 1.347e-3 ⇒ dies at ~818 s, matching the
drawn ~800 s). Since the overlay compares against *that curve*, the case uses the plot's
value: `bc.txt rp=5.25048e-4` ⇒ **`D0=1.050095 mm`**. This is a pure bookkeeping choice —
`β` (Eq. 28) is `D0`-independent and **numerically unchanged** by this choice — but it
brings the IGLOO-vs-M7 **lifetime** gap from **+10.9 %** (the caption's `D0`) **to +0.8 %**
(825.3 s vs 818.5 s), which is exactly the residual `β` difference.

**Box length `Lx=0.17 m` (68 cells, `dx = 2.5 mm`).** At `Lx=0.15` the drop leaves the
domain still at 15 % of `d0²`, truncating the curves at `t/t_life≈0.86`. The 0.17 m box lets
it evaporate to **1.2 % of `d0²`** (`t/t_life=0.989`), so the tight kernel gate validates
**98.9 % of the droplet lifetime** instead of ~85 %. Every drop therefore burns out inside
the box: a solver whose full evaporation poisons the Eulerian source with NaN fails the mass
telescoping and gate 2.

## What `check.py` gates
Three gates, tightest first:
1. **IGLOO-vs-kernel (tight).** IGLOO `d²(x)` vs the LK-corrected CEM rate law (written
   from the literature, not the code) RK4-integrated along the MEASURED `Tp(x)`; plus the
   coupled-physics identity `K = f₂·(T_G − T_wb)·8k_g/(ρ_l L_v)` (Re=0; MHB98 eqs. 3, 17–19:
   `f₂ = β/(e^β − 1)`, `β = ρ_l Pr K/(8μ_G)` with `blowing = LK`, `f₂ = 1` here), K and T_wb from the
   last half of the rows, within **5e-4** (measured 4.4e-5; the tail form's own floor, a coupled RK4 of
   the droplet, is 1.8e-5; a missing or doubled f₂ moves it by 3.1e-3, and with f₂ applied while
   blowing is off the identity fails by 3.1e-3), and the mass-source telescoping audit. At `D₀=1.1 mm` the LK correction is small (`L_K/d≈2e-4`) but the
   kernel carries it — so it is *not* a VLE-discrimination case (unlike `lk-neq`).
2. **Two-phase MASS CONSERVATION (the burnout deposit).** Every droplet here is consumed
   *inside* the domain, so all injected droplet mass must reach the gas:
   `sum(wdot) == sum(mdot_inj)`, tol 1 % (measured closure **3.1e-6**). This is the direct
   assertion that the burnout remnant is handed to the Eulerian source rather than vanishing —
   without that deposit a coupled run leaks one remnant per parcel. (Only meaningful because the
   drops fully evaporate in-domain; a case whose droplets exit alive legitimately carries mass out.)
3. **Tier-P paper reproduction — on the DIGITIZED M7 model
   line.** IGLOO runs CEM+LK = the uniform-temperature Langmuir-Knudsen model = MHB98's
   **M7** (M8 is the finite-conductivity variant; superimposed on M7 in Fig. 2, so only M7
   is digitized). Two sub-gates:
   - (a) IGLOO's integrated `β = (ρ_d Pr/8μ_G)·|dD²/dt|` (Eq. 28) vs **`β_M7`** = the LSQ
     slope of the WebPlotDigitizer M7 curve (`reference/mhb98_fig2_M7_d2.csv`), **±10%**.
     `β_M7 ≈ 6.35e-3` is guarded to lie within ±20% of the paper's stated `β ≈ 6×10⁻³`
     (p. 1038, all M1–M8 at `Re_d=0`; the fallback if the file is absent). IGLOO `β = 6.30e-3`
     → **0.8 %**. (Eq. 28's coefficient `ρ_d Pr_G/8μ_G` is built from MHB98's own `ρ_l`,
     `Pr_G`, `μ_G`; with the generic constants the two numbers read 6.44/6.51e-3. The gate is
     a slope ratio, in which the coefficient cancels; in the paper's own constants the
     digitized slope reads 6.35e-3, i.e. *closer* to their stated ≈6×10⁻³.)
   - (b) IGLOO wet-bulb vs the **M7 temperature plateau** 282.43 K (`..._Td.csv`), **±0.5 K**.
     IGLOO 282.33 K → **|ΔT| = 0.10 K**.
   Both panels are genuine IGLOO-vs-digitized-model-**curve** gates, not an analytic-`β`
   proxy.

## Comparison plot (`verify.py` → `OUTPUT/mhb98-water.svg`, non-gating)
Two panels vs the **normalized** `t/t_life` (t_life = D0²/K from the LSQ slope, extrapolated to
D²→0 rather than taken as the last sample, which lands at d²/d0²≈0.012): (a) `d²/d0²` — IGLOO, its LK kernel (gated oracle), and the
digitized **M7 model line**; (b) `T_d` — IGLOO wet-bulb history and the digitized M7 line.
Normalized (not absolute) because the paper's caption `D0=1.1 mm`⇒`D0²=1.21` disagrees with its
own plot (M7 starts at D0²≈1.10; axis tops at 1.2); since β (Eq. 28) is D0-independent and
matches to <1 %, the two normalized d²-law lines **coincide** (an honest presentation — the
coincidence is partly tautological, so the β annotation + panel (b) carry the real evidence).
Panel (a) also draws the **Ranz–Marshall experiment**, the 23 circles of Fig. 2(a) re-digitized
from the page image, normalized by the M7 line's own `D0²` and lifetime (the paper's axes: the
circles keep their place relative to M7; their own `D0` would be an extrapolation, the first
circle being at 77 s), with `β_exp` in the annotation; Fig. 2(b) has no experimental points.

**Lifetime matches too.** With the figure's `D0²` and the 0.17 m box, `t_life = D0²/K`
is **825.3 s (IGLOO) vs 818.5 s (M7) — +0.8 %**, i.e. exactly the residual `β` difference; both
curves span the full `0→1` normalized lifetime. (With the caption's `D0²=1.21` the gap is
**+10.9 %**, and in a 0.15 m box IGLOO's curve stops at `t/t_life=0.86` because the drop leaves
the box still at 15 % of `d0²` — **truncation, not a shorter lifetime**; IGLOO's lifetime is
in fact *longer*.)

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
measurement of Fig. 2b (400 dpi render, axis calibrated on the 280/282/284/286 major
ticks ⇒ **100.7 px/K**) puts the **M7/M8 marker symbols at 282.44 ± 0.03 K**. The digitized 282.43
is a faithful read.

***…but the figure cannot resolve 0.1 K between models.*** The same measurement shows the eight
plotted curves occupy **282.08 → 282.55 K (≈0.47 K)**:
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

*What the gap is NOT.* Sensitivity: a **1 % change
in any of k_g, D_v, L_v or p_sat moves the wet-bulb ≈0.076 K**, so 0.10 K ≈ 1.3 % of one property.
Every property is **MHB98's own Appendix-A value at their `T_R`** (table above); against the
generic values that moves the wet-bulb 282.30 → 282.33 K and β from 1.1 % → 0.8 %. So the residual is **not**:
- **`Le`** — MHB98 also run water at `Le=1` (Appendix A gives no `Γ_V`/`Sc` for water and states
  `Le=1` as the default in that case). IGLOO's `Le=1`, "13 % below real H₂O–air", is a
  deviation from *reality* that the paper **shares**; it cannot open an IGLOO-vs-M7 gap.
  Contrast tc-hexadecane, where `Le` (2.5, it sets `D_v`) is a genuine lever.
- **the `p_sat` closure** — their eq. (10) is the same single-anchor CC at `T_B` with constant
  `L_V` (both ~10 % below Magnus at 282 K).
- **`k_g`, `μ_g`, `c_p,g`, `L_v`, `ρ_l`** — the paper's own numbers, worth +0.025 K in total.
- **`f₂`** — M7's eq. (3) carries the evaporative heat-transfer correction `f₂ = β/(e^β−1)`
  (eq. 19) which IGLOO's CEM path lacks (`f₂≡1`), but here `β≈6.3e-3 ⇒ f₂=0.9969`: 0.31 %, and
  its sign *lowers* the wet-bulb, i.e. adding it moves IGLOO ~0.02 K **further** from M7. (It is
  a first-order term only at high evaporation rate — see the Fig. 4 decane scope note.)

*Honest limit (from the pixel measurement):* ~0.10 K remains unattributed, and
**it is not decidable against Fig. 2b** — it is smaller than the figure's own model-to-model
spacing (see above). In rate terms 0.10 K ≈ **1.3 %** of the group `ρ_g D_v L_v/k_g`, i.e. the
size of a single Table-1 *form* difference at this condition: `B_M,eq` vs `ln(1+B_M,eq)` is 0.3 %,
`(Y_s−Y_G)` vs `B_M,eq` is 0.6 %, `f₂=(1+B_T)⁻¹` (M4/M6) is 1.5 %. IGLOO's own known form
difference (`f₂≡1` instead of eq. 19) is 0.3 % and the wrong way. So the residual is a
model-**form**/bookkeeping difference of the same order as the spread MHB98 themselves call
negligible — **not a fixable IGLOO error**, and well inside the ±0.5 K gate.

*How it WOULD be decided:* only on a case where the models actually separate — **Fig. 4 decane**
(`T_G=1000 K`, `Re_d=17`, `test_mhb98_decane`), where `B_T=O(1)` so `f₂` is genuinely `≠1`, the LK
non-equilibrium term is not vacuous, and M1/M2/M5 visibly over-predict while only M7/M8
track the experiment. At Fig. 2's condition every discriminating term is ~1 %.

## Result & finding
Wet-bulb **T_p = 282.33 K vs the M7 plateau 282.43 K (|ΔT| 0.10 K)**; D²-law reproduced.
**IGLOO's integrated β = 6.30×10⁻³ matches the digitized M7 line β_M7 = 6.35×10⁻³ to 0.8%**
(K_M7 = 1.347×10⁻³ mm²/s) — IGLOO lands *on* MHB98's M7 curve, not just our re-coded kernel.
This is a genuine IGLOO-vs-paper-**model-curve** reproduction (both panels), tighter than an
analytic-β proxy against the paper's rounded 6×10⁻³ (~7% with the generic constants).

**Against the experiment (drawn, not gated).** The 23 Ranz–Marshall circles of Fig. 2(a),
re-digitized from the 1000-ppi page image (`reference/mhb98_fig2_exp_d2.csv`), give
**β_exp = 6.48×10⁻³** (their LSQ slope K = 1.3745×10⁻³ mm²/s, Eq. 28 with `check.py`'s
coefficient). They sit 0.018 mm² below the M7 line on average (about 0.7 of a marker radius)
and their slope is 2.0 % steeper than M7's; **IGLOO's β = 6.30×10⁻³ is 2.8 % below the
experiment**. The experiment is drawn in panel (a); no gate reads it.

**Scope note (out of this case):** Fig. 2 (water, Re=0) is a *weak* discriminator — all
eight MHB98 models coincide. The strong model discriminator is **Fig. 4 decane**
(`T_G=1000 K`, `Re_d=17`), where M1/M2/M5 visibly over-predict and only M7/M8 (non-eq)
track the experiment. That is a *convective* case, outside this one: `test_mhb98_decane`
(a unit test at fixed slip, see `VERIFICATION_MATRIX.md`).

## Reference data
- `reference/mhb98_fig2_M7_d2.csv`, `..._Td.csv` — **WebPlotDigitizer** exports of the M7
  model line (Fig. 2a/2b). The **gated** references: β_M7 (LSQ slope) and the T_d plateau.
- `reference/mhb98_fig2_exp_d2.csv` — the Ranz–Marshall **experiment**, the 23 filled circles of
  Fig. 2(a) re-digitized from the 1000-ppi page image (Hough transform; t ±0.9 s, D² ±1.2e-3 mm²);
  **plotted, not gated**. Method, calibration and the numbers in `reference/PROVENANCE.md`.
