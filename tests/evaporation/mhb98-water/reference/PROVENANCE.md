# MHB98 Fig. 2 reference data — provenance

Miller, R. S.; Harstad, K.; Bellan, J., "Evaluation of equilibrium and non-equilibrium
evaporation models for many-droplet gas–liquid flow simulations," Int. J. Multiphase Flow
24(6) (1998) 1025–1055, DOI 10.1016/S0301-9322(98)00028-7. Case: single water droplet,
D0=1.1 mm, T_d,0=282 K, T_G=298 K, Re_d=0 (quiescent, 1 atm), **Fig. 2**.

## `mhb98_fig2_M7_d2.csv` / `mhb98_fig2_M7_Td.csv` — MODEL line M7 (GATED reference)
**WebPlotDigitizer** exports of the **M7** curve (Langmuir-Knudsen, uniform-temperature) from
Fig. 2(a) D²[mm²] (80 pts) and Fig. 2(b) T_d[K] (31 pts). M8 (finite-conductivity LK) was
**not** digitized: it is superimposed on M7 in both panels. IGLOO runs CEM+LK = the
uniform-temperature LK model = **M7**, so M7 is the correct reference (at Fig. 2's low rate
all M1–M8 coincide anyway, so the choice does not bite here — but it justifies the mapping).

- **`d2`** extrapolates to **D0²_M7 ≈ 1.10 mm²** at t=0 — the figure's own start value. The
  caption says `D0=1.1 mm` (→ 1.21); this ~9% internal offset is a property of the paper's
  plot. The gated quantity is the **slope / β** (Eq.28, D0-independent):
  `K_M7 ≈ 1.347e-3 mm²/s → β_M7 ≈ 6.35e-3`, computed by LSQ in `check.py` at gate time (β_M7 was
  6.51e-3 before 2026-07-29, when Eq. 28's coefficient still used ad-hoc `ρ_l`/`Pr_G`/`μ_G`
  instead of MHB98's Appendix-A values; the gate is a slope ratio, so the coefficient cancels).
- **`Td`** plateau (tail mean) = **282.43 K**.

### GATE (`check.py`, Tier-P), built on these curves
- (a) IGLOO's integrated β vs **β_M7** (LSQ slope of `d2`), **±10%**. β_M7 is guarded to lie
  within ±20% of the paper's stated **β ≈ 6×10⁻³** (p.1038, all M1–M8 at Re_d=0), which is
  also the fallback if the file is absent. IGLOO β = 6.30e-3 → **rel 0.8%** vs β_M7.
- (b) IGLOO wet-bulb vs the **M7 plateau** 282.43 K, **±0.5 K**. IGLOO 282.33 → **|ΔT| 0.10 K**.

### Property provenance (2026-07-29)
MHB98 **Appendix A** publishes every correlation and its source: **air and water both Harpole
(1981)** (benzene/hexane Reid et al. 1987 + *Petroleum Refining Data Book* 1992; decane Abramzon &
Sirignano 1989; heptane Park & Aggarwal 1995), frozen once at `T_R = T_WB`(eq. 27) = 293.9676 K.
The case now runs those numbers (see `../INFO.md`). Earlier revisions of this file claimed the
property tables were "not published" — **wrong**. Water lists no `Γ_V`/`Sc`, so the Appendix's
stated `Le=1` default applies: MHB98's water runs at `Le=1`, as IGLOO does.

`verify.py` overlays the digitized M7 line in both panels, **normalized**: d²/d0² against
t/t_life, each curve by its own D0² and its lifetime t_life = D0²/K (the LSQ slope extrapolated
to D² → 0), because the caption's D0 (D0² = 1.21 mm²) and the plotted start (≈ 1.10 mm²)
disagree while β (Eq. 28) is D0-independent. β is annotated in-panel. Panel (a) also draws the
experiment below, normalized by the M7 line's own D0² and lifetime, so the circles keep their
place relative to M7.

## `mhb98_fig2_exp_d2.csv` — EXPERIMENT (Ranz & Marshall 1952b), re-digitized (plotted, NOT gated)
The 23 filled circles of Fig. 2(a) as two blank-separated columns `t_s D2_mm2` (D² in mm², t in
s; the `#` lines carry the method). Read by `verify.py` only.

- **Source image.** PDF page 14 (journal p. 1038) holds one image, a 2318×4326 CCITT bilevel
  stencil placed at 1000 ppi (`pdfimages -list -f 14 -l 14`), read losslessly at native
  resolution with `pdfimages -f 14 -l 14 -png` (no resampling).
- **Calibration.** The tick marks are the ink runs in a band 6–13 px inside the frame, each
  assigned the nearest nominal value (x pitch 40 s, y pitch 0.04 mm²) after a first map on the
  frame centrelines. Least squares on **23 x ticks** (40…960 s; 800 s lies under the curves) and
  **27 y ticks** (0.04…1.04 and 1.16 mm²; 1.08 and 1.12 lie under the start markers) gives
  `t = 0.5645567 s/px · px − 219.536` and `D² = −7.954332e-4 mm²/px · py + 1.384063`, max
  residual **0.67 s** (1.2 px) and **9.2e-4 mm²** (1.2 px). The thick frame (12–13 px) maps to
  t = 2.3 s and 1003.4 s, D² = −0.002 and 1.198 mm²: the ticks are the reference, not the frame.
- **Detection.** A fixed-radius circle Hough transform on the ink boundary, R = 31.3 px (the
  inscribed radius of the isolated first disc; 17.7 s, 0.025 mm²). A centre is kept where the
  vote fraction is ≥ 0.30, it is the local maximum within 14 px and the ink at the centre is a
  filled disc (distance transform ≥ 0.85 R, which rejects the open M7 squares and M8 diamonds);
  the legend box is masked. **23 circles, t = 77.4 … 727.8 s**; the set is the same for vote
  thresholds 0.15–0.45 and separations 14–20 px. At a 10 px separation two weak sub-peaks appear
  inside the overlapping tail (t ≈ 711 and 724 s, votes 0.29 and 0.32, 14 and 19 px from resolved
  centres): they cannot be told from the line bundle and are not listed.
- **Precision per point.** t ±0.9 s, D² ±1.2e-3 mm² (calibration residual ⊕ 1 px centroid) for
  the isolated discs; about ±3.5 s in the overlapping tail (697–728 s).
- **Calibration check, independent of this file.** The centre of the model-line bundle at 14
  marker-free columns lies **+0.0050 mm²** (rms 0.0050) from the WebPlotDigitizer M7 export,
  which is calibrated independently; the bundle is 0.016 mm² wide.
- **What the figure does not carry.** Fig. 2(b) has no experimental points: its legend lists
  M1–M8 only. The caption's D0 = 1.1 mm (D0² = 1.21 mm²) is not a plotted circle: the first
  circle is at t = 77 s, D² = 0.998 mm².

| quantity | the circles (this file) | M7 line (`mhb98_fig2_M7_d2.csv`) |
|---|---|---|
| K = −dD²/dt [mm²/s] (LSQ) | **1.3745e-3** | 1.3473e-3 |
| β (Eq. 28 with `check.py`'s `BETA_COEF` = ρ_l Pr/(8 μ_G) = 4.711318e6 s/m², Pr = 0.70826) | **6.48e-3** | 6.35e-3 |
| D²(t = 0) of the LSQ line [mm²] | 1.093 | 1.100 |
| scatter about its own line (rms / max) [mm²] | 0.0095 / 0.031 | — |

The circles sit **0.018 mm² below** the M7 line on average (rms 0.021, range −0.032…+0.007),
about 0.7 of a marker radius, and their slope is **2.0 % steeper** than M7's. IGLOO's
β = 6.295e-3 lies **2.8 % below** the experiment. Not gated: the experiment is drawn.

See `../INFO.md` for the finding.
