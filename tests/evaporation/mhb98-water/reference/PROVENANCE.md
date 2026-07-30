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

`verify.py` overlays the digitized M7 line in both panels, in **absolute** D²[mm²] (not
normalized: normalizing IGLOO and M7 by their different D0² would fake a ~10% lifetime
divergence despite β matching to ~1%). β is annotated in-panel.

## `mhb98_fig2.csv` — EXPERIMENT (eyeball, NOT plotted, NOT gated)
Ranz & Marshall (1952b) points, **eyeball-digitized** (±2–3% of axis range). Reads
systematically high (start 1.21, `K_exp ≈ 1.53e-3` → `β_exp ≈ 7.4e-3`) whereas in the figure
the exp points sit **on M7** (≈6.5e-3). So the earlier "IGLOO ~13% below experiment" is largely
an eyeball-digitization artifact, not physics. Kept for the record but no longer overlaid;
a proper WebPlotDigitizer re-digitization of the exp circles is a **follow-up** (not gated).

See `../INFO.md` for the finding.
