# Beckstead (2005) Fig. 1 reference data — provenance

M. W. Beckstead, "Correlating Aluminum Burning Times," *Combust. Explos. Shock Waves* 41(5) (2005)
533–546, DOI 10.1007/s10573-005-0067-2 (`[Beck05]`). **Fig. 1** (journal p. 538, PDF page 6): "Aluminum
burning time measured under a wide variety of conditions and test techniques" — 12 sources ([45], [41],
[37], [39], [33,36], [38], [35], [43,44], [50], [42], [47,48], [46]), log-log t [ms] against D [µm], with
the label "t = 0.003 D^1.99; r² = 0.87".

## `beckstead_fig1.csv` — the measured burning times (plotted, NOT gated)
Columns `D_um, t_ms, source, kind`; the `#` lines carry the method. Read by `verify.py` only.

- **Source image.** PDF page 6 holds one raster, a 945×827 8-bit gray JPEG placed at 300 ppi
  (`pdfimages -list -f 6 -l 6`), extracted as stored (`pdfimages -f 6 -l 6 -j`; no resampling).
- **Calibration.** The frame and the gridlines are located from dark-pixel profiles: x lines at px
  86.0, 511.5, 934.0 (D = 10, 10², 10³ µm), y lines at px 15.0, 167.5, 319.5, 471.5, 623.5, 774.0
  (t = 10⁴ … 10⁻¹ ms). Least squares in log space gives `log10 D = 2.358484e-3·px + 0.79599` and
  `log10 t = −6.585150e-3·py + 4.10224`, max residual **1.0 px in D (0.54 %)** and **0.81 px in t
  (1.24 %)**; 1 px is 0.54 % in D and 1.53 % in t. The drawn fit line (below) has the exponent 1.989
  against the label's 1.99: the relative scale of the two axes is right to 0.05 %.
- **Detection.** The 2-px gridlines and the tracked fit line are erased first (a band pixel becomes the
  AND of its two outside neighbours, so a marker stroke crossing the line survives). The templates are
  the 12 legend glyphs; a binary hit-or-miss matching pursuit with explaining-away accepts the best
  glyph placement while score = new + old − miss ≥ 0.65 and at least 35 % of the glyph's ink is new, the
  explained pixels counting as hits for the next marker, so a stack resolves top-down. **183 markers;
  4.8 % of the plot's ink is left unexplained.** Score thresholds 0.60 and 0.65 give the same set.
- **Labels (`source`, `kind`).** A plain NCC match (≥ 0.55) at the same place confirms the glyph label
  of 91 markers (`kind = isolated`). The other 92 lie in stacks (`kind = stack`, the dense columns at
  D ≈ 20, 60, 90–100, 200 and 230–400 µm), where the glyphs merge into blobs: their positions hold,
  their source labels are not recovered and `source` is empty. `source` is empty also on the `isolated`
  row at (258.66 µm, 78.74 ms): its NCC match lies on no `+` glyph — Fig. 1 has no `[38]` marker beyond
  D ≈ 80 µm, and the edges of the neighbouring open squares lie under it. 90 rows carry a label; all 183
  positions are kept. `verify.py` draws every point as one glyph.
- **Precision per point.** ±0.5 px centring: D ±0.27 %, t ±0.76 %; the triangle glyphs carry a ±1.5 px
  centring ambiguity (bounding box against centroid), ±2.3 % in t.
- **Not recovered.** The abstract (p. 533) counts "almost 400 datum points"; the page resolves 183
  markers, since markers that coincide on the page cannot be separated.

## Facts of the paper the case relies on
1. **The label rounds the fit.** The drawn line, tracked column by column over 623 columns
   (D = 10.3–981 µm, rms 1.1 px), is **t = 0.00258 D^1.989**; at the label's exponent 1.99 the
   coefficient is 0.00257. The label's 0.003 D^1.99 lies 16.5–16.8 % above the drawn line at every D.
   Fig. 9's lower line coincides with the drawn one. `verify.py` draws the line as drawn.
2. **Fig. 9's upper line is not Eq. 7.** Its two clean ends give t' = 0.358 D^1.934 (local slopes 1.95
   at D = 10–18 µm and 1.93 at the upper end); forced to n = 1.8 its coefficient would be 0.51 at 13.5 µm
   and 0.88 at 860 µm, so no constant unit factor on X_eff maps Eq. 7 (a = 0.00735, n = 1.8) onto it.
   With X_eff as a mole fraction, Eq. 7 puts the raw times inside the Fig. 1 cloud (air: 1.2–2.8× above
   the drawn fit); with X_eff in percent it lies 100× below. The case reads X_eff as a mole fraction,
   IGLOO's convention (`X-eff` in (0, 1]). Fig. 9 is not drawn.
3. **T0** is Eq. 7's "initial temperature", not defined further; the case sets the gas and the particle
   to the same 300 K, so both readings coincide.

See `../INFO.md` for the case.
