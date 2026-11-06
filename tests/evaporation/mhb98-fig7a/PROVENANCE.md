# MHB98 Fig. 7(a) reference data — provenance

Miller, R. S.; Harstad, K.; Bellan, J., "Evaluation of equilibrium and non-equilibrium evaporation models for
many-droplet gas–liquid flow simulations," Int. J. Multiphase Flow 24(6) (1998) 1025–1055,
DOI 10.1016/S0301-9322(98)00028-7, **Fig. 7(a)** (journal p. 1045). Caption: "Temporal evolution of the
non-equilibrium contribution to the surface mole fraction for Langmuir–Knudsen model M7 for the base conditions of
Fig. 4; i.e. decane (T_B = 447.7 K) with T_G = 1000 K, T_d,0 = 315 K, D0 = 2 mm and Re_d,0 = 17: calculations end
when m_d = 0.01 m_d,0 and the results are for various (a) initial droplet diameters". Text (p. 1044): "for an
initial 10 µm diameter droplet, the non-equilibrium contribution ranges from approximately 20% to 80% of the surface
mole fraction during the course of evaporation".

The panel plots r = 100 (2 L_K β/D)/χ_s,neq [%] — the depression of the surface mole fraction by the Knudsen layer
(eq. 15) relative to the non-equilibrium value — on a logarithmic axis against t/τ_d,0, τ_d,0 = ρ_l D0²/(18 μ_G),
for D0 = 10 µm (solid), 50 µm (dashed), 500 µm (dash-dot) and 2000 µm (long dash), all at Re_d,0 = 17 (the
caption's base conditions: the values at t = 0 scale as 1/D0). The traced curves span 21 to 69 % for 10 µm and
0.1 to 0.3 % for 2000 µm, Fig. 4's drop.

## `mhb98_fig7a_M7.csv` — the four M7 curves at fixed nodes (gated by `test_mhb98_decane` DC4)

- **Source image.** PDF page 21 holds one image, a 2291×6500 CCITT bilevel stencil placed at 1000 ppi
  (`pdfimages -list -f 21 -l 21`), read losslessly at native resolution with `pdfimages -png -f 21 -l 21`
  (no resampling). Panel (a) is framed by rows 232–248 and 1738–1755, columns 402–420 and 2173–2187.
- **y calibration.** The tick marks are the ink runs at most 6 px tall in a band 11–31 px inside the left frame;
  the three labelled majors (10², 10¹, 10⁰) are the 12 px thick ticks 31–41 px inside it (10⁻¹ lies under the
  2000 µm curve, 10⁻² on the frame). Each minor is assigned the nearest nominal m·10^k after a first map on the
  majors. Least squares on **40 ticks** gives log10 r = −3.32103517e-3 · row + 3.799193 (301.11 px per decade),
  max residual **1.62 px**, rms 0.72 px: 0.765 % in r per pixel.
- **x calibration.** The ticks are the ink runs in a band 20–30 px above the bottom frame, the run at the right
  frame (t/τ = 1.2) excluded, each assigned the nearest multiple of 0.04 after a first map on the majors 0.2 and 1.0.
  Least squares on the **25 ticks of [0.04, 1.00]** gives t/τ_d,0 = 6.80368068e-4 · col − 0.279201
  (1469.79 px per unit), max residual **1.70 px**, rms 0.72 px. The four ticks beyond 1.00 lie 1.55 to 3.67 px off
  that line, so above 1.00 the map is piecewise linear between consecutive ticks (1.00, 1.04, …, 1.16); it serves the
  printed nodes only.
- **Tracing.** In each column the curves are the ink runs 8–22 px tall (the strokes are 12 px; a run of two merged
  strokes is longer and is not taken), each curve the run nearest its centre in the previous column, at most 30 px
  away, from four seeds in column 527 (t/τ = 0.079, the first column right of the y ticks that crosses all four
  strokes). The legend box (right edge column 1336, bottom edge row 650) is masked. The stroke centre is the run's
  midpoint.
- **Nodes.** t/τ_d,0 = 0.10, 0.15, …, 1.10; the value at a node is the median of the stroke centres within ±0.006
  (±8.8 px) of it, kept when at least 5 columns are traced there, with the number of columns and the spread of their
  centres (0 to 6.5 px). Gated: the 57 nodes of [0.20, 1.00] (10 µm 17, 50 µm 11, 500 µm 13, 2000 µm 16). Not
  listed, a gap of a dashed curve: 50 µm at 0.30, 0.45, 0.55, 0.70, 0.80, 0.95; 500 µm at 0.40, 0.65, 0.70, 0.95;
  2000 µm at 0.80 (and outside the window 50 µm at 0.15, 1.10, 500 µm at 0.15, 1.05, 2000 µm at 0.15).
- **Precision per node.** y: the calibration's worst residual and ±1 px of stroke centre, √(1.62² + 1²) = 1.9 px =
  1.45 % in r; x: 1.70 px = 1.16e-3 in t/τ times the steepest d ln r/d(t/τ) of the window, 2.7 (10 µm at 1.0),
  0.31 %; together **≤ 1.8 %**.

## The model behind the curves, and the window

The re-coded M7 chain of `test_mhb98_decane.f90` (`m7_mdot`, `m7_f2`, fixed slip, Appendix A's decane and air at
T_R = 420.19 K, Le = 1, MHB98's Sh/Nu coefficient 0.552), integrated per D0 to m = 0.01 m0, sits within **3.24 %**
of the 57 gated nodes (10 µm at 0.20) and within 2.72 % on the others; with IGLOO's own coefficient 0.600 and no
rescale it sits 4.2 to 13 % above them. At t/τ = 0.10 the chain lies 4.7 to 7.5 % above the figure on every curve,
an offset that decays by 0.4 and that the chain does not explain; beyond 1.0 the x calibration's residuals grow
and r steepens. The window [0.20, 1.00] was chosen after these residuals were seen; the nodes outside it are
printed, not gated.

The Schmidt number: with the Appendix's vapour diffusivity (Sc = 3.018) the chain's t = 0 values fall 15 times
below the plotted ones; Le = 1 (Sc = Pr), the property set the Fig. 4 legs already use, reproduces them
(22.46, 4.485, 0.4483, 0.1121 % at t = 0, against 21.3, 4.26, 0.430, 0.106 % in the traces' first column,
t/τ = 0.079).
