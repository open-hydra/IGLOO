# Gallier et al. 2021, Figs. 7 and 11 — the deposition law's data, provenance

Gallier, S.; Braconnier, A.; Godfroy, F.; Halter, F.; Chauveau, C., "The role of thermophoresis on aluminum oxide
lobe formation," *Combust. Flame* 228 (2021) 142–153, DOI 10.1016/j.combustflame.2021.01.039 (`[Gal21]`, HAL
preprint hal-03135480; the PDF read here has md5 `34c31575e49519379999dd88fd4e21cc`). Their 0-D model of the lobe's
growth closes the retained fraction of the oxide the burn forms as `f_dep = g(d_Al) h(f_ox)`, "h ... directly taken
from results in Fig. 11", "g obtained ... from a linear regression on results presented in Fig. 7" (PDF pp. 13–14).

## The files

- `gallier21_fig11_fdep.csv` — Fig. 11 (PDF p. 13): `f_dep` [%] against the lobe's oxide mass fraction `f_ox` at
  d_Al = 70 µm, the six computed geometries f_ox = 0, 0.1, 0.3, 0.5, 0.7, 0.9.
- `gallier21_fig7_fdep.csv` — Fig. 7 (PDF p. 11), the thermophoresis + diffusiophoresis series (green squares):
  `f_dep` [%] against the aluminium diameter `d_Al` [µm] of a drop without lobe, the six computed sizes 20, 50, 70,
  100, 220, 400 µm. Fig. 11's f_ox = 0 value (6.9522 %) is this series at 70 µm (6.9633 %), not the
  thermophoresis-only series (7.5576 % at 70 µm): Fig. 11 carries the green series.

Each table holds the computed abscissae (the cases Gallier ran, read back within 0.08 µm and 0.0006 by both reads)
and the mean of two reads of the ordinate, to 4 decimals.

## Source images

Both figures are embedded JPEG rasters, no vector data (`pdfimages -list`): Fig. 7 1000 × 636 px, Fig. 11
1000 × 572 px, both at 324 ppi, read losslessly with `pdfimages -j -f 11 -l 11` and `pdfimages -j -f 13 -l 13`
(md5 `32261d00063151e888a114a8d8a429e1` and `52f0b47edba11ef4418a13c32a16a58c`), no resampling.

**Calibration.** The labelled ticks point inward from the axis lines (Fig. 7: x-axis row 559, y-axis column 99;
Fig. 11: row 488, column 114). Fig. 7: x = 0, 50, …, 400 µm at columns 99.0 … 950.5 (2.1288 px/µm), y = 0, 2, …,
10 % at rows 559.0 … 45.5 (51.35 px per %). Fig. 11: x = 0, 0.2, …, 1 at columns 122.5 … 966.5 (844 px per unit),
y = 0, 10, …, 60 % at rows 488.5 … 53.0 (7.26 px per %).

## The two reads

- **Read 1.** Every labelled tick, its centre the darkness-weighted mean over ±4 px of its run in a 6-px band inside
  the axis; a least-squares line per axis (largest residual: Fig. 7 0.067 µm, 0.0031 %; Fig. 11 0.00019, 0.088 %).
  A marker is its colour mask (red: R > 150, G < 90, B < 90; green: G > 150, R < 140, B < 120) opened with a 9 × 9
  square, which drops the 2–3 px connecting lines; its centre is the blob's area centroid.
- **Read 2.** The two end ticks of each axis only, each at the midpoint of its dark run on one row or column. A
  marker is the hole-filled colour blob; a circle's centre is the least-squares circle through its boundary pixels,
  a square's the midpoint of its extent over the rows and columns where it spans more than half its width.

Fig. 7, green squares (`f_dep` in %):

| d_Al (µm) | read 1 x | read 2 x | read 1 | read 2 | spread | table |
|---|---|---|---|---|---|---|
| 20 | 19.92 | 19.96 | 7.8612 | 7.8578 | 0.0034 | 7.8595 |
| 50 | 49.98 | 50.03 | 7.3939 | 7.3905 | 0.0034 | 7.3922 |
| 70 | 70.02 | 69.99 | 6.9646 | 6.9620 | 0.0026 | 6.9633 |
| 100 | 99.78 | 99.82 | 6.1477 | 6.1441 | 0.0036 | 6.1459 |
| 220 | 220.04 | 220.08 | 3.2659 | 3.2619 | 0.0040 | 3.2639 |
| 400 | 399.97 | 400.00 | 1.6108 | 1.6066 | 0.0042 | 1.6087 |

Fig. 11, red circles (`f_dep` in %; the f_ox = 0 marker is cut by the y-axis line on its left, which moves its x,
not its y):

| f_ox | read 1 x | read 2 x | read 1 | read 2 | spread | table |
|---|---|---|---|---|---|---|
| 0 | 0.0005 | 0.0004 | 6.9507 | 6.9538 | 0.0031 | 6.9522 |
| 0.1 | 0.1001 | 0.1000 | 10.1326 | 10.1359 | 0.0032 | 10.1342 |
| 0.3 | 0.2998 | 0.2998 | 15.0341 | 15.0369 | 0.0027 | 15.0355 |
| 0.5 | 0.4999 | 0.4999 | 19.2849 | 19.2806 | 0.0044 | 19.2828 |
| 0.7 | 0.7002 | 0.7001 | 28.0499 | 28.0395 | 0.0104 | 28.0447 |
| 0.9 | 0.8995 | 0.8994 | 52.0898 | 52.0640 | 0.0259 | 52.0769 |

The largest spread is 0.0042 % on Fig. 7 (0.22 px) and 0.026 % on Fig. 11 (0.19 px). A marker read 2 px off would
move its value by 0.04 % (Fig. 7) or 0.28 % (Fig. 11).

## The law built from them (`xi-law=Gallier21`)

Gallier gives the product form and the figures `g` and `h` come from. The normalisation of `g` at 70 µm and the
values held outside the computed points are IGLOO's readings, not stated in the paper.

- h(f_ox): Fig. 11's points over 100, linear between them; past f_ox = 0.9, where Fig. 11 has no point, held at
  0.520769 (IGLOO's reading).
- g(d_Al): the least-squares line through Fig. 7's six points, `f = 8.007948 − 0.01722580 d_Al` (% and µm), over
  its value at 70 µm, `g = f(d_Al)/f(70)` (IGLOO's reading): g(70) = 1, so the law at 70 µm is Fig. 11, whose
  f_ox = 0 value is Fig. 7's at 70 µm. Outside the computed 20–400 µm g holds its end values (1.126620 below 20 µm,
  0.164305 above 400 µm; IGLOO's reading).
- ξ = g(d_Al) h(f_ox) at the drop's current metal diameter and oxide mass fraction, clamped to [0, 1].

## Checked against Gallier's own numbers

The cap's equation in the metal mass, dm_ox/dm_Al = −ν ξ(d_Al, f_ox) with d_Al = d_Al0 (m_Al/m_Al0)^{1/3}, holds no
time, burn law or stretch: the residue's size ratio β = d_res/d_Al0, β³ = (m_ox,end/m_Al0)(ρ_Al/ρ_ox), depends on
the law alone. With the densities at the lobe's 2570 K (ρ_Al 1673.7, ρ_ox 2820.3 kg/m³) and RK4 in m_Al on 20 000
steps:

| quantity (Gallier, PDF pp. 13–14) | Gallier | the tables |
|---|---|---|
| β̃ at d_Al0 = 70 µm, f_ox0 = 0, 0.05, 0.10, 0.15, 0.20 (Table 3: the deposited oxide only) | 0.55, 0.59, 0.61, 0.63, 0.65 | 0.565, 0.595, 0.617, 0.636, 0.652 |
| β at d_Al0 = 200 µm, f_ox0 = 0 | 0.49 | 0.507 |
| β at d_Al0 = 400 µm, f_ox0 = 0 | 0.40 | 0.400 |
| β of a small drop (20 µm), f_ox0 = 0 | 0.57 | 0.581 |

The largest difference is 0.017. Table 3's β column with an initial lobe (0.64, 0.70, 0.75, 0.81) is not compared:
it depends on how Gallier sizes a drop that carries an initial lobe, which the paper leaves open.
