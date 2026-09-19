# wedge-axis-row — a wedge whose axis row sits at r = 1e-8 stays a wedge

**Purpose.** The converse of `planar-slab`: the slab/wedge classifier (`allocation.f90`) compares
the k-plane span at the outermost node with the span at an inner node (a slab's decays as 1/r, a
wedge's is constant). Its first version took "inner" as the smallest `r > 0`, and hydra's JPL
nozzle mesh (ATLAS BCB) keeps its axis row at `r = 1e-8` with roundoff `z` — a node whose azimuth
is noise. Measured on the JPL 3-solver case: the wedge was taken for a slab (no fold, no 2.5D gas
sampling), IGLOO injected 1.85× its share and the centreline exit Mach went from 2.104 to 2.449.
The inner node must be clearly off the axis (`r > 1e-3 r_max`).

**Fixture (`make_fixture.py`, output committed).** 20 × 8 cells, 1-degree wedge about the x axis,
k-planes at ∓delthe/2, j = 0 row at `r = 1e-8` with `z = 0`, uniform gas u = 5 m/s (`GAM`, `mil`,
`mit`, `kl` present so nothing is unbound); `bc.txt` in IGLOO's record order (x-min 300, x-max 400,
axis 200, outer 300, wedge planes 200). Two parcels in the meridian plane at the gas velocity
(y = 0.01 and 0.03). Phase and properties are the common ones.

**Gate (`check.py`).** Log reports `Axisymmetric wedge: delthe = 0.01745329` and no `Planar slab`;
no give-up; both parcels leave through the outlet (x ≥ 0.0999) with y and z at their injection
values to 1e-9 (a meridian-plane parcel without swirl stays there).

**Red/green.** One-pass classifier: `Planar slab ... 0.00000000 rad at the innermost` (FAIL).
Two-pass classifier with the radius floor: wedge, 21 rows per parcel, both exits at x = 0.1.
