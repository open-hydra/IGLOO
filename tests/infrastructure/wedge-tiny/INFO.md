# wedge-tiny — a single-layer sector of 5e-7 rad is a planar single layer

**Purpose.** The wedge detector (`allocation.f90`) takes a single-layer mesh for an axisymmetric
wedge when the azimuthal span of its two k-planes at the outermost node exceeds 1e-6 rad. At or below
that span the mesh takes the planar path, with one notice: there the fold's relative slack nears
roundoff, and below 1.46e-9 rad `nint(theta/delthe)` overflows a default integer at `|theta|` ≈ π.

**Fixture.** `tools/make_wedge_case.py INPUT --delthe-deg 2.8647890e-05 --omega 0` writes
`INPUT/solfile.tec`: `swirl-wedge`'s annular mesh (x ∈ [0, 2], r ∈ [0.19, 0.99], 200 × 40 cells, one
cell in the azimuth) with the two k-planes 5e-7 rad apart and no swirl (U = 1, V = W = 0). Its
`bc.txt` is byte-identical to `swirl-wedge`'s and is a symlink to it, as are `phase.txt` and
`properties.dat`. `input.ini` is `swirl-wedge`'s with `wp = 0`: three DB parcels at cell centres
(x = 0.105; y = 0.5, 0.5, 0.3) at the gas velocity.

**Gate (`check.py`).**
- (a) the log carries `Sector of  5.000E-07 rad at the outermost node, below 1e-6 rad: planar single
  layer` (the span within 1e-3 of the generator's) and `2D path: planar single layer`, and neither an
  `Axisymmetric wedge` nor a `Planar slab` line; no give-up;
- (b) three exit records, read by the header's column names, each at x = 2 with y and z at the
  injection values to 1e-6 (no force acts across the gas direction).

**RED.** With the threshold at 1e-9 rad the log reports `Axisymmetric wedge: delthe = 0.00000050 rad`
and the wedge's 2D path, so (a) fails three ways (no notice, no planar line, a wedge line); the exits
of (b) hold on both paths.
