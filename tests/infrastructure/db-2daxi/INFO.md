# db-2daxi — 2Daxi wedge + DB injection + euler output

**Purpose.** Exercises, together: the axisymmetric-wedge mesh path (delthe fold,
flattened gas dual), a REAL MOSE flow solution whose header carries particle vars
(`rho_p`), DB injection (`[IGLOO-BC] x/y/diam/mdot`), and euler-only output
(`out-file = e`). The `rho_p` header exercises the species count: `allocate_blocks` counts
as gas species only the variables `import_gas` fills (`Roi(` / `rho(`, with the bracket), so
a particle variable allocates no species row that is never imported. The euler-only output
exercises the accumulator sizing: both accumulators take their dims from the source blocks,
which `allocate_blocks` sizes whatever `out-file` selects. `axis-200` and `wedge-fold` run on
the same MOSE solution.

**Fixtures.** `INPUT/solfile.tec → tests/common/solfile_mose.tec` (the real MOSE
nozzle solution, 21 MB); `bc.txt`/`properties.dat` case-local (differ from the common
box set); `phase.txt` common. The solver reads **geometry and gas from the solfile** —
MOSE's output — so nothing here specifies a mesh: IGLOO parses no `[GRIB-*]` section (no
reference anywhere in `src/`), and the case carries none and no `MESH/`.

**Gate (check.py) — deliberately md5-free.** Trajectory bytes drift by ~1 ULP across
compiler AND cmake-configure generations (measured: a plain reconfigure of the
verification tree moves this case's trajectory-multiset md5 and another DB-injection case's,
while a third case's holds). Byte regression stays with the manual A/B protocol (both
builds in one session); this gate asserts the
behavioral contract: 2 particles, ≥50 rows each, both exit at the outlet (x > 2.0),
all trajectory rows finite with T∈[200,3700] K and dp∈(0,1.2e-4], and the ord2 eulerian
deposit, audited as in axis-200: **E1** mass retention
Σρ_p V / Σṁ t, total AND per parcel (the two-mass decomposition of every cell is unique:
w₁ = n(q−m₂)/(m₁−m₂), w₂ = n(m₁−q)/(m₁−m₂)), 1.001438 / 1.000821 (ID 1) / 1.001482 (ID 2) ± 1e-4
(≈1 by physics; volumes recomputed from the tec nodes, residence by quadrature over the crossing
rows in FILE order — ID 2 moves backward in x for its first 75 rows; the floor is ZERO —
`euler1.tec` bit-identical over 3 runs × OMP {1, 5}); **E3** the per-drop mass identity with TWO
admissible masses, since ID 1 carries d = 1e-4 and ID 2 d = 2e-5: every depositing cell's ρ_p/n_p
equals ρπd³/6 for one of them to 1e-12 or lies strictly between (a cell both smeared deposits
reach); measured 4691 / 5013 / 232 mixed / 0 outside; **E4** the outlet-band share (deposit in the
last 3 geo columns over the total) 2.514026e-3 ± 1e-3 relative — the sharp witness of the boundary
clip, which also pins the smoother (re-measure if `mollify-passes` changes). No E2 (both parcels
off-axis). The helpers are imported from `axis-200/check.py`. The ~1e-3 offsets of E1 from 1 are
NOT attributed: the same run at `gas-order = 1` reads 0.996428 / 1.001233 / 0.996097, so they are a
property of the projection (ord2 dual→geo remap on a curvilinear mesh) and/or the quadrature; the
gate pins the measured ord2 values.

**Proven RED.** (1) Unclipped boundary dual cells (`insideFrac` returning 1, in a
throwaway build): the output differs from GREEN only in
`euler1.tec` (trajectories multiset-identical, outloc byte-identical): E1 = 1.001230 / 1.000453 /
1.001285 — 2.1e-4, 3.7e-4 and 2.0e-4 away, small because both parcels touch only the outlet
boundary — hence the tight 1e-4 tolerance (2× the weakest signature; half-ULP noise on every
printed row moves E1 by ≤ 1.4e-6); E4 reads 2.331509e-3, −7.3 % (the deficit sits in the last
columns: −2.2e-3 … −9.5e-2 on columns 194–199, 8 passes + 1 of binomial reach). (2) A factor-2
field, which a gate without the eulerian audit passes intact: E1 = 2.0029 and E3 fails in 4757
cells. Limitation stated: E3's
"strictly between" rule accepts any error factor in (1, 125) on a mixed cell — the per-parcel E1
is the stronger statement and carries the magnitude.
