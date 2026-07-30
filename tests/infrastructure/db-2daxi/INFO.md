# db-2daxi — promoted legacy `test/assigned-pos` (2Daxi wedge + DB injection + euler)

**Purpose.** The only case exercising, together: the axisymmetric-wedge mesh path
(delthe fold, flattened gas dual), a REAL MOSE flow solution whose header carries
particle vars (`rho_p` — the B5 phantom-species trigger), DB injection
(`[IGLOO-BC] x/y/diam/mdot`), and euler-only output (`out-file = e`, the B6
accumulator path). Fixed by B5+B6 (85e8739) and green since; this promotion makes it
regression-protected in ctest.

**Fixtures.** `INPUT/solfile.tec → tests/common/solfile_mose.tec` (the real MOSE
nozzle solution, 21 MB); `bc.txt`/`properties.dat` case-local (differ from the common
box set); `phase.txt` common. The solver reads **geometry and gas from the solfile** —
MOSE's output — so nothing here specifies a mesh. The old `[GRIB-tecgen]` section and
the `MESH → tests/common/MESH` symlink were removed 2026-07-30: IGLOO never parsed
`[GRIB-*]` (no reference anywhere in `src/`), and the `mesh.dat` that section named
never existed in the first place.

**Gate (check.py) — deliberately md5-free.** Trajectory bytes drift by ~1 ULP across
compiler AND cmake-configure generations (observed 2026-07-15: a plain build_verif
reconfigure moved assigned-inj/assigned-pos multiset md5s while sym-bc held). Byte
regression stays with the manual same-session A/B protocol; this gate asserts the
behavioral contract: 2 particles, ≥50 rows each, both exit at the outlet (x > 2.0),
all trajectory rows finite with T∈[200,3700] K and dp∈(0,1.2e-4], euler field finite.

**History.** Bugs B5/B6, fixed 2026-07-15.
