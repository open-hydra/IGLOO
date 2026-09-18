# wedge-fold — the axisymmetric 200-face FOLD with a swirling parcel (ledger O22)

**Purpose.** The wedge fold — `axisymFold` rotating a parcel's position AND velocity back
into the ±delthe/2 sector about the symmetry axis — had never rotated anything on the suite:
every wedge fixture (`db-2daxi`, `axis-200`) has z ≡ 0, so the fold returned at its first
test, and `test_axis_dispatch` D1–D7 check only the scalar sign of the fold. This case gives
the near-axis parcel an azimuthal injection velocity, so it leaves the 1° sector within its
first ODE segment and the fold is exercised for real.

**The defect (found 2026-09-16 by the `af307dd` wp witness, not by reading).**
`rotateVector` (`src/lib/vectors.f90`) listed the Rodrigues matrix ROW-major and fed it to
`reshape(..., [3,3])`, which fills COLUMN-major, so it returned R^T = R(−θ): rotateVector(ŷ,
x̂, +1°) had z = −0.017452. `axisymFold` therefore rotated a parcel AWAY from the sector and
its 1000-iteration guard loop (no convergence check) walked it to ±180° and returned
silently. On this fixture the parcel was put at y = −9.6e-5 with its velocity turned by 178°
after its first segment, then 2-cycled with zero advance until outer maxIter (500000),
flagged gone. Fix: `order=[2,1]` on the reshape; the third caller (`initialization.f90`,
ds-injection directions at kd·π/3, live on 29 of 35 fixtures) negates its angle so the
kd↔direction pairing — hence ID and RNG-stream order — is unchanged (byte-inert, proven by
A/B); the fold is now one rotation by −n·|delthe| with n = nint(θ/|delthe|), uses |delthe|
(the loop tested against the SIGNED value, wrong for a reversed k ordering), and asserts its
post-condition.

**Fixture.** axis-200's `INPUT/` (all symlinks, including its `bc.txt` with face 3 tagged
`axisymmetric`) and its two DB parcels on the same axial station; ID 1 (y0 = 1e-4)
additionally carries `wp = 0.5` (≈ 5000 rad/s about the axis at injection), ID 2 (y0 = 0.55,
no wp) is the control and must reproduce axis-200's ID 2 exactly. `out-file = e` is kept so
the ODE system is axis-200's. INI traps as in axis-200 (single-space values, no `=` in
comment lines — this file tripped the second one once while being written).

**Gates** (`check.py`; all proven RED on the pre-fix binary `798e037e`, scratch
`wplan/o22/e2e_RED.txt`, and GREEN on `fc9612eb`):

| id | assertion | RED signature (pre-fix) |
|---|---|---|
| G1 | none of the give-up markers in `run_out.txt` (`no net progress`, `stuck in cell`, `Inner loop`, `outer maxIter`, `non-finite state`) | `Particle 1 hit outer maxIter (500000); flagging gone` |
| G2 | both parcels' last trajectory row has x > 2.0 (the outlet; the last row is the exit point) | ID 1 last x = −0.489651 |
| G3 | every trajectory row is inside the sector, print-aware: \|z\| ≤ \|y\|·tan(delthe/2) + 1e-6 and y > −1e-6 (rows are F12.6 and the sector half-width in z at r = 1e-4 is 8.7e-7, below one print unit, so `atan2` on printed values is ill-conditioned there) | 500000 rows with y = −9.6e-5 |
| G4 | vacuity: ID 1's first row has W = 0.5, and the sign of z flips between consecutive ID-1 rows ≥ 9 times | (passes on the broken run too: 499999 flips) |
| G5 | ID 2's last x equals axis-200's (2.051541) to 1e-6 | (passes: ID 2 is byte-identical to axis-200's) |

`Run_ODESolver err=`/`=> marking gone` (`Lib_Integration.f90`) are not in G1's list; a parcel
lost that way fails G2 (never reaches the outlet).

**Measured on the fixed binary** (OMP 1 and 5, three runs each): every output identical as
a sorted multiset (record order is the only OMP effect); ID 1: 347 rows, 18 z-sign flips,
min y = 4.5e-5 (its closest approach to the axis, at |p × v̂| = 1e-4·0.5/√1.25 = 4.47e-5),
exits at x = 2.058712 (axis-200's swirl-free ID 1 exits at 2.058713); ID 2: 260 rows,
byte-identical to axis-200's. `N_FOLD_MIN = 9` is half the measured count; it is
deterministic per binary (one thread per parcel) and only 1-ULP generation drift could move
a near-zero landing.

**What this case does NOT pin.** The fold's cadence is not "every 1° crossing": under
`mesh2D` the containment test is z-blind (`geometry.f90`, quad in x–y), so a k-face crossing
never interrupts a segment and the unfolded azimuth at a segment end can reach several
sectors (the closed-form fold handles any n). The gas sample between folds is exact since
W-plan item E (`Lib_Equations::sampleGas2D`: evaluated at the parcel's (x, r), velocity
rotated to its azimuth — `standard/swirl-wedge` gates it at the print floor); this gate is
behavioral (in-sector, no stall, exits) and measures nothing quantitative. What remains open
for an INTERIOR swirling parcel far from the sector (no k-face cap on its segment, the
unfolded azimuth bounded only by the dual's y-range through y = r cos θ): the dual cell is
still LOCATED by (x, y = r cos θ), so past a few sectors the interpolation extrapolates from
a cell inward of the parcel's true radius, and the source/euler deposits land in that cell —
the many-sector sub-case of ledger O23 — **reproduced 2026-09-18** on the swirl-wedge fixture with one
over-spun parcel (`y = 0.2`, `up = 1`, `wp = 10`, or `up = 0.05`, `wp = 2`): the first segment sweeps
≈ 27° / 34° (27 sectors folded at once), the true radius runs 0.200 → 0.224 / 0.234 while the dual
cell is located by y = r cos θ = 0.200, so 8 % / 23 % of that segment's deposit lands a row inward;
the sweep then decays (8°, 6°, 5° … < 1° as r grows with angular momentum). A transient here; a
parcel swinging past the axis (r_min ~ dx) would sustain it. Mechanism: the segment is ended by
the z-blind (x, y) containment, never by a k-face, and the fold runs only at the segment end; the
coherent fix is to end the segment at the sector edge (a `sectorOut` interrupt in `solout`), which
no plan has scoped — recorded, not built. The deposits'
*frame* (Cartesian at the parcel's azimuth, not the meridian plane) was O23's measurable
deposit half: fixed by `Lib_Equations::toMeridian` and gated by `standard/swirl-wedge-deposit`;
the sampling half is closed by E. This case's `euler1.tec` moved with that fix (its ID 1
leaves the plane) and its trajectory rows flipped by one F12.6 ULP (the rotated euler
states change the SDIRK4 step sequence: 349 rows for ID 1, 347 before).
