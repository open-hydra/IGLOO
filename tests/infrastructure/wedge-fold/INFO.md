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
a sorted multiset (record order is the only OMP effect); ID 1: 337 rows (347 before the
meridian-frame deposits, 349 after them, 337 with sector-edge segments), 18 z-sign flips (10 with
sector-edge segments),
min y = 4.5e-5 (its closest approach to the axis, at |p × v̂| = 1e-4·0.5/√1.25 = 4.47e-5),
exits at x = 2.058712 (axis-200's swirl-free ID 1 exits at 2.058713); ID 2: 260 rows,
byte-identical to axis-200's. `N_FOLD_MIN = 9` is half the measured count; it is
deterministic per binary (one thread per parcel) and only 1-ULP generation drift could move
a near-zero landing.

**What this case does NOT pin.** The fold's cadence IS now "every 1° crossing": since
2026-09-18 the containment test carries the azimuth band (`geometry.f90::isPointInsideCell`,
optional `sectorOut`), so a k-plane crossing ends the ODE segment like any other face, the
refinement pins it to the plane (`IGLOO_bcBox::sectorDs`) and `axisymFold(force=.true.)` rotates
by exactly one sector before the cell logic runs; the run prints `wedge sector folds: N
(multi-sector: M)` and M must be 0. This gate is still behavioral (in-sector, no stall, exits) and
measures nothing quantitative about it — `standard/swirl-wedge-spin` does (fold count against the
oracle sweep, eulerian mass per cell against the projected residence). History: the many-sector
sub-case of ledger O23 was reproduced 2026-09-18 on the swirl-wedge fixture with one over-spun
parcel (`y = 0.2`, `up = 0.05`, `wp = 2`): the first segment swept ≈ 34° (34 sectors folded at
once) while the dual cell stayed the one located by `y = r cos θ = 0.200`, so 23 % of that
segment's deposit landed a row inward (+33 % / −37 % on the first two rows of the injection
column) — the sector-edge segments closed it the same day (`swirl-wedge-spin` RED-proves it on
the previous binary). The deposits' *frame* (Cartesian at the parcel's azimuth, not the meridian
plane) was O23's other half: fixed by `Lib_Equations::toMeridian` and gated by
`standard/swirl-wedge-deposit`; the sampling half is closed by E. This case's `euler1.tec` moved
with both fixes (its ID 1 leaves the plane; ≤ 3.6e-6 of scale for the sector-edge segments) and
its trajectory rows changed with each (347 → 349 → 337 for ID 1: the segment sequence sets the
SDIRK4 steps).
