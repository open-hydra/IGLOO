# axis-200 — the axisymmetric AXIS face tagged `axisymmetric` (bcdef 200)

**Purpose.** The axis of an axisymmetric wedge case, when ATLAS tags it `axisymmetric`
rather than `sym`. No other case covered this: `db-2daxi` — the only other 2Daxi case —
tags face3 `sym` (bcdef 300) and injects both particles far off-axis, so neither the
tagging nor the geometry ever exercised the axis path.

**The defect (diagnosed 2026-09-03, from `test/JPL-Lagrangian-20micron`).** `bcDef`'s 200
branch was written for the wedge k-faces and rotates the particle by ±`delthe` about x.
ATLAS emits 200 for the axis face too. Rotation about x is an isometry, so it leaves
`hypot(y,z)` unchanged and cannot bring a particle that reached the axis back inside the
domain. The result was an exact period-2 cycle — `bcDef` rotated `+delthe`, `axisymFold`
rotated it back — with x, u and Tp byte-identical across the pair and zero net
displacement, ended only by the `nStall` displacement guard discarding the particle.

Confirmed on the JPL nozzle (635 particles/sweep, 190 axis-face cells tagged 200): the
six innermost particles were lost every sweep, i.e. exactly the ones nearest the
centreline, where that case takes its Mach-number comparison against `JPL_chang.dat`.

**Fixture.** `db-2daxi`'s mesh, solfile and properties (symlinked), with a case-local
`bc.txt` that is db-2daxi's with the 200 face-3 entries retagged 300 → 200. Two DB-injected
particles on the same axial station: ID 1 at y = 1e-4 (inside the first radial cell, whose
outer edge is r ≈ 5.6e-3 there) which the converging flow carries onto the axis, and ID 2
at y = 0.55 as an off-axis control, so a failure is identifiably axis-specific.

**Gate (`check.py`) — behavioral, md5-free**, following db-2daxi's reasoning (trajectory
bytes drift by ~1 ULP across compiler/configure generations). It asserts: no give-up
message in `run_out.txt`; **coverage** — ID 1's minimum radius over its trajectory is
< 1e-5, so the case cannot pass vacuously if a future flow or mesh change stops carrying
it inward; both particles exit at x > 2.0; all trajectory rows finite with physical T and
dp.

That radius is read from the trajectory file's F12.6 columns, so it floors at ~5e-7 m and
reads as exactly 0 once the particle is inside that — it answers *did it reach the axis*,
not *how close*. Measured 0.0 both before and after the fix, and stable across repeated
5-thread ctest runs, so the 1e-5 margin is set by the print format rather than by
`grazeStandoff` and will not drift if that constant is retuned.

**Verified RED at HEAD** (ID 1: `no net progress ==> marking gone`, exit x = −0.2609
instead of 2.0587) **and GREEN with the fix** (ID 1 reaches the axis, is projected back to
r ≈ 1.01e-6 by the grazing branch, and exits at x = 2.0587).

**Behaviour after the fix.** Every axis impact takes the *grazing* sub-branch of the
reflection — measured v_n/|v| ≈ 3e-4 against `grazeFrac` = 2e-2 — so the radial velocity is
projected out and the particle slides along the axis rather than bouncing. 31 such events
over the full traverse in this case; no reflection storm.

**Note on `grazeStandoff`.** The graze places the particle `grazeStandoff` = 1e-6 m inside
the face (measured: r 1.0e-8 → 1.01e-6). That constant is absolute, and the axis is where
meshes get radially thin. Here the first cell reaches r ≈ 5.6e-3 and on the JPL mesh
r ≈ 5.7e-4, so the standoff is ≤ 0.2 % of the first cell — but on a mesh whose first radial
cell is thinner than ~1e-6 m it would place the particle beyond cell j=1. Unguarded by
design (it is pre-existing behaviour for every 300 face); revisit if such a mesh appears.
