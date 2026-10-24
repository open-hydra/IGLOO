# axis-200 — the axisymmetric AXIS face tagged `axisymmetric` (bcdef 200)

**Purpose.** The axis of an axisymmetric wedge case, when ATLAS tags it `axisymmetric`
rather than `sym`. `db-2daxi` tags face3 `sym` (bcdef 300) and injects both particles far
off-axis, so neither its tagging nor its geometry exercises the axis path.

**What fails without the axis branch (seen on hydra's `test/JPL-Lagrangian-20micron`).**
ATLAS emits 200 for the wedge k-faces and for the axis face. A 200 branch that rotates the
particle by ±`delthe` about x on every 200 face cannot handle the axis: rotation about x is
an isometry, so it leaves `hypot(y,z)` unchanged and cannot bring a particle that reached
the axis back inside the domain. The result is an exact period-2 cycle — `bcDef` rotates
`+delthe`, `axisymFold` rotates it back — with x, u and Tp byte-identical across the pair
and zero net displacement, ended only by the `nStall` displacement guard discarding the
particle.

On the JPL nozzle (635 particles/sweep, 190 axis-face cells tagged 200) such a solver loses
the six innermost particles every sweep, i.e. exactly the ones nearest the centreline,
where that case takes its Mach-number comparison against `JPL_chang.dat`.

**Fixture.** `db-2daxi`'s mesh, solfile and properties (symlinked), with a case-local
`bc.txt` that is db-2daxi's with the 200 face-3 entries retagged 300 → 200. Two DB-injected
particles on the same axial station: ID 1 at y = 1e-4 (inside the first radial cell, whose
outer edge is r ≈ 5.6e-3 there) with an inward radial velocity `vp = −1.0`, and ID 2 at
y = 0.55 as an off-axis control, so a failure is identifiably axis-specific. ID 1's `vp` is
load-bearing: the gas radial velocity vanishes on the axis, so a parcel in equilibrium with
the gas only approaches it (measured r_min ≈ 0.24·y0 at every injection radius tried) and
never touches the face: without `vp` the coverage assertion below is the one that fails, and
the others pass while exercising none of the axis path. ID 2's `up`/`vp` are its gas
velocity at injection (an absent component is filled with 0, not with the gas value).

**Gate (`check.py`) — behavioral, md5-free**, following db-2daxi's reasoning (trajectory
bytes drift by ~1 ULP across compiler/configure generations). It asserts: no give-up
message in `run_out.txt`; **coverage** — ID 1's minimum radius over its trajectory is
< 1e-5, so the case cannot pass vacuously if a future flow or mesh change stops carrying
it inward; both particles exit at x > 2.0; all trajectory rows finite with physical T and
dp; and the ord2 eulerian deposit (`out-file = e`, `euler1.tec`): **E1** total mass
retention Σρ_p V / Σṁ t = 1.000787 ± 1e-3, **E2** the near-axis band's deposit (y < 0.05)
over ID 1's own ṁ t = 0.999005 ± 1e-3 (a dual ghost reflected across the axis reads
0.931625 / 0.499505), **E3** the per-drop mass ρ_p/n_p = ρπd³/6 in every depositing cell to
1e-12 (worst measured 9.1e-15). Their floor is zero: `euler1.tec` is bit-identical across
repeats and across 1, 2 and 5 threads.

The coverage radius is read from the trajectory file's F12.6 columns, so it floors at ~5e-7 m and
reads as exactly 0 once the particle is inside that — it answers *did it reach the axis*,
not *how close*. Measured 0.0, also on a solver without the axis branch, and stable across
repeated 5-thread ctest runs, so the 1e-5 margin is set by the print format rather than by
`grazeStandoff` and will not drift if that constant is retuned.

**RED on a solver whose 200 branch only rotates** (ID 1: `no net progress ==> marking gone`,
exit x = −0.2609 instead of 2.0587; measured on a variant without ID 1's `vp` and with the
dual ghost reflected across the axis, whose spurious radial gas velocity carries ID 1 onto
it) **and GREEN on this one** (ID 1 reaches the axis, is projected back to r ≈ 1.01e-6 by
the grazing branch, and exits at x = 2.0587).

**Behaviour at the axis.** ID 1 meets the axis face once (one reflection event, stable
across `vp` = −0.5 … −5.0), and the impact takes the *grazing* sub-branch of the
reflection (v_n/|v| below `grazeFrac` = 2e-2), so the radial velocity is projected out and
the particle slides along the axis rather than bouncing; no reflection storm.

**Note on the graze standoff.** The graze places the particle `grazeOffset` inside the face:
`grazeStandoff` = 1e-6 m, capped at `grazeCellFrac` = 1e-2 of the cell's extent along the
face normal (measured: r 1.0e-8 → 1.01e-6). The axis is where meshes get radially thin.
Here the first cell reaches r ≈ 5.6e-3 and on the JPL mesh r ≈ 5.7e-4, so the cap is
inactive and the standoff is ≤ 0.2 % of the first cell; on a first radial cell thinner than
1e-4 m the cap sets it ([graze_standoff](../graze_standoff/INFO.md)).
