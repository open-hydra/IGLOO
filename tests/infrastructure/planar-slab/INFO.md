# planar-slab — a single-layer (Nk = 1) box is a planar 2D case, not a wedge

**Purpose.** The axisymmetric-wedge detector (`allocation.f90`) must not fire on a planar slab
of finite thickness. No other case covered this: every 2D case in the suite is a genuine wedge
(`db-2daxi`, `axis-200`, `wedge-fold`, `swirl-wedge*`), and the standard boxes have Nk > 1.

**The defect.** The detector measured the angular span of the two k-planes about the axis at
the OUTERMOST node only. A slab of thickness t also spans `atan(t/r)` there (0.04995840 rad on
this box), so every finite-thickness slab was taken for a wedge. Consequences, measured on
hydra's `test/evap-box-*` (this mesh): the parcels were folded about the x axis on their first
step and the fold fought the reflecting k-faces (dozens of `stuck in cell` events per run);
after the 2026-09-17 sector-centring check the run refused the mesh at setup (`wedge sector
must be centred on the azimuth origin`), because a slab's k-planes sit at 0 and +span, never
at ±span/2.

**The fix.** A wedge spans the same angle at every radius; a slab's span decays as 1/r. The
detector now also measures the innermost node clearly off the axis (r > 1e-3 r_max, see `wedge-axis-row`) and takes the planar 2D path when the
two spans disagree (relative 1e-6). On this box: 0.04995840 rad outermost, 1.57079633 rad
innermost.

**Fixture.** `INPUT/solfile.tec`: uniform gas (u = 5 m/s, v = w = 0, T = 350 K, p = 1 atm) on
the evap-box-ta box, 20 × 20 × 1 cells over 0.1 × 0.1 × 0.005 m, written by `make_fixture.py`
(kept for provenance; note `import_gas` binds `GAM`, not MOSE's `g`). `INPUT/bc.txt`:
evap-box-ta's ATLAS `part-bc.txt` (six symmetry faces, 300) with face 2 retagged 400 by
`retag_bc.py`, so the parcels leave the box. Phase and properties are the common ones. Two
parcels injected at the gas velocity in the mid-plane z = 0.0025 at y = 0.025 and 0.075: two
different radii from the x axis, which is exactly what tells a slab from a wedge.

**Gate (`check.py`).** Log reports `Planar slab (parallel k-planes)` and no wedge; no give-up;
both parcels exit at x ≥ 0.0999; y and z stay at their injection values to 1e-9 on every
trajectory row (a fold rotates z; a mis-oriented k-face reflection moves y).

**Red/green.** Pre-fix binary: exit 128, `Axisymmetric wedge: delthe = 0.04995840`, refused.
Post-fix: 21 rows per parcel, y and z bit-identical to injection, both exits at x = 0.1.
