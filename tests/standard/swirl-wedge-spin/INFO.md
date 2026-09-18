# swirl-wedge-spin — sector-edge segments on the wedge (ledger O23, many-sector case)

**Purpose.** The third wedge twin. `swirl-wedge` pins the parcel dynamics and
`swirl-wedge-deposit` the deposits' frame; this one pins that an ODE segment on the
axisymmetric wedge never sweeps more than one sector, and therefore that a deposit lands
in the cell the parcel is actually in. One DB parcel at r = 0.2 (the first cell row, a
hair off the inner wall at 0.19) with `up = 0.05` and `wp = 2.0`, ten times the
co-rotating speed: its azimuth advances at ~10 rad/s at injection, so inside its first
(x, r) cell it would sweep ~34 one-degree sectors. Same fixture as the parent (`INPUT` is a
symlink), `out-file = ALL`, `mollify = off`. Registered `swirl-wedge-spin` (labels
`e2e;standard;drag;gas;axisym;euler`); runtime ~2 s.

**What was wrong.** Under `mesh2D` the point-in-cell test was azimuth-blind (a quad in
x–y), so a k-plane crossing never ended a segment: the segment ran until the (x, y) quad
was left, the fold at its end rotated back by nδ at once, and the whole segment's source
and euler deposit went to the dual cell located by `(x, y = r cos θ)` at the segment
start. For this parcel that cell sits a row inward of where it flies (true r 0.200 →
0.234 over the first segment while the locator stays at y = 0.200), and 23 % of the
segment's deposit landed there. **Fix (2026-09-18):** `isPointInsideCell` takes an
optional `sectorOut` (two dot products against the k-plane normals, `IGLOO_variables::
sectorNorm`, set at wedge detection); `solout` passes it, the interrupt fires like a
face crossing, the refinement pins the segment to the plane with the exact distance
(`IGLOO_bcBox::sectorDs`), and the outer loop calls `axisymFold(force=.true.)` before
any cell logic, so every fold rotates by exactly one sector. A sector crossing is neither
a cell crossing (no trajectory row, `nCross` untouched) nor residency (`Ncell` untouched,
so a spinning parcel does not trip `nMaxCell`). The run prints
`wedge sector folds: N (multi-sector: M)` at the end of `solve` (MPI-reduced).

**Oracle** (`check.py`, all derived). The parent's 6-state Cartesian Stokes ODE in the true
field `(U0, −Ωz, +Ωy)`, RK4 from the exact injection row, integrated to the outer wall
(r = 0.99, bisected) — 0.5540 s, unwrapped sweep 1.3725 rad = 79 sectors. Rows are matched
to the oracle at the **closest approach in (x, r)** rather than by x: rows sit at x- and
r-crossings, `u = 0.05` here, and matching by x alone turns the F12.6 floor of x into a
1e-5 s clock error (7.5e-5 in v_r at 20 m/s²). The projected residence uses the solver's
own dual-to-geo rule (`obj_block.f90::finalizeEUL`): mass in geo cell (i, j) =
Σ over its four corner nodes n of `subVol(i,j;n)/V_dual(n) · ṁ T_dual(n)`, with the
sub-quadrant volumes `(dx/2)·δ·(r_b² − r_a²)/2` from the tec nodes and `V_dual` the
tiling sum of the adjacent sub-quadrants; `T_dual(n)` bins the fine path by nearest node.

**Gates.**

| id | assertion | tolerance |
|---|---|---|
| S0 | `run_out.txt` carries the witness; M = 0; N ≥ 60; \|N − oracle sectors\| ≤ 2 | witness (measured 79 = 79) |
| S1 | every row inside the sector, \|θ\| ≤ δ/2 + 1e-6/r | in-sector |
| S2 | per row `w_p`, `v_r`, `u` vs the oracle at closest approach (r is absorbed by the matching) | 2e-5 / 2e-5 / 3e-6 — the floor `a·Δx/\|v\| = 20·5e-7/2 = 5e-6` (measured 3.8e-6 / 7.0e-6 / 5.8e-7) |
| S3 | mass per geo cell `ρ_p V` vs the projected oracle residence, every cell ≥ 1 % of the peak (97 cells; ≥ 30 required); `Σρ_pV = ṁ T_wall` | 2 % (measured 4.0e-3); 5e-3 (1.000044) |

**Proven RED (2026-09-18, the pre-fix binary `fc6c5f4`, same INI).** S0: no witness line.
S3: 11 of 97 cells off — the injection column's first rows at **+33 %** (r = 0.20) and
**−37 %** (r = 0.24): the first segment's mass sits a row inward. The mass total reads
1.000044 on both binaries: the tiling invariant (`axis-200` E1, `swirl-wedge-deposit` D5)
cannot see *where* the mass went, which is why this gate is per cell. S2 was already green
there (3.8e-6 / 4.5e-6 / 5.5e-7): the trajectory was never the problem — O23 was where the
sample was taken (closed by W-plan E) and where the deposit landed (closed here).

**Floor.** One parcel, one thread: deterministic per binary.

**Suite A/B of the fix** (`fc6c5f4` → this binary, full serial suite as sorted multisets):
byte-inert everywhere except the three off-plane fixtures — `swirl-wedge` (trajectories
identical; `source.tec` 3.4e-14 of scale; scatter cloud steps differ), `swirl-wedge-deposit`
(`euler1.tec` ρ_p ≤ 3e-7, v_p/w_p ≤ 1e-7 of scale; `source.tec` 2.3e-13) and `wedge-fold`
(ID 1 rows 349 → 337, `euler1.tec` ≤ 3.6e-6 of scale; gate green) — those parcels sweep
≤ 0.115 of a sector per segment, so bounding the segment to the sector barely moves them.
`db-2daxi` and `axis-200` (parcels on z = 0) are byte-identical: the sector test is inert
on the meridian plane.

**Cost.** Two dot products per accepted ODE step on a wedge, and one segment restart per
sector swept. The number of sectors a parcel sweeps is bounded by its total azimuth
(≈ π for one pass by the axis, whatever its closest approach), so the segment count is
≈ sweep/δ plus the cell crossings — 79 + 51 here; `wedge-fold`'s near-axis parcel went
from 349 to 337 rows in the same wall time. If a case ever spends its runtime in sector
restarts, the next step is capping `deltat` against the sector plane in `computeDeltat`.

**Not pinned here.** The `ord1` path (`gas-order = 1`): the geo containment call passes
`sectorOut` too, but no wedge case runs ord1. A parcel *injected* outside the sector is not
re-sectored at injection (the location searches stay azimuth-blind on purpose — an
off-sector point would match no cell); its first accepted step trips `sectorOut`, the
refinement finds no plane ahead (`sectorDs` = 0, the interior /nStep fallback) and the
forced fold rotates it home by `nint(θ/δ)` sectors at once — a multi-sector fold, so M ≠ 0
in the witness tells you an injection station was off-sector.
