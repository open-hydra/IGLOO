# periodic-2d — translational periodic faces (201) on a one-layer box

**Purpose.** `periodic-y` on a one-layer (2D) mesh. On a mesh of one k-layer ATLAS writes a 201 (and
101) connection line with nine integers — block, i, j, k (= 1), the partner face and the four
orientation integers — and `read_cdp_bc_file` must take the partner face from the fifth of them: the
periodic transport moves a particle by `T = partner face centre − exit face centre`, and the partner
face selects that centre. A line of eight integers (block, i, j, partner face, orientation) is read too.

**Fixture** (`python3 ../../tools/make_box_case.py INPUT --nz 1 --ny 25 --periodic-y --kv 1.0`, output
committed; `phase.txt` and `properties.dat` from `common/`, symlinked): `drag-stokes`'s uniform gas on
a planar single-layer box of 60 × 25 × 1 cells (0.15 × 0.05 × 0.05 m; the 2D path, no wedge), face 1
the inlet (25 cells), faces 3/4 (y = 0 and y = 0.05) a 201 pair whose lines carry nine integers, every
other face 400. `input.ini` is `periodic-y`'s: g_y = −8000 m/s², so each particle drifts about 2.8 box
heights while it crosses the box.

**Gates** (`check.py` runs `periodic-y`'s oracle, unchanged, on this case's output): v_y(x) on the
anchored closed form across every wrap, y(x) on the unwrapped closed form modulo L_y, at least two wraps
per particle (20 of the 25 must verify), every particle out at x = 0.15, and the u, w, T decoupling
guards; and S-P, the source smoothed across the periodic pair (2D sweeps in x and y): the case against
its unsmoothed sibling smoothed by `check.py` (measured 1.5e-14 of the field scale, tolerance 1e-12), and
the same pair with only face 1's first row injecting (measured 8.6e-15; periodic and zero-flux smoothing
differ there by 0.29).

**RED.** A reader that takes four index integers on a one-layer mesh stores the partner's k (1) as the
partner face, so the transport takes the centre of a face-1 facelet (or of no facelet at all): no
particle wraps as it should, 18 of the 25 cycle until they hit the outer iteration cap (500000
segments) near the inlet, and the run still had not ended after six minutes (the ctest timeout is 300 s).
The stored face itself is gated by `test_bc_families` F7. S-P is RED with zero flux at the 201 faces:
6.1e-8 of the scale on the case (its deposits are uniform in y), 0.29 on the first-row run.
