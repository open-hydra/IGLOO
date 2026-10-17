# periodic-2d-ord2 — translational periodic faces (201) at gas-order 2

**Purpose.** At gas-order 2 the gas is sampled on the dual mesh, whose nodes are the cell centres plus a ghost
ring. The ring is first filled by linear extrapolation; at a 201 face `fillGhostPartners` then overwrites each
ghost with the gas of the partner side, at the ghost's position translated by the periodic offset (partner face
centre − own face centre). The sampled field is then the periodic extension of the cell-centre field, and a
parcel crossing the periodic pair sees a continuous gas. `periodic-y` and `periodic-2d` run gas-order 1 in a
uniform gas, where the ring is never read; this case is the end-to-end gate of the periodic ghosts.

**Fixture** (`python3 ../../tools/make_box_case.py INPUT --nz 1 --ny 25 --periodic-y --kv 1.0 --kt 0.5
--wave-y 100`, output committed; `phase.txt` and `properties.dat` from `common/`, symlinked): `periodic-2d`'s
one-layer box (60 × 25 × 1 cells, 0.15 × 0.05 × 0.05 m, faces 3/4 a 201 pair, every other face 400 but the
inlet, face 1) whose gas temperature carries a wave across the periodic direction, T_c(j) = 600 + 100 cos(2π
y_c(j)/L_y) at the cell centres (uniform in x; every other gas variable uniform). Parcels (d = 11.89 µm,
ρ_p = 2950, c_p = 1250, one per inlet cell) start at rest in y with u = u_g and T = 0.5 T_c(j). `input.ini` is
`periodic-2d`'s (Stokes, Ranz-Marshall, g_y = −8000 m/s², about 2.8 wraps of L_y per parcel) with
`gas-order = 2` and `out-time = on`.

**Oracle** (`check.py`; it never reads a solver quantity to build the reference). The path is `periodic-y`'s
closed form (the drag is Stokes, so v_y and y do not depend on T). T(t) is a classical RK4 of
dT/dt = K(t)(T_g(y(t)) − T), K = 6 Nu k_g/(ρ_p c_p d²), Nu = 2 + 0.6 Re^½ Pr^⅓, Re = ρ_g |v_y(t)| d/μ,
Pr = μγR/((γ − 1)k_g), from the injection state, with T_g(y) piecewise linear through the node lines y_c(j),
j = 0 … N_y + 1, and the periodic extension T(0) = T(N_y), T(N_y + 1) = T(1). The steps are split where y(t)
crosses a node line (the kinks of T_g), graded near t = 0 (Re^½ grows as t^½ there), and the integral is checked
against a run at half the steps (≤ 1e-9 K).

**Gates.**
- P0 each parcel's first row is the injection state (x = 0, y = y_c(j), v = 0, T = 0.5 T_c(j)) to the print.
- P1 `periodic-y`'s gates with its tolerances: v_y on the anchored closed form, y on the unwrapped closed form
  modulo L_y, u = u_g, w = 0.
- P2 every row with |T − T_g(y)| > 1 K: |T − T_or(t)| ≤ tol_T = δ_T + |dT/dt| δ_t + ε_int, δ_T = 5e-7 K (F12.6),
  δ_t = 5e-9 t (the ES16.8 time), ε_int = 1e-6 K. The worst residual over every row is printed as well.
- P3 at least two wraps per verified parcel, 20 of the 25 verified.
- P4 25 exits, every one at x = 0.15.
- P5 non-vacuity: for every verified parcel some row separates T_or from the same RK4 through the ring's linear
  fill (T(0) = 2T(1) − T(2), T(N_y + 1) = 2T(N_y) − T(N_y − 1)) by more than 50 tol_T.

**Measured.** P2 worst residual 4.33e-6 K against tol_T = 5.59e-6 K (0.78 of it; the time print dominates the
tolerance there, |dT/dt| ≈ 8e4 K/s); the same over every row. Rows whose time-print term is below 1e-8 K: worst
residual 4.87e-7 K, at `relative-tol = absolute-tol` 1e-11 and at 1e-13 alike — the integration error is below
the T print, so ε_int = 1e-6 K is not binding. The oracle against itself at half the steps: 4.4e-10 K. P5: the
two fills separate by at least 1.0e5 tol_T on every parcel (the ghost line at y = −L_y/50 holds 699.2115 K
periodic, 705.4453 K by the linear fill). The parcels injected next to the corners of the ring (rows 1 and 25)
pass from their first row: the corner ghosts take the partner block's x-ghosts, periodic here because the
field is uniform in x.

**RED.**
- 201 dropped from the ghost links (`buildGhostLinks`, `case (101, 103, 201)`): the 201 ghosts keep the
  linear fill; P2 fails on all 25 parcels, the per-parcel worst residual 0.290 K to 0.660 K (up to 3.1e5
  tol_T), while P1 (v and y), the wrap count and P4 stay as on the correct code.
- `handleGeoEvent` dropped after the exit's `updateCell` (`Lib_Integration`, gas-order 2): the gas stencil
  is not re-sought after the periodic transport; every parcel stops at its first wrap (`no net progress`),
  no parcel reaches x = 0.15: P3 (at most one wrap) and P4 (0 of 25 exits at x = 0.15) fail.
