# vie-plait-dense — 250 strands at gas-order 2 in vie-plait's compressive field

**Purpose.** `vie-plait` follows six assigned-position parcels through the 3D trilinear gas sample. This case
samples the same field with 250 strands injected from the inlet face, one per inlet cell, so that the sampler
runs over the whole cross-section (the boundary cells next to the side faces included) and the strands cross each
other at the axis; it also holds the 3D gas-order-2 source field — the deposits reduced from the dual cells to the
geometry cells and smoothed — to the momentum and kinetic energy the strands hand to the gas.

**Fixture.** `INPUT/solfile.tec` and `INPUT/properties.dat` link `vie-plait`'s (u_g = 0.2 m/s,
v_g = −ε(y − 1), ε = 1 s⁻¹, 120 × 50 × 5 cells over 6 × 2 × 0.05 m; ρ_p = 1620 kg/m³), `INPUT/phase.txt` links
`common/`. `INPUT/bc.txt` is `python3 ../../tools/make_vie_case.py INPUT --kv 1.0 --rp 5.0e-4` (the generated
`solfile.tec` is byte-identical to `vie-plait`'s and is replaced by the link): the inlet face's 250 cells
carry κ_v = 1 along the cell's gas velocity, κ_T = 1 and r_p = 0.5 mm (d = 1 mm, τ = 5 s, St = 5).
`input.ini` is `vie-plait`'s without the assigned-position keys, with `ds = 10` (cm, above the 4 cm × 1 cm
cells: one parcel at each cell's centre).

Each parcel starts at its cell's gas velocity: U = u_g, so x = u_g t is the clock, and V = −εY_a (Y = y − 1) —
not the V = 0 start of `vie-plait`'s parcels. The damped oscillator Y'' + Y'/τ + (ε/τ)Y = 0 then swings each
strand to −1.48 Y_a: the 80 strands with |Y_a| ≥ 0.70 (rows 1–8 and 43–50 of the face) reach the side face
y = 0 or 2 (code 400) before the outlet, the other 170 leave at x = 6.

**Gates** (`check.py`).
- O1 `vie-plait`'s anchored closed form, with its tolerance, on every row of every strand: 250 of 250 strands,
  41 963 rows, worst 0.712 of the tolerance; u = u_g, z and T = T_g unchanged on every row; one strand from the
  centre of each inlet cell, starting at the cell's gas velocity; the strands that leave through a side face are
  exactly those whose closed form reaches |Y| = 1 before x = 6 (every peak |Y| at least 0.023 from 1), the others
  leave at the outlet, every exit with code 400 and every last row at its exit point; the 190 strands with
  |Y_a| ≥ 0.25 cross the axis.
- O2 the source totals against the strands' telescoped flows. `source.tec` holds each cell's extensive totals
  (N, W; the smoother works on densities and gives the totals back, so their sum is conserved): Σ_cells F =
  Σ_p ṁ_p(v_first − v_last) per component and Σ_cells E = Σ_p ṁ_p(|v_first|² − |v_last|²)/2 (T stays at T_g:
  no enthalpy term), within the print floor of the rows (F12.6 velocities, E13.6 mass flows). F_x is zero (no
  axial slip) and F_y cancels between mirror strands, so the energy total carries the check: measured
  1.779984E-03 W against the strands' 1.779982E-03 W, 0.22 of the 8.1E-09 W floor.

**RED.**
- The reduction of the source to the geometry cells dropped (`finalizeSRC`'s gas-order-2 branch): the deposits
  stay on the dual cells and the run stops in the smoother (`binomial_smooth: work scratch shape mismatch`,
  exit 128) without writing `source.tec`; O2 fails.
- Each dual cell shared among its geometry cells in equal eighths instead of by octant volume
  (`finalizeSRC`'s weight): the boundary dual cells lose part of their deposits and the energy total falls by
  10.7 % (1.589E-03 W, 2.3E+04 floors); F_x and F_y stay at zero by symmetry and pass.
- The sampler's Newton stopped at an absolute residual of 1e-5 m moves this case's rows by less than the
  print: O1's worst stays 0.712. The stopping rule is gated by `test_gas_reconstruction` (E7, E8) and
  `test_drag` (E6a, E6b).

All 250 parcels leave the domain (170 at the outlet, 80 through the side faces).
