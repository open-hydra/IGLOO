# temp-relax-varcp-euler — the Eulerian energy moment of a varying-cp material

**Purpose.** Under a varying cp the particle state is the enthalpy, and the equivalent-Eulerian field must turn the
energy moment back into a temperature through the table. Model 1 deposits per cell ∫|v| dt and ∫h|v| dt;
`normalizeEUL` divides them and inverts the table (`comp_TfromTab`). With |v| constant along the path the field's
temperature of cell c is the table's T of the time-mean enthalpy over the parcel's residence:
T_p(c) = T(h̄_c), h̄_c = (1/(t_out − t_in)) ∫ h(T(t)) dt.

**Fixture.** `INPUT/` is `temp-relax-varcp`'s (symlink): the shared box gas (600 K, 10 m/s), a table on 250..800 K
with cp = 1000 + (T − 250), κ_v = 1 (each parcel coasts at u_g along its own column of cells, Re = 0, Nu = 2),
d = 16 µm, four inlet groups at 210, 300, 750 and 840 K (below, inside and above the table). `input.ini` is
`temp-relax-varcp`'s with `out-file = ALL` and `out-time = on`; `mollify = off`.

**Oracle** (`check.py`). T(t) is `temp-relax-varcp`'s closed form (its module: the table rebuilt from its recipe, the
segments marched from the injection state). On a table segment of slope s, h(T(t)) is an exponential in t, so
∫h dt is exact segment by segment; the window [t_in, t_out] of each cell is read from the trajectory rows (the row
entering the cell and the next one, the exit record for the last cell). T(h) is the check's own inverse of the
piecewise-linear table, the end segments extended.

**Gates.**
- W the witnesses: `Solving enthalpy equation`, `Field mollification OFF`; every first row the injection state.
- E1 |T_p(c) − T(h̄_c)| ≤ 1e-5 K on every crossed cell (the moments and the trajectory are ODE states under
  rtol 1e-11; the trajectory alone carries up to 8e-7 K, `temp-relax-varcp`'s V1).
- E2 T_p(c) inside [T_in, T_out] of the two rows bounding the cell, within 1e-6 K (the F12.6 print).
- E3 coverage: all 1500 cells of the box (the 25 inlet columns) have ρ_p > 0 and are checked.

**Measured.** E1 worst 2.04e-6 K (0.20 of the tolerance, a cell of the 840 K group); E2 0; E3 1500/1500, 25 parcels.

**RED.**
- The inversion skipped (`if (cpVariable)` → `if (.false.)` in `normalizeEUL`): T_p is the mean enthalpy itself,
  2.3e5–1.0e6 "K"; E1 and E2 fail on every cell. `temp-relax-varcp` stays green (it writes no Eulerian field).
- A constant-cp inversion (h/1250 in place of `comp_TfromTab`): E1 fails by up to 75.0 K (7.5e6 tolerances) and
  E2 on every cell.
