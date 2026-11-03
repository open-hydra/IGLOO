# inlet-402-403 — the particle inlets 402 and 403

**Purpose.** Besides the 401 inlet (number density, speed and temperature as factors on the gas), the boundary
file has two particle inlets whose payload sets the mass flux: 402 gives the particle speed and temperature
absolutely, 403 the speed as a factor on the local gas speed and the temperature absolutely. In both the payload's
first column is the particle mass flux gp [kg/(m² s)], so a cell injects gp·A, A the cell's face area, shared by
the N_inj parcels of the cell (`initializePart`, `obj_particles.f90`; the reader adds gp·A to the cell's
particle flow, `IO.f90`). No other fixture carries a 402 or 403 record.

**Fixture** (`python3 make_fixture.py`, `bc.txt` committed; `solfile.tec` is `drag-stokes`'s, symlinked;
`phase.txt` and `properties.dat` from `common/`): `drag-stokes`'s box and gas (uniform, u_g = 10 m/s along +x,
T_g = 600 K) with the 25 inlet records recoded: rows j = 1–2 are 402 (gp = 5.0, |v| = 4.0 m/s, Tp = 450 K),
rows j = 3–5 are 403 (gp = 2.0, kV = 0.5, Tp = 350 K), direction `normal` (the gas direction), r_p =
5.945e-6 m, Dirac; `ds = 10` cm gives one parcel per cell (N_inj = 1). `make_fixture.py`'s `recode_inlet` is the
writer `stuck-cell` imports. `input.ini` is `drag-stokes`'s with `out-time = on`.

**Gates** (`check.py`; the oracle is the inputs).
- I1 each parcel's first row: (u, v, w) = (4, 0, 0) m/s and T = 450 K on the 402 rows, (5, 0, 0) m/s and
  T = 350 K on the 403 rows, x = 0, to the print.
- I2 each exit record's mass flow (read by its header name) = gp·A, A = d_y d_z = 1e-4 m²: 5.0e-4 kg/s on the
  402 rows, 2.0e-4 kg/s on the 403 rows, to the E13.6 print.
- I3 `drag-stokes`'s Stokes closed form x(u), anchored at each parcel's first row (u₀ = 4 or 5 m/s against
  u_g = 10 m/s), with its tolerance; v = w = 0 on every row. Measured worst 0.012 of the tolerance.
- I4 25 exits, every one at x = 0.15.

**RED.**
- The velocity statements of 402 and 403 exchanged: 402 parcels start at 40 m/s, 403 parcels at 0.5 m/s; I1
  fails on 25 of 25 parcels.
- 402's temperature taken as a factor on the gas (`kT*T_g`): the parcels start at 2.7e5 K, past the F12.6
  print; the run's output is refused as non-finite (`check_finite.py`) and I1 cannot read the rows.
- The face area dropped from 402's mass flow: the 10 parcels of the 402 rows carry 5.0 kg/s; I2 fails.

The division by N_inj is not seen here (one parcel per cell).
