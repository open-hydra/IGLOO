# INFO — infrastructure/db-varrho (database streams sized at their injection cell's gas, e2e)

## Reference
Infrastructure, no external reference: the oracle is the gas file and the property table. The scatter cloud's
quantum is `dNscat = Ṅ τ_ref/(fsample-traj · nStreams)`, `τ_ref = L_ref/V_ref` (`docs/user/using.md`), sized
by a serial pre-pass over throwaway copies of the parcels before they are integrated.

## Setup
`make_box_case.py`'s box (60 × 5 × 5 cells, 0.15 × 0.05 × 0.05 m) with `--gradient 1000`: the gas temperature
600 + 1000 x at the cell centres, 10 m/s along x (`INPUT/solfile.tec`; `bc.txt` and `phase.txt` linked to
`db-injection`'s). `INPUT/properties.dat` is `temp-relax-varrho`'s table: ρ = 2950 − 0.5 (T − 250) kg/m³ and a
varying c_p on 250–800 K. Four database streams (`[IGLOO-BC] x = 0.01125 0.05125 0.09125 0.13125`, y = z = 0.025,
the centres of cells 5, 21, 37 and 53 along x), d = 11.89 µm, ṁ = 1e-4 kg/s, no `temp0` and no `up`/`vp`/`wp`:
each is injected at the gas temperature and velocity of the cell holding it (611.25, 651.25, 691.25, 731.25 K;
10 m/s). Stokes drag, Ranz–Marshall heat, `gas-order = 1`.

## What it verifies
`check.py`:
- **G-DB1** the logged `dNscat` against `Ṅ (L_ref/V_ref)/(100 × 4)`: Ṅ = Σ ṁ/(π/6 · ρ(Tᵢ) · d³), Tᵢ the gas
  temperature of the cell holding stream i and V_ref the mean gas speed of those cells, both read from
  `INPUT/solfile.tec` (the cell found from the mesh's nodes), ρ from `INPUT/properties.dat`, L_ref the mesh's
  bounding-box diagonal (0.165831 m). Tolerance: the `ES10.3` half-ULP plus 1e-6 relative (0.507 on 6878.59).
- **G-DB2** each parcel's first trajectory row: m_p = π/6 · ρ(Tᵢ) · d³ within its `E13.6` half-ULP plus 1e-9
  relative, T = Tᵢ within the `F12.6` half-ULP. It holds `integrate`'s injection state to the state G-DB1 assumes.
- **G-DB3** the scatter file holds markers of each of the four parcels.

## Pass criterion
G-DB1–G-DB3 hold. Measured: `dNscat = 6.879E+03` (expected 6878.593); first-row masses 2.437400E-12,
2.419790E-12, 2.402190E-12, 2.384590E-12; 189 markers, of all four parcels.

## Falsification
- A pre-pass that leaves a database stream without `temp0` at T = 0 on the density table (no number flow) and
  without a velocity out of V_ref: `dNscat = 0.000E+00` (1.4e4 times the tolerance), no marker (G-DB1, G-DB3).
- The temperature taken from the injection cell and the velocity not (V_ref = 1 m/s): `6.879E+04` (1.2e5 times the
  tolerance, G-DB1).
- Every stream given cell (1,1,1)'s gas: `6.791E+03` (173 times the tolerance, G-DB1).
G-DB2 passes under all three: `integrate`'s injection is not what they change.
