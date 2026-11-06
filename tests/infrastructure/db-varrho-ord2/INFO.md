# INFO — infrastructure/db-varrho-ord2 (database streams sized at their injection gas under gas-order 2, e2e)

## Reference
Infrastructure, no external reference: the oracle is the gas file, the property table and the injection state
the run prints. The scatter cloud's quantum is `dNscat = Ṅ τ_ref/(fsample-traj · nStreams)`, `τ_ref = L_ref/V_ref`
(`docs/user/using.md`), sized by a serial pre-pass over throwaway copies of the parcels.

## Setup
`db-varrho`'s case (`INPUT/` linked to its files: the box with the gas temperature 600 + 1000 x at the cell centres,
10 m/s along x, `temp-relax-varrho`'s density table, four database streams at the centres of cells 5, 21, 37 and
53 along x without `temp0` and without a velocity) with `gas-order = 2`. The gas is then held at the cell centres,
the nodes of the dual cells; a stream at a cell centre sits on a node shared by eight dual cells, and `integrate`
injects it at the gas of the first node of the one it finds: here the centre of the previous cell along x
(608.75, 648.75, 688.75, 728.75 K; 10 m/s).

## What it verifies
`check.py`:
- **G-DB1** the logged `dNscat` against `Ṅ (L_ref/V_ref)/(100 × 4)`: Ṅ = Σ ṁ/(π/6 · ρ(Tᵢ) · d³), Tᵢ and the speed
  the injection state `integrate` gave stream i (its first trajectory row), V_ref their mean, ρ from
  `INPUT/properties.dat`, L_ref the mesh's bounding-box diagonal (0.165831 m). Tolerance: the `ES10.3` half-ULP
  plus 1e-6 relative (0.507 on 6875.455).
- **G-DB2** each first row's T is the temperature of a cell centre of the dual cell around its stream (its own
  cell or a neighbour along x, from `INPUT/solfile.tec`) within the `F12.6` half-ULP, its speed 10 m/s, its m_p =
  π/6 · ρ(T) · d³ within the `E13.6` half-ULP plus 1e-9 relative.
- **G-DB3** the scatter file holds markers of each of the four parcels.

## Pass criterion
G-DB1–G-DB3 hold. Measured: `dNscat = 6.875E+03` (expected 6875.455); first rows at 608.75, 648.75, 688.75,
728.75 K.

## Falsification
- A pre-pass that leaves a database stream without `temp0` at T = 0 on the density table: `dNscat = 0.000E+00`, no
  marker (G-DB1, G-DB3).
- A pre-pass that takes the dual cell's last node under `gas-order = 2` (`integrate` takes its first):
  `6.879E+03`, 7 times the tolerance (G-DB1).
G-DB2 passes under both: `integrate`'s injection is not what they change.
