# inj-edge-payload — a parcel carries the payload of the inlet cell that pinned it

**Purpose.** An injection point on an edge two inlet cells share lies in both cells. The parcel must start
in the cell that pinned it and carry that cell's loading, velocity and temperature, whatever cell a search
would return first. The case gives neighbouring inlet rows different payloads and
puts lattice points on the edges between them.

**Fixture.** `INPUT/solfile.tec` symlinks `standard/drag-stokes/INPUT/solfile.tec` (the shared box);
`INPUT/bc.txt` is written by `make_bc.py` from drag-stokes's: the inlet rows m = 1 and 3 (y 0–1 cm and
2–3 cm) carry κ_ρ = 0.68 and κ_v = 0.2, rows 2, 4 and 5 κ_ρ = 0.34 and κ_v = 0.1, and row 2 a per-cell
spacing of 2.5 mm (column 9), so that its advancing front reaches the edges y = 1 and 2 cm it shares with
rows 1 and 3. The global `ds = 10` (cm) puts one parcel at the centre of each cell of the other rows.
112 parcels. `phase.txt` and `properties.dat` symlink `common/`.

**Gates** (`check.py`):

| id | kind | assertion |
|---|---|---|
| P1 | target | every parcel's first trajectory row moves at its own row's κ_v·\|U\| (U = 10 m/s): a parcel starting at y in [1, 2] cm belongs to row 2 (1 m/s), one at y = 0.5 or 2.5 cm to row 1 or 3 (2 m/s) |
| P2 | target | the exit records carry the inlet's mass flow, Σ κ_ρ/(1 − κ_ρ) ρ U A over the inlet cells (ρ = 1.2 kg/m³, A = 1 cm²), to 1e-5 |
| W1 | control | every parcel leaves at x = 0.15 |
| W2 | control | some row-2 parcel starts on a row edge (12 do): the gate reaches the case it is for |

**RED.** When the search locates an injected parcel, the 9 row-2 parcels on the edge y = 2 cm start in row
3's cell: P1 9 of 112 at 2 m/s, P2 exit mass flow +2.778e-2 of the inlet's. **GREEN** P1 112 of 112,
P2 1.0e-7.
