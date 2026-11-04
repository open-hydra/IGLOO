# inj-edge — injection points on the edges of the inlet cells

**Purpose.** With a spacing below the cell size, the 3D advancing front puts injection points on the edges
of the inlet cells, the edges of the inlet face included, where rounding leaves a point a few 1e-18 m
outside the box. The pinning test accepts such a point in its cell, within the cell's containment
tolerance. The case gates that every parcel starts in the cell that pinned it, is found again at each cell
crossing, and leaves at the outflow.

**Fixture.** `INPUT/` symlinks `standard/drag-stokes/INPUT` (the shared box: 60 × 5 × 5 cells,
x 0–0.15 m, y and z 0–0.05 m, uniform gas at U = 10 m/s); `input.ini` is drag-stokes's with
`ds = 0.25` (cm) in place of 10: 472 parcels from the 25 inlet cells.

**Gates** (`check.py`):

| id | kind | assertion |
|---|---|---|
| E1 | target | the exit records number the run log's particle count and every one is at x = 0.15 |
| E2 | target | no parcel is reported in no cell, out of the domain at injection, stuck or without progress, and `run_err.txt` is empty |
| E3 | control | more parcels than the 25 inlet cells: the front placed a lattice, not the cell centres |

**RED.** Without the rule (the search alone locates an injected parcel, with exact bounding boxes): 13 of the
points lie outside every cell's bounding box, the parcel keeps block index 0 and the gas lookup reads
block 0: SIGSEGV, 5 of 472 exit records written (ctest: `solver died on signal`). With the no-cell guard
alone: no crash, `[WARNING] 13 parcel(s) whose injection point no cell holds: not integrated`, 459 exits.
With the parcel in its own cell but the exact boxes only: the 13 parcels start and are lost at their first
cell crossing (12 exit records at inf/nan through the wall y = 0.05, one parcel stuck). **GREEN** 472 of
472 exits at x = 0.15, no warning.

**Census** (the rule's decisions over this run, counted by a scratch probe outside the suite): 48 parcels
start in a cell the search would not have returned, 13 that it found in no cell and 35 that it put in the neighbouring inlet cell across a shared edge
(the inlet payload is the same in every cell here: `inj-edge-payload` gates that case); the search with
widened boxes finds the cell 767 times.
