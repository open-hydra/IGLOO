# inj-pack — several parcels per inlet cell, packed at the spacing ds

**Purpose.** With `[IGLOO-BC] ds` below the size of an inlet cell, the boundary placement (`pin_particles_bc_ds`,
`initialization.f90`) packs each cell: on a 3D face an advancing front builds a hexagonal lattice at the spacing
R = ds, seeded at the first inlet cell's centre, and refuses a candidate closer than R/2 to a placed point
(`check_distance`). A cell's mass flow is shared by its N_inj parcels — for a 402 cell gp·A/N_inj. At this spacing
lattice points fall on the face-cell edges, the face's outer edges included: a point on an edge two cells share lies
in both, is pinned by one of them and must start in that cell with that cell's payload. The 402 row borders the 401
rows, so the payloads differ across the edge y = 4 cm.

**Fixture** (`python3 make_fixture.py`, `bc.txt` committed; `solfile.tec` is `drag-stokes`'s, symlinked; `phase.txt`
and `properties.dat` from `common/`): `drag-stokes`'s box and gas with inlet rows j = 1–4 left 401 (κ_ρ = 0.34)
and row 5 recoded 402 (gp = 5.0 kg/(m² s), |v| = 4.0 m/s, Tp = 450 K; `inlet-402-403`'s `recode_inlet`).
`input.ini` is `drag-stokes`'s with `ds = 0.25` (cm) and `out-time = on`: 472 parcels, 16 to 23 per cell; 72 of the
pinned points lie on face-cell edges, 48 of them on edges two cells share.

**Gates** (`check.py`; its `run()` serves `inj-dscell`, `inj-dsdegen` and `inj-mdotmax` too; the oracle is the
inputs and conservation).
- J0 every parcel is integrated to an exit at x = 0.15, and the pinned points (a parcel's first row) on the
  face-cell edges exist and each carries the payload of a cell that holds it: the velocity and temperature of the
  cell's code and its share F_c/N_c of the cell's flow, N_c read from the cell's interior parcels; the edge parcels
  are matched to the cells so that every cell receives exactly N_c. A case whose points all lie off the edges
  (`inj-dscell`, `inj-dsdegen`) asserts that instead, by more than 2e-6 m, and a parcel's cell is the cell holding
  its first row.
- J1 every first row on the inlet plane x = 0, at t = 0.
- J2 the first rows within a cell at least R/2 apart, J2b at least R (the lattice, on this uniform face), less
  1.5e-6 m for two printed points; measured 2.499945e-3 m against R = 2.5e-3 m.
- J3 per inlet cell the parcels' mass flows (exit records, by header name) sum to the cell's flow —
  κ_ρ/(1 − κ_ρ)·ρ_g U A for 401, gp·A for 402 — to the E13.6 print, whatever N_inj is; and the 402 cells carry
  16 to 20 parcels each.

**RED.**
- 402's mass flow not divided among the cell's parcels (`initializePart`, `/self%Ninj` dropped): the five 402
  cells carry 16 to 20 times their flow; J3 fails on the five cells and J0 on their counts.
- The injected parcel located by the search instead of starting in the cell that pinned it (`integrate`,
  `updateCell` called whether or not the pinning cell holds the point): 6 of the 12 row-4 parcels on the edge
  y = 4 cm start in a row-5 cell, at its velocity and temperature (4 m/s, 450 K) and a share of its flow; J0 fails
  (no matching fills the cells' counts), J3 on three row-4 cells (6.00e-4 to 6.12e-4 kg/s against 6.181818e-4),
  the exit mass flow 0.26 % low.
- Without that rule, the guard against a point no cell holds and the widened containment search: a point on the
  face's outer edge lies in no cell's bounding box, the parcel keeps block index 0 and the run stops with a
  segmentation fault, 5 of 472 exit records written; J0 and J3 fail.
