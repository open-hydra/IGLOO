# inj-pack — several parcels per inlet cell, packed at the spacing ds

**Purpose.** With `[IGLOO-BC] ds` below the size of an inlet cell, the boundary placement (`pin_particles_bc_ds`,
`initialization.f90`) packs each cell: on a 3D face an advancing front builds a hexagonal lattice at the spacing
R = ds, seeded at the first inlet cell's centre, and refuses a candidate closer than R/2 to a placed point
(`check_distance`). A cell's mass flow is shared by its N_inj parcels — for a 402 cell gp·A/N_inj. No other
case injects more than one parcel per cell from a boundary file.

**Fixture** (`python3 make_fixture.py`, `bc.txt` committed; `solfile.tec` is `drag-stokes`'s, symlinked; `phase.txt`
and `properties.dat` from `common/`): `drag-stokes`'s box and gas with inlet rows j = 1–4 left 401 (κ_ρ = 0.34)
and row 5 recoded 402 (gp = 5.0 kg/(m² s), |v| = 4.0 m/s, Tp = 450 K; `inlet-402-403`'s `recode_inlet`).
`input.ini` is `drag-stokes`'s with `ds = 0.27` (cm) and `out-time = on`: 407 parcels, 14 to 20 per cell.

**Gates** (`check.py`; its `run()` serves `inj-dscell` and `inj-dsdegen` too; the oracle is the inputs and
conservation).
- J0 every pinned point (a parcel's first row) lies more than 2e-6 m off the face-cell edges, so the cell of a
  parcel is the cell holding its first row, and every parcel is integrated to an exit at x = 0.15.
- J1 every first row on the inlet plane x = 0, at t = 0.
- J2 the first rows within a cell at least R/2 apart, J2b at least R (the lattice, on this uniform face), less
  1.5e-6 m for two printed points; measured 2.699767e-3 m against R = 2.7e-3 m.
- J3 per inlet cell the parcels' mass flows (exit records, by header name) sum to the cell's flow —
  κ_ρ/(1 − κ_ρ)·ρ_g U A for 401, gp·A for 402 — to the E13.6 print, whatever N_inj is; and the 402 cells carry
  16 to 20 parcels each.

**RED.** 402's mass flow not divided among the cell's parcels (`initializePart`, `/self%Ninj` dropped): the five
402 cells carry 16 to 20 times their flow; J3 fails.

**The spacing is 0.27 cm, not 0.25 cm.** A lattice point that lands on a face-cell edge is placed by one
containment test (`isPointInsideHexahedron12`, a face tolerance of 1e-15 m) and then sought at the start of the
integration by another (`searchInBlock`, whose bounding-box pre-check has no tolerance); the two disagree on such a
point. At ds = 0.25 cm the front reaches the face's outer edge y = 0 (rounding leaves the point at y = −2.2e-18 m):
no cell holds it, the parcel's cell index stays 0 and the run stops with a segmentation fault
(`Lib_Integration.f90`, `gasblock(b)` with b = 0). At an interior edge the parcel is moved to the neighbouring cell
and takes that cell's payload (none, next to a non-inlet cell: a segmentation fault). ds = 0.2, 0.25, 0.3 and
0.5 cm stop that way on this box, 0.27 and 0.33 cm do not; `mdot-max` sets R = d_cell/N, which reaches the edges
for every N. Reproduce from this directory:

    mkdir -p /tmp/inj-pack-025/OUTPUT && ln -sfn "$PWD/INPUT" /tmp/inj-pack-025/INPUT
    sed 's/^ds = 0.27$/ds = 0.25/' input.ini > /tmp/inj-pack-025/input.ini
    (cd /tmp/inj-pack-025 && OMP_NUM_THREADS=1 "$IGLOO")    # $IGLOO the executable: forrtl severe (174), SIGSEGV
