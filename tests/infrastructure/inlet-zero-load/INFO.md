# inlet-zero-load — an inlet cell that injects nothing pins no parcel

**Purpose.** ATLAS writes a 401 of krho 0 for a phase without inlet keys at an inflow grain (its mixture rule over the
dispersed phases; `grain-bc = inflow`), and a cell outside a krho profile carries 0 too. Such a cell injects nothing,
yet the pinning places parcels by geometry alone (`isInj` reads the code): without a filter each one would carry zero
mass from the inlet to its exit. `dropZeroFlow` (`obj_IGLOO.f90`), called once every parcel of a family is pinned,
removes the parcels of the cells whose column 1 is zero for that family (401: krho/(1 − krho_tot); 402/403: gp),
renumbers the rest 1…n and prints `[IGLOO] family <f>: <n> parcel(s) of inlet cells with a zero loading (field 1 = 0)
not injected`. The draws of the parcels kept are those of the pinning without the filter.

**Fixture.** `two-fam-bc`'s case (drag-stokes's box and gas, materials A and B from `two-mat`, one centred parcel per
inlet cell at `ds = 10`; `out-file = S`, `out-time = on`) with its own `INPUT/bc.txt` (`make_fixture.py`, from
`two-fam-bc`'s): family A (copy 1) field 1 0.34/(1 − 0.34) on rows j = 3–5 and 0 on rows 1–2, family B (copy 2) 0 on
every record, as ATLAS writes A's krho 0.34 on rows 3–5 beside a B without keys; G 0.

**Gates** (`check.py`): C1 the fixture (A 15 records loaded, 10 zero; B 25 zero); Z1 the log: family 1, 10 parcels
not injected, family 2, 25; Z2 family A: 15 parcels, IDs 1–15, all pinned at y > 2 cm, each leaving at x = 0.15 with
ṁ = 0.34/0.66 × 1.2e-3 = 6.181818e-4 kg/s (1e-5); Z3 family B: no trajectory row and no exit record, and no exit of
either family at a zero mass flow. B's group runs empty (0 particles) to the stop condition.

**RED.** Without the call to `dropZeroFlow`: Z1 (no line), Z2 (25 parcels, IDs 1–25, 10 pinned at y ≤ 2 cm leaving
at ṁ = 0) and Z3 (B: 1525 trajectory rows and 25 exits; 35 exits at a zero mass flow) fail.
