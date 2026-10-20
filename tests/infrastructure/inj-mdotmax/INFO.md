# inj-mdotmax — injection placement from a cap on the parcel's mass flow

**Purpose.** `[IGLOO-General] mdot-max` (g/s) caps the mass flow a parcel carries. The placement (`compute_ds`,
`initialization.f90`) asks each inlet cell for N_pmin = ⌈F/ṁ_max⌉ parcels — F the cell's mass flow,
κ_ρ/(1 − κ_ρ,tot)·(|ṁ_gas| + ṁ_part) for 401 — and packs the cell at the spacing R = d_cell/N_pmin (d_cell the
cell's mid-edge size; a positive `ds` below that wins), by `inj-pack`'s advancing front. No other case sets
`mdot-max`.

**Fixture.** `INPUT/` symlinks `standard/drag-stokes/INPUT` (the shared box, its uniform gas ρ_g U = 12 kg/(m² s)
and its 401 inlet, κ_ρ = 0.34); `input.ini` is `drag-stokes`'s with `ds = 0`, `mdot-max = 0.25` (g/s) and
`out-time = on`. Each 1 cm² inlet cell carries F = 0.34/0.66 × 12 × 1e-4 = 6.181818e-4 kg/s, so N_pmin = 3 and
R = 1/3 cm: 264 parcels, 9 to 14 per cell; 54 of the pinned points lie on face-cell edges, 36 of them on edges two
cells share.

**Gates** (`check.py`): a control (N_pmin = 3 from the inputs); `inj-pack`'s J0–J3 (`run()`: every parcel
integrated to an exit, the points on the edges each with the payload and share of a cell that holds them, the
inlet plane, the spacing R/2 and R within a cell — measured 3.332349e-3 m against R = 3.333333e-3 m — and each
cell's mass flow); and J4: every parcel's mass flow (exit records, by header name) at most ṁ_max — the largest
6.868690e-5 kg/s — and F/ṁ_p ≥ 2 in every cell (9 to 14).

**RED.**
- The `mdot-max` branch of `compute_ds` skipped (the spacing set as without the cap): one centred parcel per cell
  carrying 6.181820e-4 kg/s; J4 fails (above the cap, one parcel per cell), J0 (no point on an edge).
- Without the start in the pinning cell, the guard against a point no cell holds and the widened containment
  search: the run stops with a segmentation fault, 68 of 264 exit records written; J0 and J3 fail.
- The start in the pinning cell alone removed changes nothing here: every inlet cell carries the same payload and
  the widened search finds the points on the outer edges. `inj-pack` (a 402 row beside the 401 rows) and
  `inj-edge-payload` carry that RED.
