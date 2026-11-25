# Infrastructure tests

The machinery every run depends on: gas sampling and cell tracking, boundary faces, wedges, multi-block and chimera grids, injection placement, the coupling outputs, the input and file readers, solver failures and the setup refusals. The results are discussed in the documentation, [Geometry & Boundaries](../../docs/vv/geometry.md), [Injection](../../docs/vv/injection.md), [Coupling & Outputs](../../docs/vv/coupling.md) and [Contracts & Robustness](../../docs/vv/contracts.md); this file holds what is needed
to run the tests and find a case.

## Run

```bash
./tests/test.sh infrastructure                                 # the whole folder (builds build/verif/ first)
./tests/test.sh <name>                                # one test by its CTest name, e.g. periodic-y
ctest --test-dir build/verif -L infrastructure --output-on-failure
```

The setup refusals live in [`refusals/`](refusals/INFO.md): each is registered as `refuse-<folder>` and checked by `../tools/check_refusal.py` (exit 128, the `error stop` message, nothing integrated, no output). Several refusal inputs are written by the `make_*.py` scripts in that folder.

An end-to-end case can also be run by hand from its folder, after one `./tests/test.sh` build:

```bash
cd tests/infrastructure/periodic-y
../../../bin/IGLOO > run_out.txt 2> run_err.txt && python3 -B check.py
```

The CTest name of a case is its folder name; of a unit family, the `test_*` names in the table.
Each folder's `INFO.md` holds the reference, the oracle and its tolerance, and what a deliberately
broken build gives. General instructions (requirements, reading the report, diagnosing a failure,
adding a test) are in [`../README.md`](../README.md).

## Unit families

| Folder | CTest | What it checks |
|---|---|---|
| [`axis_dispatch/`](axis_dispatch/INFO.md) | `test_axis_dispatch` | Bcdef-200 face classification |
| `bc_families/` | `test_bc_families` | The boundary file's per-family copies and chimera (102) records |
| [`cell_volume/`](cell_volume/INFO.md) | `test_cell_volume` | Cell and octant volumes: the trilinear volume |
| [`containment/`](containment/INFO.md) | `test_containment` | The 12-plane point-in-hexahedron test and its tolerance |
| [`dual_clip/`](dual_clip/INFO.md) | `test_dual_clip` | The ord2 dual mesh tiles the domain |
| [`fmt/`](fmt/INFO.md) | `test_fmt` | The particle files' edit descriptors without the runtime |
| [`gas_reconstruction/`](gas_reconstruction/INFO.md) | `test_gas_reconstruction` | In-cell gas sampling |
| [`geometry_rules/`](geometry_rules/INFO.md) | `test_geometry_rules` | The boundary search's face |
| `ghost_bc/` | `test_ghost_bc` | The gas-order-2 ghost nodes of connected faces |
| [`graze_standoff/`](graze_standoff/INFO.md) | `test_graze_standoff` | Reflection standoff ≤ cell |
| [`hop_crossing/`](hop_crossing/INFO.md) | `test_hop_crossing` | A reflection's crossing from the step's two end states |
| [`host_assoc/`](host_assoc/INFO.md) | `test_host_assoc` | Internal procedures address the main program's arrays |
| [`ini_pipeline/`](ini_pipeline/INFO.md) | `test_ini_pipeline` | INI → module variables |
| [`material_tokens/`](material_tokens/INFO.md) | `test_material_tokens` | The material-line tokens of the phase file |
| [`mollify_blocks/`](mollify_blocks/INFO.md) | `test_mollify_blocks` | The field mollifier across 101 block interfaces, 201 periodic faces and chimera faces; the chimera overlap mirror |
| [`pin_grid/`](pin_grid/INFO.md) | `test_pin_grid` | The bucket grid of the 3D injection placement |
| [`properties_reader/`](properties_reader/INFO.md) | | The `properties.dat` reader |
| [`psd_table/`](psd_table/INFO.md) | `test_psd_table` | Injection size distributions |
| [`rng_stream/`](rng_stream/INFO.md) | `test_rng_stream` | Per-particle RNG stream contract |
| [`slot_order/`](slot_order/INFO.md) | `test_slot_order` | The packed property slots of the right-hand sides |
| [`source_reduction/`](source_reduction/INFO.md) | `test_source_reduction` | Conservative ord2 source reduction |
| [`tec_read/`](tec_read/INFO.md) | | The gas file's one-pass reader |
| [`tec_write/`](tec_write/INFO.md) | | The grid files' stream writer |

## End-to-end cases

| Folder | What it checks |
|---|---|
| [`axis-200/`](axis-200/INFO.md) | The axisymmetric AXIS face tagged `axisymmetric` (bcdef 200) |
| [`bc-center-2grp/`](bc-center-2grp/INFO.md) | Multi-group `bc_center` pinning |
| [`body-burn-m5/`](body-burn-m5/INFO.md) | The body-force accumulators of a burning parcel (ODE model 5) |
| [`body-evap-m2/`](body-evap-m2/INFO.md) | The body-force accumulators of an evaporating parcel (ODE model 2) |
| [`body-evap-m4/`](body-evap-m4/INFO.md) | The body-force accumulators of an evaporating parcel under ODE model 4 |
| [`chimera-2blk/`](chimera-2blk/INFO.md) | Two overlapping blocks joined by chimera (102) faces |
| [`chimera-2blk-mollify/`](chimera-2blk-mollify/INFO.md) | The field mollifier across chimera (102) faces |
| [`chimera-2blk-rim/`](chimera-2blk-rim/INFO.md) | A partly covered chimera face: rim exits |
| [`chimera-2blk-shift-mollify/`](chimera-2blk-shift-mollify/INFO.md) | The overlap mirror and the field mollifier between misaligned cells |
| [`chimera-2blk-shift-ord2/`](chimera-2blk-shift-ord2/INFO.md) | Chimera faces between misaligned cells at gas order 2 |
| [`chimera-2blk-shift3/`](chimera-2blk-shift3/INFO.md) | The overlap mirror between cells shifted by 0.3 |
| [`conn-103/`](conn-103/INFO.md) | A 103 face written header only |
| [`coupled-body/`](coupled-body/INFO.md) | Euler + source + body-force combined output path |
| [`db-2daxi/`](db-2daxi/INFO.md) | 2Daxi wedge + DB injection + euler output |
| [`db-injection/`](db-injection/INFO.md) | Assigned-position injection |
| [`db-varrho/`](db-varrho/INFO.md) | Database streams sized at their injection cell's gas |
| [`db-varrho-ord2/`](db-varrho-ord2/INFO.md) | Database streams sized at their injection gas under gas-order 2 |
| [`exit-census/`](exit-census/INFO.md) | Every parcel whose integration stops leaves one exit record |
| [`foreign-keys/`](foreign-keys/INFO.md) | Presence tests ask IGLOO's own section |
| [`foreign-keys-db/`](foreign-keys-db/INFO.md) | The optional keys of the assigned-position path |
| [`ini-comment-eq/`](ini-comment-eq/INFO.md) | `=` inside an in-section `;` comment |
| [`inj-dscell/`](inj-dscell/INFO.md) | A cell's own injection spacing |
| [`inj-dsdegen/`](inj-dsdegen/INFO.md) | The injection degeneracy floor |
| [`inj-edge/`](inj-edge/INFO.md) | Injection points on the edges of the inlet cells |
| [`inj-edge-payload/`](inj-edge-payload/INFO.md) | A parcel carries the payload of the inlet cell that pinned it |
| [`inj-mdotmax/`](inj-mdotmax/INFO.md) | Injection placement from a cap on the parcel's mass flow |
| [`inj-pack/`](inj-pack/INFO.md) | Several parcels per inlet cell, packed at the spacing ds |
| [`inj-pack-polar/`](inj-pack-polar/INFO.md) | One hexagonal lattice over an annular-sector inlet |
| [`inlet-401-g/`](inlet-401-g/INFO.md) | The 401 inlet's field G |
| [`inlet-402-403/`](inlet-402-403/INFO.md) | The particle inlets 402 and 403 |
| [`inlet-zero-load/`](inlet-zero-load/INFO.md) | An inlet cell that injects nothing pins no parcel |
| [`inner-loop-cap/`](inner-loop-cap/INFO.md) | The inner loop's cap, reached by every parcel |
| [`periodic-2d/`](periodic-2d/INFO.md) | Translational periodic faces (201) on a one-layer box |
| [`periodic-2d-ord2/`](periodic-2d-ord2/INFO.md) | Translational periodic faces (201) at gas-order 2 |
| [`periodic-2d-ord2-dual/`](periodic-2d-ord2-dual/INFO.md) | The gas-order-2 dual across a periodic (201) face |
| [`periodic-y/`](periodic-y/INFO.md) | Translational periodic BC (bcdef 201) |
| [`planar-slab/`](planar-slab/INFO.md) | A single-layer (Nk = 1) box is a planar 2D case, not a wedge |
| [`plane-back/`](plane-back/INFO.md) | Plane crossings in both directions |
| [`psd-bare-token/`](psd-bare-token/INFO.md) | A size-distribution file named without '.' or '/' |
| [`rc-gate/`](rc-gate/INFO.md) | The e2e harness fails a non-zero exit |
| [`rot-2blk/`](rot-2blk/INFO.md) | An ATLAS two-block boundary file with a rotated interface |
| [`solver-fail-consumed/`](solver-fail-consumed/INFO.md) | `err < 0` exit hands the remnant to the gas |
| [`solver-fail-consumed-m4/`](solver-fail-consumed-m4/INFO.md) | `err < 0` exit hands the remnant to the gas, model 4 |
| [`solver-fail-midway/`](solver-fail-midway/INFO.md) | The Eulerian deposit of a segment the solver gives up on |
| [`solver-fail-refine/`](solver-fail-refine/INFO.md) | A solver failure on a refinement call ends the parcel at its last accepted step |
| [`stuck-cell/`](stuck-cell/INFO.md) | The give-up guard of a parcel that does not move |
| [`swirl-wedge-2blk-ord2-dual/`](swirl-wedge-2blk-ord2-dual/INFO.md) | The gas-order-2 dual across a block interface on a wedge |
| [`sym-plane-bounce/`](sym-plane-bounce/INFO.md) | Slow parcels pressed onto a symmetry plane |
| [`sym-plane-bounce-ord2/`](sym-plane-bounce-ord2/INFO.md) | Slow parcels pressed onto a symmetry plane at gas order 2 |
| [`sym-plane-creep/`](sym-plane-creep/INFO.md) | Parcels creeping on a symmetry plane |
| [`sym-plane-creep-ord2/`](sym-plane-creep-ord2/INFO.md) | Parcels creeping on a symmetry plane at gas order 2 |
| [`sym-plane-hop/`](sym-plane-hop/INFO.md) | Parcels bouncing in place on a symmetry plane: the reflection's crossing |
| [`sym-plane-skate/`](sym-plane-skate/INFO.md) | Parcels launched along a symmetry plane |
| [`sym-plane-skate-ord2/`](sym-plane-skate-ord2/INFO.md) | Parcels launched along a symmetry plane at gas order 2 |
| [`timing-report/`](timing-report/INFO.md) | `[IGLOO-General] timing = on` |
| [`two-fam-bc/`](two-fam-bc/INFO.md) | A particle boundary file with one copy of the inlet table per family |
| [`two-fam-bc-2blk/`](two-fam-bc-2blk/INFO.md) | Two families on a box of two blocks joined by a 101 connection |
| [`two-fam-bc-2blk-mollify/`](two-fam-bc-2blk-mollify/INFO.md) | The field mollifier across a 101 block interface |
| [`two-fam-bc-2blk-ord2-dual/`](two-fam-bc-2blk-ord2-dual/INFO.md) | The gas-order-2 dual across a block interface |
| [`two-fam-bc-2blk-ord2-dual-1l/`](two-fam-bc-2blk-ord2-dual-1l/INFO.md) | The gas-order-2 dual across a block interface, one layer |
| [`two-fam-bc-2blk-ord2-dual-default/`](two-fam-bc-2blk-ord2-dual-default/INFO.md) | The gas-order-2 dual across a block interface, smoothed |
| [`two-fam-bc-2blk-stretch-ord2/`](two-fam-bc-2blk-stretch-ord2/INFO.md) | Ghost nodes of a block interface between unequal cells |
| [`two-mat/`](two-mat/INFO.md) | Two materials (nm = 2) that differ only by density |
| [`varcp-tmin/`](varcp-tmin/INFO.md) | A varying-cp table from 280 K, the parcels above its end |
| [`wall-approach/`](wall-approach/INFO.md) | The ord2 gas ghost ring must know its boundary conditions |
| [`wedge-axis-row/`](wedge-axis-row/INFO.md) | A wedge whose axis row sits at r = 1e-8 stays a wedge |
| [`wedge-fold/`](wedge-fold/INFO.md) | The axisymmetric 200-face FOLD with a swirling parcel |
| [`wedge-tiny/`](wedge-tiny/INFO.md) | A single-layer sector of 5e-7 rad is a planar single layer |
| [`zero-group-line/`](zero-group-line/INFO.md) | An unnamed 0-group material line is inert and named in the log |
