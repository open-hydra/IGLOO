# INFO — solidification/solid-box-2mat (a model-6 material followed by a model-1 one, e2e)

## What it verifies
Two materials on `solid-box`'s box and inlet (`bc.txt`, `solfile.tec` symlinked; `properties.dat` is
`infrastructure/two-mat`'s two-zone table): `A` carries the solidification tokens (ODE model 6), `B` is a plain
constant-mass material of density 1000 (model 1). The reader copies each inlet line to every family, so both
inject the same 25 parcels. Materials are set up and integrated one after the other through module-level
layouts (`setupRHS`), so `B`'s group runs with whatever `A`'s left behind unless every index is reset; this is
also the suite's mixed model-6/model-1 run.

| id | gate | measured |
|---|---|---|
| A | `solid-box`'s regime oracle (E1, E2, E2c, E3, E5) on `trajectories-A.dat`, ≥ 20 parcels in all three regimes | 25 parcels; worst E1 0.010, E3 0.006 of tol |
| B1 | `T(x) = T_g + (T_a − T_g) e^{−(x−x_a)/L_B}`, `L_B = u_g c_p ρ_B d²/(12 k_g) = 3.606e-2 m`, over `|T − T_g| > 1 K`, ≥ 20 parcels with ≥ 3 window points | 25 parcels; worst 0.010 of tol |
| B2 | `u = u_g`, `v = w = 0`, `d` and `m = ρ_B π d³/6` constant, the injection row at `T_0` | pass |
| W | the log reports model 6 for `A` and integrates `B` | pass |

## RED (measured)
With the reset of the phase slot index removed from the auxiliary-state layout, `B` (no auxiliary state) packs
its phase into a zero-length buffer. The release build runs green (the slot lands in the one-element scratch
array); ifx 2023.1 `-check bounds` and `-check all` do not flag an index past the end of a zero-size
explicit-shape dummy; a gfortran build with `-fcheck=bounds` stops with `Index '1' of dimension 1 of array
'auxst' above upper bound of 0` in `packAuxState`, exit 2, while the correct code under the same build passes
this oracle. Before model 6 existed, setup refused `solidification=on`, exit 128.
