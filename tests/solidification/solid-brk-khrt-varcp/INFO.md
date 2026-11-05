# INFO — solidification/solid-brk-khrt-varcp (KH-RT on a solidifying material with a varying heat capacity and density; e2e)

## What it verifies
[solid-brk-khrt](../solid-brk-khrt/INFO.md) (solid-brk-tab's box, inlet and material, symlinked; KH-RT with
`B0 = 0.2`; `out-file = ALL`) with its own `INPUT/properties.dat`: `Cp = 1250 + 0.5 (T − 2400)` J/(kg K), the
`Enthalpy` column its exact integral, and `Density = 2950 + 0.5 (2400 − T)` kg/m³. ODE model 7 keeps the temperature
as its state; the sites that read it there and that the box case cannot reach run here: the KH shed's enthalpy, a
child's first state, the density at a breakup event, the breakup's properties in the step test, and the Eulerian
temperature. The drops slip (Stokes drag, Ranz-Marshall), so the gates are solid-brk-khrt's books with the table's
enthalpy and density, read from the table file itself; the injected droplet mass is `ρ(T_0) π d_0³/6`
(4.170464e-11 and 4.550331e-11 kg for the two families). Family 1 sheds 175 children, as on the constant table.

| id | gate | tolerance | measured |
|---|---|---|---|
| K1 | solid-brk-khrt's: family 2 sheds nothing, family 1's logged children are its child exit records, every child with two rows | exact | 175 logged, 175 records; 0 in family 2 |
| K2 | from a parcel's first frozen row on every row carries that row's `m` (`d` follows `ρ(T)` at a constant mass); every family-2 parcel's `m` decreased on that row; on every trajectory row `\|m − ρ(T) π d³/6\|` | exact; the E13.6 half-ULPs of `d` and `m` | 13425 rows; 25 of 25 stripped; 13675 rows, worst 0.963 of the bound |
| K3 | solid-brk-khrt's: every exit at the outlet, each family's exit flows sum to 25 injected flows | 1e-5 | 225 records, worst 1.1e-6 |
| K4 | `Σ E` of `source.tec` = `Σ ṁ (h(T_0) + v_0²/2)` injected minus `Σ ṁ (h_exit + \|u\|²/2)` over the exit records, `h` from the table by the record's phase | 1e-5 | 1.3e-7 |
| K5 | solid-brk-khrt's Eulerian mass per family | 1e-5 | 1.5e-7, 1.1e-6 |
| K6 | every scatter row: `0 < m ≤` the larger injected droplet mass, `\|m − ρ(T) π d³/6\|` within the row's print bound, `0 ≤ f ≤ 1`; rows with `0 < f < 1` and with `f = 1` | first-order print bound | 87455 rows, worst 0.988; 342 and 87010 |
| K7 | every Eulerian cell holding particles: `T_p` between the gas temperature and the hotter injection temperature | 1e-9 relative | 600.000–2042.794 K and 871.005–2327.000 K |

K2 states solid-brk-khrt's frozen-size gate in the mass: with a density table the diameter of a frozen drop follows
its temperature. K7 has no counterpart in solid-brk-khrt. The scatter cloud holds 87455 rows against
solid-brk-khrt's 953832: it is saturated at one marker per accepted step during the RT cascade (solid-brk-khrt's
INFO), so its size follows the integrator's steps there, which the table's properties change.

## RED (each measured by one mutation of the code)
| mutation | what fails |
|---|---|
| a solver that refuses a varying `Cp` under solidification | setup refuses the material, exit 128, no output |
| the KH shed's enthalpy left as the state (`T` for a table) | K4 1.0e-1 |
| a child's first state an enthalpy | K4 rel 2.6e2, K7 12615.7 K, the scatter file unreadable |
| the density at a breakup event from `T` read as an enthalpy | K3 family 1 7.2 (relative), K5 7.6e-3, K6 a droplet mass above the injected one |
| the breakup's properties in the step test at `T` read as an enthalpy | K2 a family-1 parcel with no frozen row; K3 25 exits not at the outlet, family 1 9.0e-1; K5 family 1 9.6e-3; 29 violations |
| the Eulerian temperature read from its moment as an enthalpy | K7 11.354–34.800 K and 16.120–38.953 K |
| `updatePart` reading `T` as an enthalpy | K3 7.6e-2, K4 2.9e-1, K5 2.3e-4 |
| the injection state an enthalpy | the trajectory file unreadable (a value overflows its field) |
