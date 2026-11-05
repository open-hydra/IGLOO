# INFO — solidification/solid-box-varcp (a varying liquid heat capacity under solidification; e2e)

## What it verifies
[solid-box](../solid-box/INFO.md) (its `bc.txt` and phase file, symlinked) with its own `INPUT/properties.dat`:
`Cp = 1250 + 0.5 (T − 2400)` J/(kg K) on 1..5000 K, the `Enthalpy` column its exact integral `50 T + 0.25 T²` (the
table check's trapezoid is exact), `Density` 2950. The state of a solidifying material stays the temperature: in the
liquid phases `m c_i dT/dt = Q̇`, `c_i = h(i+1) − h(i)` the chord of the kelvin cell holding `T` (the exact derivative
of the enthalpy the source term books); the recalescence takes `f0 = (h(T_m) − h(T_n))/h_fus`; the solid keeps
`cp-solid`. 25 particles at `T_0 = 2400 K` move with the 600 K gas (`Re = 0`, `Nu = 2`, `x = u_g t`).

`check.py` integrates the liquid cell by cell from the table file, `t = Σ (m c_i/G) ln((T_a − T_g)/(T_b − T_g))`,
`G = π d k_g Nu`: the exact solution of the discrete model. Values: `x_n = 0.033493 m` (`c_l(T_n) = 980.75`),
`f0 = 0.477209`, `s = 18.9670 /m`, `x_s = 0.061056 m`, 11 plateau faces, `L_s = 0.051058 m`, the exit at 902.511 K,
`h(T_0) − h_solid(T_exit) = 2014610.94 J/kg`.

| id | gate | tolerance | measured |
|---|---|---|---|
| V0 | the cell-by-cell oracle against the continuous closed form of the linear `Cp` (oracle self-check) | 1e-4 K | 5.09e-5 K |
| E1 | the liquid rows (`x < x_n − dx`): `T` against the oracle from `(0, T_0)` | solid-box E1's budget with the local length `L(T) = u_g m c(T)/G`, times 4 (the chord-seam term, below) | 325 rows, worst 0.190 of tol (1.17e-2 K at 30 mm) |
| E2, E2c | the plateau rows at `T_m`: 11 ± 1 per parcel, contiguous, the first in `(x_n, x_n + dx]` | exact | 275 rows |
| E3 | the solid rows (`x > x_s + dx`): anchored exponential with `L_s` | solid-box E3's budget | 875 rows, worst 0.007 of tol |
| E4, E4b | `Σ E` of `source.tec` = `Σ ṁ [h(T_0) − h_solid(T_exit)]` on the table; the plateau cells; no cell with `E < 0` | 1e-5 | 2.9e-7; 250 cells, 2.9e-7 |
| E5 | `u = u_g`, `v = w = 0`, `d`, `m` constant; the mollify-off and model-6 witnesses; no enthalpy-state notice; ≥ 20 parcels in all three regimes | as solid-box | 25 parcels |
| E6 | `f`: 0 liquid, 1 solid and at every exit; `f − s x` one constant per parcel equal to `f0 − s x_n` within solid-box's event bound at `c_l(T_n)` | spread 2.0e-5, offset 1.875e-4 | 5.2e-7, 1.65e-5 |
| S7 | the plane records: ten before `x_n` on the oracle (`f = 0`), three on the plateau at `f0 + s (x − x_n)` | solid-box S7's, with `L(T)` | 325 records, worst 0.090 and 0.085 of tol |
| S8 | every scatter marker's `T` between the `T` of its parcel's two bracketing rows (the bracket holding `x_n` excepted) | the two F12.6 half-ULPs | 2226 markers, 0 excess |

**E1's tolerance.** solid-box's budget (the F12.6 half-ULPs of `x` and `T` over the local length, the integrator floor
1e-8 K) times 4 for the chord seams: the liquid's heat capacity steps by 0.5 J/(kg K) at every integer kelvin, and
SDIRK4 stepping over these jumps leaves a residual of one sign that grows with `x` (the solver cools about 2.5e-5 faster
in relative time than the oracle). Against the bare budget, over the ODE tolerance: worst 0.874 at 1e-9, 1.12 at
1e-10, 0.761 at 1e-11 (the case's), 0.673 at 1e-12, 0.428 at 1e-13 under ifx; 0.860, 1.16, 0.735, 0.598, 0.465 under
gfortran. The factor is 3.6 times the worst of the ifx scan, 3.4 times gfortran's. With `c_l = 1250` held in the
liquid arm E1 fails by 691 times the tolerance at 30 mm and 5.3 times at 2.5 mm. The
energy deposit closes at the print floor.

The re-melt through the table (a plateau ending at `f = 0`) and the whole freeze at the recalescence
(`h(T_m) − h(T_n) ≥ h_fus`) do not occur here: the unit family's SG14 runs both against closed forms.

## RED (each measured by one mutation of the code)
| mutation | what fails |
|---|---|
| a solver that refuses a varying `Cp` under solidification | setup refuses the material, exit 128, no output |
| `c_l = 1250` in the liquid arm (the table's value at 2400 K) | E1 +42.86 K at 30 mm (691 × tol; 5.3 × at 2.5 mm); E2c, the first plateau row at 40 mm; E6 offset 8.18e-2 (436 × the bound); S7 52.1 K; 483 violations |
| `f0 = c_l(T_m)(T_m − T_n)/h_fus` | E6 offset 5.07e-2 (270 × the bound), the phase rule at 60 mm; 127 violations |
| the source term on the constant-`cp` enthalpy (`cp` unset for a table) | E4 4.46e-2 |
| the plane records reading `T` as an enthalpy | S7 32.1 K against 1863.4 K, 325 violations |
| the scatter markers reading `T` as an enthalpy | S8 33.2 K against [1914.8, 1953.3] K, 2226 violations |
| the injection state an enthalpy | the trajectory file unreadable (a value overflows its field) |
| `updatePart` reading `T` as an enthalpy | E1 2318 K, no plateau row; 2291 violations |
| the enthalpy-equation notice for a solidifying material | E5's witness |
