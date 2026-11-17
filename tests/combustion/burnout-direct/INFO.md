# INFO — combustion/burnout-direct (drops of `combustion=direct` burned out at `d-burnout`, e2e)

## What it verifies
A drop of `combustion=direct` that heats into its burning temperature burns on the d²-law at a constant temperature
and is removed at the removal size `d-burnout`, inside the box, what it still holds handed to the gas. The burning
branch (`rhsAlCombustion`, `src/lib/Lib_RHS.f90`) and the burn-out removal (`solout`,
`src/lib/Lib_Integration.f90`) test the same temperature, `directOnset` × `T-burn` (`directOnset` = 1 − 1e-10,
`src/lib/Lib_Combustion.f90`): a drop burns from the step that ends at or above it, and its temperature then stays
where that step left it, since every heat it receives feeds the burning.

## Reference
The d²-law of a drop whose vapour burns in the gas (`[GG09]`, eqs. 4-5), at Nu = 2 (Ranz-Marshall at
Re = 0: no slip):

    d²(t) = d_a² − K (t − t_a),   K = 4 k_g Nu ln(1 + B)/(ρ_p c_p,g),   B = (c_p,g (T_g − T_b) + q-comb)/Lv,

c_p,g = γR/(γ − 1) = 1004.5 J/(kg K). At T_b = `T-burn` = 302 K: B = 2.898087, K = 9.549606e-8 m²/s.

## Case
burn-box's uniform gas (600 K, 10 m/s, k_g 0.026 W/(m K), γ 1.4, R 287 J/(kg K); `solfile.tec` symlinked),
`common/properties.dat` (ρ_p 2950, c_p 1250, constant: the temperature is the state), its own `bc.txt`:

    python3 tests/tools/make_box_case.py INPUT --kv 1 --kt 0.5 --rp-sweep 5.00e-06,5.40e-06,…,1.46e-05

25 inlet cells, one radius each, 5.0 to 14.6 µm by 0.4 µm (the reader's cell order, n outer, m inner): 25 drops of
10 to 29.2 µm enter at 300 K with the gas velocity. `phase.txt` `A 1 combustion=direct`; `models.txt`

    A combustion direct T-burn=302 q-comb=3.1e7 d-burnout=3.0e-6

and `[IGLOO-Properties] Lv = 1.08e7`. Each drop reaches 302 K within 7.9 µs (10 µm) to 67 µs (29.2 µm) of its
injection, before its first cell crossing, and burns out between x = 9.6 mm and 89.0 mm. SDIRK4 at
`relative-tol = absolute-tol = 1e-11`, `out-time = on`.

`check.py` runs the same input a second time into `ref-varcp/` with `standard/temp-relax-varcp`'s table (c_p = 1000 +
(T − 250) J/(kg K), ρ_p 2950, relative enthalpy): there the enthalpy is the state and the burning branch and the
removal read the temperature back from the table (`comp_TfromTab`). Both runs pass D1-D6; `ref-varcp/` also prints
the enthalpy banner.

## Pass criterion (`check.py`, the inputs and the printed records only)

| id | gate | tolerance | measured (this run; `ref-varcp/`) |
|---|---|---|---|
| D1 | no `Run_ODESolver err`, `non-finite state`, `maxIter` or `EXIT OF` line in `run_out.txt` | exact | 0; 0 |
| D2 | 25 exit records, IDs 1-25, every `bc` = 0 and X < 0.15 m: each drop removed inside the box | exact | 25; 25 |
| D3 | each exit diameter at the removal size, 0.5 `d-burnout` ≤ d_e ≤ `d-burnout` | the crossing step's burn | d_e/`d-burnout` 0.992713 … 0.999980; 0.992823 … 0.999137 |
| D4 | the injection row at 300 K; one temperature T_b on every later row and on the exit record, `T-burn` ≤ T_b ≤ `T-burn` + 1 K | F12.6 print; 1 K | T_b − `T-burn` 0 … 1.96e-2 K; 0 … 4.24e-3 K |
| D5 | every row after the drop's first burning row (t_a, d_a), and its exit record, on the d²-law with K at the drop's own T_b | the half-ULPs of d, d_a (E13.6), t, t_a (ES, 8 decimals), T_b (F12.6); 1e-9 of d² for the integrator | 417 rows, 0.966 of the tolerance at most; 417 rows, 0.935 (the rounding of a drop's anchor, shared by its rows) |
| D6 | Σ `wdot(A)` of `OUTPUT/source.tec` = the flow injected, Σ col1 ρ_g u_g A over the 25 inlet records of `bc.txt` (A = 1e-4 m², G = 0): every stream burns out inside the box and its remnant is handed | 1e-12 rel. | 8.98e-16; 7.86e-16 |

D3 and D4's bounds are not closed forms: how far below the removal mass the crossing step ends and how far past
`T-burn` the step that reaches it ends depend on the step sizes. A drop that is not removed goes on burning to
d_e < 1e-11 m or NaN, far outside D3's band. D5 with K off by 2e-5 fails every drop (3.1 times the tolerance at most).

## RED
The removal tested at `T-burn` itself, the burning branch's margin left out: a drop whose crossing step ends in
[`T-burn`(1 − 1e-10), `T-burn`) burns but is never removed. In this run 18 of the 25 drops end there, in
`ref-varcp/` 24: they burn past `d-burnout` down to d → 0, the solver stops (`EXIT OF SDIRK4 AT X= 0.1055069D-02
H= 0.9409357D-18` and `[IGLOO] Run_ODESolver err=-1 part=<n> ... => marking gone`) and the drop leaves on its last
accepted state. `tools/check_finite.py` fails first (5 exit records print d_p = NaN); `check.py` gives D1 36 and 48
failure lines, D3 18 and 24 of 25 outside (d_e/`d-burnout` from 0 to 0.998700 and 0.998363), D5 5 and 3 drops off
(NaN). D2, D4 and D6 pass: the last accepted state is handed to the gas (`revertToLastAccepted`), the books stay
closed.

Full test matrix in the family overview [../INFO.md](../INFO.md).
