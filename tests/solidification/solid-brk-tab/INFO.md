# INFO — solidification/solid-brk-tab (TAB breakup on a solidifying material; e2e)

## What it verifies
Breakup acts on the liquid of a solidifying material only (phases 0 and 1). TAB is an event model: the
material stays ODE model 6 and its breakup event is gated by the phase in force over the step.

Two inlet families of one material enter the shared box (`tools/make_box_case.py --kv 3.5 --families 2
--family-kt 4.0,3.1043333 --family-rp 15e-6,15e-6`; the material line of `INPUT/phase.txt` is
`A 2 solidification=on h-fus=1.07e6 cp-solid=600`, `T-melt = 2327 K` and `T-nuc = 1861.6 K` at their defaults;
`ρ = 2950`, `c_l = 1250` from `common/properties.dat`) at 35 m/s against the 10 m/s, 600 K gas: `d = 30 µm`,
`Re = 50`, `We_r = 11.25` with the fixture liquid's `σ = 1e-3 N/m`, `μ = 1e-3 Pa s`; `drag = Stokes`,
`heat = Ranz-Marshall`, `breakup = TAB` (method 2, `n = 3.5`), `mollify = off`, `out-time = on`. Design estimates
at the injection state: TAB's first breakup from rest about 64 µs after the injection; family 1 (`T_0 = 2400 K`)
nucleates about 1.3 ms after it (`x ≈ 4.6 cm`), so it breaks up while liquid; family 2 (`T_0 = 1862.598 K`, 1 K
above `T_nuc`) nucleates within 3 µs (`x ≈ 0.1 mm`, inside the first cell), about 20 times sooner than its first
breakup would come, so it must never break up.

A row's phase comes from its `(f, T)` (`tools/datcols.py` `solid_phase`); a parcel's first frozen row `k` is its
first row of phase 2 (plateau) or 3 (solid). A breakup inside the nucleation cell, before the nucleation, lands
on row `k` itself, so T1 counts a decrease on `k` as liquid and T3 starts at `k`.

| id | gate | tolerance | measured |
|---|---|---|---|
| T1 | every family-1 parcel's `d` decreases at least once on a row up to `k` | exact | 25 of 25 (none on `k`); exit `d` 2.9–13.1 µm |
| T2 | every family-2 row carries the injected `d` | exact (the printed value) | 1525 rows |
| T3 | every row from `k` on carries the `d` of row `k`, both families | exact (the printed value) | 2937 rows |
| T4 | every exit record ends at the outlet (code 400) with its family's injected `ṁ = krho/(1 − Σ krho) ρ_g u_g A` (8.326531e-4, 4.163265e-4 kg/s), and each family's total is 25 of them | 1e-5 (E13.6 half-ULP ≤ 5e-6 per record) | 50 records at the outlet; worst 1.1e-6 |
| T5 | `Σ E` of `source.tec` = `Σ ṁ_inj (c_l T_0 + v_0²/2) − Σ ṁ_exit (h_exit + \|u_exit\|²/2)`, `h` by the phase of the exit record; `ẇ = 0` | 1e-5 | 1.6e-7 |
| T6 | ≥ 20 parcels per family with the injection row, a frozen row and one exit record; the model-6 witness | — | 25 + 25 |

## RED (each measured by one mutation of the code)
| mutation | what fails |
|---|---|
| the code before breakup was allowed on a solidifying material | setup refuses `solidification=on` with breakup, exit 128 |
| the breakup event not gated by the phase | T2: 1500 of the 1525 family-2 rows, the frozen drops breaking up before their first face; T3 cannot see it (the breakup lands on row `k`) |
