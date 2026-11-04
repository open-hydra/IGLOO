# INFO — breakup/pe-sigma-table (e2e)

## Purpose
Gates the `Sigma` and `Viscosity` columns of `INPUT/properties.dat` end to end: a breaking
material that takes its surface tension and its liquid viscosity from the table, with no
`[IGLOO-Properties]` section, read at the drop temperature by the right-hand side of a
continuous-rate breakup model. The physics is `pilch-erdman-e2e`'s; only the route of σ and
μ_p into it changes.

## Case construction
- `INPUT/solfile.tec` and `INPUT/bc.txt` are `pilch-erdman-e2e`'s (symlinked): the uniform gas
  box (0.15 × 0.05 × 0.05 m, 60 × 5 × 5 cells, U = 200 m/s, ρ_g = 1.2, T = 300 K) and its 25
  face-1 inlet cells, one diameter each (0.12–6.0 mm), κ_v = 0.5 (slip 100 m/s), κ_t = 1.0.
- `INPUT/phase.txt`: one material `A`. `INPUT/properties.dat`: water, ρ_p = 1000, c_p = 4182,
  `h = c_p T`, 101 rows on 270..370 K, with `Sigma = 0.036 + 1e-3 (T − 300)` N/m (six decimals)
  and `Viscosity = 0.05 exp(−(T − 300)/10)` Pa s (ten significant digits).
- `input.ini`: `pilch-erdman-e2e`'s models (`Schiller-Naumann`, `Pilch-Erdman`, `Cd = 0.5`,
  `B = 0.0758`) with `heat = NoHeat`, so every drop stays at its injection temperature,
  T0 = 300 K, and the solver's σ and μ_p are the columns at 300 K: 0.036 N/m (half
  `pilch-erdman-e2e`'s constant, so every Weber number doubles: We 40…2000) and 0.05 Pa s
  (fifty times its constant).

μ_p enters the model through the Ohnesorge number: the critical Weber number
`12 (1 + 1.077 Oh^1.6)` at every drop and, for Oh > 0.1 and We < 228, the Ohnesorge branch of
T*, which the drops with We 80–220 take (Oh 0.32–0.54 at μ_p = 0.05 Pa s).

## Oracle and gate (`check.py`)
- Premise: the drop temperature is T0 on every trajectory row (else the oracle's σ and μ_p
  would not apply).
- Oracle: `pilch-erdman-e2e`'s PE87 kernel, imported from its `check.py` unchanged, with its
  module constants `SIGMA` and `MUP` set to the columns' values at T0, read from the file by
  `tests/tools/proptab.py` (`table_value`, the solver's end-value rule) — never through the
  solver's reader.
- Gate (`pilch-erdman-e2e`'s first gate): IGLOO's initial breakup rate of every drop with
  We ≥ 45 within 1 % of the kernel's, the same leading-window least-squares slope on both
  sides; at least 20 drops gated.

## Result
24 drops gated, 0 violations, worst relative difference 6.5e-5. Drop 2 (We 60) sits at its
stable diameter on both sides (its rate is round-off on both).

## RED
- The reader without the two columns: setup stops (`Breakup requires sigma`, exit 128).
- σ and μ_p read at T + 1 K in the right-hand side: 6 violations, worst 0.29 — the We 80–220
  drops move 3.4–29 % (μ_p falls 9.5 % per kelvin), the drops with We ≥ 260 0.7–0.8 % (inside
  the tolerance).
- The σ table filled from the `Viscosity` column: 21 violations, worst 0.35.
