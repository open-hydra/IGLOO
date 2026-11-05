# INFO — evaporation/burnout-threshold (the burnout mass on both sides of 1e-14 kg, e2e)

## The rule
A parcel of a consuming model (ODE models 2, 4, 5 and 8) is consumed at the first accepted step whose end has its
mass at or below

    m_b = min(1e-6 m₀, 1e-20 kg),

m₀ its reference mass (the injection mass here). The two bounds cross at m₀ = 1e-14 kg: lighter drops burn out at
1e-6 of their mass (d_b = 0.01 d₀ at a fixed density), heavier ones at 1e-20 kg. The oracles' copy of the rule is
`tools/burnout.py`, written from the rule's statement.

## Reference
The classical d²-law at Nu = Sh = 2 (Re = 0), `d²(t) = d₀² − K t`, `K = 8 k_g ln(1 + B_T)/(c_p,g ρ_p)`,
`B_T = c_p,g (T_g − T_p)/L_v` (Godsave 1953 `[God53]`, Spalding 1954 `[Spa53]`), with `K0` at the injection
temperature: `B_T0 = 1.50675`, `K0 = 2.378662e-8 m²/s` from the inputs alone. The drop warms, `B_T` and `K` fall,
so the exact life is never shorter than `(d₀² − d_b²)/K0`.

## Case construction
`d2law-line`'s case with one diameter per inlet cell:

    python3 tests/tools/make_box_case.py INPUT --kv 1 --kt 0.5 --rp-sweep 5.00000e-08,6.05764e-08,…,5.00000e-06

radii `0.05 µm · 10^{k/12}`, k = 0 … 24 (D₀ = 0.1 … 10 µm, log-spaced; the reader's cell order, n outer, m inner).
The gas file is the shared box's (600 K, 10 m/s; `d2law-line`'s `solfile.tec`, symlinked), `phase.txt` the common
one, `properties.dat` `d2law-line`'s (ρ_p 8000, c_p 125000 J/(kg K), constant; symlinked): the 100-fold c_p keeps the
drops within 20 K of their 300 K injection (313.8–319.7 K at burnout). `kv = 1`: zero slip, so `x = u_g t` and a parcel's exit (Y, Z) names
its inlet cell. `out-time = on` puts `t` on every exit record; SDIRK4 at `relative-tol = absolute-tol = 1e-11`,
`d2law-line`'s.

| D₀ [µm] | m₀ [kg] | m_b [kg] | d(m_b) |
|---|---|---|---|
| 0.100 | 4.19e-18 | 4.19e-24 (1e-6 m₀) | 1.000 nm |
| 1.000 | 4.19e-15 | 4.19e-21 (1e-6 m₀) | 10.00 nm |
| 1.468 | 1.32e-14 | 1e-20 | 13.37 nm |
| 10.00 | 4.19e-12 | 1e-20 | 13.37 nm |

## What it verifies (`check.py`, inputs and the printed records only)

| id | gate | tolerance | measured |
|---|---|---|---|
| B0 | 25 exit records, one per inlet cell, every `bc` = 0, `t > 0`, `X < Lx − 1e-6`: every parcel burns out inside | exact | 25 of 25 |
| B1 | the remnant: every exit record's `d_p ≤ d(m_b)`(1 + 1e-5), m₀ = π/6 ρ_p D₀³ (ρ_p from the table at 300 K, D₀ = 2 rp of the cell's `bc.txt` record), `d(m) = (6m/(π ρ_p))^{1/3}` | twice the E13.6 print of `d_p` | `d_p/d(m_b)` 0.99168 … 0.99989 |
| B2 | the full life, one-sided: `t_exit ≥ (D₀² − d(m_b)²)/K0 (1 − 1e-3)`; the lag `t_exit/t_b − 1`, `t_b = D₀²/K0`, printed | 1e-3 (`d2law-line`'s whole residual is 3.1e-4 of d₀²) | lag +3.261e-3 … +3.362e-3 |

The burnout step overshoots `m_b` by little: `d_p ≥ 0.9917 d(m_b)`, the mass at least 0.975 `m_b`. B1 is the rule's
witness; the lag (+3.3e-3) is the slow heating's, the same at every size up to 1 µm and growing to +3.36e-3 at 10 µm.

## RED (scratch builds, `check.py` on each)

| solver | what fails |
|---|---|
| a 1e-15 kg threshold | B0: the 10 parcels below 0.62 µm (m₀ ≤ 1e-15 kg) consumed at injection, `t = 0`; B1 on 25, `d_p/d(m_b)` 46 to 100; B2 on 24, all but the 10 µm parcel |
| the absolute bound alone (`m ≤ 1e-20 kg`) | B1 on the 14 parcels below 1.34 µm, `d_p` up to 13.3 `d(m_b)`; B2 on the 4 below 0.2 µm, up to 1.4e-2 short |
| the relative bound alone (`m ≤ 1e-6 m₀`) | B1 on the 11 parcels above 1.34 µm, up to 7.4 |
| the larger of the two | B1 on all 25; B2 on 4 |
