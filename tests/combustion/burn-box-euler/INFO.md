# INFO — combustion/burn-box-euler (equivalent-Eulerian field of ODE model 5, e2e)

## Reference
No literature reference: the burn law is [`combustion/burn-box`](../burn-box/INFO.md)'s and is gated there.
This case gates the **equivalent-Eulerian deposit of the burning material (ODE model 5)**.

## Why it exists
`Lib_Integration::computeEulField` deposits, for models 2, 5 and 8, the number density `npdot·Tstay/V` and the
density `npdot·intE(5)/V`, where `intE(5) = ∫ m dt` is the mass integral of **one** droplet over its stay in
the cell (`F(14) = m` in `rhsAlCombustion`); `finalizeEUL` divides the velocity and temperature moments by
that density. `burn-box` writes the source field only, `infrastructure/body-burn-m5` writes the Eulerian field
without reading it, `evaporation/tc-box-euler` covers model 2 and `combustion/burn-box-cap` model 8: no other
fixture gates model 5's Eulerian field.

## Fixture
`INPUT/` is four per-file symlinks to `../../burn-box/INPUT/`; `input.ini` is burn-box's with two changes:

| change | why |
|---|---|
| `out-file = S` → **`ALL`** | turns `eulerSwitch` on, so the Eulerian moments are integrated and `OUTPUT/euler1.tec` is written |
| **`mollify = off`** | the binomial smoother mixes `ρ_p` and `n_p` with the same kernel; their per-cell ratio is exact only unsmoothed |

25 parcels (d₀ = 30 µm, `kv = 1`) ride their own `(j,k)` column of 60 cells at u_g = 10 m/s, ignited from
injection, and lose 65 % of their mass by the last checked row. Constants come from `burn-box/check.py`
(ρ_p = 2950 kg/m³ from `properties.dat`, d₀, u_g) plus the loading `krho = 0.34` and the gas density 1.2 kg/m³:
`m₀ = 4.170464e-11 kg`, `ṁ = 0.34/0.66·1.2·10·1e-4 = 6.181818e-4 kg/s`, `npdot = 1.482285e7 s⁻¹`.

## What it verifies (`check.py`)
Per cell of every column, with both bounding trajectory rows present (all 1500 cells: `print-dcell = 1`
writes one row per cell entry and the exit row at the outlet):

| id | gate | tolerance | measured |
|---|---|---|---|
| W1, W2 | `out-file = ALL` and `mollify = off` are read (log witnesses) | — | pass |
| W3 | the droplets burn: mass lost by the last checked row | > 0.2 of m₀ | 0.650 |
| N | cells checked | ≥ 1000 | 1500 |
| E3 | `ρ_p/n_p` (the time-mean droplet mass in the cell) against `(m_in + m_out)/2` of the bounding rows | 1e-4 relative | 2.11e-5 |
| E4 | `u_p = u_g`; `|v_p|`, `|w_p|` | 1e-11 relative; 1e-18 m/s | 5.5e-14; 7.1e-23 m/s |
| E5 | `T_p` inside `[T_in, T_out]` of the bounding rows | 0.05 K | 0 |
| E1 | `Σ ρ_p V` against `Σ npdot·∫m dt` over the rows (coarse) | in [0.5, 2] | 1.000 |

Floor: `euler1.tec`, `source.tec` and the check printout are byte-identical at 1, 2 and 5 threads (twice
at 5), the trajectories identical as sorted multisets — no two threads deposit into one cell. The E3 bound
sits five times above the trapezoid's own error and four decades below the defects below.

## RED (each one line of `computeEulField`, model 5 only)

| mutation | what fails |
|---|---|
| model 5 deposits `intE(5)·factor`, without the number rate | E3 1.000; E4 `u_p` off by 1.48e7 relative; E5 1.05e10 K; E1 0.0000 |
| model 5 takes the default branch (the injected `ṁ·Tstay`) | E3 1.754 (the injected mass against the burning mean); E4 1.0; E5 707 K |

Under the first mutation `evaporation/tc-box-euler` (model 2, the same `case`) stays green: this case is the
only one that sees model 5's deposit.
