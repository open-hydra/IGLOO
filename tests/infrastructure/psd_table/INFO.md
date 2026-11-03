# INFO — psd_table (injection size distributions, unit)

`test_psd_table` (ctest `test_psd_table`, labels `unit;infrastructure;injection;t1`) samples every
injection law through the production `IGLOO_Lib_Statistics::sampleDiameter` and tests the draws
against the law's cumulative distribution, typed in the test and never read back from the solver.
The end-to-end counterpart, on 20 000 injected parcels, is [`standard/psd-table`](../../standard/psd-table/INFO.md).

## The laws

| law (bc.txt column 8) | parameters | CDF of the injected diameter d |
|---|---|---|
| a size-distribution file | per line the radius r_k [m] and its mass fraction m_k (any scale) | histogram: class k holds w_k = m_k / Σm uniformly between the class edges; edges at the mid-points between neighbouring radii, the end classes extended by half the adjacent spacing (the first edge floored at r_1/2), all doubled to diameters; F is piecewise linear through (e_k, Σ_{j≤k} w_j) |
| `Normal` | mean 2 r_p, standard deviation σ_p | Normal truncated to d > 0 (the sampler rejects d ≤ 0): (Φ((d−μ)/σ) − Φ(−μ/σ)) / (1 − Φ(−μ/σ)) |
| `LogNormal` | mean 2 r_p, standard deviation σ_p | Φ((ln d − m)/s), s² = ln(1 + (σ/μ)²), m = ln μ − s²/2 |
| `RosinRammler` | x0 = 2 r_p, spread n = σ_p | 1 − exp(−(d/x0)^n) |

## Tests

| id | what | oracle | tolerance |
|---|---|---|---|
| T1 | the knots of `table6` (radii 2.5, 4, 6, 15, 22.5, 30 µm; mass fractions 0.05, 0.10, 0.15, 0.35, 0.25, 0.10): class edges 3.5, 6.5, 10, 21, 37.5, 52.5, 67.5 µm (diameter), cumulative mass fraction 0, 0.05, 0.15, 0.30, 0.65, 0.90, 1 | hand values | 1e-15 (edges relative) |
| T2 | the same token loaded twice: one registry entry | — | exact |
| T3 | `nowhere/psdtest-fallback.dat` resolves to `INPUT/psdtest-fallback.dat` | — | exact |
| T4 | a one-row table: d = 2 r_1 = 24 µm, and the intrinsic random stream is not advanced | `random_seed(get=)` before and after | bitwise |
| T5 | Kolmogorov-Smirnov D_n and Anderson-Darling A² (case 0) of n = 10⁵ draws, seed 42, for `table6`, Normal (20 µm, 10 µm), LogNormal (20 µm, 10 µm) and Rosin-Rammler (30 µm, 2.5) | the CDFs above | α = 10⁻³: D_n < 1.94947/√n = 6.165e-3 and A² < 5.969 |

D_n = max_i max(i/n − F(d_(i)), F(d_(i)) − (i−1)/n); A² = −n − (1/n) Σ (2i−1)[ln F(d_(i)) + ln(1 − F(d_(n+1−i)))],
with F clamped to [1e-15, 1 − 1e-15] (a draw u = 0 maps onto the lower end of the support). The
Kolmogorov points are 1.358 / 1.628 / 1.949 and the case-0 Anderson-Darling points 2.492 / 3.878 / 5.969 at
α = 0.05 / 0.01 / 0.001 [Ste74], [AD54], [MM04]; `tools/ks_ad.py` recomputes both sets. A² weights the
empirical-CDF error by 1/(F(1 − F)), so a defect confined to a tail — a truncated or clipped lower tail, a
misread end class — moves A² by orders of magnitude where D_n moves little.

The intrinsic `random_number` is processor dependent, so the draws, and the statistics, change with the
compiler: the gate holds for any seed with probability 1 − α per statistic. Seeds 43..61 are drawn and
printed (the largest D_n and A² over the twenty seeds, and how many exceed a threshold), not gated.

## Measured (ifx 2023.1)

| law | D_n seed 42 | A² seed 42 | max D_n, seeds 42..61 | max A², seeds 42..61 | above a threshold |
|---|---|---|---|---|---|
| table6 | 3.619e-3 | 1.206 | 4.342e-3 | 2.946 | 0 of 20 |
| Normal | 3.002e-3 | 1.196 | 4.755e-3 | 2.272 | 0 of 20 |
| LogNormal | 3.570e-3 | 1.479 | 4.366e-3 | 2.267 | 0 of 20 |
| Rosin-Rammler | 3.619e-3 | 1.206 | 4.342e-3 | 2.946 | 0 of 20 |

`table6` and Rosin-Rammler share their statistics: both invert their CDF with one draw u, so their
probability-integral transforms are u and 1 − u.

Sensitivity — each a one-line change of the sampler (or, last row, of the test's oracle), seed 42:

| change | D_n | A² |
|---|---|---|
| table: every draw at its class centre (discrete reading) | 0.1749 | 2710 |
| table: number fractions m_k/r_k³ taken as the weights | 0.7022 (T1 fails too) | 1.48e5 |
| table: the draw mapped onto the next class | 0.3499 | 2.89e4 |
| table: the last class edge shortened to the last radius | 0.0500 (T1 fails too) | 202.0 |
| Normal: rejection removed, a draw below 1 µm clipped to 1 µm | 0.0221 | 200.7 |
| LogNormal: σ_p taken as the log-space standard deviation | 0.5933 | 4.22e4 |
| Rosin-Rammler: x0 and n exchanged | 0.6310 | 1.75e6 |
| Rosin-Rammler: n 10 % high | 0.0246 | 199.1 |
| oracle: the Normal untruncated | 0.0231 | 184.4 |

Each change fails its law at all twenty seeds and leaves the other laws green.
