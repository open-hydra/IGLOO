# khrt-varcp — KHRT breakup with a varying heat capacity (the enthalpy state of ODE model 3)

`khrt-e2e`'s case (KHRT breakup, ODE model 3: 25 water drops on the 25 inlet cells, a uniform 200 m/s, 300 K gas
along x) with `khrt-varrho`'s `bc.txt` (`kT` 0.9: the drops enter at 270 K and 100 m/s) and its own
`INPUT/properties.dat`, written by `make_table.py`: 1..1000 K, cp = 4182 + 10 (T − 270) J/(kg K), density 1000,
the enthalpy h(1 K) = cp(1 K) and then its trapezoid sum (exact for a linear cp). `out-time = on` writes the parcel
time on every row. The 25 parents shed 241 KH children, which heat by up to 17 K toward the gas.

## What it verifies
With a tabulated, varying cp the particle state is the enthalpy: the heat enters dh/dt with no 1/c_p factor, the
temperature comes back through the table, a KH shed hands its child the parent's enthalpy, and the gas receives each
segment's energy at the table's enthalpy. The model-3 right-hand side (`rhsBreakupOnly`) carries its own copy of the
h-state lines; `evaporation/evap-breakup-box-varcp` covers model 4's.

    dh/dt = π d k_g Nu (T_g − T)/m = 6 k_g Nu (T_g − T)/(ρ_p d²),  Nu = 2 + 0.6 Re^½ Pr^⅓,  Re = ρ_g (u_g − U) d/μ_g

| id | gate | tolerance | measured |
|---|---|---|---|
| P | premises: 25 injection rows at 270 K and 100 m/s; at least 200 parcels with 5 rows; the median parcel heats by at least 1 K | exact | 25 rows; 263 parcels; median 13.707 K |
| E | Σ E of `source.tec` = Σ_inj ṁ (h(T₀) + ½\|v₀\|²) − Σ_exit ṁ (h(T) + ½\|u\|²), h from the table | the exit records' prints (ṁ E13.6, T and \|u\| F12.6) + 1e-12 of the gross | 0.268 of 0.58 W (gross 4.86·10⁵ W) |
| M | `khrt-e2e`'s exit mass-flow total; its mass and momentum budget of `source.tec` | 1e-5; 1e-12 and 5e-6 | 6.6e-7; 0 and 6.4e-7 |
| K, runs | on every run of rows at one printed d (no breakup between them; 8373 rows): an RK4 of dh/dt from the run's first row, the velocity between two rows the Schiller-Naumann solution from the earlier row's U, against the measured T and U of every row | the T print at the row and at the anchor, the d, U and time prints carried through the rates (the march repeated at the prints' extremes), 5e-6 K and 5e-6 m/s for the integrator | worst 0.915 of the tolerance; no row exceeds the print terms by more than 7e-7 |
| K, intervals | across every interval in which d changes (an RT event, a KH shed, the stripping; 1977 intervals): the measured T at its end inside the bracket of two marches from its start at constant d and slip, the larger d with the smaller slip and the smaller d with the larger slip | the prints and the integrator floor | 1977 of 1977 inside |

The tight check of K lives on the runs at one d: the children and the parents between events. A parent under
continuous KH stripping changes d on every interval and is checked by the bracket alone. The bracket is a bound by
comparison: the heating rate falls with d and rises with the slip, and both only fall along a path. K's velocity is
the oracle's premise, not the property under test: a failure of its U part alone means the Schiller-Naumann solution
no longer reproduces the rows, not that the enthalpy state is wrong.

## Falsification (each a mutation of the code, measured)
| mutation | what fails |
|---|---|
| the KH shed's enthalpy in the constant-cp form, cp T + h_off | E: Σ E 4.92·10⁴ W against −9.73·10³ W (1.0·10⁵ tolerances); P, M, K hold |
| the child's enthalpy slot set to its temperature | P (median heating 0), E (1.1·10⁵ tolerances), K: 116 rows (a child at 0.18 K against 270 K) and 240 intervals |
| model 3's heat rate with a T-state factor, 1/4182 | P (median heating 0.004 K), K: 8373 of 8373 rows (14.3 K, 1.9·10⁵ tolerances) and 1977 of 1977 intervals; E holds (the deposit telescopes whatever the rate) |
| model 3's heat rate divided by the material's scalar cp slot, which a tabulated cp does not have | the drops do not heat: P, K 8373 of 8373 rows |

## Pass criterion
All gates.
