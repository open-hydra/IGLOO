# evap-breakup-box-varcp — ODE model 4 with a varying heat capacity

`evap-breakup-box`'s case (ODE model 4: Reitz-Diwakar breakup with d²-law evaporation; 25 identical 30 µm drops at
40 m/s and 300 K in a uniform 50 m/s, 600 K gas, every drop consumed inside the box; `bc.txt`, `solfile.tec` and
`phase.txt` symlinked) with its own `INPUT/properties.dat`, written by `make_table.py`: 1..1000 K,
cp = 1250 + (T − 300) J/(kg K), density 2950, the enthalpy h(1 K) = cp(1 K) and then its trapezoid sum.
`out-time = on` writes the parcel time on every row. The drops heat by 291 K before they burn out at x = 0.448 m.

## What it verifies
Model 4's right-hand side (`rhsEvapBreakup`) carries its own copy of the h-state lines: with a varying cp the state
is the enthalpy and the latent term enters its rate with no heat-capacity factor,

    dh/dt = (Q + ṁ L_v)/m,  Q = π d k_g Nu (T_g − T),  ṁ = −2π d (k_g/c_p,g) ln(1 + B_T),  m = ρ_p π d³/6,

the temperature coming back through the table and the gas receiving each segment's energy at the table's enthalpy.

| id | gate | tolerance | measured |
|---|---|---|---|
| P | premises: 25 injection rows at 300 K, 40 m/s and 30 µm; every parcel with at least 100 rows, heating by at least 100 K | exact | 180 rows each; 291.4 K |
| E | Σ E of `source.tec` = Σ_inj ṁ (h(T₀) + ½\|v₀\|²) − the exits through a face, h from the table | the exit records' prints + 1e-12 of the gross | 0.002 of 2.6e-8 W |
| M | `evap-breakup-box`'s mass audit: Σ ẇ = the injected flow, no cell takes mass from the gas, the exit flows | 1e-12; −1e-6 of the largest cell; 1e-3 of the injected | 2.3e-15; 0; 1.9e-4 |
| K, runs | on every run of rows on one Reitz-Diwakar branch (4450 rows): an RK4 of dh/dt from the run's first row along the measured d (d² the cubic through the four nearest rows), the velocity between two rows the Schiller-Naumann solution from the earlier row's U, against the measured T and U of every row | the T print at the row and at the anchor, the d, U and time prints carried through the rates, the row sampling (Richardson, every row against every other row), 5e-6 K and 5e-6 m/s for the integrator | worst 0.802 of the tolerance |
| K, switches | across the interval in which the breakup branch switches (one per parcel, at x ≈ 0.0375 m): the measured T inside the bracket of the marches at constant d and slip (the larger d with the smaller slip and the reverse) | the prints and the integrator floor | 25 of 25 inside (the bracket is 0.4 K wide by construction) |

Every drop is consumed in the box, so no exit goes through a face: E then states that the gas receives the injected
energy at the table's enthalpy, which is the valuation at injection. K's velocity is the oracle's premise, as in
`khrt-varcp`.

## Falsification (a mutation of the code, measured)
| mutation | what fails |
|---|---|
| the latent term with a T-state factor (1/1250 in place of the h state's 1) | K: 7975 of 7975 rows (61 K, 1.2·10⁵ tolerances) and 25 of 25 switches; M (the drops leave through the outlet: Σ ẇ 0.0527 of 0.0773 kg/s, the exit flows 0.32 of the injected); E holds |

## Pass criterion
All gates.
