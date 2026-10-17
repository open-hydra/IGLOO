# d2law-varcp — evaporation under a varying cp (the enthalpy state of ODE model 2)

**Purpose.** A tabulated, varying cp makes the particle state the enthalpy: the latent term enters
dh/dt = (Q + ṁ L_v)/m with no 1/c_p factor, the temperature comes back through the table, and the gas receives each
segment's mass at the table's enthalpy. `temp-relax-varcp` covers the enthalpy state without mass exchange; this case
covers it under evaporation.

**Fixture.** `d2law`'s case (`bc.txt`, `solfile.tec` and `phase.txt` symlinked): the shared box, 25 drops of 30 µm,
κ_v = 1 (Re = 0, x = u_g t), T₀ = 300 K in a 600 K gas, d2-law evaporation (L_v = 2·10⁵), with its own
`INPUT/properties.dat` written by `make_table.py`: 250..900 K, cp = 1000 + 2(T − 250) J/(kg K), density 2950, the
enthalpy its trapezoid sum from h(250 K) = 2.5·10⁵ J/kg (h = 1000 T + (T − 250)²).

**Oracle** (`check.py`). With Nu = Sh = 2:

    dh/dt = (Q + ṁ L_v)/m,  Q = 2π d k_g (T_g − T),  ṁ = −2π d (k_g/c_p,g) ln(1 + B_T),  B_T = c_p,g (T_g − T)/L_v

T = T(h) through the table (linear between its nodes); every constant is a case input.

**Gates.**
- V1 `d2law`'s own mass-rate oracle along the measured T (its `check.py`, unchanged), with its mass-source audit.
- V2 the enthalpy history of every parcel: an RK4 of dh/dt along its measured d (d²(t) the cubic through the four
  nearest rows), T = T(h) by the check's own inversion of the table, against the measured T on every row, within
  5e-7 K (the T print) + the row sampling (Richardson, all rows against every other row) + the d print carried
  through the rate + 1e-5 K (the integrator: the parcels' histories differ by up to 6e-6 K, the table's nodes being
  kinks of T(h)).
- V3 the energy handed to the gas: ΣE = Σ_p ṅ[m(h + |v|²/2)]_injection − [..]_exit with h from the table. The 25
  parcels share one history (asserted within 5e-5 K), so the exit flow is Σ_p ṅ m₀ − Σẇ (V1's audited mass source)
  at the exit records' T and u_g; within the T print and spread of the exit records and the sums' 15 digits.

**Measured.** V1 worst 0.34 of its tolerance; V2 worst 0.67 (1.4e-5 K at the first row: the d print through the rate);
V3 1.6e-8 W (0.004 of its tolerance, 3e-12 of the gross 5.24·10³ W).

**RED.**
- The latent term with a T-state factor (1/1250 in place of the h-state's 1): V2 fails on 1475 rows, by up to 22.4 K
  (3.2e5 tolerances); V1 holds (it reads the measured T). The factor of the material's own scalar cp,
  1/`aux(ind_cp)`, divides by the zero cp of a tabulated material: every parcel fails at injection and V1 fails.
- The deposit's enthalpy forced to cp·T + hOff (the constant-cp form): V3 fails by 4.1·10³ W (1.1e9 tolerances);
  `temp-relax-varcp`'s energy gate fails too, the constant-cp cases hold.
