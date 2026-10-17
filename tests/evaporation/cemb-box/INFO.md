# cemb-box — the CEM-B evaporation model end to end

**Purpose.** `evaporation = CEM-B` is the classical model (Spalding transfer number, Ranz-Marshall Sherwood number)
with the gas properties of the film at the 1/3-rule temperature T_f = T + (T_g − T)/3. Its film chain is pinned at
Re = 0 by `test_evaporation` (EV3); this case runs it end to end, with a convective Sherwood number.

**Fixture.** `lk-neq`'s box, properties and droplet (`solfile.tec`, `properties.dat`, `phase.txt` symlinked): 20 µm,
ρ_p 2950, T₀ = 300 K in a 600 K gas, L_v 2·10⁵, M_v 100, T_b 900, Le 1, Y_∞ 0; an equilibrium interface
(`interface = VLE`); its own `bc.txt` (`make_box_case.py INPUT --kv 0.5 --kt 0.5 --rp 1e-5`): the drops enter at half
the gas speed, so the slip starts at 5 m/s (Re 6.7) and relaxes under Stokes drag. `out-time = on`.

**Oracle** (`check.py` → `film_check.py`, shared with `asm-box`). The surface state is the equilibrium one (Clausius-
Clapeyron, p = ρ_g R_g T_g); the rate law is CEM-B as the theory page documents it (`docs/theory/evaporation.md`:
film density by the ideal gas, viscosity and conductivity by a 0.7 power law, c_p unchanged, Re_f corrected for the
film's density and viscosity, Sh_f = 2 + 0.6 Re_f^½ Sc_f^⅓, ṁ = −π d ρ_f D_f Sh_f ln(1 + B_M)). Hubbard, Denny and
Mills' primary text is not in hand and the page follows the code, so this leg pins the documented film rule.

**Gates.**
- G0 every first row is the injection state.
- K1 d² on every row against an RK4 of d(d²)/dt = 4ṁ/(πρ_p d) along the measured T and slip (linear in t between
  rows, Re from the slip and the current d), within `lk-neq`'s error model: the d print at the row and the anchor,
  the Richardson estimate of the row sampling, the T and u prints through the rate, 1e-7 d₀².
- N1 non-vacuity: the CEM curve (far-field properties) along the same measured state misses every parcel by more
  than 1.5 tolerances.

**Measured.** K1 worst 0.33 of its tolerance; N1: CEM misses by ≥ 43 tolerances; 25/25 parcels.

**RED.** The CEM-B selector calling the classical model: K1 fails by 44.5 tolerances (and `test_evaporation`'s EV3).
