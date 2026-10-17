# asm-box — the Abramzon-Sirignano evaporation model end to end

**Purpose.** `evaporation = ASM` is the Abramzon-Sirignano model: Frössling's 0.552 in Sh₀ and Nu₀, the film
corrections F(B_M) and F(B_T), the φ–B_T loop, and its own gas-side heat in place of the drag-side convective heat.
Its sign and the heat override are pinned by `test_evap_probes` (XE3); this case runs it end to end, with a
convective Sherwood number, on the mass and on the temperature.

**Fixture.** `cemb-box`'s (`lk-neq`'s box, properties and droplet, an equilibrium interface, the drops at half the gas
speed: Re 6.7 at injection) with its own `bc.txt` (`make_box_case.py INPUT --kv 0.5 --kt 0.535 --rp 1e-5`): the
injection temperature (321 K) is chosen so every state stays at B_M ≥ 1e-2, outside the near-equilibrium branch.
`out-time = on`.

**Oracle** (`check.py` → `../cemb-box/film_check.py`). The surface state is the equilibrium one (Clausius-Clapeyron,
p = ρ_g R_g T_g). The rate law is typed from `[MHB98]`'s statement of the model (eqs. 4, 6, 11–13 and Table 1):
Sh₀ = 2 + 0.552 Re^½ Sc^⅓, Nu₀ = 2 + 0.552 Re^½ Pr^⅓, F(B) = (1 + B)^0.7 ln(1 + B)/B, Sh* = 2 + (Sh₀ − 2)/F(B_M),
ṁ = −π d ρ_g D_v Sh* ln(1 + B_M), B_T = (1 + B_M)^φ − 1 with φ = (c_p,v/c_p,g)(Sh*/Nu*)/Le and
Nu* = 2 + (Nu₀ − 2)/F(B_T), iterated to convergence. The gas-side heat into the drop is Abramzon and Sirignano's
own form, Q = −ṁ c_p,v (T_g − T)/B_T with the vapour's c_p (their eq. 20, as the code and the theory page take it;
MHB98's restatement writes the gas c_p there, a factor 2 with these inputs), and m c_L dT/dt = Q + ṁ L_v.

**Gates.**
- G0 every first row is the injection state.
- K1 d² on every row against an RK4 of d(d²)/dt along the measured T and slip, as in `cemb-box`.
- K2 T on every row against an RK4 of m c_L dT/dt = Q + ṁ L_v along the measured d² and slip (each the cubic in t
  through the four nearest rows), within the T print, the Richardson estimate of the row sampling, the d print
  through the rate and 1e-6 K.
- B the minimum B_M over every row of every parcel ≥ 1e-2 (printed), and B_T > 1e-3: the loop's fallback never fires.
- N1 non-vacuity: the CEM curve misses every parcel by more than 1.5 tolerances.

**Measured.** B_M ≥ 1.37e-2 (its injection value) and B_T ≥ 2.7e-2 on every row; K1 worst 0.33, K2 worst 0.08 of the
tolerances; N1: CEM misses by ≥ 15.8 tolerances; the loop's five iterations against its converged B_T: 1.1e-9; 25/25.

**RED.**
- Sh₀'s 0.552 → 0.6: K1 fails by 11 tolerances, K2 by 177 (0.75 K).
- Nu₀'s 0.552 → 0.6: K2 fails by 682 tolerances (2.66 K); K1 holds (the heat path alone).
