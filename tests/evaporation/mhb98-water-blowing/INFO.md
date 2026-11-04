# mhb98-water-blowing — the blowing correction in the right-hand side, and the wet-bulb identity with f₂

**Purpose.** `[IGLOO-Models] blowing = LK` multiplies the droplet's convective heat by MHB98's Stefan-blowing factor
f₂ = β/(e^β − 1), β = −(3/2) Pr τ_d ṁ/m, τ_d = ρ_l d²/(18 μ_G) (`[MHB98]` eqs. 17–19), in models 2 and 4 (model 4's call
site is pinned by `test_evap_breakup` U4). No other case sets it. On the wet-bulb plateau of a quiescent droplet
(Re = 0, Nu = Sh = 2, dT/dt ≈ 0) the heat balance f₂·2πd k (T_G − T_wb) = −ṁ L_v gives

    K = −d(d²)/dt = f₂ (T_G − T_wb) 8 k / (ρ_l L_v),

and on the d² line ṁ/m = −(3/2) K/d², so β = Pr ρ_l K/(8 μ_G) — MHB98's eq.-28 group, computed from the fitted K,
never from the code.

**Fixture.** `INPUT/` is `mhb98-water`'s (symlink): the MHB98 Fig. 2 water droplet (D₀ 1.05 mm, T_d0 282 K in dry air
at 298 K, Re_d = 0, CEM + LK, MHB98's Appendix-A properties). `input.ini` is `mhb98-water`'s with `blowing = LK`.

**Gates** (`check.py`; it runs `mhb98-water`'s `check.py` on this run, which reads `blowing = LK` from this
`input.ini`).
- W1 `mhb98-water`'s gates — the CEM + LK d² kernel along the measured T, the mass telescoping and closure, the M7
  reproduction — and its wet-bulb identity with f₂, K and T_wb from the last half of the rows, within 5e-4.
- W2 non-vacuity: the identity without f₂ misses by more than 2e-3 on every parcel (the blowing acted).
- W3 a sibling `ref-noblow/` runs `mhb98-water`'s input; a coupled (d², T) RK4 of the droplet (CEM + LK rate of
  `mhb98-water`'s kernel, m c_p dT/dt = f₂ Q + ṁ L_v, f₂ typed from eqs. 17–19), integrated at the rows' times, gives
  the shifts blowing causes. IGLOO's K_blow/K_none − 1 must match the kernel's within 2e-5 and T_wb,blow − T_wb,none
  within 1e-4 K. On `mhb98-water` the kernel reproduces IGLOO's K to 5.5e-7 and T_wb to 1.4e-7 K; the bands are 20×
  and ~700× those.
- W4 the d² kernel of W1 takes the temperature from the rows, so it cannot see f₂: a control.

**Measured.** W1 4.1e-5 (25/25); W2 3.1e-3; W3 K_blow/K_none − 1 = −1.612492e-3 against the kernel's −1.612584e-3
(9.2e-8) and T_wb shift −0.024036 K against −0.024036 K (5.9e-9 K).

**RED.**
- The model-2 site never applying f₂ (`blowSelect == 2` in `rhsEvaporation`): W1 3.2e-3 (0/25), W2 4.4e-5, W3 off by
  1.6e-3 and 0.024 K.
- b doubled in `blowingFactor`: W1 3.1e-3, W3 −3.22e-3 and −0.048 K (`test_evap_breakup` U4 fails too).
- f₂ applied while blowing is off: this case fails W1–W3 and `mhb98-water`'s identity fails by 3.1e-3, which a 2 %
  tolerance would let through.
