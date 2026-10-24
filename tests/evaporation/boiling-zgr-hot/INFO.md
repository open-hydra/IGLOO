# INFO — evaporation/boiling-zgr-hot (boiling of a superheated drop in a hotter gas, e2e)

## Reference
- **Model:** the Zuo–Gomes–Rutland boiling branch (`boiling = ZGR`) `[ZGR00]`, `[OFLEB]`: the flash fed by the superheat
  and the boiling fed by the heat the gas conducts through the blown film, the film conductivity at the 1/3-rule film
  temperature, the vapour leaving at the boiling temperature.
- **Tags:** `[ZGR00]`, `[OFLEB]` — see [../../REFERENCES.md](../../REFERENCES.md).

## What it verifies
`boiling-zgr-flash` with T₀ = κ_t T_g = 490.0002 K (`bc.txt` carries κ_t = 0.816667) and `boiling-temperature =
422.3191381`: the boiling temperature at the gas pressure is 449.999999987 K, 150 K below the gas (A = 0.066670). At the
outlet the drop has shed 35.16 K of superheat and kept m/m₀ = 0.829109; the flash alone would leave 0.9343.
`INPUT/bc.txt` comes from `make_fixture.py`; `solfile.tec` is `drag-stokes`'s, `properties.dat` `boiling-zgr-flash`'s.

`check.py` (inputs only, the trajectory read by name, scipy needed): T(x) and m(x) on every row against the boiling
branch integrated by `tools/zgr_oracle.py` (scipy `solve_ivp`, LSODA, rtol 1e-12) at t = x/u_g; T within 5e-7 K +
|dT/dx| 5e-7 m + 1e-6 K — worst 0.627 of it; m within the E13.6 half-ULP of m + |dm/dx| 5e-7 m + 1e-8 m — worst 0.598
of it; 25 parcels, no row below T_B, every exit record at the outlet, no `non-finite state` or `Run_ODESolver err` line.

The step of α_s at a superheat of 5 K (−0.59 %, the published laws kept) leaves the integrated temperature 5.2e-6 K
below the oracle's after the drop crosses it (x ≈ 0.135 m; 1.5e-6 K or less before), which the row tolerance holds by
its x-print term.

## Fails on
- a solver that refuses `boiling = ZGR` at setup: no trajectory file;
- the regime test bypassed (the boiling clamp, CEM): T 381.46 K at the outlet against 454.84 K (8.9e6 times its
  tolerance), m 1.5e5 times;
- the vapour leaving at the mean temperature (the carry dropped): T 4.2e4, m 2.9e3 times;
- the film conductivity replaced by the gas's (k_f = k_g): T 4.5e3, m 1.2e4 times;
- the blowing of the film dropped (G = B ln(1 + A)): T 9.0e3, m 2.0e2 times.
