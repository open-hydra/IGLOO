# INFO — evaporation/boiling-zgr-flash (the adiabatic flash of a superheated drop, e2e)

## Reference
- **Model:** the Zuo–Gomes–Rutland boiling branch (`boiling = ZGR`) `[ZGR00]`, its boiling heat-transfer coefficient
  as `[OFLEB]` codes it; the vapour leaves at the boiling temperature.
- **Tags:** `[ZGR00]`, `[OFLEB]` — see [../../REFERENCES.md](../../REFERENCES.md).

## What it verifies
Water drops (ρ 958, c 4186, L_v 2.26e6 J/kg, M_v 18, d₀ 50 µm) enter the shared box at the gas velocity (κ_v 1: Re 0,
Nu 2) at T₀ = κ_t T_g = 640.002 K (`bc.txt` carries κ_t = 1.06667) in the 600 K, 2.0664 bar gas, with `evaporation =
CEM` and `boiling = ZGR`. `boiling-temperature = 551.7783324` puts the boiling temperature at the gas pressure at
599.999999959 K, the gas temperature to 4e-8 K (A = 1.8e-11): the branch draws no heat from the gas, the superheat
flashes the liquid through α_s, and the drop cools toward T_B, which it does not reach inside the box.
`INPUT/bc.txt` and `properties.dat` come from `make_fixture.py`; `INPUT/solfile.tec` is `drag-stokes`'s.

`check.py` takes T₀ and d₀ from `bc.txt`, the rest from `input.ini` and `tools/make_box_case.py`'s gas, and reads the
trajectory by its header names (`tests/tools/datcols.py`). It needs scipy.

- (i) on every row the flash identity with the vapour leaving at T_B, `m/m₀ = (L_v − c (T₀ − T_B))/(L_v − c (T − T_B))`
  (the heat-fed mass over the box is 2e-11 of m₀), within the E13.6 half-ULP of m + m₀ |r′(T)| 5e-7 K (T's F12.6
  print) + 1e-8 m — worst 0.967 of it;
- (ii) T(x) against the boiling branch integrated by `tools/zgr_oracle.py` (scipy `solve_ivp`, LSODA, rtol 1e-12) at
  t = x/u_g, within 5e-7 K + |dT/dx| 5e-7 m + 1e-6 K (3.2e-3 K at injection, 9.5e-6 K at the outlet) — worst 0.122 of it;
- (iii) 25 parcels, 1525 rows, every exit record at the outlet (x = 0.15, code 400), no `non-finite state` or
  `Run_ODESolver err` line in `run_out.txt`.

## Fails on
- a solver that refuses `boiling = ZGR` at setup: no trajectory file;
- the regime test bypassed (the boiling clamp, CEM): T at the outlet 468.07 K against 604.58 K, (ii) 1.5e7 times its
  tolerance, (i) 1.2e5 times;
- the vapour leaving at the mean temperature (the carry dropped): (i) 3.3e3 times, (ii) 3.5e3 times its tolerance;
- the exponent 2.33 of α_s typed 2.3: (ii) 1.3e4 times its tolerance.
