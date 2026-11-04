# INFO — evaporation/boiling-zgr-switch (the switch of a boiling drop into evaporation at the boiling point, e2e)

## Reference
- **Model:** the Zuo–Gomes–Rutland boiling branch (`boiling = ZGR`) `[ZGR00]`, `[OFLEB]` above the boiling temperature
  at the gas pressure, the classical evaporation model (CEM) below it.
- **Tags:** `[ZGR00]`, `[OFLEB]` — see [../../REFERENCES.md](../../REFERENCES.md).

## What it verifies
`boiling-zgr-flash` with the gas at 1 m/s (`make_box_case.py --ug 1.0`) and d₀ = 20 µm: the drop (T₀ = 640.002 K)
flashes in the boiling branch until its mean temperature reaches the boiling temperature at the gas pressure
(599.999999959 K, x = 0.061238 m, m/m₀ = 0.925908), then cools in CEM to its wet bulb (463.83 K) and is consumed at
x = 0.086235 m. `INPUT/bc.txt` and `solfile.tec` come from `make_fixture.py`; `properties.dat` is `boiling-zgr-flash`'s.
The absolute ODE tolerance is 1e-22, below the mass (4e-12 kg at injection, 1e-15 kg at burnout): with 1e-11 the
integrator leaves the mass of the last rows uncontrolled (up to 1.4 % off at m = 4e-14 kg), and at 1e-26 or below the
step across the switch at T_B falls under the integrator's floor (`Run_ODESolver err` at x = 0.06 m).

`check.py` (inputs only, the trajectory read by name, scipy needed): T(x) and m(x) on every row against
`tools/zgr_oracle.py` — the boiling branch integrated by scipy `solve_ivp` (LSODA, rtol 1e-12) to T = T_B, then CEM
(Radau, rtol 1e-10) to m = 1e-15 kg — at t = x/u_g; T within 5e-7 K + |dT/dx| 5e-7 m + 1e-6 K — worst 0.310 of it;
m within the E13.6 half-ULP of m + |dm/dx| 5e-7 m + 1e-8 m — worst 0.967 of it. Witnesses: every parcel has rows above
and below T_B; every parcel burns out inside the box (an exit record of code 0 within 2.5e-4 m of the oracle's
burnout, 2e-6 m measured, its diameter that of m ≤ 1e-15 kg); no `non-finite state` or `Run_ODESolver err` line.

## Fails on
- a solver that refuses `boiling = ZGR` at setup: no trajectory file;
- the regime test bypassed (the boiling clamp): every parcel consumed by x = 0.025 m, T 3.7e6 times its tolerance;
- the 0.5 K floor of the superheat dropped: the drop approaches T_B without reaching it, no row below T_B, every
  parcel leaves at the outlet;
- the vapour leaving at the mean temperature (the carry dropped): T 2.9e2, m 2.1e3 times its tolerance.
