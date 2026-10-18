# etab-e2e-dopri5 — the ETAB product kick under the explicit solver

**Purpose.** At an ETAB breakup event the product drops receive a velocity kick normal to the path
(Tanner 1997, eqs. 8–10). The event is evaluated in the solver's output callback, which receives the
solver's state array: SDIRK4 passes a copy, DOPRI5 its own array, which is the state the segment then
continues from. The kick is computed on a copy of the velocity and added once, by the segment that ends at
the event; this case gates that under `ode-solver = H-dopri5`, the solver whose output callback sees the
live state.

**Fixture.** [`etab-e2e`](../etab-e2e/INFO.md)'s case (`INPUT/` symlinked: 25 drops across the onset
We_r = 6 and the bag-to-stripping transition, gas at 50 m/s, slip 25 m/s) with `ode-solver = H-dopri5`;
every other key is `etab-e2e`'s.

**Gates** (`check.py`).
- E `etab-e2e`'s oracles, unmodified, on this run: no breakup below the onset, the first-breakup time of the
  ORA87 oscillator, the Tanner cascade size, the kick reaching the parcel state (gate D).
- K the same input under `H-sdirk4`, run by `check.py` into `ref-sdirk4/` (the `IGLOO_EXE` executable): the
  same drops carry a kick (15 of the 20 that break), and each one's largest lateral speed √(V² + W²) over its
  trajectory rows equals the SDIRK4 run's within 1e-4 (about twenty times the F12.6 print quantum relative to
  the smallest kicked lateral speed, 1e-6 on 0.188 m/s; measured 3.2e-7).

**RED.** The kick applied to the solver's state as well as by the segment: E passes (gate D's envelope
admits a doubled kick, 0.034 of the axial speed against 0.10) and K fails, the ratio between 1.99918 and
2.00099 on all 15 kicked drops.
