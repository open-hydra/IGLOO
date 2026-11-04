# INFO — combustion/burn-box (Beckstead d^n aluminium burn, e2e)

## Reference paper
- **Title:** Correlating Aluminum Burning Times
- **Author(s):** M. W. Beckstead
- **Year:** 2005
- **Venue:** *Combustion, Explosions, and Shock Waves*, 41(5), pp. 533–546
- **DOI:** 10.1007/s10573-005-0067-2
- **Tag:** `[Beck05]` (see [../../REFERENCES.md](../../REFERENCES.md))
- **Reproduced:** the final burn-time correlation `t_b·X_eff·p^0.1·T0^0.2 = 0.00735·D^1.8`
  (X_eff **linear**). Beckstead's Fig. 1, the experimental `t_b`-vs-`D` scatter across
  many datasets, is digitized and drawn by [../beckstead-air/](../beckstead-air/INFO.md) beside its
  burnout times; this case's `K-burn` is a verification constant, 15× below Eq. 7's at 1 atm and 300 K
  (4.5e-7 against 6.747e-6 m^1.8/s), so its parcels leave the box half burned.

## What it verifies
Uniform-gas box, `kv=1 kt=1`, `T-ign < Tp0` ⇒ burning from injection. `check.py` builds an
independent oracle from the LITERATURE correlation (not the code):

    d^n(x) = d_a^n − Keff·(x − x_a)/u_g,   Keff = K-burn·X_eff   (X_eff linear, [Beck05])
    mdot   = −(rho_p·pi·Keff/(2n))·d^(3−n),   n = 1.8

plus an independent RK4 of the lumped heat-release energy balance (Nu=2 stagnant film,
`beta·q-comb·|mdot|` release) and a mass-telescoping audit `Σwdot = Σnpdot·Δm`.

## Comparison plot
`verify.py` writes `OUTPUT/burn-box.svg` (non-gating; no-op without matplotlib): for a sample
particle it overlays the **IGLOO output** (markers) against the **Beckstead correlation curve**
from `reference/beckstead.txt` (line) for `d(x)` and the RK4 energy-balance oracle for `Tp(x)`.
A match visualises that IGLOO integrates the published `d^1.8` law and the heat-release balance.
The **experimental Fig. 1** (`t_b` vs `D`) is drawn by `beckstead-air`, beside the burnout times of
25 particles at a stated condition (air, 1 atm, 300 K) against Eq. 7.

## Pass criterion
Truncation-budget tol on `d^n`; 0.1 K on `Tp`; 1% mass telescoping; ≥ 20 verified burners.
Provenance: `make_box_case.py --kv 1.0 --kt 1.0 --rp 15e-6`; `input.ini` X-eff=0.5.
Full test matrix (CB1–CB4 + E2E) in the family overview [../INFO.md](../INFO.md).
