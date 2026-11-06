# INFO — evaporation/temp-relax-p2t-mucol (P2T under slip with the liquid viscosity a properties.dat column, e2e)

## Reference
- **Model:** the parabolic temperature profile inside the drop `[DS03]` (`liquid-conduction = P2T`) with the internal
  circulation of Abramzon & Sirignano `[AS89]` (the surface velocity, eqs. 30 and 32, and the conductivity factor
  χ(Pe_L), eq. 39), the liquid viscosity read from the property table's `Viscosity` column at the surface
  temperature; the gas side is Ranz-Marshall `[RM52]`. The surface balance is typed in `check.py`, so the oracle is
  the root of that balance at each printed state.
- **Tags:** `[DS03]`, `[AS89]`, `[RM52]` — see [../../REFERENCES.md](../../REFERENCES.md).

## What it verifies
The shared box (u_g = 10 m/s, T_g = 600 K, μ_g 1.8e-5, k_g 0.026, γ 1.4, R 287) and an inlet at κ_v = 0.1 and
κ_t = 0.5 (`make_box_case.py --kv 0.1 --kt 0.5 --rp 15e-6`): 30 µm drops enter at 1 m/s (9 m/s slip, Re 18, Stokes
drag) and 300 K; ρ 2950, c 1250, the evaporating right-hand side (model 2) with `Yinf = 1` (no mass transfer,
B_M = 0) and `k-liq = 0.026` on the material line, no `mu-liq`. The own `INPUT/properties.dat` (1..1000 K) carries
`Viscosity` μ_l(T) = 0.1 exp((300 − T)/60) Pa s (17 significant digits), so that Pe_L spans 20.8 to 78.2 across
χ's centre at 30. On every trajectory row (T̄, u), slip = u_g − u, Nu = 2 + 0.6 Re^1/2 Pr^1/3, and T_s is the root of

    G(T_s) = T_s − T̄ − Nu k_g (T_g − T_s)/(10 χ(Pe_L) k_l),   Pe_L = U_s d ρ_l c_l/k_l,
    U_s = (12.69/32) slip (μ_g/μ_l(T_s)) Re^1/3,   χ = 1.86 + 0.86 tanh(2.245 log10(Pe_L/30))

with μ_l(T_s) the column's linear interpolation between its kelvin nodes, read from the file. At the injection state
T_s = 364.075565 K; μ_l at the mean temperature gives 383.721503 K, μ_l = 0 (χ = 2.72) 340.584238 K.

| id | gate | tolerance | measured |
|---|---|---|---|
| M0 | the column against its closed form at mid-cell (oracle self-check) | its chord error, μ/28800 relative (+0.1 %) | 3.472e-5 (bound 3.476e-5) |
| M1 | the surface-temperature column `T<sub>s` on every row | 5e-7 K + \|∂T_s/∂T̄\| 5e-7 K + \|∂T_s/∂u\| 5e-7 m/s (the F12.6 prints, derivatives by central differences of the oracle) + 4ε T_s | worst 0.922 of it (8.8e-7 K; gfortran 0.874); Re 1.25 to 18 |
| M2 | `m = ρ π d³/6` (E13.6 half-ULP); `d = d₀`; `v = w = 0` | 5.05e-17 kg | 4.2e-17 kg |
| M3 | non-vacuity: parcels with ≥ 3 rows at Re > 1 where μ_l at T̄ and μ_l = 0 both put T_s ≥ 100 tolerances away | ≥ 20 | 25; the nearer at ≥ 1.99e4 tolerances |
| M4 | 25 exit records at the outlet (`bc` 400, x = 0.15 m); the `mu-liq tabulated from properties.dat` notice; no guard or solver-failure line; the run ends | — | 25 |

## Fails on
| mutation | what fails |
|---|---|
| the P2T branch of a group's property flags deleted (`propFlags(4)` of a non-breaking P2T material) | μ_l = 0 reaches the surface solve, which has no root: 25 `Run_ODESolver err=`, every parcel ends at its injection (bc 0, x = 0), T_s NaN on the 25 injection rows: `check_finite.py`, M1, M3, M4 |
| the material-to-group copy of the column's flag deleted | the same |
| either of the two with the surface solve's refusal of μ_l ≤ 0 removed | χ = 2.72 on every row: T_s 23.49 K low at injection, M1 on all 1525 rows, M3 |
| μ_l at the mean temperature (the column not passed to the surface solve) | T_s 19.65 K high at injection, M1 on all 1525 rows, M3 |

## Not seen here
T̄(x) is not gated: the rows are tested as states, so a right-hand side that printed the right T_s but heated the drop
with another would pass M1 (`temp-relax-p2t` and `temp-relax-p2t-kcol` gate T̄(x) at Re = 0). The breakup
viscosity as μ_l (model 4) is `test_liquid_p2t`'s leg (l).
