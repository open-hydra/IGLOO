# INFO — combustion/burn-box-cap (the oxide cap on a burning droplet, ODE model 8, e2e)

## Reference
`[Beck05]` (the burn law), `[Gal21]` (the exposed-surface factor of an oxide lobe on the drop, eqs. 9–14 and the
property laws of its p. 5). The oracle is written from the model equations, never from IGLOO's table.

## Model
A burning parcel carries the metal `m_Al` (Z(8)) and its oxide cap `m_ox` (Z(9)). With `xi-cap = ξ`, the cap
retains `ξν` of the metal burned (`ν = M(Al₂O₃)/(2 M(Al)) = 1.889441`); the Beckstead d^n law acts on the metal's
diameter, scaled by Gallier's exposed-surface factor `F(f_ox)` (`f_ox = m_ox/m`) and by the stretch
`c = ∫₀¹ ds/F` that keeps Beckstead's burn time; drag and heat act on the equal-volume sphere of metal and cap;
the heat capacity is `m_Al c_Al + m_ox c_res`:

    ds/dt = −K_eff c F(f_ox(s))/d₀ⁿ,  s = (d_Al/d₀)ⁿ,  m_Al = m₀ s^{3/n},  m_ox = ξν (m₀ − m_Al)
    (m_Al c_Al + m_ox c_res) dT/dt = Nu π d k_g (T_g − T) + β q_comb |ṁ_Al|,   d = (6/π (m_Al/ρ_Al + m_ox/ρ_res))^{1/3}

The gas receives `(1 − ξν)` of the metal burned.

## Fixture
`burn-box`'s case (25 parcels, d₀ = 30 µm, `kv = 1`, `kt = 1`; `bc.txt` and `solfile.tec` symlinked to
`../../burn-box/INPUT/`) with `burn-box-euler`'s output settings (`out-file = ALL`, `mollify = off`, one
trajectory row per cell entry). Its own `INPUT/phase.txt`:

    A 1 combustion=Beckstead K-burn=4.5e-07 n-burn=1.8 X-eff=0.5 T-ign=550 beta-part=0.3 q-comb=1000000 xi-cap=0.2 product=OX
    OX 0

and a two-zone absolute table (`make_table.py`, committed): zone A the metal (Cp 1250, Density 2950,
`Enthalpy_abs` 1250 T − 372687.5, zero at 298.15 K), zone OX the product, liquid alumina (Cp 1887.6, Density
3950, `Enthalpy_abs` 1887.6 T − 1.723894e7, JANAF's liquid line at 2327 K). `K_eff = 2.25e-7 m^1.8/s`,
`t_b = d₀ⁿ/K_eff = 3.210966e-2 s`; the 15 ms residence is 0.467 t_b: no parcel burns out. The log prints
`stretch c = 1.397685` for the line.

## What it verifies (`check.py`)
The oracle evaluates `F` from Gallier's lobe geometry at 2570 K by bisection in the lobe angle and takes
`c = 1.397686` (`ξ = 0.2`, `n = 1.8`, on that exact `F`); RK4 in time, `x = u_g t`.

| id | gate | tolerance | measured |
|---|---|---|---|
| W1–W3 | `out-file = ALL`, `mollify = off` read; the capped group's banner | — | pass |
| C0 | 1525 trajectory rows of 25 parcels | ≥ 1000 rows | 1525 |
| C1 | every row's `m = m_Al + m_ox` against the oracle at its x (`m₀ = ρ_Al π d₀³/6` from the inputs) | 4e-6 relative | 2.12e-6 |
| C2 | every row's `d` (the equal-volume sphere) | 4e-6 relative | 1.85e-6 |
| C3 | every row's `T` | 0.1 K | 2.8e-5 K |
| C4 | the coast: `u = u_g`, no transverse velocity | 1e-4, 1e-6 m/s | 0 |
| E3 | `euler1.tec`: per cell `ρ_p/n_p` against `(m_in + m_out)/2` of the bounding rows, `m = m_Al + m_ox` | 1e-4 relative | 1.71e-5 |
| E4, E5 | `u_p = u_g`; `T_p` inside the rows' band | 1e-11; 0.05 K | 6.4e-14; 0 |
| M1 | `Σ wdot(A) = Σ npdot (1 − ξν)(m₀ − m_Al,exit)`, `m_Al,exit` from the outlet row's `m` | 1e-5 relative | 1.27e-6 |
| M2 | `wdot(OX) = 0` in every cell; no parcel of the product (no residue is asked) | exact | pass |

At the outlet: `m/m₀ = 0.530926`, `d/d₀ = 0.771267`, `T = 704.805826 K`. The tolerances on `m` and `d` hold the
E13.6 print of `m` and `d` (2.3e-6) with the cell-crossing tolerance; the oracle's `F` and IGLOO's 513-node table
differ by at most 1.6e-6. Each parcel rides its own `(j, k)` column: the fields and the printout do not depend on
the thread count.

## RED (each one line of `rhsAlCap`)

| mutation | what fails |
|---|---|
| `F(9) = 0` (no cap) | C1 67.71 %, C2 27.98 %, C3 50.86 K |
| the burn law without `F` and `c` | C1 12.21 %, C3 18.23 K |
| `F` without the stretch (`c = 1`) | C1 18.27 %, C3 24.44 K |
| the heat capacity `m c_Al` for the whole parcel | C3 3.859 K |
| the Eulerian tail weighted by `m_Al` alone | E3 53.1 % (the outlet cell) |

A solver that does not know the `product` token stops on `IGLOO: unknown per-material key in the phase file` and
writes no output.
