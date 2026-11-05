# INFO — combustion/burn-box-cap-gallier (the oxide cap under Gallier's deposition law, ODE model 8, e2e)

## Reference
`[Gal21]`: the deposition law `f_dep = g(d_Al) h(f_ox)` (PDF pp. 13–14), its data Figs. 7 and 11 digitized in
`reference/` (`PROVENANCE.md`: the two reads, their spread, the check of the law against Gallier's residue sizes);
the exposed-surface factor of the lobe (eqs. 9–14) as in `burn-box-cap`; `[Beck05]` (the burn law).

## Model
`burn-box-cap`'s parcels with `xi-law=Gallier21`: the cap retains `ξ(d_Al, f_ox) = g(d_Al) h(f_ox)` (clamped to
[0, 1]) of the oxide the burn forms, at the metal's current diameter and the cap's mass fraction; `h` Fig. 11's
points over 100, linear, held past `f_ox = 0.9`; `g` the least-squares line through Fig. 7's points over its value at
70 µm, held at its 20 and 400 µm values outside. The cap's mass is no function of the metal's alone, so the state
carries it (`y = m_ox/m₀`, `x = m_Al/m₀ = s^{3/n}`):

    ds/dt = −K_eff c F(f_ox)/d₀ⁿ,   dy/dt = −ν ξ(d₀ s^{1/n}, f_ox) dx/dt,   f_ox = y/(x + y)
    (m_Al c_Al + m_ox c_res) dT/dt = Nu π d k_g (T_g − T) + β q_comb |ṁ_Al|

Each parcel takes its own stretch `c = ∫₀¹ ds/F(f_ox(s))` on the law's path from its injected size, so its burn time
stays `d₀ⁿ/K_eff`.

## Fixture
`burn-box-cap`'s case (`burn-box`'s `bc.txt` and `solfile.tec`, `burn-box-cap`'s absolute two-zone table, all
symlinked), `out-file = S`, `mollify = off`. Its own `INPUT/phase.txt`:

    A 1 combustion=Beckstead K-burn=4.5e-07 n-burn=1.8 X-eff=0.5 T-ign=550 beta-part=0.3 q-comb=1000000 xi-law=Gallier21 product=OX
    OX 0

25 parcels of 30 µm; the 15 ms residence is 0.467 of Beckstead's time: no parcel burns out. Over it `d_Al` falls to
18.8 µm (`g` from 1.101 to its held 20 µm value) and `f_ox` rises to 0.407 (`ξ` 0.195 at the outlet).

## What it verifies (`check.py`)
`law_path.py` builds the law from the committed tables (never from IGLOO) and the path's stretch,
`c = 1.340400` at 30 µm (Gallier's `F` by bisection, 4096 cells; IGLOO's 513-node table gives the same within
6.2e-7, `test_combustion` CB15). The oracle integrates `(s, y, T)` by RK4 in time at 1e-6 s, `x = u_g t`.

| id | gate | tolerance | measured |
|---|---|---|---|
| W1, W2 | the capped group's banner and the cap's notice name the law; the notice names the stretch per parcel | — | pass |
| C0 | 1525 trajectory rows of 25 parcels; the oracle at 1e-6 s against 2e-6 s, `m` and `d` | ≥ 1000 rows; 1e-8 | 1525; 3.5e-10 |
| C1 | every row's `m = m_Al + m_ox` against the oracle at its x | 4e-6 relative | 2.34e-6 |
| C2 | every row's `d` (the equal-volume sphere) | 4e-6 relative | 2.34e-6 |
| C3 | every row's `T` | 0.1 K | 2.2e-5 K |
| C4 | the coast: `u = u_g`, no transverse velocity | 1e-4, 1e-6 m/s | 0 |
| M1 | `Σ wdot(A) = Σ npdot (m₀ − m_exit)`, `m_exit` the outlet row's `m` | 1e-5 relative | 1.82e-6 |
| M2 | `wdot(OX) = 0` in every cell; no parcel of the product (no residue is asked) | exact | pass |

At the outlet: `m/m₀ = 0.413454`, `d/d₀ = 0.718455`, `T = 714.657 K`. The tolerances on `m` and `d` are
`burn-box-cap`'s: the E13.6 prints with the cell-crossing tolerance and the lobe table; the law adds no term (the
code and `law_path.py` interpolate the same nodes linearly).

## RED

| mutation | what fails |
|---|---|
| the stretch left at 1 (no per-parcel stretch) | C1 22.7 %, C2 8.7 %, C3 25.5 K |
| `g` read at the injected diameter for the whole burn | C1 6.3e-3, C2 1.8e-3, C3 0.18 K |
| the law read at the drop's diameter in place of the metal's | C1 1.1e-3, C2 3.2e-4 |

A solver that does not know the `xi-law` token stops on `IGLOO: unknown per-material key in the phase file`.
