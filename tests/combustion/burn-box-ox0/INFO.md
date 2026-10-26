# INFO — combustion/burn-box-ox0 (an oxide fraction at injection: the initial oxide cap, e2e)

## Reference
The burn time is Beckstead's for the injected diameter, `t_b = d_p0^n/K_eff` (`[Beck05]`); the exposed-surface factor
`F(f_ox)` is Gallier et al.'s lobe (`[Gal21]`), evaluated by bisection in the lobe angle (`burn-box-cap`'s oracle).
The state at injection, the cap's growth, the residue and the books follow from the model's definitions and the
inputs.

## Model
`xi-ox0 = x0` on the material line of a burning, capped material (ODE model 8) is the oxide's mass fraction of the
injected particle, `m_ox/(m_Al + m_ox)`. The inlet's diameter `d_p0` and mass flow describe the whole particle, metal
and oxide at the equal volume:

    m_p0  = (π/6) d_p0³ / ((1 − x0)/ρ_Al + x0/ρ_res),   m_Al0 = (1 − x0) m_p0,   m_ox0 = x0 m_p0,   ṅ = ṁ/m_p0
    m_ox  = m_ox0 + ξν (m_Al0 − m_Al),   f_ox = m_ox/(m_Al + m_ox)
    ds/dt = −K_eff c F(f_ox(s))/d_Al0^n,   s = (d_Al/d_Al0)^n
    c     = (d_Al0/d_p0)^n ∫₀¹ ds/F(f_ox(s)) = 1.388383   (ξ 0.2, n 1.8, x0 0.1)

so the particle burns out in `t_b = d_p0^n/K_eff` whatever `x0`, and the residue carries the initial oxide:
`m_res = m_ox0 + ξν (m_Al0 − m_Al,b)`. The parcel's first source segment brings `ṁ (h_mix(T₀) + ½u²)` into the books,
`h_mix = (1 − x0)(c_Al T + hOff_Al) + x0 (c_res T + hOff_res)`.

## Fixture
`burn-box-residue`'s case with `xi-ox0=0.1` on the burning line and an inert product (no model on its line, so the
files carry no `f` column): `burn-box`'s `bc.txt` and `solfile.tec` and `burn-box-cap`'s absolute two-zone table
(symlinked); `out-file = S`, `out-time = on`, `mollify = off`. Its own `INPUT/phase.txt`:

    A 1 combustion=Beckstead K-burn=2.0e-6 n-burn=1.8 X-eff=0.5 T-ign=550 beta-part=0.3 q-comb=1000000 xi-cap=0.2 product=OX residue=on xi-ox0=0.1
    OX 0

`m_p0 = 4.278788e-11 kg` (`ρ_mix` 3026.62 kg/m³), `m_Al0 = 3.850909e-11 kg`, `d_Al0 = 2.921332e-5 m`. The 25
parcels burn out near x = 0.072 m at 765.25 K; each residue (`m = 1.883055e-11 kg`, `d = 2.0881e-5 m`, 0.696 of
`d_p0`) relaxes to 678.16 K at the outlet.

## What it verifies (`check.py`)

| id | gate | tolerance | measured |
|---|---|---|---|
| W | the log: `xi-ox0= 1.000E-01`; `stretch c =  1.388382 from xi-ox0 = 0.1000` read back against 1.388383; the capped banner; `combustion product (residues: 25)` | 3e-6 relative | 7.2e-7 |
| I | each burning parcel's injection row: `m = m_p0`, `d = d_p0` | 2e-6 relative; exact | 4.7e-7; exact |
| C | every trajectory row of A before its parcel's exit record against the oracle at the row's `t` (RK4 on `(s, T)`, Gallier's exact `F`, `c` typed): `m = m_Al + m_ox`, the equal-volume `d`, `T`; the coast `u = u_g`, `v = w = 0` (725 rows) | 4e-6, 4e-6 relative; 0.1 K | 2.6e-6, 2.1e-6; 3.5e-5 K |
| B | each parent's exit record (bc 0) at a time inside `[t(s_b), t_b] = [7.1974437e-3, 7.2246741e-3] s`, `s_b = (1e-15 kg/m_Al0)^{n/3}`: the burnout fires on the first accepted step with `m_Al ≤ 1e-15 kg`, and no step with `m_Al < 0` is accepted | 3e-6 relative at each end | 7.1985146e-3 to 7.2002811e-3 s |
| R1 | one residue per parent burned out inside the box: IDs 1–25 = the parents' exit records with `X < Lx − 1e-6`; the product's banner 25 | exact | 25 |
| R2 | each residue's birth row: `m = m_ox0 + ξν(m_Al0 − 1e-15 kg) = 1.883055e-11 kg`; `d = (6m/(π ρ_res))^{1/3}` of the printed `m` | 3e-5; 5e-6 relative | 2.9e-6; 1.6e-6 |
| R3 | `Σ wdot(A) = Σ ṁ(parents) − Σ ṁ(residues)` (exit records); `wdot(OX) = 0` | 1e-5 relative | 5.2e-7 |
| R3e | `Σ E = Σ ṁ_p (h_mix(600 K) + u²/2) − Σ ṁ_r (c_res T_exit + hOff_res + u²/2)`, `h_mix(600 K) = −1.271057e6 J/kg` | 1e-5 relative | 6.2e-8 |
| R5 | each residue's `ṁ = (m_res/m_p0)` × its parent's | 1e-5 relative | 3.3e-6 |
| S | the scatter quanta of A and OX: `dNscat = (Σ_p ṁ_p/m_p0)(L_ref/u_g)/(100 × 25)`, `L_ref` the box diagonal | 1e-3 relative | 5.7e-5 |

Each parcel and its residue ride one `(j, k)` column: the printout is the same at one and at five threads.

## RED

| mutation | what fails |
|---|---|
| the stretch not re-evaluated (`capStretch(ξ, n)`, 1.397685) | W (6.7e-3); C1 3.2e-3, C3 0.74 K; B: the burnouts at 7.1518e-3 s, 0.63 % before the window |
| the factor `(d_Al0/d_p0)^n` dropped (the metal burning in its own Beckstead time) | W (4.9e-2); C1 2.2e-2, C3 5.3 K; B: 6.8636e-3 s, 4.6 % before the window |
| the split on `(π/6) ρ_Al d_p0³` | I 2.5e-2; C1 3.2e-2; B 1.5 % early; R2 2.5e-2; R5 2.5e-2 |
| the split skipped (the whole `m_p0` as metal, the cap from zero) | C1 0.14, C3 32.9 K; B 1.9 % late; R2 0.14; R3e 0.35 (I holds: the injection row is unchanged) |
| the number rate left at `ṁ/((π/6) ρ_Al d_p0³)` | R3 4.7e-2; R3e 5.6e-3; R5 2.5e-2 |
| the scatter pre-pass's number rate on `(π/6) ρ_Al d_p0³` | S 2.6e-2 on both clouds |
| the `xi-ox0` line missing from the log | W |

A solver that does not know the `xi-ox0` token stops on `IGLOO: unknown per-material key in the phase file`.
