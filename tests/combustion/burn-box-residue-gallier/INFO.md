# INFO — combustion/burn-box-residue-gallier (the burnout residue under Gallier's deposition law, e2e)

## Reference
`[Gal21]`: the deposition law (PDF pp. 13–14, its data in `../burn-box-cap-gallier/reference/`) and the residue
sizes it gives: `β = d_res/d_Al0` 0.55 at 70 µm (Table 3, no initial lobe), 0.49 at 200 µm, 0.40 at 400 µm and 0.57
for a small drop (PDF p. 14). The burnout residue's model is `burn-box-residue`'s.

## Model
Under the law the cap grows in the metal mass as `dm_ox/dm_Al = −ν ξ(d_Al, f_ox)`, `d_Al = d_Al0 (m_Al/m_Al0)^{1/3}`:
no time, burn law or stretch enters, so the residue's mass over the injected one, `z = m_ox/m₀` at burnout, depends
on the law and the injected size alone, and so does `β³ = z ρ_Al/ρ_ox`. A parcel burns out at the first accepted
step with `m_Al ≤ m_b = min(1e-6 m₀, 1e-20 kg)` (`tools/burnout.py`; 1e-20 kg at every size here): its residue lies
between the path's cap at `m_Al = m_b` and at `m_Al = 0`, its burnout time between the path's time to `m_b` and Beckstead's `d₀ⁿ/K_eff` (the stretch per parcel keeps every burn
time on it).

## Fixture
`burn-box-residue`'s case with the law and a faster burn. Its own `INPUT/phase.txt`:

    A 1 combustion=Beckstead K-burn=1.5e-4 n-burn=1.8 X-eff=0.5 T-ign=550 beta-part=0.3 q-comb=1000000 xi-law=Gallier21 product=OX residue=on
    OX 0 solidification=on T-melt=2327 h-fus=1.08949e6 cp-solid=1360

and its own `INPUT/bc.txt`: `burn-box`'s with the inlet radius set by the inlet row `j` (10, 35, 100, 200 and 50 µm:
diameters 20, 70, 200, 400 and 100 µm, five parcels each); `burn-box`'s `solfile.tec` and `burn-box-cap`'s table
(symlinked); `out-file = S`, `out-time = on`, `mollify = off`. `K_eff = 7.5e-5 m^1.8/s`: the parcels burn out at
x = 0.46, 4.4, 8.4, 29.3 and 102 mm (20 to 400 µm), inside the 0.15 m box.

## What it verifies (`check.py`)
`../burn-box-cap-gallier/law_path.py` builds the law from the committed tables (never from IGLOO) and integrates the
path in `s = (d_Al/d₀)ⁿ` (RK4, 4096 cells): `z` at `m_b` and at 0, and the fraction of `d₀ⁿ/K_eff` spent to `m_b`.
The books are `burn-box-residue`'s.

| id | gate | tolerance | measured |
|---|---|---|---|
| R1 | one residue per parcel burned out inside the box (`X < Lx − 1e-6`); the product's banner counts 25 | exact | 25 |
| R2 | each residue's mass at birth over `m₀ = ρ_Al π d₀³/6` (`d₀` the inlet size) inside `[z(m_b), z(0)]` of the path, each end widened | 1e-5 relative (the E13.6 print of `m`, the integrator, the last step's overshoot past `m_Al = 0`: its metal clamped to 0, its cap kept) | 0 of 25 outside; the farthest 3.4e-6 above `z(0)` and 2.8e-6 below `z(m_b)` (70 µm) |
| T1 | each parent's burnout time (its exit record's `t`) inside `[t(m_b), d₀ⁿ/K_eff]`, each end widened | 1e-5 relative | 0 of 25 outside |
| R3 | `Σ wdot(A) = Σ ṁ(parents) − Σ ṁ(residues)`; `wdot(OX) = 0` | 1e-5 relative | 5.2e-7 |
| R3e | the energy books with the residues' outflow at the product's datum (4.285015e6 J per kg injected) | 1e-5 relative | 3.2e-7 |
| B | `β = (z ρ_Al/ρ_ox)^{1/3}` with Gallier's densities at 2570 K (1673.7 and 2820.3 kg/m³: `z` does not depend on the densities, `β` does) from each size's residue mass, against Gallier's 0.57 (20 µm), 0.55 (70 µm), 0.49 (200 µm), 0.40 (400 µm) | 0.03 absolute (their two digits, the digitization) | 0.581, 0.565, 0.507, 0.400: worst 0.017 |

`z` at burnout: 0.330888 (20 µm), 0.303117 (70 µm), 0.283405 (100 µm), 0.220175 (200 µm), 0.107929 (400 µm). The
400 µm size is the one B needs: with `g` read at the injected diameter the 200 µm residues give `β` = 0.460, inside
Gallier's 0.49 ± 0.03, and only 400 µm (0.256) fails B (R2 and T1 fail at every size but 20 µm). Each
parcel rides its own `(j, k)` column: the fields and the printout do not depend on the thread count.

## RED

| mutation | what fails |
|---|---|
| the constant law at `h(0)` (`xi-cap=0.069522` in place of the law) | R2 25 of 25 (`z` 0.131 at every size), B: `β` 0.427 at every size |
| `g` read at the injected diameter for the whole burn | R2 and T1 20 of 25 (the 20 µm parcels hold `g(20)` either way), B: `β` 0.256 at 400 µm |
| the law read at the drop's diameter in place of the metal's | R2 and T1 20 of 25 |
| the stretch left at 1 | T1 25 of 25: every burnout at 1.20–1.34 times `d₀ⁿ/K_eff` (R2 holds: the path holds no time) |

A solver that does not know the `xi-law` token stops on `IGLOO: unknown per-material key in the phase file`.
