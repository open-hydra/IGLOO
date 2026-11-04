# INFO — combustion/burn-box-residue (the burnout residue: a parcel of the combustion product, e2e)

## Reference
No literature value: the residue's birth, its books and its freezing follow from the model's definitions and the
inputs. The cap's growth and burn law are `burn-box-cap`'s (`[Beck05]`, `[Gal21]`); the freezing is
`solid-box`'s (ODE model 6).

## Model
A capped parcel (ODE model 8) with `residue=on` ends at burnout (`m_Al ≤ 1e-15 kg`), at a solver failure or at a
non-finite state with its cap `m_ox > 0`; the cap becomes **one parcel of the combustion product** (`product=OX`, a
0-group line of the phase file with its own table zone and solidification data), born at the parent's position,
velocity, temperature and time — the state its exit record prints — with the cap's mass (`d-residue = 0`: one
residue per droplet) or at `d-residue` with the mass flow kept. The gas books lose the residue's birth flux at that
point (`ṁ_res (c_res T + hOff_res + ½|v|²)` in energy). A freezing product's residue is born liquid at or above
`T-melt`, else undercooled, with `f = 0`; its first segment applies the nucleation jump. Residue IDs are the
group's offset (the parcels of the earlier groups naming the same product) plus the parent's ID.

## Fixture
`burn-box`'s case (`bc.txt`, `solfile.tec` symlinked) and `burn-box-cap`'s absolute two-zone table (symlinked);
`out-file = S`, `out-time = on`, `mollify = off`. Its own `INPUT/phase.txt`:

    A 1 combustion=Beckstead K-burn=2.0e-6 n-burn=1.8 X-eff=0.5 T-ign=550 beta-part=0.3 q-comb=1000000 xi-cap=0.2 product=OX residue=on
    OX 0 solidification=on T-melt=2327 h-fus=1.08949e6 cp-solid=1360

(`T-nuc` 0.8 T-melt = 1861.6 K.) The 25 parcels burn out near x = 0.072 m at 792.81 K; each residue
(`d = 1.9678e-5 m`) is born undercooled, freezes whole (the jump to 998.8 K) and relaxes to the 600 K gas with
`τ = ρ_res c_s d²/(12 k_g) = 6.667e-3 s`. `check.py` runs two siblings with the same executable (`IGLOO_EXE`):

- `ref-dres/`: `d-residue=5.0e-6` and the burning line in two groups (`A 2 …`, `bc.txt`'s face table written twice:
  one copy per family);
- `ref-fail/`: `[IGLOO-ODE] fail-after = 3.6e-3` — the solver failure of the cap model forced once a parcel's last
  accepted time reaches it; the parcel reverts to its last accepted step, and its residue is born there.

## What it verifies (`check.py`)

| id | gate | tolerance | measured |
|---|---|---|---|
| R1 | one residue per parcel burned out inside the box: the product's IDs (1–25) = the parents' exit records with `X < Lx − 1e-6`; the product's banner counts 25 | exact | 25 |
| R2 | each residue's birth row = its parent's exit record (X, Y, Z, T printed alike; the speed within its print); `m = ξν(m₀ − 1e-15 kg)`; `d = (6m/(π ρ_res))^{1/3}` of the printed `m` | 3e-5; 5e-6 relative | 5.2e-6; 8.4e-7 |
| R3 | `Σ wdot(A) = Σ ṁ(parents) − Σ ṁ(residues)` (exit records); `wdot(OX) = 0` | 1e-5 relative | 8.5e-7 |
| R3e | `Σ E = Σ ṁ_p (h_A(600 K) + u²/2) − Σ ṁ_r (h_OX(T_exit) + u²/2)`, `h_OX` on the product's datum by phase (`f` of the exit record): 6.467406e6 J per kg injected | 1e-5 relative | 5.8e-7 |
| R4 | each residue's rows: born undercooled (`f = 0`), solid after its first segment at the enthalpy-conserving jump's temperature, then `T = T_g + (T_a − T_g) e^{−(t−t_b)/τ}`; per row the F12.6 budget of `temp-relax` with the jump's anchor (`c_l/c_s`) and the print of `d` in `τ` | the budget | worst 0.096 of it (800 rows) |
| R5 | each residue's `ṁ = (m_res/m₀)` × its parent's | 1e-5 relative | 1.5e-6 |
| D | `ref-dres/`: 50 residues with IDs 1–50 (the second group's at 25 + its parents' IDs); `d = 5e-6 m`; `ṅ_res/ṅ_parent = 60.958139`; R3, R3e (6.531021e6 J/kg); `T = 600 K` at the outlet | exact; 4e-5; 1e-5; 1e-4 K | 3.7e-6; 3.0e-7, 4.7e-7; 5e-6 K |
| F | `ref-fail/`: every parent ends at `t* ≥ 3.6e-3 s` (3.7278e-3 to 3.7491e-3 s, its last accepted step); its exit record is its residue's birth row; `m = ξν(m₀ − m_Al(t*))`, `m_Al(t*)` from the burn law (Gallier's `F` by bisection, `c = 1.397686`) at the printed `t*`; R3, R3e | 1e-5 relative | 3.8e-6 |

Each parcel and its residue ride one `(j, k)` column: the fields and the printout do not depend on the thread count.

## RED

| mutation | what fails |
|---|---|
| the residue pushed at a clean burnout only | F: no residue after the solver failure (banner 0) |
| the residue sized with the burning material's density (2950) | R2: d 2.1689e-5 m against 1.9678e-5 m (9.3 % of the printed d); R4 |
| the birth flux left in the gas books | R3: 0.61 of the expected flow (main, D), 0.44 (F); R3e |
| the product's datum left out of the birth flux | R3e: 100.7 % (main), 99.7 % (D), 98.8 % (F) |
| the residue born solid at its parent's temperature | R4: exit 659.854 K against 723.779 K; R3e 1.6 % |
| `d-residue` ignored | D: d 1.97e-5 m, the number-rate ratio 0.98 off |
| the group's offset dropped from the residue ID | D: 25 distinct IDs among the 50 residues |
| no revert on the failure exit (the segment-start state) | F: `t* = 3.5e-3 s < fail-after` on every parent |

A solver that does not know the `residue` token stops on `IGLOO: unknown per-material key in the phase file`.
