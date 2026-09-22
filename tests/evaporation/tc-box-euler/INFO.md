# INFO — evaporation/tc-box-euler (Eulerian field of a mass-losing material, e2e)

## Reference
No literature reference and no oracle curve: the evaporation physics is
[`evaporation/tc-box`](../tc-box/INFO.md)'s and is gated there. This case gates the
**equivalent-Eulerian deposit for ODE models 2 and 5** — see [BUGS.md](../../BUGS.md) §E
row **O31**.

## Why it exists
`Lib_Integration::computeEulField` branches on the ODE model. For models 2 and 5 it
deposited

    rho = particle%intE(5) * factor

where `intE(5)` is `∫ m dt` for **one** droplet (`F(14) = m` in the model-2 and model-5
RHS): units kg·s/m³, the number rate missing. The sibling branches are all consistent —
model 4's integrand is already `m·npdot` (`F(15) = mTraj`), and models 1 and 3 deposit
`mdot·Tstay` with `mdot = npdot·m`. Even the body-force source correction multiplies the
same slot by `npdot`.

Because `finalizeEUL` normalizes velocity and temperature **by** `rho_p`, the missing factor
does not merely scale the density: `u_p` and `T_p` come out multiplied by `npdot` (5.0e7
here), and `rho_p/n_p` — the mean droplet mass — reads `1/npdot` of the truth.

**It is a regression, not a design gap.** The deposit read
`rho = particle%npdot*particle%intE(5)*factor` at `b2ed810` (2025-11-20) and lost the factor
in `047b11e` (2026-03-24, "refactoring - phase 1") while the integrand stayed `F(14) = m`.
Its parent `842c180` still carries the correct line. Those commits are reachable only on the
local `archive/dima-*` branches — the public history is squashed at `c7195f1`, so `git blame`
of the current line points at the squash. The fix restores the Nov-2025 line verbatim.

**Nothing in the suite could see it.** Every Eulerian writer here is model 1 (no
evaporation/combustion key), and every model-2/5 fixture writes `out-file = S`. The
intersection was empty.

## Fixture
`INPUT/` is four per-file symlinks to `../../tc-box/INPUT/` (the `datum-abs` precedent);
`input.ini` is tc-box's with **two** changes, both load-bearing:

| change | why |
|---|---|
| `out-file = S` → **`ALL`** | turns `eulerSwitch` on, so `OUTPUT/euler1.tec` is written at all |
| **`mollify = off`** | the binomial smoother mixes `rho_p` and `n_p` across cells with the same kernel; the per-cell *ratio* is only exact unsmoothed, and that ratio is the oracle |

`gas-order = 1` and `print-dcell = 1` are kept. Model 2 comes from tc-box's
`evaporation = TC`. 25 parcels inject at the inlet cell centres with `kV = 1` ⇒ zero slip ⇒
`u ≡ 10 m/s`, so each rides its own `(j,k)` column of 60 cells and never changes column.

Constants the oracle derives itself, never from production output:
`ρ_p = 2950 kg/m³` (`INPUT/properties.dat` column 3 — ⚠ the INI's `[GPB-Phase1] rho = 8000`
is **not read by IGLOO**; do not "correct" the oracle to it), `d₀ = 20 µm` ⇒
`m₀ = 1.235693e-11 kg`; `ṁ_inj = 0.34/0.66·1.2·10·1e-4 = 6.181818e-4 kg/s` ⇒
`npdot = 5.002713e7 s⁻¹`; box 60×5×5 over `x∈[0,0.15]`, `y,z∈[0,0.05]`.

## What it verifies
`euler1.tec` is BLOCK-packed: 3·nn nodal values then 6·nc cell values in the order
`rho_p u_p v_p w_p T_p n_p`, `i` fastest. A parcel's column is its first trajectory row's
`(y,z)`; for cell `i` the bounding rows are at `x = (i−1)Δx` and `x = iΔx`. `print-dcell = 1`
writes one row per cell **entry**, so the outlet layer `i = 60` has no `m_out` row and is
skipped — **1475** of 1500 cells are checked and the count is asserted.

| gate | assertion | RED (pre-fix) | GREEN |
|---|---|---|---|
| W1/W2 | both INI changes were read — the two log witnesses | present | present |
| N | ≥ 1000 cells checked (non-vacuity) | 1475 | 1475 |
| E3 | `abs(rho_p/n_p − (m_in+m_out)/2) / m̄ ≤ 1e-4` | **1.000** | 2.930e-05 |
| E4 | `abs(u_p − 10)/10 ≤ 1e-11` | **5.003e+07** | 7.19e-14 |
| E4 | `abs(v_p), abs(w_p) ≤ 1e-18` | 2.17e-14 | 4.34e-22 |
| E5 | `T_p` inside `[T_in, T_out] ± 0.05 K` | **2.591e+10 K** | **0.0 exactly** |
| E1 | `Σ rho_p·V / Σ npdot·∫m dt ∈ [0.5, 2]` (coarse conservation) | **0.0000** | 1.0148 |

W1/W2 are **positive** witnesses. Without them a run that never read the INI would pass the
absences, and `out-file`/`mollify` are exactly the two keys this fixture exists to set.

E3 is the sharp gate: `rho_p/n_p` is an *intensive* quantity, so the defect's `1/npdot`
shows as a deviation of ~1 whatever the cell — four decades above the bound. E1 is the
independent extensive half, deliberately coarse (its 1.5 % offset is the trapezoid rule plus
the uncounted exit segment, not a physics residual).

## Pass criterion and tolerances
All eight PASS. Every tolerance comes from measurement on this fixture and every oracle
quantity is **bit-stable at `OMP_NUM_THREADS` 1, 2 and 5** — each parcel owns its cell
column, so no two threads accumulate into one cell and the `!$OMP ATOMIC` never
re-associates. Margins: E3 3.4×, E4 `u` 133×, E4 transverse 2300×, E5 exact. Do not loosen
E3 past ~1e-3; the residual it measures is the second-order error of the two-row mean
against the ODE's own time integral, and it is what makes the gate an oracle rather than a
self-comparison.
