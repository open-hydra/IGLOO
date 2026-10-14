# INFO — standard/psd-table, psd-normal, psd-lognormal, psd-rosinrammler (injection size distributions, e2e)

## Reference
The prescribed size distributions themselves (see [`psd_table`](../../infrastructure/psd_table/INFO.md) for
the four CDFs), tested on the parcels the solver injects with the Kolmogorov-Smirnov and the case-0
Anderson-Darling statistics [Kol33], [Ste74], [AD54], [MM04] (`tools/ks_ad.py`).

## What a sampled law means
Inside an inlet cell every parcel carries the same mass flow, ṁ_cell/N_inj, whatever its diameter: the
injected **mass-flow** distribution over d is the sampled law, and the number distribution is ∝ f(d)/d³.
A size-distribution file therefore lists mass fractions, sampled directly. These cases use the cell-centre
injection mode (`fsample = 1`, no `ds`, no `mdot-max`: one parcel per inlet cell, N_inj = 1) on a uniform
inlet, so every parcel carries the same ṁ and the unweighted distribution of the parcel diameters is the
mass-flow distribution; the check asserts both halves of that statement.

## Fixture
One box from `tools/make_box_case.py --nx 1 --ny 200 --nz 100 --kv 1.0`: x ∈ [0, 0.15] m in one cell, a
200 × 100 inlet on face 1 (cells 0.25 × 0.5 mm), uniform gas (ρ 1.2, U 10 m/s, T 600 K), parcels injected
at the gas velocity and temperature, so each crosses the box in a straight line. `psd-table/INPUT/` holds
the gas file; the other three link it, and link `psd-table/check.py`. 20 000 parcels per case, about 4 s
at 5 threads.

| case | bc.txt columns 6-8 (`--rp`, `--law`) | law |
|---|---|---|
| `psd-table` | `1.0e-4 0.0 INPUT/table6.dat` | `INPUT/table6.dat`: radii 2.5, 4, 6, 15, 22.5, 30 µm, mass fractions 0.05, 0.10, 0.15, 0.35, 0.25, 0.10 (bimodal), histogram reading |
| `psd-normal` | `1.0e-5 1.0e-5 Normal` | mean 20 µm, σ 10 µm, truncated to d > 0 |
| `psd-lognormal` | `1.0e-5 1.0e-5 LogNormal` | mean 20 µm, standard deviation 10 µm |
| `psd-rosinrammler` | `1.5e-5 2.5 RosinRammler` | x0 30 µm, n 2.5 |

In `psd-table` column 7 is 0, which must not turn the table into a Dirac (the `sigmap ≤ 0` rule applies to
the law names only), and column 6 (r_p = 100 µm) lies outside the table on purpose: the table sets the
size, and IGLOO prints one notice.

## What it verifies (`check.py`)
Each parcel's diameter and position come from its **first trajectory row**, written at injection (the
exit record carries no diameter), its mass flow from its exit record, joined by ID.

| gate | criterion |
|---|---|
| W | the cell-centre mode ran (`Placing particles at the bounday every  1 cells`, `number of particles =  20000`) |
| N | 20 000 parcels, IDs 1..20000 once each in both files; first rows at x = 0, one per inlet cell centre |
| M | every parcel's mass flow equal (relative spread ≤ 1e-12), and the mass-weighted D_n equal to the unweighted one |
| L | `INPUT/bc.txt` columns 6-8 (and `table6.dat`) equal to the constants typed in `check.py`: the oracle and the fixture describe the same law |
| KS | D_n < 1.94947/√n = 1.3785e-2 (α = 10⁻³) |
| AD | A² < 5.969 (α = 10⁻³) |
| W | `psd-table` only: exactly one notice that column 6 lies outside the table |

The seed is fixed (`[IGLOO-General] seed = 42`), pinning is serial and replicated, so the diameters do not
depend on the thread or rank count (the printout is the same at 1 and 5 threads); the intrinsic random
stream does depend on the compiler, which is why both statistics use α = 10⁻³.

## Measured (ifx 2023.1)

| case | D_n | p | A² | p |
|---|---|---|---|---|
| `psd-table` | 5.135e-3 | 0.667 | 0.457 | 0.791 |
| `psd-normal` | 7.510e-3 | 0.209 | 1.008 | 0.353 |
| `psd-lognormal` | 7.395e-3 | 0.224 | 0.917 | 0.404 |
| `psd-rosinrammler` | 5.134e-3 | 0.668 | 0.457 | 0.791 |

Mass flow per parcel 7.72727e-7 kg/s in every case, spread 0. Smallest Normal diameter 1.3 nm.

Sensitivity — each a one-line change of the reader or the sampler:

| change | case | D_n | A² |
|---|---|---|---|
| the table read with a list-directed read, which ends a path at its slash | `psd-table` | refused: `unsupported distribution law: INPUT`, no output | — |
| table: every draw at its class centre | `psd-table` | 0.1753 | 542.6 |
| table: number fractions m_k/r_k³ taken as the weights | `psd-table` | 0.7033 | 2.95e4 |
| Normal: rejection removed, a draw below 1 µm clipped to 1 µm (a lower-tail defect) | `psd-normal` | 0.02334 | 40.9 |
| LogNormal: σ_p taken as the log-space standard deviation | `psd-lognormal` | 0.5933 | 8435 |
| Rosin-Rammler: n 10 % high | `psd-rosinrammler` | 0.02711 | 41.3 |
| Rosin-Rammler: x0 and n exchanged | `psd-rosinrammler` | 0.6315 (diameters 0 and ∞, 276 parcels lost, N fails) | 3.4e5 |
