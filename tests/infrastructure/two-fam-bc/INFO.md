# two-fam-bc — a particle boundary file with one copy of the inlet table per family

**Purpose.** ATLAS BCB writes `<name>-bc.txt` with one copy of every mesh block's boundary table
per (material, population) pair of the phase — block-outermost, material-major, population-minor
(`write_dp_bc` in ATLAS's `src/BCB/io_output.f90`). IGLOO numbers its families the same way
(`famID` counts the groups of `phase.txt` material by material), so copy `c` of a block is family
`c`. This case is the end-to-end gate on that contract: two families whose inlet payloads differ
must each inject with their own.

**Reader contract** (`read_cdp_bc_file`, `src/lib/IO.f90`). A first pass counts the records:

| records in the file | layout | families get |
|---|---|---|
| one per boundary face of every block | one copy | every family the same payload (every single-population phase, every box fixture) |
| `totFam` × that, `totFam > 1` | one copy per family, ATLAS order | copy `c` → family `c`; copies 2.. must repeat copy 1's header integers and codes face by face |
| anything else | — | refused at setup (`refuse-bc-copies`); a copy out of order is refused too (`refuse-bc-copy-order`) |

`krhoTot` is then the sum of the families' own `krho` and `mdotPart` of their own `gp · area`; the
per-cell spacing (column 9) and the diameter law are per family. Only copy 1's headers tag the face
cells; for 101/201 cells the connection line of the last copy is the one kept (ATLAS writes the same
line in every copy).

**Fixture.** `two-mat`'s case with its own `bc.txt`: the `drag-stokes` box (60 × 5 × 5 cells,
0.15 × 0.05 × 0.05 m, uniform gas ρ 1.2, U 10, T 600, μ 1.8e-5; face 1 the inlet, 25 cells of
0.01 × 0.01 m, every other face 100), `phase.txt` `A 1` / `B 1` and the two-zone `properties.dat`
(ρ_A 2950, ρ_B 1000) symlinked from `two-mat`, `drag = Stokes`, `[IGLOO-BC] ds = 10.0` (one
centred parcel per inlet cell). `INPUT/bc.txt` is written by `make_fixture.py` (committed; run once
from this directory, `python3 make_fixture.py`): 2 × 1250 records, family 1 (A) = `drag-stokes`'s
inlet line, family 2 (B) half the loading and twice the radius:

```
A:   3.40000E-01   1.00000E-01   normal,   normal,   1.00000E+00   5.94500E-06   0.00000E+00   Dirac   0.00000E+00
B:   1.70000E-01   1.00000E-01   normal,   normal,   1.00000E+00   1.18900E-05   0.00000E+00   Dirac   0.00000E+00
```

The headers use ATLAS's `b i j k f code` (6I8, `(i, j, k)` from `obj_block%fmn2ijk`'s face rule), so
the file is a faithful layout witness; the box generators (`tools/make_box_case.py`,
`wall-approach/make_fixture.py`) write `b f m n 1 code`, which the reader also accepts — it never
interprets columns 1–5 against the mesh, only compares the copies with each other.

Two materials rather than two populations of one material because each family then has its own
output files and its own closed-form oracle; the reader sees only the family count
(`test_bc_families` F2 runs both shapes). ATLAS reads one key set for every material of a phase, so
two materials with different inlet payloads cannot come from it today; two populations of one
material can (`<name>-krho = a b`, one entry per population).

**Numbers** (all from the inputs; ṁ_gas per inlet cell = ρ U A = 1.2 · 10 · 1e-4 = 1.2e-3 kg/s,
exact on the uniform field):

| quantity | read per family | copy 1 fed to both families |
|---|---|---|
| `krhoTot` | 0.34 + 0.17 = 0.51 | 0.68 |
| ṁ_A = krho_A/(1 − krhoTot) · ṁ_gas | 8.326531e-4 kg/s | 1.275e-3 kg/s |
| ṁ_B | 4.163265e-4 kg/s | 1.275e-3 kg/s |
| d_B | 23.78 µm | 11.89 µm |
| τ_B = ρ_B d_B²/(18 μ) | 1.745e-3 s | 4.36e-4 s |
| m_B/m_A at injection | (1000/2950) · 8 = 2.711864 | 0.338983 |
| ρ_p/n_p in `euler2.tec` = ρ_B π d_B³/6 | 7.0410e-12 | 8.8012e-13 |

**Gates** (`check.py`; the Stokes oracle is imported from `standard/drag-stokes/check.py`, the
parsers from `two-mat/check.py`):

| id | kind | assertion |
|---|---|---|
| T1 | target | B's parcels pass the Stokes closed form at d = 23.78 µm, τ = 1.745e-3 s (the oracle checks every row's diameter, so the import is the diameter gate) |
| T2 | target | every `outloc-A.dat` record carries ṁ = 8.326531e-4 and every `outloc-B.dat` record 4.163265e-4 (1e-5 relative) |
| T3 | target | m_B/m_A = 2.711864 at the injection row of every parcel (1e-5 relative) |
| T4 | target | `euler2.tec`: ρ_p/n_p = 7.0410e-12 in every deposited cell (1e-9 relative); `euler1.tec`: A's 2.596e-12 |
| C1 | control | A's parcels pass the Stokes closed form at d = 11.89 µm; A and B inject from the same 25 stations with the same (x, y, z, U, V, W, T) |
| C2 | control | 25 exits per material, non-empty `scatter-<mat>.dat` |
| C3 | control | `source.tec` carries `wdot(A)`, `wdot(B)`, both zero; the shared Fx and E slots balance both materials' parcels with each parcel's own ṁ to 1e-6 |

**RED / GREEN.** Run on a reader that feeds copy 1 to every family (the frozen binary before the
per-family read), T1–T4 fail with exactly the right-hand column: T1 `dp=1.1890e-05 != 2.3780e-05`
on all 25 parcels of B (0 verifiable relaxers), T2 ṁ = 1.275000e-03 on all 25 exits of both
materials, T3 m_B/m_A = 0.338983 on 25 of 25 parcels, T4 ρ_p/n_p = 8.8012e-13 in 1500 of 1500
cells of `euler2.tec`; C1–C3 pass (Σ Fx −0.573749, Σ E −3.155615 — two-mat's values, since that run
is two-mat's). With the per-family read every gate passes with the left-hand column: τ_B =
1.745335e-3 s (25/25 verified, worst resid/tol 0.021), ṁ 8.326531e-4 / 4.163265e-4, m_B/m_A
2.711864, ρ_p/n_p 7.0410e-12 (B) and 2.5964e-12 (A) in 1500 cells each, Σ Fx −0.2810128,
Σ E −1.545536.

**Not covered here.** One mesh block: no gas file in the suite has two zones, so the two-block
order (block 1's copies, then block 2's) is gated in-process by `test_bc_families` F3 only.
