# two-mat — the first case with more than one material (nm = 2)

**Purpose.** Every `tests/**/INPUT/phase.txt` carried exactly one material line, so the
per-material code had only ever executed at arity 1 — in every serial, repeatability and MPI
gate (`plan-bucket/mpi-residuals.md` R2, coverage audit). This case runs `standard/drag-stokes`
for two materials that differ only by density, so each material keeps a closed-form oracle.
Paths that execute at arity 2 here for the first time: `obj_IGLOO%solve`'s material loop
(shard opens, dNscat pre-pass, per-group child-ID bands), the material index of `sourceMass`
and its finalize/mollify, `nfam = 2` (`euler(:, famID)`, `euler1.tec`/`euler2.tec`), the
per-material output files, `setupRHS`'s sequential per-material setup, `reset_state`'s per-material restore, and — under MPI — six
shard merges per sweep and the per-family euler reduction.

**Fixture.** `drag-stokes`'s box, gas, `bc.txt` and `solfile.tec` (symlinked). `phase.txt`
lists `A 1` and `B 1`; `properties.dat` has two zones, A = `tests/common/properties.dat` (ρ =
2950, cp = 1250) and B the same table at ρ = 1000. `bc.txt` carries ONE property line per inlet
cell and the reader copies it to every family, so both materials inject the same 25 parcels
(same d, v₀, krho) and differ only through their density — hence their Stokes time
τ = ρ_p d²/(18μ). The fan-out also sums krho into `krhoTot` for every family (0.68), so each
parcel of either material carries ṁ = krho/(1 − krhoTot)·ṁ_gas = 0.34/0.32 · 1.2e-3 = 1.275e-3 kg/s
(drag-stokes: 0.618e-3) — a 68 % condensed mixture; internally consistent, oracle-blind (Stokes
relaxation is ṁ-independent), pinned by M7. `out-file` is omitted on purpose: the default enables both the source and
the eulerian output (with a "not recognized or not given" warning — ledger O21), so
`source.tec` carries two `wdot` slots and one euler file per family is written.

**Gates** (`check.py`; the Stokes oracle is *imported* from `standard/drag-stokes/check.py`
and re-parametrised per material — one oracle, not two copies):

| id | assertion |
|---|---|
| M1 | the drag-stokes oracle passes on `trajectories-A.dat` at ρ = 2950 and on `trajectories-B.dat` at ρ = 1000 |
| M2 | identical injection rows (x, y, z, U, V, W, T, d) for every ID in A and B; parcel mass ratio m_B/m_A = 1000/2950 |
| M3 | at the first interior row U_B > U_A for every parcel (the lighter material has relaxed further); the files differ |
| M4 | `source.tec` variables include `wdot(A)` and `wdot(B)`, both finite and identically zero (by construction — no phase change — so this pins only the slot names); the SHARED slots balance both materials' parcels: Σ Fx = −Σ_p ṁ_p (u_exit − u₀) and Σ E = −Σ_p ṁ_p (u²_exit − u₀²)/2 to 1e-6 (measured −0.573749 and −3.155615; an nm-fold double count or a dropped material in `sourceMom` misses by O(1)) |
| M5 | `euler1.tec` is material A and `euler2.tec` is material B: ρ_p/n_p = ρ_mat π d³/6 in every deposited cell (2.596e-12 / 8.801e-13, uniform to 7e-15) — a family↔material swap or a cross-deposit changes the ratio |
| M6 | 25 exits in `outloc-<mat>.dat` and a non-empty `scatter-<mat>.dat` for each material |
| M7 | every exit record of either material carries ṁ = 1.275e-3 kg/s — the krho fan-out and `krhoTot` at arity 2 |

**Falsification (2026-09-16, throwaways in the session scratchpad).** (1) `phase.txt` lines
swapped (`B 1` / `A 1`): material "A" then integrates on zone 2's ρ = 1000 and M1 fails on every
parcel (|x − x_pred| 4.9e-3..1.5e-2 vs tol 1e-6) — `properties.dat` zones are bound to `phase.txt`
lines by **order**; the zone title is not read (ORION's points reader keeps no name), so a
reordered file is a silent density swap that only a per-material oracle can see (ledger O24).
The R2 plan expected the per-material files to swap and both oracles to pass — that presumed
name-keyed properties, which the code does not have; the gate sees the swap as a failure instead.
(2) One zone for two materials: used to SIGSEGV at `block(2)` with no message; `read_cdp_properties`
now error-stops with the counts (this commit), and also refuses a zone whose row count or start
temperature differs from zone 1's (zone 1 sets `Tmin/Tmax` for every zone — a shorter zone 2 was
an unchecked out-of-bounds read; falsified with a 3000-row zone B). The plan's other falsification
— delete a `[GPB-Phase2]` property — is inapplicable: IGLOO reads neither ρ nor cp from the INI. The four registrations are all green on the same
fixture: `two-mat` (serial), `repeat-two-mat` (sweep 1 ≡ sweep 0 on six `.dat` and three
`.tec`), `mpi-two-mat` (n = 4 × 2 threads, rank witness) and `mpi-consistency-two-mat`
(n = 1 vs 4: six `.dat` multisets and three `.tec` fields identical).

**Measured.** Serial: 25 parcels per material, 1500 trajectory rows each, 296 (A) / 386 (B)
scatter rows; A's rows differ from `drag-stokes`'s in U on 94 rows — 92 by 1 ULP of F12.6 and 2 by
2 ULP (euler on → neq = 12 vs 7 changes the step sequence) — inside the oracle's theoretical
tolerance. Byte-inert on the 62 pre-existing entries (A/B capture, floor-only).

**What it does not pin.** Both materials share one injection line, so a per-material
injection contract (different d or krho per material) cannot be expressed in `bc.txt` today —
hydra's MI2 memo records the single ATLAS DP slot as structural. The `[IGLOO-Properties]`
vector-size checks (`size == nm`) are exercised only for the vectors present, i.e. none here
(no evaporation/breakup). The famID-keyed RNG streams are *seeded* for both materials but never
*consumed* here: no breakup and a Dirac diameter (bc.txt col 7 = 0) mean zero draws, so
cross-material stream distinctness (R2's fifth bullet) is not observable on this fixture —
`test_rng_stream` is the unit-level pin.
