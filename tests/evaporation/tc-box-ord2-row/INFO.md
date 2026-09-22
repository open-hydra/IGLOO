# INFO — evaporation/tc-box-ord2-row (ord2 source reduction conserves, e2e)

## Reference
No literature reference: this is a conservation identity end to end. The evaporation physics
is [`evaporation/tc-box`](../tc-box/INFO.md)'s and is gated there. See
[BUGS.md](../../BUGS.md) §E row **O28**; the unit half is
[`infrastructure/source_reduction`](../../infrastructure/source_reduction/INFO.md).

## Why it exists
Under `gas-order = 2` the source accumulators live on the **dual** mesh and
`obj_sourceblock%finalize` reduces them to geo cells. Before O28 that reduction volume-
**averaged** extensive rates, which hands each dual cell on with weight
`Σ_{c∋d} subVol/V_geo(c)` — 1 in the interior, but **½** on a face row, ¼ on an edge, ⅛ on
a corner. Every boundary dual row therefore leaked most of its deposit, and hydra-MI2
consumes exactly this field.

**Nothing in the suite could see it.** Every telescoping gate ran at `gas-order = 1`; the
only ord2 source gate, `standard/swirl-wedge-deposit`, never enters a boundary row except
the outlet and is blind by magnitude — measured, its band sums were **0.14 %** off before
the fix, well inside its own 2 % tolerance.

This case is also the **only** gate on the *wiring*: `test_source_reduction` calls
`precomputeDualWeights` itself, so it would stay green if the call vanished from
`allocate_blocks`.

## Fixture
`INPUT/{solfile.tec,phase.txt,properties.dat}` symlink `../../tc-box/INPUT/`;
`INPUT/bc.txt` is written by the committed `make_fixture.py` (`python3 make_fixture.py`,
1250 lines, every face `bcdef = 100`). The shared `tests/tools/make_box_case.py` cannot be
reused: its `write_bc` hardcodes `401` on face 1, and this case injects by coordinate
instead of through an inlet face.

`input.ini` is tc-box's with four changes:

| change | why |
|---|---|
| `gas-order = 1` → **`2`** | the dual mesh and `finalizeSRC`'s reduction are what is under test |
| **`mollify = off`** | the smoother is itself conservative, but it spreads the boundary row 60+ cells inward, which would hide *where* the mass went |
| `[IGLOO-BC] ds` → **DB injection** (`x/y/z/diam/mdot/up/vp/wp`) | places parcels by coordinate, so they can be put in the leaking row on purpose |
| `INPUT/bc.txt` | no `401` face, so nothing is injected by the boundary |

**Placement** (box 60×5×5 over `x∈[0,0.15]`, `y,z∈[0,0.05]`: `Δx = 2.5 mm`,
`Δy = Δz = 10 mm`; dual nodes sit at the geo cell centres 5, 15, 25, 35, 45 mm):

- `y = 2 mm` → dual row **`j = 1` = [0, 5] mm, a FACE row** — the row that leaked;
- `z = 12/22/32/42 mm` → dual `k = 2..5`, interior, so exactly one index is a boundary;
- `x = 5.5 mm` → dual `i = 3` and inside geo cell 3. **Not 5.0 mm**: that is a geo node
  plane, and `isPointInsideHexahedron12` uses an absolute `eps = 1e-15 m`, so a point on a
  shared face plane is inside *both* neighbours and which one wins is decided by the
  locator's sweep order, not by geometry.

`up = u_g = 10`, `v = w = 0`, no body force ⇒ zero slip ⇒ the whole path stays in row
`j = 1` and exits through face 2. `mdot = 6.18e-4 kg/s` per stream with
`ρ_p = 2950` (`properties.dat` col 3 — the INI's `[GPB-Phase1] rho = 8000` is **not read**)
and `d = 20 µm` ⇒ `m = 1.2357e-11 kg` ⇒ `npdot = 5.0e7 s⁻¹`, the same number flow as one
tc-box inlet cell.

⚠ **This is not tc-box's thermal regime.** There is no `temp0` key, so the injection
temperature is 0, which `Lib_Integration` replaces with the gas value: the drops start at
**600 K** and cool to the ~536 K wet bulb, where tc-box injects at `kt = 0.5` ⇒ 300 K and
heats. (Under ord2 that replacement value is dual node 1 — for row `j = 1` a ghost — which
equals 600 K here only because this box's gas field is spatially uniform. Do not restate it
as a general "starts at the local gas temperature" rule.)

## What it verifies
`source.tec` keeps its geo-shaped cell-centred layout under ord2 (`finalizeSRC`
`move_alloc`s the accumulators down before the write), `wdot(A)` first of five cell
variables.

| gate | assertion |
|---|---|
| W1/W2 | both INI changes were read — no ord1 fallback notice, and the mollification-OFF witness |
| P | all 4 parcels injected and exited |
| X | the last trajectory row reaches the outlet (`L − x_last ≤ Δx/2`) |
| T | `Σ wdot` over every cell equals `Σ_p npdot·(m_inject − m_exit)` within `1e-3` |

The oracle reads trajectory rows and the injected mass flow only — never the Eulerian field
it is checking. Nothing is dropped from the sum: the whole path is inside the domain.

Under ord2 the **exit row is printed** (an `IamOut` crossing zeroes `Ncell` before the print
test, and the position is snapped to the face intersection), so `x_last = L` exactly and the
closure term vanishes — measured gap `0.000e+00`. It is kept as a guard in case the last row
ever lands at the final dual face instead; there its error is second order on the plateau
(< 1e-5). This differs from ord1, where the last row sits at the final *interior* crossing —
which is why tc-box's own oracle drops the exit layer and this one must not.

## Red–green evidence (2026-09-22, release build, sprop2)

| | `bin/IGLOO.pre-A` | with O28 |
|---|---|---|
| `Σ wdot` | 4.336423e-04 kg/s | **8.699755e-04 kg/s** |
| evaporated (oracle) | 8.699722e-04 kg/s | 8.699722e-04 kg/s |
| **resid** | **5.015e-01** | **3.779e-06** |

`rc = 0` in both cases: half the source silently vanished.

## Pass criterion
All five PASS. Tolerance `1e-3` sits 265× above the measured floor (3.779e-06, identical at
`OMP_NUM_THREADS` 1, 2 and 5 — the four parcels ride separate `(j,k)` columns, so no two
threads accumulate into one cell) and 500× below the RED signature. Do not loosen past
~1e-2, where the gate stops separating half a deposit from a whole one.
