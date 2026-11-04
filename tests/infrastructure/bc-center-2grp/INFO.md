# INFO — infrastructure/bc-center-2grp (multi-group `bc_center` pinning, e2e)

## Reference
Not a physics case: there is **no** literature reference and no oracle curve. The physics
is `standard/drag-stokes`' (uniform-gas box, Stokes drag) and is gated there. This case
gates **particle pinning** only.

## What it pins
`pin_particles` dispatches on

    method=='FB' .and. ds==0 .and. mdotMax==0 .and. .not.dsSwitch

so `pin_particles_bc_center` runs only on a case that sets neither `ds` nor `mdotMax`
(`ds` ⇒ `pin_particles_bc_ds`; `[IGLOO-BC] x/y/z` ⇒ assigned/DB), and `setup` calls
`pin_particles` once per (material, group). The routine's particle index `p` and cell
counter `nn` must therefore start from 0 on every call: Fortran initializes a local *with an
initializer* once at program start, not per call (implicit SAVE), so the routine assigns
both at run time. This case drives it with two groups.

## Fixture
`INPUT/` symlinks `standard/drag-stokes/INPUT/{bc.txt,solfile.tec}` and
`common/properties.dat`; only `phase.txt` is a real file. `input.ini` is drag-stokes'
with two deliberate changes, **both load-bearing**:

| change | half of the carried-counter defect it triggers |
|---|---|
| `ds = 10.0` **deleted** | selects `bc_center` over `bc_ds` — nothing runs without it |
| `phase.txt` `A 1` → **`A 2`** | ≥2 groups ⇒ the particle index `p` carries from group 1 and group 2's writes land past the end of its own freshly allocated array |
| **`fsample = 2`** | `fsample > 1` ⇒ the cell counter `nn` carries into the **counting** pass, `mod(nn,fsample)` shifts phase, and the array is sized from a different cell subset than the fill pass walks |

⚠ At `fsample = 1` the second half is **invisible** (`mod(nn,1) == 0` always) and both
groups return 25 particles either way — a population gate would be vacuous. Do not
"simplify" this back to the default.

5×5 = 25 inlet cells on face 1, every 2nd taken ⇒ **12 particles per group, 24 total**.

## What it verifies
`check.py` builds every expectation from known inputs (mesh size, `fsample`, group count),
never from production output. Four gates:

1. **population** — both groups hold exactly 12 particles.
2. **identity** — each group's ID set is exactly `{1..12}`, no duplicates.
3. **clone** — group 2's exit records equal group 1's particle for particle (same cells,
   same properties, same physics ⇒ group 2 must be a faithful clone). The strongest
   available statement that group 2 was *pinned* correctly, not merely *counted* correctly.
4. **integration** — every particle writes ≥5 trajectory records and exits on the outflow
   plane. A pinning bug that injects everything out of domain still exits 0, so "the
   solver ran" is never the criterion.

## RED signature (counters carried between calls)
Measured on a release build whose `p` and `nn` carry between calls: `rc = 0`, **stderr
empty** — the defect announces itself in no way at all.

| gate | result with carried counters |
|---|---|
| population | group 1 = 12, **group 2 = 13** |
| identity | group 1 = `{1..12}`; **group 2 = twelve `0`s and a stray `13`**: of its 13 entries (the population row) only index 13 is written — `p` resumes at 12, so writes land at indices 13…24 of a `1:13` array (one in bounds, **eleven past the end**) — and indices 1…12, group 2's real IDs, are never assigned and keep their default `0` |
| clone | 13 rows vs 12 |
| integration | **1 of 13** particles reach the outflow plane — group 2's `pInj` is overwritten out of bounds too |

## Pass criterion
All four gates PASS. No tolerance tuning is involved: counts and IDs are exact integers,
and the clone comparison is bit-identical in practice (`max |Δ| = 0.0`), with a `1e-9`
slack only to absorb a formatting ULP.
