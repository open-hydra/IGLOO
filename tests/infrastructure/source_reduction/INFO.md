# INFO — infrastructure/source_reduction (conservative ord2 source reduction, unit)

## Reference
No literature reference: this is a conservation identity, not a model.

## What it catches
Under `gas-order = 2` the source accumulators live on the **dual** mesh (`1..N+1`) and
`obj_sourceblock%finalize` reduces them to geo cells. A volume-**average** reduction,

    S_geo(c) = Σ_oct subVol(oct,c)·S_dual(c+oct) / Σ_oct subVol(oct,c)

is wrong, because `sourceMass/Mom/En` are **extensive rates** (kg/s, N, W):
`computeSource` deposits a segment's whole contribution, it is not a density. Summing the
averaged field over geo cells hands each dual cell on with total weight
`Σ_{c∋d} subVol/V_geo(c)`, which is 1 only for an interior dual cell. On a uniform
Cartesian mesh a face row gets **½**, an edge **¼**, a corner **⅛**; on a wedge the axis
row's inner half-annulus is ¼ of the cell area, so the axis dual row hands on ¼; on a
general curvilinear cell it is off by the octant inequality.

So every boundary dual row leaks most of its deposit, and the loss scales with how much
of the source sits near a boundary — which for an injected spray is *all of it at the
inlet*. hydra-MI2 consumes exactly this field.

## The reduction, and why it is exact
A second array, `dualW(d)`, assembled once per block beside `subVol`:

    dualW(d) = Σ over every (cell, octant) pair touching d of subVol

    S_geo(c) = Σ_oct subVol(oct,c)/dualW(c+oct) · S_dual(c+oct)

Then

    Σ_c S_geo = Σ_d S_dual · [Σ_{(c,oct): c+oct=d} subVol(oct,c)] / dualW(d) = Σ_d S_dual

because the bracket **is** `dualW(d)` by construction. Exact on any mesh, no tolerance.

The clipped physical dual volume is deliberately **not** used: a conservation identity
enforced through an approximate geometry can mask a geometry defect, such as the dual ghost
reflected across the axis that `axis-200`'s E1/E2 catch.

## What it verifies

| case | mesh | deposit | RED (averaging reduction) `Σ_geo/Σ_dual` | GREEN |
|---|---|---|---|---|
| S1 | 3D 4×3×2 uniform | every dual cell | **0.400000** (= 24/60) | 1 |
| S2 | 3D 2×2×2 uniform | every dual cell | **0.296296** (= 8/27) | 1 |
| S3 | 2D 4×3 uniform | every dual cell | **0.600000** (= 12/20) | 1 |
| S4 | 2D 2×2 uniform | every dual cell | **0.444444** (= 4/9) | 1 |
| S5 | 3D 4×3×2 uniform | one dual cell: interior / face / edge / corner | 1 / **½** / **¼** / **⅛** | 1 each |
| S6 | 2D 4×3 uniform | one dual cell: interior / face / corner | 1 / **½** / **¼** | 1 each |
| S7 | 3D 4×3×2 **curvilinear** | all / interior / face | 0.400000 / **1.000052** / **0.499861** | 1 each |
| U2 | two 3D boxes joined by 101 along x (3 cells of h_A, 2 of h_B; h_B = h_A, then 1.5 h_A) | one interface dual cell of block A, source and euler | the class sum and the combine skipped: source A **1.0** / B **0.0**; euler total of A + B **0.5** / **0.4** | source A h_A/(h_A + h_B) = 0.5 / 0.4, B the rest; every total 1 (err ≤ 4.4e-16) |

S7 is the case that matters most: on a warped mesh the averaging reduction's error is the
octant inequality rather than a clean rational, and the `dualW` identity must **still** be exact.
It is — every assertion measures `0.00E+00`, so the gate needs no tolerance at all (the
declared `1e-12` is four decades of headroom over the worst double-precision residual seen
anywhere, 2.2e-16).

Each case also asserts the partition-of-unity invariant `Σ_d dualW == Σ_c Σ_oct subVol`,
and the mass, momentum and energy slots are checked separately so a reduction that is right
on only one of the three cannot pass.

⚠ **`Σ dualW` is NOT asserted against `Σ cellVol` on S7.** `cellVol` is `hexVolume` of the
parent hex, and `hexVolume` is a 5-tetrahedron split taking `abs()` per tet — **not additive
under octant subdivision** when the hex has non-planar faces. Measured on S7's node formula
the two differ by **5.1e-3** relative in total and **1.1e-2** worst per cell. An assertion
against `cellVol` there would be RED on the **correct** code. It is asserted on the
Cartesian cases S1–S6, where all faces are planar and the identity is exact.

## The 101 leg (U2)
Two boxes share a 101 face. Each block's dual cell at an interface vertex holds only its own
octants, so the setup pass (`assembleInterfaceDual`) joins the two dual cells of every interface
vertex into one class and writes the class's summed octant volume to both `dualW` entries. After
the cross-rank reduction `combineInterfaceDual` gives every member the class's source sum and,
for the euler field, the class's extensive total over that summed volume. `finalizeSRC` and
`reduceEUL` then hand each block its own octants' share: h_A/(h_A + h_B) to block A.

The euler deposit is written as the solver writes it, the extensive amount over the dual cell's
volume (`gas%cellVol`); the euler total of a block is Σ_c ρ_geo(c)·Σ_oct subVol(oct,c).

## Fixture notes
The mesh is built through the **production** geometry path (`compute_geometry` +
`precomputeMetric`), so the octant volumes under test are the shipped ones. Two things the
`dual_clip` template does **not** do and this test must:

- `use IGLOO_variables` and set **`ord2 = .true.` before `compute_geometry`** — `subVol` is
  allocated only under `ord2`, so a verbatim copy of the template leaves it unallocated and
  `precomputeDualWeights` crashes;
- set the module `mesh2D` per case (not the `is2D` dummy the template passes), plus
  `nm = 1` (`finalizeSRC` sizes `accMass(nm)`) and `mollifyPasses = 0` (the reduction alone
  is under test).

## Pass criterion
Every assertion exact — measured `0.00E+00` on all fourteen cases. U2 within 1e-14 relative
(measured ≤ 4.4e-16).
