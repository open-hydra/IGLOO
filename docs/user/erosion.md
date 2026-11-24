# Wall Erosion

Add an `[IGLOO-Erosion]` section and every parcel that hits a wall erodes it. IGLOO writes the erosion of each wall
face to `OUTPUT/wall.tec`. Without the section nothing erodes and no file is written. The laws:
[Wall Erosion](../theory/erosion.md).

Which faces are walls: the dispersed-phase wall (code 301), a chimera face toward a block of another phase (104),
and a connection toward another phase with no partner cell (103). See
[Boundary & Injection](boundary-conditions.md).

---

## The keys

Each constant takes **one value per material**, in the order of the phase file. A material uses the values of its
own model; give any number for a key its model does not use. The keys are case-sensitive: `K` is Oka's, `k`
Neilson–Gilchrist's.

| Key | Model | Meaning | Required |
| :--- | :--- | :--- | :--- |
| `model` | | one word per material: `Oka`, `Neilson-Gilchrist` or `none` (any case). `none` erodes nothing, but its flux still counts | yes |
| `wall-density` | | $\rho_w$ [kg/m³], one value for every wall | yes, > 0 |
| `Hv`, `a`, `b` | Oka | the wall's Vickers hardness [GPa] and load relaxation; one value each | with an Oka material; $H_v$ > 0, $a$ > 0, $b \ge 0$ |
| `K` | Oka | the particle factor [mm³/kg] | > 0 |
| `k1`, `k3` | Oka | the exponents of hardness and diameter | |
| `v-ref`, `D-ref` | Oka | the reference speed [m/s] and diameter [m] | > 0 |
| `k2`, `n1`, `n2` | Oka | the speed and angle exponents: each as a value, or as `<key>-s` and `<key>-q` with `<key>` = s $H_v^q$ | one form each; $k_2$, $n_1$ > 0 |
| `phi`, `epsilon` | Neilson–Gilchrist | the energy of cutting and of deformation per unit eroded mass [J/kg] | > 0 |
| `n` | Neilson–Gilchrist | the angle constant | > 0 |
| `k` | Neilson–Gilchrist | the normal speed below which no deformation erosion occurs [m/s] | ≥ 0 |

## Example

Material A erodes under Oka's law with SiO₂'s constants; material B under Neilson–Gilchrist's.

```ini
[IGLOO-Erosion]
model        = Oka Neilson-Gilchrist     ; material A, material B
wall-density = 1800.0
Hv           = 1.43                      ; the wall, for Oka
a            = 0.05
b            = 1.2
K            = 65.0 0                    ; A: Oka's K [mm^3/kg]
k1           = -0.12 0
k2-s         = 2.3 0                     ; k2 = 2.3 Hv^0.038
k2-q         = 0.038 0
k3           = 0.19 0
n1-s         = 0.71 0                    ; n1 = 0.71 Hv^0.14
n1-q         = 0.14 0
n2-s         = 2.4 0                     ; n2 = 2.4 Hv^-0.94
n2-q         = -0.94 0
v-ref        = 104.0 0
D-ref        = 326.0e-6 0
phi          = 0 4841625.9               ; B
epsilon      = 0 11354975.3
n            = 0 5.0
k            = 0 20.0                    ; B: the threshold speed [m/s]
```

Setup stops the run on an unknown model word, a key whose number of values is not the number of materials, a
missing constant of a model in use, `k2`, `n1` or `n2` given in neither or both forms, a value out of range, a word
where a number is expected, and a wall constant with several values. The log names each material's model and its
constants.

## The wall file

`OUTPUT/wall.tec` has one zone per block face holding a wall cell, the values at the cell centres. For each material
`<mat>`:

| Variable | Meaning | Unit |
| :--- | :--- | :---: |
| `er-<mat>` | the erosion ratio: eroded mass per unit mass of impacting particles, weighted by their flow (0 where none hit) | kg/kg |
| `g-<mat>` | the incident mass flux | kg/(m² s) |
| `rdot` | the recession rate of the wall, all materials together | m/s |
