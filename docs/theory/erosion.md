# Wall Erosion

With an `[IGLOO-Erosion]` section in `input.ini`, every parcel that leaves the domain through a **wall cell**
erodes it: IGLOO evaluates the erosion law of the parcel's material at the parcel's impact speed, angle and
diameter, sums per wall face cell and per material the parcels' flow and eroded flow, and writes
`OUTPUT/<prefix>wall<sweep>.tec` ([Output](#output-walltec)). The parcel is removed at its first impact, as at
any wall; a rebound is not modelled. Without the section nothing is summed and no file is written. The code is
`src/lib/Lib_Erosion.f90` (laws, input checks, wall sums), `depositErosion` in `src/lib/Lib_Integration.f90`
(the deposit at the exit) and `write_wallfield` in `src/lib/IO.f90`.

---

## Wall cells

A wall cell is a boundary cell of code **301** (the dispersed-phase wall), **104** (a chimera interface with a block
of another phase) or **103** (a connection with a block of another phase, which ATLAS writes header only in a
dispersed-phase file: the walls a coupled solver gives). A 103 that names a partner cell in the file is a
connection: the parcel crosses it and nothing erodes. Every other code — outlets (400), symmetries (300),
inlets, connections — is no wall. ATLAS writes a manifold (501) and every gas-surface interaction (502–506) as a
symmetry (300) for a dispersed phase, so they do not erode either; a solid-propellant grain (502) whose section
gives the phase's inlet keys is an inlet ([Boundary Conditions](../user/boundary-conditions.md)).

## The impact

At the exit step, `bcDef` finds where the step crosses the face and interpolates the parcel's velocity linearly
along the step to that point, as it does for a reflection ([Geometry & Cell Tracking](geometry.md)); the parcel
keeps that velocity $\mathbf v_\mathrm{imp}$ and the face cell it crossed. The law takes

- the impact speed $V = |\mathbf v_\mathrm{imp}|$;
- the impact angle **from the surface**, $\alpha = \arcsin\left(|\mathbf v_\mathrm{imp}\cdot\hat{\mathbf n}|/V\right)
  \in [0, \pi/2]$, $\hat{\mathbf n}$ the face cell's normal (its orientation does not matter): $\alpha = 90°$ is a
  normal impact, $\alpha \to 0$ a grazing one;
- the parcel's diameter $d$ at the exit;

and returns $E$, the wall mass eroded per unit mass of impacting particles [kg/kg]. The parcel's flow is
$\dot m_p = \dot n_p\, m_p$, its eroded flow $\dot m_p E$.

## Output: `wall.tec`

One ORION/Tecplot zone `B<b>F<f>` per block face that holds a wall cell, on the face's nodes, the values
cell-centred; per face cell and for every material of the phase file (its name in the variable):

| variable | definition | unit |
| :--- | :--- | :---: |
| `er-<mat>` | $\sum_p \dot m_p E_p \,/\, \sum_p \dot m_p$ over the material's parcels that hit the cell (0 where none did): the mass-weighted erosion ratio | kg/kg |
| `g-<mat>` | $\sum_p \dot m_p \,/\, A$: the material's incident mass flux, $A$ the area of the cell's face | kg/(m² s) |
| `rdot` | $\sum_\mathrm{mat} g_\mathrm{mat}\, er_\mathrm{mat} \,/\, \rho_w$: the recession rate of the wall | m/s |

$\rho_w$ is `wall-density`, one value for every wall. The variables come in the order of the table: every `er`,
every `g`, then `rdot` ($2 n_\mathrm{mat} + 1$ of them). The sums of a run are reduced across MPI ranks in rank
order with the grid fields ([Parallelization](../development/parallelization.md)).

---

## Oka et al. (2005)

Oka, Okamura and Yoshida (Part 1) and Oka and Yoshida (Part 2) write the erosion damage as the wall **volume**
removed per unit mass of particles, at any impact angle as a fraction of the damage at normal impact:

$$
E(\alpha) = g(\alpha)\, E_{90}, \qquad
g(\alpha) = (\sin\alpha)^{n_1}\,\bigl(1 + H_v(1 - \sin\alpha)\bigr)^{n_2}
\qquad \text{(Part 1, eqs. 1–2)}
$$

$$
E_{90} = K\,(a H_v)^{k_1 b}\left(\frac{v}{v'}\right)^{k_2}\left(\frac{D}{D'}\right)^{k_3}
\qquad \text{(Part 2, eq. 4)}
$$

$E$ and $E_{90}$ in mm³/kg; $H_v$ the wall's (initial) Vickers hardness in **GPa**; $v$ the impact speed and $v'$
its reference [m/s]; $D$ the particle diameter and $D'$ its reference (any one unit: IGLOO takes both in m).
The first factor of $g$ carries the repeated deformation (or the brittle fracture of a ceramic, Oka et al. 2009),
the second the cutting, which acts most at shallow angles; $g(90°) = 1$ and $g(0) = 0$.

- $a$ and $b$ are the wall's load relaxation in a quasi-static indentation, $(F_0 - F_1)/F_0 = a\,(d/D)^{-b}$
  (Part 2, eq. 1), within $0 \le aH_v \le 1$, $b \ge 0$. With $a = b = 1$ the law is Part 1's,
  $E_{90} = K H_v^{k_1}(v/v')^{k_2}(D/D')^{k_3}$ (Part 1, eq. 4).
- $K$, $k_1$, $k_3$ are particle properties; $k_2$, $n_1$, $n_2$ depend on the particle and the wall. The papers
  correlate them with the hardness, $X = s\,H_v^{\,q}$ (Part 1, Table 3 for $n_1$, $n_2$; Part 2, Table 1 for $k_2$);
  IGLOO takes each of $k_2$, $n_1$, $n_2$ as a value or as its $s$ and $q$.

IGLOO's ratio is $E[\mathrm{kg/kg}] = 10^{-9}\,\rho_w\,E[\mathrm{mm^3/kg}]$, so $\rho_w$ cancels in `rdot`.

For reference, the papers' constants (they are calibrated for metallic targets; the user calibrates a pair of
interest and gives the constants in the input):

| particle | $K$ | $k_1$ | $k_2$ | $k_3$ | $v'$ [m/s] | $D'$ [µm] | $n_1 = s H_v^q$ | $n_2 = s H_v^q$ |
| :--- | ---: | ---: | :--- | ---: | ---: | ---: | :--- | :--- |
| SiO₂-1 | 65 | −0.12 | $2.3\,H_v^{0.038}$ | 0.19 | 104 | 326 | $0.71\,H_v^{0.14}$ | $2.4\,H_v^{-0.94}$ |
| SiO₂-2 | 50 | −0.12 | $2.3\,H_v^{0.038}$ | 0.19 | 104 | 326 | as SiO₂-1 | as SiO₂-1 |
| SiC | 45 | −0.05 | $3.0\,H_v^{0.085}$ | 0.19 | 99 | 326 | $0.71\,H_v^{0.14}$ | $2.8\,H_v^{-1.00}$ |
| glass beads | 27 | −0.16 | 2.1 | 0.19 | 100 | 200 | $2.8\,H_v^{0.41}$ | $2.6\,H_v^{-1.46}$ |

The 2009 paper applies the same $g(\alpha)$ to ceramic targets (alumina, zirconia, magnesia, silicon carbide), with
$n_1$, $n_2$ fitted per pair: the maximum then sits at 90°, as for a brittle material.

## Neilson and Gilchrist (1968)

Neilson and Gilchrist split the erosion of a mass $M$ of particles into **cutting** wear, which takes $\varphi$ of
the parallel kinetic energy per unit eroded mass, and **deformation** wear, which takes $\varepsilon$ of the normal
one (eqs. 3–5):

$$
\frac{W}{M} =
\begin{cases}
\dfrac{V^2\cos^2\alpha\,\sin n\alpha}{2\varphi} + \dfrac{(V\sin\alpha - K)^2}{2\varepsilon}, & \alpha < \alpha_0 \\[2ex]
\dfrac{V^2\cos^2\alpha}{2\varphi} + \dfrac{(V\sin\alpha - K)^2}{2\varepsilon}, & \alpha \ge \alpha_0
\end{cases}
\qquad \alpha_0 = \frac{\pi}{2n}
$$

$W/M$ is already a mass ratio [kg/kg]. Below $\alpha_0$ a particle leaves the surface again with part of its
parallel velocity, $1 - \sin n\alpha$ of its parallel energy; at $\alpha_0$ the two branches meet. $K$ (the input
key `k`) is "the velocity component normal to the surface below which no erosion takes place": IGLOO sets the
deformation term to zero where $V\sin\alpha \le K$ (the equations as printed would give a positive square there,
growing as the normal velocity drops below $K$, against the definition; with $K$ the paper's glass spheres on glass
do not erode below 35°, as its Fig. 11 shows). A large $\varphi$ leaves deformation wear alone, the paper's brittle
case.

The paper's constants are in ft·lbf per gram of eroded material and ft/s: 1 ft·lbf/g = 1355.818 J/kg,
1 ft/s = 0.3048 m/s. For aluminium and 210 µm angular alumina at 424 ft/s (Table I) $\varphi$ = 3571 ft·lbf/g
(4.842 MJ/kg), $\varepsilon$ = 8375 ft·lbf/g (11.35 MJ/kg), $n$ = 5; the paper reports that $\varphi$ and
$\varepsilon$ depend on the particle velocity, their ratio little.

---

## Inputs

`[IGLOO-Erosion]` in `input.ini` ([Input reference](../user/input.md)). The constants are model-specific and carry
no prefix: each key takes one value per material of the phase file, in its order, and each value is read for the
material's own model (the value of a material whose model does not use the key is ignored; give any number). Every
key belongs to one model, and the keys are case-sensitive: `K` is Oka's particle factor, `k` Neilson-Gilchrist's
threshold velocity (the paper's $K$).

| Key | Meaning | Required |
| :--- | :--- | :--- |
| `model` | one token per material: `Oka`, `Neilson-Gilchrist` or `none` (any case); `none` erodes nothing, its incident flow still counts | yes |
| `wall-density` | $\rho_w$ [kg/m³], every wall | yes, > 0 |
| `Hv`, `a`, `b` | Oka: the wall's $H_v$ [GPa], $a$, $b$ — one value each | with an Oka material; $H_v > 0$, $a > 0$, $b \ge 0$ |
| `K`, `k1`, `k3`, `v-ref`, `D-ref` | Oka: $K$ [mm³/kg], $k_1$, $k_3$, $v'$ [m/s], $D'$ [m] | with an Oka material; $K$, $v'$, $D' > 0$ |
| `k2` or `k2-s`, `k2-q` | Oka: $k_2$, or $k_2 = s H_v^q$ | one of the two forms; $k_2 > 0$ |
| `n1` or `n1-s`, `n1-q` | Oka: $n_1$, or $n_1 = s H_v^q$ | one of the two forms; $n_1 > 0$ |
| `n2` or `n2-s`, `n2-q` | Oka: $n_2$, or $n_2 = s H_v^q$ | one of the two forms |
| `phi`, `epsilon`, `n`, `k` | Neilson-Gilchrist: $\varphi$, $\varepsilon$ [J/kg], $n$, $K$ [m/s] | with a Neilson-Gilchrist material; $\varphi, \varepsilon, n > 0$, $K \ge 0$ |

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
k            = 0 20.0                    ; B: the threshold velocity K [m/s]
```

Setup refuses a model token it does not know, a key whose number of values is not the number of materials, a
missing constant of a model in use, $k_2$, $n_1$ or $n_2$ in neither or both forms, a value out of its range, a
token that is not a number and more than one value of a wall constant; the log names each material's model and the
constants it uses ($n_1$, $n_2$, $k_2$ evaluated when given as $s$, $q$).

---

## V&V

Unit family `tests/erosion/laws/` (Oka's $g(\alpha)$ against the maxima the papers state, $E_{90}$ against its
closed form and Part 1's limit, the conversion to kg/kg; Neilson-Gilchrist's branches, their meeting at $\alpha_0$,
the maximum against eq. 7 for Table I's three speeds, the 35° threshold of Fig. 11), the e2e cases `wall-flux` (the
wall file and the incident flux with `model = none`) and `wall-laws` (both laws in the wall file against an
independent oracle, eight impacts from 10° to 90°, the two forms of $n_1$, $n_2$, $k_2$), and fifteen setup
refusals (`tests/erosion/refusals/`). See `tests/erosion/*/INFO.md`.

## References

- Y.I. Oka, K. Okamura, T. Yoshida, Practical estimation of erosion damage caused by solid particle impact. Part 1:
  Effects of impact parameters on a predictive equation, *Wear* 259 (2005) 95–101.
- Y.I. Oka, T. Yoshida, Practical estimation of erosion damage caused by solid particle impact. Part 2: Mechanical
  properties of materials directly associated with erosion damage, *Wear* 259 (2005) 102–109.
- Y.I. Oka, S. Mihara, T. Yoshida, Impact-angle dependence and estimation of erosion damage to ceramic materials
  caused by solid particle impact, *Wear* 267 (2009) 129–135.
- J.H. Neilson, A. Gilchrist, Erosion by a stream of solid particles, *Wear* 11 (1968) 111–122.
