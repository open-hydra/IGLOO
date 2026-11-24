# Wall Erosion

When erosion is on, every parcel that leaves the domain through a **wall cell** erodes it. The erosion law of the
parcel's material is evaluated at the parcel's impact speed, angle and diameter, and summed per wall face cell and
per material, weighted by the parcels' flow. The parcel is removed at its first impact, as at any wall: a rebound is
not modelled. The input and the wall file: [Wall Erosion](../user/erosion.md).

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

At the exit step, IGLOO finds where the step crosses the face and interpolates the parcel's velocity linearly
along the step to that point, as it does for a reflection ([Geometry & Cell Tracking](geometry.md)); the parcel
keeps that velocity $\mathbf v_\mathrm{imp}$ and the face cell it crossed. The law takes

- the impact speed $V = |\mathbf v_\mathrm{imp}|$;
- the impact angle **from the surface**, $\alpha = \arcsin\left(|\mathbf v_\mathrm{imp}\cdot\hat{\mathbf n}|/V\right)
  \in [0, \pi/2]$, $\hat{\mathbf n}$ the face cell's normal (its orientation does not matter): $\alpha = 90°$ is a
  normal impact, $\alpha \to 0$ a grazing one;
- the parcel's diameter $d$ at the exit;

and returns $E$, the wall mass eroded per unit mass of impacting particles [kg/kg]. The parcel's flow is
$\dot m_p = \dot n_p\, m_p$, its eroded flow $\dot m_p E$.

## What is summed

Per wall face cell and material, over the material's parcels $p$ that hit it:

$$
er = \frac{\sum_p \dot m_p E_p}{\sum_p \dot m_p}, \qquad
g = \frac{\sum_p \dot m_p}{A}, \qquad
\dot r = \frac{1}{\rho_w}\sum_\mathrm{materials} g\,er,
$$

the mass-weighted erosion ratio, the incident mass flux ($A$ the face area of the cell) and the recession rate of the
wall, $\rho_w$ the wall density.

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

IGLOO's ratio is $E[\mathrm{kg/kg}] = 10^{-9}\,\rho_w\,E[\mathrm{mm^3/kg}]$, so $\rho_w$ cancels in the recession
rate.

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
parallel velocity, $1 - \sin n\alpha$ of its parallel energy; at $\alpha_0$ the two branches meet. $K$ is "the velocity component normal to the surface below which no erosion takes place": IGLOO sets the
deformation term to zero where $V\sin\alpha \le K$ (the equations as printed would give a positive square there,
growing as the normal velocity drops below $K$, against the definition; with $K$ the paper's glass spheres on glass
do not erode below 35°, as its Fig. 11 shows). A large $\varphi$ leaves deformation wear alone, the paper's brittle
case.

The paper's constants are in ft·lbf per gram of eroded material and ft/s: 1 ft·lbf/g = 1355.818 J/kg,
1 ft/s = 0.3048 m/s. For aluminium and 210 µm angular alumina at 424 ft/s (Table I) $\varphi$ = 3571 ft·lbf/g
(4.842 MJ/kg), $\varepsilon$ = 8375 ft·lbf/g (11.35 MJ/kg), $n$ = 5; the paper reports that $\varphi$ and
$\varepsilon$ depend on the particle velocity, their ratio little.

---

## V&V

Oka's $g(\alpha)$ against the maxima the papers state and $E_{90}$ against its closed form, Neilson–Gilchrist's
branches and their meeting at $\alpha_0$, the 35° threshold of their Fig. 11, and both laws in the wall file against an
independent oracle, eight impacts from 10° to 90°: the erosion tests of the suite ([Testing](../development/testing.md)).

## References

- Y.I. Oka, K. Okamura, T. Yoshida, Practical estimation of erosion damage caused by solid particle impact. Part 1:
  Effects of impact parameters on a predictive equation, *Wear* 259 (2005) 95–101.
- Y.I. Oka, T. Yoshida, Practical estimation of erosion damage caused by solid particle impact. Part 2: Mechanical
  properties of materials directly associated with erosion damage, *Wear* 259 (2005) 102–109.
- Y.I. Oka, S. Mihara, T. Yoshida, Impact-angle dependence and estimation of erosion damage to ceramic materials
  caused by solid particle impact, *Wear* 267 (2009) 129–135.
- J.H. Neilson, A. Gilchrist, Erosion by a stream of solid particles, *Wear* 11 (1968) 111–122.
