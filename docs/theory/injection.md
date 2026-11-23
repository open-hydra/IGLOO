# Injection & Packing

Parcels are placed once, at setup, for every injection group. Two ways are available
([Placing the parcels](../user/boundary-conditions.md#placing-the-parcels-igloo-bc)):

- **On the inlet faces** (the default): parcels are seeded on the boundary faces whose code injects (inlets 401–420,
  and the solid-propellant grains 501–502).
- **At given points**: the input lists each parcel's position, velocity, temperature, mass flow and diameter.

A cell that injects nothing for a family (a zero loading: a 401 of $k_\rho$ = 0, as ATLAS writes it for a phase without
inlet keys at a grain, or a 402/403 of zero flux) keeps no parcel of that family. Every parcel is placed first, so the
diameters drawn for the others are those of a run without the filter; then the parcels of such cells are removed and
the rest renumbered, and the log counts them. A family left without parcels runs empty.

---

## Placing parcels on the inlet faces

With no spacing and no mass-flow cap, each inlet cell gets **one parcel at its centre** (or every $n$-th cell, with
subsampling). Otherwise the parcels are placed at a **spacing** $R$ chosen per cell:

1. a spacing given for the cell in the inlet table overrides the global one;
2. a cell no larger than the spacing (or with no spacing) gets exactly one centred parcel;
3. with a cap on the mass flow of one parcel, the cell needs at least $N_\mathrm{min}$ parcels, from its mass flux and
   area, and $R = \min(R, d_\mathrm{cell}/N_\mathrm{min})$.

$d_\mathrm{cell}$ is the cell's tangential size along the face. Cells smaller than a set floor are skipped: they are
usually collapsed cells at edges, and would start parcels in ill-conditioned places.

### On a 2D face

A front advances along the face: the next point is where a circle of radius $R$ around the previous one meets the cell
boundary. Single cells get their centred point first; a last pass gives one centred parcel to every inlet cell the
front missed.

### On a 3D face: hexagonal packing

A front spreads over the face from the seeded single cells:

1. Each placed point proposes six neighbours at distance $R$, by turning the point's tangent frame about its cell's
   normal in steps of $\pi/3$. The first point's frame is its cell's first edge; each new point inherits its parent's
   frame projected onto its own cell's plane (or its own cell's first edge where that projection vanishes, at a 90°
   fold between two inlet cells). On a planar face the frame is then the same everywhere, and the points form **one
   hexagonal lattice over the whole face**, whatever the cells' edges do.
2. Each candidate is projected onto its cell's plane and kept if the cell holds it and no placed point lies within
   $R/2$.
3. A last pass gives one centred parcel to every inlet cell the front missed.

### The cell a parcel starts in

A parcel starts in the inlet cell that placed it whenever that cell holds its point (with the cell's tolerance,
[Geometry](geometry.md#is-a-point-inside-a-cell)): a point on an edge two inlet cells share lies in both, and one on the
domain boundary may lie a few rounding units outside. Otherwise the cell search locates it. The parcel takes the
loading, velocity, direction and temperature of the cell it starts in (its diameter was drawn when it was placed). A
parcel whose point no cell holds is not integrated, and the log counts such parcels at every sweep.

---

## Stochastic diameter sampling

Once every position is placed, each parcel draws its diameter from three properties of its inlet cell: a radius $r_p$,
a width $\sigma_p$ and a law ([Injection patches](../user/boundary-conditions.md#injection-patches)).

| Law | Parameters | Cumulative distribution of the injected diameter $d$ |
| :--- | :--- | :--- |
| Dirac | $2r_p$ (also every law with $\sigma_p \le 0$) | a step at $2r_p$; no draw |
| Normal | $\mu = 2r_p$, $\sigma = \sigma_p$ | truncated to $d > 0$ by rejection: $\dfrac{\Phi((d-\mu)/\sigma) - \Phi(-\mu/\sigma)}{1 - \Phi(-\mu/\sigma)}$ |
| LogNormal | mean $\mu = 2r_p$, standard deviation $\sigma_p$ | $\Phi\big((\ln d - m)/s\big)$, $s^2 = \ln(1 + \sigma_p^2/\mu^2)$, $m = \ln\mu - s^2/2$ |
| Rosin–Rammler | $x_0 = 2r_p$, $n = \sigma_p$ | $1 - \exp\big(-(d/x_0)^n\big)$ |
| size-distribution file | radii $r_k$, mass fractions $w_k$ (normalised) | piecewise linear through the class edges $e_k$ (mid-points between radii, the end classes extended by half a spacing, doubled to diameters): $F(e_k) = \sum_{j \le k} w_j$; one row is a Dirac at $2r_1$ |

**What a law describes.** In one inlet cell every parcel carries the same mass flow, $\dot m_\mathrm{cell}/N_\mathrm{inj}$,
whatever diameter it drew, and its number rate is $\dot m/m$ with $m = \pi\rho_p d^3/6$. The sampled law $f(d)$ is
therefore the distribution of the injected **mass flow** over the diameter, and the number of particles goes as
$f(d)/d^3$ (normalised). That is why a size-distribution file lists mass fractions. (TAB's Rosin–Rammler product size
is a number distribution, a different object.)

**Random stream.** The draws come from the compiler's random generator, seeded from the run's seed (42 by default) just
before the parcels are placed, and consumed in placing order. Placing is serial and repeated identically on every
process, so the diameters do not depend on the thread or process count; they do depend on the compiler. The draw is
made once: every sweep injects the same diameters. Draws made during the integration (TAB and ETAB) come from each
parcel's own stream ([New parcels](breakup.md#new-parcels)).

The laws are checked with $10^5$ draws each, and end to end with 20 000 injected parcels, against the
Kolmogorov–Smirnov and Anderson–Darling statistics at $\alpha = 10^{-3}$: [Injection V&V](../vv/injection.md).

---

## Parcels at given points

The listed positions, velocities, temperatures, mass flows and diameters become the parcels as they are; the cell
search locates each one when its integration starts.
