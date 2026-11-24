# Field Mollification

After every parcel is integrated, the deposited source and Euler fields can be smoothed by a dimension-split binomial
smoother, on by default with 8 passes ([Run Settings](../user/general.md#smoothing-the-fields)).

---

## Rationale

Lagrangian deposition concentrates particle contributions in cells visited by particle
trajectories.  With few particles per cell the resulting field has strong Nyquist-scale
(2-cell) noise.  A standard diffusion smoother with a fixed physical width amplifies
near-wall mesh non-uniformities (face-distance $d_{CN}\to 0$ on boundary-layer grids).

The binomial smoother avoids both problems:

- **Annihilates Nyquist modes exactly**: the per-sweep amplification factor for a mode of
  wavenumber $k$ is $g(k) = 1 - 2\theta(1 - \cos k)$, so $g(\pi) = 1 - 4\theta = 0$ for
  $\theta = 1/4$.  A single sweep kills the 2-cell mode completely in each axis.
- **Mesh-local width**: the conservative interior-face flux
  $\phi_f = \theta\,\min(V_c, V_n)\,(q_n - q_c)$ uses only cell volumes, not face
  distances.  The diffusion number $\nu_c = \theta\sum_f \min(V_c,V_n)/V_c \le \theta \le 1/2$
  is unconditionally stable.
- **Conservative**: $\sum q V$ is preserved to machine precision at every pass.

---

## Algorithm

Each pass applies one 1-D sweep per axis (x, then y, then z in 3D) in a dimension-split
sequence, each sweep reading the values of the previous one:

$$
\phi_f = \theta\,\min(V_c,\,V_n)\,(q_n - q_c) \quad \text{[interior faces only]}
$$

$$
q_c^\mathrm{new} = q_c + \frac{\phi_{c+} - \phi_{c-}}{V_c}
$$

Domain boundaries are Neumann (zero flux): no ghost cells, no artificial influx from
outside the domain.

The result after $n$ passes is a Gaussian-equivalent filter with width

$$
\sigma \approx \sqrt{n/2} \text{ cells}
$$

(per-pass variance $= 1/2$ cell² for the $\theta = 1/4$ binomial kernel; independence
from mesh cell size is by design).

### Block interfaces and periodic faces

A block face joined to another block by a conformal connection (code 101) is an interior
face of the domain, and the smoother treats it as one: before each sweep direction of each pass,
every face cell of such an interface takes a one-cell halo, the current value $q_p$ and the volume
$V_p$ of its partner cell in the other block, and the flux across the face is

$$
\phi_f = \theta\,\min(V_c,\,V_p)\,(q_p - q_c),
$$

the interior-face formula with the same operands. The two sides of the interface use the same two
values in the same stage, so their fluxes are equal and opposite and $\sum q V$ over all blocks is
conserved; a mesh split into blocks along whole index planes is smoothed as the unsplit mesh (the
arithmetic is the same; the default optimised build may differ in the last bits). The halo needs the
two faces to be normal to the same index direction (an i-face joined to an i-face, and so on); an
interface whose faces are normal to different index directions keeps zero flux, and the setup prints
one warning line per such face pair. A translational periodic pair (code 201) is smoothed the same
way: each face cell's partner is the cell its connection names on the paired face, in the same block
or another, so the cells next to the two faces are neighbours, as in one unrolled period (periodicity
is translational, so the momentum source needs no rotation across the face). Every other block face —
walls, inlets and outlets — keeps zero flux. The conservation guard sums $q V$ over every block.

A chimera face (code 102) lies inside another block, whose cells there are a second copy of the
overlap. Before each sweep along the face's index direction, each facelet with layer-1 donors takes a
halo blended from them, $q_h = \sum_d w_d\,q_d$ and $V_h = \sum_d w_d\,V_d$ (ATLAS's volume weights, the
donors' current values), and its cell takes the flux $\theta\,\min(V_c, V_h)\,(q_h - q_c)$; the donors
are cells of their own block, smoothed there, and are not debited. Each block's total therefore changes
by exactly what its chimera cells took, $F_{102}$ (summed link by link in a fixed order), and the guard
checks $|\sum qV_\mathrm{after} - \sum qV_\mathrm{before} - F_{102}| \le 10^{-10} \sum |q|V$ ($F_{102} = 0$
without chimera faces). For the source, each block's exchange is added to its mirrored totals (see
[Boundary conditions](../user/boundary-conditions.md#chimera-faces-code-102)), so the sum over the blocks
minus the mirrored totals is still the once-counted total. Where the overlap's cells coincide (one donor
of weight 1 per facelet), the overlap's deposits mirrored into both blocks and the halo give each block
the one-block smoothed field on its own cells; with misaligned cells (two donors of weight 1/2 at a
half-cell shift) the halo is the average of the two donor cells and the result meets the one-block field
only to the interpolation error of the halo and of the overlap's volume split: on the shifted two-block test
each block's smoothed fields lie within 2(max\|ΔR\|/16 + 2 max\|Δ²R\|/8) of the one-block ones, R the unsmoothed
one-block cells around the seam (measured at 1/13 to 1/40 of that bound). A facelet with an empty layer 1
keeps zero flux.

---

## Default parameters

| Parameter | Value |
| :--- | :---: |
| $\theta$ | 0.25 |
| passes | 8 |
| effective width $\sigma$ | 2 cells ($\sqrt{8/2}$) |

With 8 passes:

- Modes of wavelength $\le 4$ cells are attenuated to $< 0.4\%$ (Nyquist killed on pass 1).
- Modes of wavelength $\ge 16$ cells are preserved at $\ge 73\%$.

This provides effective deposition-noise suppression while leaving resolved physical
structure intact; more passes give the same noise kill with more diffusion of the
$\lambda \approx 8$-cell band.

---

## Applied fields

The smoother runs after the fields are reduced to the geometry cells and before the Euler fields are averaged. It
smooths the conserved quantities:

- of the Euler fields, the numerators: density ($\rho_p\,T_\mathrm{stay}/V$), the three momentum components, the energy,
  and the number density;
- of the source fields, the per-cell totals: every row of the mass source, the three momentum components and the
  energy.

Each field is divided by the cell volume before smoothing and multiplied back after, which keeps $\sum q V$.
