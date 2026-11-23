# Eulerian Feedback

The parcels leave two kinds of fields on the gas mesh:

- **Source fields**: the mass, momentum and energy the parcels exchange with the gas in each cell. They are the
  source terms of a two-way coupling.
- **Euler fields**: the averaged state of the particles in each cell: density, velocity, temperature and number
  density, built from moments integrated along each trajectory.

Which fields are written: [Run Settings](../user/general.md#what-igloo-writes). A coupled solver reads them through
IGLOO's exchange functions ([Coupling with other solvers](../development/coupling.md)).

---

## Source fields

Each cell holds the mass source of every material [kg/s], the momentum source [N] and the energy source [W]: rates
over the whole cell, not densities.

For each passage of a parcel through a cell (a **residence**), the cell receives what the parcel carried in minus what
it carried out:

$$
S_\mathrm{mass}(m) \mathrel{+}= \dot{m}_\mathrm{in} - \dot{m}_\mathrm{out}, \qquad
\mathbf{S}_\mathrm{mom} \mathrel{+}= \mathbf{P}_\mathrm{in} - \mathbf{P}_\mathrm{out}, \qquad
S_\mathrm{en} \mathrel{+}= E_\mathrm{in} - E_\mathrm{out},
$$

with $\dot m$ the parcel's mass flow, $\mathbf P = \dot m\,\mathbf v_p$ and $E = \dot m\,(h + |\mathbf{v}_p|^2/2)$, $h$ the
material's specific enthalpy at the parcel's state. For a freezing material $h$ includes the latent heat of fusion
([Solidification](solidification.md)), so a freezing parcel deposits the heat it loses on the plateau.

- A material whose parcels keep their mass deposits no mass.
- A drop consumed inside a cell (burnout, or complete evaporation) carries nothing out, so everything it brought in
  is deposited there.
- A KH-RT shed lowers the parent's number rate and starts a child at the shed point. The child's birth flux,
  $\dot m_c\,(1, \mathbf{v}_p, h + |\mathbf{v}_p|^2/2)$, counts as part of what the parent carried out, and the child
  deposits it along its own path: the shed itself gives the gas nothing.
- On an axisymmetric wedge the momentum source is stored in the **meridian frame**: each residence's momentum is
  rotated back to the meridian plane, so every cell holds (axial, radial, azimuthal) components whatever the parcel's
  azimuth ([Boundary conditions](../user/boundary-conditions.md#symmetry-handling-code-300)).

### Rows of the mass source

The mass source holds one **net** row per material, and further rows only when something needs them. With $\dot n$ a
parcel's number rate:

| Row | Present when | Holds, per cell [kg/s] |
| :--- | :--- | :--- |
| net | always | the mass the material's parcels give the gas, $\dot m_\mathrm{in} - \dot m_\mathrm{out}$ |
| carried in | a material carries another | the carried material deposited on the material's parcels, $\dot n\,(m_\mathrm{ox,out} - m_\mathrm{ox,in})$ per residence |
| carried out | same | the carried material a parcel consumed in the cell gives the gas when it forms no residue, $\dot n\,m_\mathrm{ox}$ at its end |
| processed | a coupler asks for it | the mass $W$ the process (burning, evaporation) takes off the parcels, and $W\mathbf u$, $WT$ |

Only drops that carry a material fill the carried rows. The net row of a carrying material counts both masses: the
metal burnt less the carried material deposited; for a parcel consumed in the cell, everything it carried at cell
entry, less its cap when that becomes a residue ([The burnout residue](combustion.md#the-burnout-residue)). In every
case

$$
S_\mathrm{net} + S_\mathrm{carried\ in} - S_\mathrm{carried\ out}
= \sum_\mathrm{residences} \dot n\,(m_\mathrm{Al,in} - m_\mathrm{Al,out}),
$$

the metal that left the material's parcels in the cell, whatever the share kept. A coupler reads the burn rate from
these three rows ([Energy](combustion.md#energy-the-absolute-datum-and-the-streams)).

The processed rows weight the velocity and temperature by the mass processed: over each step, the mass it processed
times the step's mean velocity and temperature. The processed mass is valued at the state of the process, not at the
end of the residence.

What the parcels **absorb** from a donor is kept apart from these rows, per geometry cell and never smoothed, so the
coupler takes from the donor exactly what was absorbed where it offered it
([Coupling](../development/coupling.md#receiving-mass)).

### The body force

With a body acceleration $\mathbf g$ (gravity), the parcels' momentum and energy change also by the body force's work,
which the gas must not receive. The deposit is corrected so that the gas gets the drag reaction only:

- for parcels of constant mass flow (inert, breaking, freezing), in closed form:
  $\mathbf{P}_\mathrm{in} \mathrel{+}= \dot m\,T_\mathrm{stay}\,\mathbf{g}$ and
  $E_\mathrm{in} \mathrel{+}= \dot m\,\mathbf{g}\cdot(\mathbf{x} - \mathbf{x}_\mathrm{entry})$, with the residence time
  $T_\mathrm{stay}$ and the entry position;
- for parcels whose mass changes, by integrating $J = \int \dot m\,dt$ and $W = \int \dot m\,\mathbf{g}\cdot\mathbf{v}\,dt$
  along the trajectory with the state.

### Second order: from the dual mesh to the cells

At the second gas order the parcels deposit on the gas dual mesh, whose cells surround the mesh nodes. The deposits are
then reduced to the geometry cells **conservatively**, since they are rates (kg/s, N, W), not densities. Each geometry
cell $c$ collects, from every dual cell $d$ its octants touch,

$$
S_\mathrm{geo}(c) \;=\; \sum_\mathrm{oct} \frac{V_\mathrm{oct}(c)}{W(d)}\, S_\mathrm{dual}(d),
\qquad
W(d) \;=\; \sum_{\text{(cell, octant) touching } d} V_\mathrm{oct},
$$

so the weights of each dual cell sum to one over the geometry cells that share it, and
$\sum_c S_\mathrm{geo} = \sum_d S_\mathrm{dual}$ **exactly, on any mesh**. Dividing by the cell's own octant sum instead
would be a volume average: right for an interior dual cell, but it would keep only ½ (face), ¼ (edge) or ⅛ (corner) of
a boundary dual cell on a Cartesian mesh, and so discard most of the deposit exactly where an injected spray puts it.
The assembled weight is used rather than the clipped dual volume, so that a conservation identity does not depend on
an approximate geometry.

### Block interfaces

The dual cell of a node on a conformal interface (101) straddles the face, and each block holds only the octants on its
side; the dual cell of a node on a periodic face (201) is split the same way, its other half one period away. At the
second gas order the halves are joined: the dual cells of every face node, matched through the face's orientation (for
a periodic face, modulo the period's offset $\mathbf t$), form one **node class**. The face's ghost nodes sit at the
partner cells' centres (brought back by $\mathbf t$ across a period), the dual metric is assembled without a clip
normal to the face, and every member takes the class's summed octant volume $W$, so that the weights sum to one over
both sides. Before any field is reduced, every member takes the class's source sum and, for the Euler fields, the
class's amount over that volume, $\sum_m V_m q_m / W$. Each block then reduces the share of its own octants, and the
two blocks together reduce what one block would.

The integration still ends a step at the block face, and at a periodic face, where the parcel is moved
([Geometry](geometry.md#crossing-a-block-interface)). The source and the Euler densities add over the cut. The velocity
and temperature numerators, arc-length means over a step times its residence, differ from those of the uncut step at
second order: a cut at the fraction $f$ of a step over which the speed and the transported quantity change by the
relative amounts $\varepsilon$ and $\eta$ changes the numerator by $-\varepsilon\eta\,f(1-f)/4$ of itself, to leading
order. On a curved mesh the octant sum $W$ differs from the dual cell's metric volume $V$ (a curved hexahedron's volume
is not the sum of its octants'): next to the face the two blocks then differ from one block by up to that ratio
($7.6\times10^{-5}$ on the wedge of the two-block test).

### Chimera overlaps

A parcel deposits into the block it is tracked in. On a mesh with chimera (overset) faces the gas solver keeps a copy
of the overlap in each block, so every deposit cell $d'$ of a block $t$ takes the parts of the other blocks' deposits
made inside $t$, split by the exact intersection of the cells:

$$
S'(d') = S(d') + \sum_{d} W(d \to d')\, I(d, t), \qquad W(d \to d') = \frac{|d \cap d'|}{|d \cap t|},
$$

$\sum_{d'} W = 1$, and $I(d, t)$ the part of $S(d)$ made inside $t$: all of it when $d$ lies inside $t$, none when
outside, and when $t$'s boundary cuts $d$ the sum over the steps in $d$ of each step's deposit times the fraction of
its chord inside $t$ (the deposit taken uniform along the chord). The Euler densities move with the amounts,
$q'(d') = q(d') + \sum_d I_q(d, t)\,\big(V(d) W\big)/V(d')$. The mirrored part of each block is reported per sweep, so
the total counted once is the sum of the final source over every block minus the mirrored parts, exactly
([Boundary conditions](../user/boundary-conditions.md#chimera-faces-code-102)).

---

## Euler fields

Per family and cell, the Euler fields accumulate, over the residences, with $T_\mathrm{stay}$ the residence time,
$\delta L$ the arc length and $V$ the deposit cell's volume:

| Field | Numerator |
| :--- | :--- |
| density | $\sum \rho_p \cdot T_\mathrm{stay} / V$ |
| number density | $\sum \dot{n}_p \cdot T_\mathrm{stay} / V$ |
| velocity | $\sum \rho_p\,\mathbf{v}_p \cdot \delta L^{-1} \cdot T_\mathrm{stay} / V$ |
| temperature | $\sum \rho_p\,T_p \cdot \delta L^{-1} \cdot T_\mathrm{stay} / V$ |

The density row is the **constant-mass** form, with $\rho_p\,T_\mathrm{stay} = \dot m\,T_\mathrm{stay}$, $\dot m$ the stream
flow. For a material that loses mass the drop's mass is inside the integral, and the **number rate multiplies it**: the
parcel carries $\dot n_p$ drops per second, not one,

$$
\rho_p^{\,\mathrm{cell}} \;=\; \frac{1}{V}\sum_p \dot{n}_p \int_{\mathrm{cell}} m_p\,\mathrm{d}t .
$$

This matters beyond the density: the velocity and temperature are their numerators **divided by** the density, so a
density missing $\dot n_p$ would scale them by it too.

The integrals come from moments integrated with the parcel's state: the arc length, the velocity and temperature
along it, and the mass or number moment. At the second gas order $V$ is the dual cell's. On a wedge the momentum moment
is in the meridian frame, as for the source. A step that ends on a numerical failure is taken back to its last accepted
state and deposited up to it.

### From the numerators to the averages

After every parcel is integrated, the fields are reduced to the geometry cells (the node classes of an interface
joined first), mirrored into every block that covers a chimera overlap (preserving the amount), smoothed
([Field Mollification](mollification.md)), and the velocity and temperature numerators are divided by the density:

$$
\langle\mathbf{v}_p\rangle = \frac{\sum \rho_p\,\mathbf{v}_p\,\delta L^{-1}\,T_\mathrm{stay}/V}{\sum \rho_p\,T_\mathrm{stay}/V}.
$$

When the state is the enthalpy, the averaged enthalpy is turned back into a temperature through the material's table.
Summed over the geometry cells, $\rho_p V$ equals the injected mass flow times the residence time.
