# Field Mollification

After all particles are integrated, the deposited source and Euler fields can be smoothed
by a dimension-split binomial smoother controlled by `mollify-passes` in
`[IGLOO-General]`.  Setting `mollify-passes = 0` (or omitting `src-switch`/`eul-switch`)
disables smoothing entirely.

The smoother is implemented in `src/lib/Lib_Mollify.f90::binomial_smooth`, called by
`mollifyEuler` and `mollifySource` in `src/lib/obj_block.f90`.

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
sequence, using a ping-pong buffer (GATHER form, race-free):

$$
\phi_f = \theta\,\min(V_c,\,V_n)\,(q_n - q_c) \quad \text{[interior faces only]}
$$

$$
q_c^\mathrm{new} = q_c + \frac{\phi_{c+} - \phi_{c-}}{V_c}
$$

Block boundaries are Neumann (zero flux): no ghost cells, no artificial influx from
outside the block.

The result after $n$ passes is a Gaussian-equivalent filter with width

$$
\sigma \approx \sqrt{n/2} \text{ cells}
$$

(per-pass variance $= 1/2$ cell² for the $\theta = 1/4$ binomial kernel; independence
from mesh cell size is by design).

---

## Default parameters

| Parameter | Value | Source |
| :--- | :---: | :--- |
| $\theta$ | 0.25 | `IGLOO_Lib_Mollify`, parameter |
| Default passes | 8 | `DEFAULT_MOLLIFY_PASSES` |
| Effective width $\sigma$ | 2 cells | $\sqrt{8/2}$ |

With 8 passes:

- Modes of wavelength $\le 4$ cells are attenuated to $< 0.4\%$ (Nyquist killed on pass 1).
- Modes of wavelength $\ge 16$ cells are preserved at $\ge 73\%$.

This provides effective deposition-noise suppression while leaving resolved physical
structure intact.  10 passes would give the same noise kill with more diffusion of the
$\lambda \approx 8$-cell band; 8 was preferred.

---

## Applied fields

`mollifyEuler` (`obj_eulerblock::finalize`, Phase 4) smooths the conserved density
numerators before Favre averaging:

- `density` ($\rho_p\,T_\mathrm{stay}/V$)
- `velocity(1:3)` (momentum numerators, each component separately)
- `temperature` (energy numerator)
- `np` (parcel-number density)

`mollifySource` (`obj_sourceblock::finalize`, Phase 4) smooths the source extensive
per-cell totals:

- `sourceMass(1:nm)` (per material)
- `sourceMom(1:3)`
- `sourceEn`

Each scalar field is divided by cell volume before smoothing and multiplied back out
after, satisfying the $\sum q V$ conservation invariant.

An optional pre-allocated scratch buffer (`work`) is reused across all 6 (Euler) or
$n_m+4$ (source) calls per finalize pass, avoiding repeated allocations in the inner
loop.

---

## Configuration

```ini
[IGLOO-General]
mollify-passes = 8   ; override default (0 = off)
```

`mollifyPasses = 0` when `mollify-on = F` (or `eul-switch = F` / `src-switch = F`).
Full parameter registry: [../user/registry.md](../user/registry.md).
