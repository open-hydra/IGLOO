# inj-pack-polar — one hexagonal lattice over an annular-sector inlet

**Purpose.** On a 3D inlet face the advancing front of `pin_particles_bc_ds` (`initialization.f90`) spawns six
candidates around every placed point at the spacing R = ds, turned about the normal of the point's cell from the
tangent frame the point carries: the seed takes its cell's first m-edge in the cell's plane, a child its parent's
frame projected onto the child cell's plane. On a planar face the frame is the same at every point, so the front
builds one hexagonal lattice over the whole face, whatever the cells' edges do. On this face the first edge of a
cell (its radial edge) turns by 5 degrees from one column of cells to the next.

**Fixture** (`python3 make_fixture.py`: `make_box_case.py --polar 0.04,0.07,1.0471975511965976 --ny 6 --nz 12`;
`solfile.tec` and `bc.txt` committed, `phase.txt` and `properties.dat` from `common/`): `drag-stokes`'s box and gas
(U = 10 m/s along x, 60 cells along x) with its cross-section mapped to an annular sector, x = 0 inlet, r 4–7 cm over
6 radial cells, θ 0–60° over 12 azimuthal cells (5 mm × 3.49–6.11 mm), the whole inlet 401 (κ_ρ = 0.34).
`input.ini` is `drag-stokes`'s with `ds = 0.25` (cm), below every cell's size, and `out-time = on`.

**Gates** (`check.py`; the oracle is the face's nodes, read from `INPUT/solfile.tec`, the inputs and conservation;
the pinned points are the parcels' first trajectory rows, F12.6).
- U0 non-vacuity: at least 300 parcels, each of the 72 cells holding at least 3 first rows (a row counting for every
  cell that holds it within 2e-6 m of its four edges), the seed (parcel 1's first row) at the vertex mean of face
  cell (1, 1) to 1e-6 m.
- U1 every pair of first rows at least R − 1.5e-6 m apart (two printed points).
- U4 every first row within 2e-6 m of the lattice {s + i a1 + j a2}: s the seed from the nodes, a1 = R ê with ê the
  unit vector of face cell (1, 1)'s first radial edge, a2 = a1 turned by 60° about x.
- U5 every lattice node at least 2e-6 m inside a face cell (from each of its four straight edges) has a first row
  within 2e-6 m.
- J3 (control) per inlet cell the parcels' mass flows (exit records, by header name) sum to
  κ_ρ/(1 − κ_ρ)·ρ_g U A_c, A_c the planar quad's area from the nodes, to the E13.6 print; a first row on an edge two
  cells share belongs to the holding cell whose interior parcels carry its mass flow (`inj-pack`'s rule).
- W (control) every parcel has an exit record at x = 0.15.

GREEN: 312 parcels, 3 to 6 first rows per cell; the closest pair 2.499945e-3 m (R at the print: offsets R/2 and
R√3/2 printed as 1.25 and 2.165 mm); every first row within 5.1e-7 m of the seed's lattice; all 304 lattice nodes
inside a cell pinned; 8 first rows on face-cell edges, 6 of them on edges two cells share, each filling the cell
whose share it carries.

**RED** — the tangent taken per cell (each point expanded with its own cell's first m-edge): 391 parcels, each
column of cells packing its own rotated lattice; U1 fails on 308 first rows (closest pair 1.258223e-3 m =
0.5033 R), U4 on 357 (up to 1.392e-3 m off the seed's lattice), U5 on 277 of the 304 lattice nodes; J3 and W stay
green (the mass split and the exits hold under either frame).
