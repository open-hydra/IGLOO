# inlet-401-g — the 401 inlet's field G

**Purpose.** ATLAS writes a 401 record under its mixture rule over the dispersed phases: field 1 is
κ_ρ/(1 − κ_ρ,tot) and field 10 is G, the fixed mass flux of the cell (the sum of `gp` over its 402/403 records in the
files of every dispersed phase, the gas not in it). IGLOO injects F = field 1 · (|ṁ_gas| + G·A) per cell
(`obj_particles.f90`) and asks the placement for N_pmin = ⌈F/ṁ_max⌉ parcels with the same F (`compute_ds`,
`initialization.f90`). Every other fixture carries G = 0; this case gives G a value.

**Fixture.** `inj-mdotmax`'s case (`drag-stokes`'s box, uniform gas ρ_g U = 12 kg/(m² s) and 401 inlet, κ_ρ = 0.34;
`ds = 0`, `mdot-max = 0.25` g/s, `out-time = on`) with its own `INPUT/bc.txt`: `drag-stokes`'s with G = ρ_g U =
12 kg/(m² s) on every 401 line (`make_fixture.py`, `tools/bc401.py`'s `set_field`; field 1 left 0.34/0.66), and
`solfile.tec`, `phase.txt`, `properties.dat` symlinked. Each 1 cm² inlet cell carries F = 0.34/0.66 × (1.2e-3 +
12 × 1e-4) = 1.236364e-3 kg/s, twice `drag-stokes`'s 6.181818e-4, so N_pmin = 5 (3 without G) and R = d_cell/5 =
2 mm: 714 parcels, 25 to 33 per cell; 84 of the pinned points on face-cell edges, 56 of them on edges two cells share.

**Gates** (`check.py`): C1 the fixture (25 lines of 401, field 1 0.34/0.66, G 12); `inj-pack`'s J0–J3 (`run()`, as
`inj-mdotmax`) with the cell flow F; G1 the exit mass flow, the sum over the 25 cells of F = 3.090909e-2 kg/s, to
1e-5 (measured 5.1e-7); G2 the lattice at R = d_cell/5: the closest two parcels of a cell at most R apart, plus the
print (measured 1.999956e-3 m).

**RED.**
- The mass flow without its G term (`self%mdot = krho/Ninj · |ṁ_gas|`): G1 (1.545457e-2 kg/s, half), J3 (each cell
  6.18e-4 kg/s of its 1.236364e-3) and J0 fail; the placement is unchanged (G2 passes).
- N_pmin without its G term (`compute_ds`): the lattice at d_cell/3, 264 parcels; G2 fails (the closest two
  parcels of a cell 3.332349e-3 m apart), the mass flow gates pass.
