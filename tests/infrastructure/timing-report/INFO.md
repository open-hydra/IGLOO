# INFO — infrastructure/timing-report (`[IGLOO-General] timing = on`, e2e)

## Reference
Not a physics case. The physics is `standard/vie-plait`'s and is gated there. This case gates the
`[IGLOO-timing]` report that `timing = on` prints after each output write.

## Fixture
`INPUT/` symlinks `standard/vie-plait/INPUT` whole (gas-order 2, 120 x 50 x 5 cells, one gas
species, six DB streams, one material). `input.ini` is vie-plait's own file with `timing = on`.
`tests/mpi/timing-report` symlinks this directory's `INPUT/`, `input.ini` and `check.py` and runs
the same oracle at four ranks.

## What `check.py` asserts
Times vary run to run and are never compared to a number; everything else in the report is.

| # | assertion | oracle |
|---|---|---|
| 1 | one header, then exactly one line per phase in the fixed order, each time finite and >= 0, inside its rank min/max | the phase list |
| 2 | one line per thread of rank 0, one line per rank | the header's thread count; the solver's `MPI ranks = n` banner (1 when absent) |
| 3 | rank r integrated the parcels `ip = r+1, r+1+n, ...`; rank 0's threads sum to rank 0's count | the strided ownership and the solver's `number of particles` line |
| 4 | the rank max/mean of the parcel-loop time is >= 1, and 1 at one rank | its definition |
| 5 | allocated gas, geometry and accumulator values per rank | the closed form of the mesh (below) |

Closed form (gas-order 2, a 3D block of Ni x Nj x Nk cells, N^3 = Ni Nj Nk, (N+k)^3 =
(Ni+k)(Nj+k)(Nk+k), nsc gas species, nm materials, source output only):

- gas fields: (nsc + 9) (N+2)^3 (density per species, velocity, temperature, five property arrays on
  the ghost-ringed dual mesh);
- accumulators: (nm + 4) (N+1)^3 (mass per material, momentum, energy on the dual mesh);
- geometry: geometry block 4 (N+1)^3 + 16 N^3 (nodes, dual weights; centres, edge lengths, minimum
  edge, octant volumes, cell volumes) and gas block 3 (N+2)^3 + (N+1)^3 (dual nodes, dual volumes);
  the source block holds no geometry.

## RED
A binary without the report prints no `[IGLOO-timing]` line (FAIL on the first check). With the
per-rank gather removed, ranks 1..3 report no parcels and the four-rank variant fails assertion 3.
