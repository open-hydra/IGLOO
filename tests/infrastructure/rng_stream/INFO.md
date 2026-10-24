# INFO — rng_stream (per-particle RNG stream contract)

Pins the contract of the per-particle random-number stream
(`IGLOO_Lib_Statistics::rngSeedFor` / `rngNext`, splitmix64) that the samplers reachable
from `!$OMP PARALLEL DO` draw from — TAB's child-size draw, ETAB's azimuth. On the intrinsic
`random_number`, whose state is per-thread, which draw a parcel receives depends on thread
scheduling, and `tab-e2e` is not reproducible run to run (measured as a trajectory MULTISET).

The three properties below are ones **no e2e gate can see**: a wrong-but-deterministic stream
(every parcel sharing one seed) reproduces perfectly and passes every e2e oracle. `repeat-tab`
proves reproducibility at the solver level; this test proves the stream itself.

## Tests (`test_rng_stream.f90`, ctest `test_rng_stream`, labels `unit;infrastructure;rng;t1`)

| id | property | compared | status |
|---|---|---|---|
| RS1 | distinctness | seeds unique over a realistic (famID, ID) grid — 8 families × 4000 IDs, 0 collisions. `ID` is unique only within a group, so `famID` must be mixed in or parcel 1 of every family shares a stream | PASS |
| RS2 | reproducibility | re-seeding replays the identical sequence — what makes a sweep, a thread count and an MPI rank count invisible in the output | PASS |
| RS3 | range + uniformity | draws strictly in (0,1) (`RosinRammler` takes `log(x)`); moment check against the uniform distribution | PASS |

## Notes

- Oracle is analytic (count collisions, replay, moments); no reference data.
- `repeat-tab` is `tab-e2e`'s case (`INPUT/` symlinked) in the repeatability suite.
