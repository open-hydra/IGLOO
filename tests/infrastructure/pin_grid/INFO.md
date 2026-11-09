# INFO — pin_grid (the bucket grid of the 3D injection placement)

On a 3D inlet face the advancing front of `pin_particles_bc_ds` (`src/lib/initialization.f90`) refuses a candidate
that lies within R/2 of a point the call has already placed, R the spacing of the cell that spawned it. The placed
points are filed in a uniform bucket grid (`pinGrid`) over the box of the call's injecting face cells, its edge h the
smallest spacing of a packed cell (doubled while the grid would hold more than `NBMAX` = 2²⁴ buckets).
`pinGridFar` walks the buckets of the interval [p − R/2 − m, p + R/2 + m] per direction, the margin
m = max(10⁻⁶·R/2, 64·ε·L) (L the largest coordinate magnitude of the box) above the rounding of that interval, and
evaluates on every point it finds the predicate `norm2(p − q) < R/2`; a coordinate outside the box maps to the edge
bucket. The grid therefore returns the boolean of a scan over every placed point, and the same candidates are
accepted in the same order.

## Tests (`test_pin_grid.f90`, ctest `test_pin_grid`, labels `unit;infrastructure;injection`)

The reference is the test's own linear scan of the stored points (`scanFar`); every query's two booleans are compared,
and a leg with no query near a stored point fails as vacuous. P3 runs first, and its P3b before any point is filed.

| id | compared | status |
|---|---|---|
| P1 | 20000 points uniform in a 0.05 × 0.05 × 2·10⁻⁸ m slab (a planar face, h = 3.3·10⁻⁴ m) and in a 0.05 m cube (a curved face's cloud, h = 2.0·10⁻³ m), h set so that about half the queries at R = h find a point; 20000 uniform queries per cloud at R = h, 1.7 h and 3.1 h: 0 mismatches (near: 9873, 17245, 19968 in the slab; 9663, 19057, 20000 in the cube) | PASS |
| P2 | near-ties: from 10000 stored points q, queries p = q + s·e (e a random unit vector) at s = (R/2)(1 + k·10⁻¹²), k = −4…4, R = h, both clouds: 0 mismatches of 90000 each | PASS |
| P3 | 3000 points whose coordinates lie on the bucket planes lo + n·h or 10⁻¹⁸ m (one ulp at the far side) outside the box, queried at 0, R/4 and (R/2)(1 ∓ 10⁻¹²) along the six axis directions: 0 mismatches of 72000 | PASS |
| P3b | the bucket of every coordinate of P3 (and of a point 10 h beyond each side of the box) lies in 1…nb, checked before any point is filed | PASS |
| P4 | the cube cloud inside a 1 m box at the requested edge 2·10⁻³ m (1.25·10⁸ buckets): the edge is doubled once, to 4·10⁻³ m (1.5625·10⁷ buckets ≤ `NBMAX`), and P1–P3 repeated on that grid at R = 2·10⁻³ m: 0 mismatches | PASS |

## RED

- **The query range cut to the point's own bucket**: every leg of P1–P4 fails — the grid misses 3091 of the 9873
  near queries of the slab at R = h and 3797 of the 9663 of the cube, 19786 and 27608 of the 90000 near-tie queries,
  11902 of the 72000 edge queries; P3b and the cap assertion stay green.
- **The bucket index unclamped** (`floor((x − lo)/h) + 1` without the clamp to 1…nb): P3b fails, 55403 coordinates
  mapped outside 1…nb; P3, filing and querying through those indices, then reaches outside the grid and the program
  stops with a segmentation fault (exit 174 under ifx) before P3's verdict.
