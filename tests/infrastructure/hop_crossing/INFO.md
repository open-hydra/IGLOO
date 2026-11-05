# INFO — hop_crossing (a reflection's crossing from the step's two end states)

Pins `hopCrossing(pa, va, pb, vb, h, pI, nn, s, xc, vc, found)`: the first outward crossing of the plane through `pI`
(normal `nn`, either sense) on the cubic Hermite of a step's two end states, (pa, va) at s = 0 and (pb, vb) at s = 1,
h the step's duration. The signed distance along the outward normal (the sense of the chord pb − pa) is a cubic in s;
[0, 1] is split at the roots of its derivative, the first piece that goes from negative to non-negative holds the
crossing and a bisection runs there to 4 ULP of s; s is the piece's non-negative side, `xc` and `vc` the Hermite's
position and velocity there (vc · n_out ≥ 0 by construction). `found = .false.` for h ≤ 0, a chord parallel to the
plane, or no outward crossing in (0, 1]. Every path below is of degree ≤ 3 in time: the Hermite is the path itself, so
the oracle is exact.

## Tests (`test_hop_crossing.f90`, ctest `test_hop_crossing`, labels `unit;infrastructure;bc`)

| id | step | compared | measured |
|---|---|---|---|
| U1 | a whole hop inside one step: z = z_a + w_a t − g t²/2, x = u t (z_a = 1e-12 m, w_a = 2e-4 m/s, u = 0.1 m/s), h = 1.5 hops | crossing time s h against the root (w_a + √(w_a² + 2 g z_a))/g within 16 ULP of h; velocity (u, 0, w_a − g t_c) within 1e-12 relative; v_c · n_out > 0; the point on the plane within 1e-18 m | 9.1e-21 s (1.3 ULP of h); 4.2e-17 m/s; 3.6e-24 m |
| U2 | three crossings in one step, g(s) = K (s − 0.1)(s − 0.3)(s − 0.9) (out, in, out) | s = 0.1 within 8 ULP | 0 |
| U3 | the path inside the domain; h = 0; a chord parallel to the plane | `found = .false.` | — |
| U4 | a straight step at constant velocity | s = the chord's fraction within 8 ULP, the point = the chord's intersection with the plane and the velocity = the step's within 64 eps | 1 ULP; 1.0e-24 m; 0 |

Tolerances: the bisection's 4-ULP bracket is a stopping rule; the error adds the cubic's rounding near the root and, in
U1, the product s h.

## RED

- The chord's root (the distance taken linear between the two end points, the crossing a straight chord gives): U1
  fails on all four assertions (the crossing 4.08e-5 s early, near the step's start; the velocity off by 4.0e-4 m/s,
  the start velocity, pointing into the domain; the point 3.0e-12 m off the plane) and U2 (s = 0.3); U4 passes, where
  the chord is exact.
- A plain bisection on [0, 1] (no split at the cubic's turning points): U2, s = 0.9, the last crossing.
- The pieces searched from the last: U2, s = 0.9.
