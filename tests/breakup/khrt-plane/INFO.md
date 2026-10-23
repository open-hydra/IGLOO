# INFO — breakup/khrt-plane (plane-crossing records and the birth row of every breakup child; e2e)

## What it verifies
[khrt-e2e](../khrt-e2e/INFO.md)'s run (`INPUT/` symlinked) with `plane-x = 0.03875 0.075 0.11125 0.15`: mid-cell,
a cell face, mid-cell and the outlet face. 25 water drops enter on the 25 inlet cells into a uniform 200 m/s gas
along x with V = W = 0 everywhere, so every parcel keeps its inlet (Y, Z) line and leaves through the outlet; the
25 parents shed 241 KH children. `OUTPUT/planes-A.dat` holds one record per crossing: `plane X Y Z U V W T d m mdot
npdot f ID parent t dir` (see the user guide's *Plane crossings*).

| id | gate | tolerance | measured |
|---|---|---|---|
| P0 | every record names a listed plane and sits at its position | exact (F12.6) | 792 of 792 |
| P1 | every plane: Σ `mdot` of the forward records − Σ of the backward ones = the injected flow 25 × 0.34/0.66 × 1.2 × 200 × 1e-4 = 17/55 kg/s (a shed moves flow from the parent to the child, an RT shatter keeps it) | 1e-8 relative (each ES16.8 value within 5e-9 of itself) | 1.25e-9, 6.96e-10, 6.53e-10, 1.59e-10 over 119, 180, 227, 266 forward records, 0 backward |
| P2 | each injected parcel crosses each plane once, forward, with parent 0; each child each plane downstream of its birth row once, forward, and none upstream; an ID's records in increasing `t` | exact | 1330 of 1330 |
| P3 | a child's records name one parent: another parcel of the run, on the child's (Y, Z) line | exact | 241 of 241 |
| P4 | each record on its ID's path: Y, Z on its line; U, T and d between two consecutive rows of the path (trajectory rows by x, closed by the exit record) that bracket the plane; m = ρπd³/6 and `mdot` = `npdot`·m | the print half-ULPs; m 2e-5 relative | 792 of 792 |
| P5 | outlet plane: each record's `mdot` = its exit record's plus those of its exit-step children (a child shed on its parent's exit step is born past the outlet: no record, its flow in the parent's record, the flow over the step). Premises: every birth row has V = W = 0 on an injected parcel's line; each exit-step child has one candidate parent on its line; a child without an outlet record is born at x ≥ 0.15 − 1e-6 | the E13.6 half-ULPs of the exit records plus the ES16.8 one | premises 482 of 482; 266 of 266 records |
| P6 | every child's first trajectory row is its birth state on its parent's path: within 2e-6 of a parent row in X, its U and T that row's; otherwise strictly between two consecutive parent rows, U, V, W, T, Y, Z inside their values; its d below the parent's at the parent's last row before the birth; m = ρπd³/6 | 1e-6 (two F12.6 half-ULPs) + 6e-6 × the larger secant slope beside the row; m 2e-5 relative | 241 of 241; 2 births within 2e-6 of a parent row |

The d leg of P6 bounds the child by the parent's last row before the birth: a KH shed needs `dStable` below the
drop's d, and the drop's d only falls along its path. The parent's first row after the birth is no bound: 7 of the
241 births lie between two parent rows across which an RT shatter cuts the parent's d by more than 1.5×, and in 3
the parent ends below that child. A child handed its parent's d would fail P1 as well: its flow `npdot`·m (the shed's
number rate times the mass at that d) would be wrong at every plane downstream.

## RED (each measured by one mutation of the code)
| mutation | what fails |
|---|---|
| `plane-x` not read (an unknown `[IGLOO-General]` key is ignored) | no `planes-A.dat` |
| no crossing check at the segment finalize (a face crossing ends in a short step that is never latched) | 346 records instead of 792; P1 rel 1.00 at x = 0.075 and 0.15 (no record), P2 845 of 1291, the P5 premises 482 of 723, P6 202 of 241 (39 children with no record) |
| records for the injected parcels only | P1 rel 1.02e-1, 1.70e-1, 2.14e-1, 2.49e-1 (the children's share), P2 397 of 1089, P6 0 of 241 |
| the birth row never written | P6 0 of 241 (each child's first row a cadence row on a parent row, U off it by up to 4.9 m/s), P2 1326 of 1330, P4 786 of 792 |
| the model-3 diameter from the auxiliary `d` instead of the stream's m | P4 719 of 792 (m against ρπd³/6) |

## Not covered
No parcel moves toward −x here (the gas drives every parcel forward): the backward branch of the crossing rule is
not exercised. No child is shed on its parent's exit step in this run (P5's exit-step term is empty). Model 3 does
not solidify: the record at a segment end that a phase change shortens is not exercised.
