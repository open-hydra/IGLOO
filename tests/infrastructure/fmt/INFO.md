# INFO — fmt (the particle files' edit descriptors without the runtime)

`IGLOO_fmt` (`src/lib/fmt.f90`) writes the edit descriptors of the particle files into a character buffer, byte for
byte what the compiler runtime's internal write produces, without calling the runtime: `F12.6`, `E13.6E2`, `I8`,
`ES16.8E2`, `F10.6`, `I3`, `I5`, `I10` and the `I0` of zone headers. `fmt_compile` turns a format such as
`'(7F12.6,2E13.6E2,I8)'` into a record layout (`[r]Fw.d`, `[r]Ew.dEe`, `[r]ESw.dEe`, `[r]Iw` and `I0`; any other item
stops the run); `fmt_row` appends one record of a layout from a row of reals, in which an integer column holds an exact
integer below 2^53; `fmt_int` appends one integer. The routines are `pure` and the module holds no variable: the
number of fields the runtime wrote is returned to the caller.

The fast path rounds the scaled value |x|·10^j, with 10^j an exact power (0 ≤ j ≤ 22), to the nearest integer. The
runtime's internal write formats every field the fast path cannot decide exactly:

- a scaled value within 4 ulp of a half-integer (a tie or a near-tie: the rounded product cannot tell its side);
- a value outside the descriptor's fast range: `Fw.d` from |x| = 10^(w−d−2) − 1/2 on (9999.5 for `F12.6`, 99.5 for
  `F10.6`); `Ew.dEe` and `ESw.dEe` when the significand's scale leaves 0 ≤ j ≤ 22, that is outside
  1e−17 ≤ |x| < 1e6 for `E13.6E2` and 1e−14 ≤ |x| < 1e9 for `ES16.8E2` (zero has a fast path of its own); every
  exponent the fast path writes is within ±22, so a wider exponent reaches the runtime through the range;
- a subnormal, NaN or ±Inf, classified from the bits (a build that flushes subnormals to zero compares a subnormal
  as zero, while the runtime formats its bits);
- an integer wider than its field (the runtime's asterisks).

Signs follow the runtime: both compilers sign a negative value that rounds to zero (`-0.000000`); gfortran signs
−0.0 and the Intel runtime does not (`__GFORTRAN__` selects the rule). An integer column that holds NaN, ±Inf or a
value of 2^53 or more prints asterisks.

## Tests (`test_fmt.f90`, ctest `test_fmt`, labels `unit;infrastructure;io`)

Every record is formatted by `fmt_row` and by the runtime's internal write with the format and the item list of the
production writer of its layout; the bytes must agree. Mismatched fields are counted per descriptor (the surface
temperature column apart), the fixed parts apart from the random ones. The first argument sets the random records
per layout (default 300 000, the registered size); P1 and the constructed parts do not depend on it.

| id | compared | status |
|---|---|---|
| P1 | fixed: the scatter layout `(7F12.6,2E13.6E2,I8)` on 2·10^6 plausible records (positions, velocities, temperature, diameter, mass, ID; every 7th with a negative that rounds to zero, every 11th with an F value k/256 (k odd), every 13th with an E tie between 10^5 and 10^6, every 17th with −0.0, every 19th with a zero mass) and 10^6 decimal half-way points of `F12.6` with their ±1 and ±2 ulp neighbours: 0 of 3 000 000 records differ; 1 151 726 of them have a field the runtime writes | PASS |
| L1–L11 C | fixed: every field runs through its descriptor's constructed values (below); a layout's C part has as many records as its longest list (5 946), field k of record i holding value i + 37(k − 1) of its list, cyclically — the trajectory layouts L1–L6 (plain, `+ t`, `+ f`, `+ t + f`, `+ t + T_s`, `+ T_s`; the plain and `+ f` ones are the scatter's, `+ t` and `+ t + f` the snapshot's), the outloc layouts L7–L10 (`+ t + f`, `+ t`, `+ f`, plain), the planes layout L11 `(I3,7F12.6,2E13.6E2,2ES16.8E2,F10.6,2I10,ES16.8E2,I3)` | PASS |
| L1–L11 R | 300 000 records per layout, every field drawn at random: every decade of the double range, the decades around the field, the fast range, and the constructed classes at random points | PASS |
| L12 | the zone headers' integers, `(I0,I0,I0,I0)`, through `fmt_row` and through `fmt_int`: C (107 values up to 2^53) and R (300 000 records) | PASS |
| S | an `I8` column holding NaN, +Inf, 2^53 or −1e300 prints `********`; `fmt_int` of −2^63 and of 123456789 in `I8` is the runtime's; each layout's record width | PASS |

The constructed values, both signs; the ties, half-way points, powers of ten, the values below them and the edges
with their neighbours up to 8 ulps away:

| descriptor | ±0 | NaN, −NaN, ±Inf | subnormals, smallest normals | rounding to zero | exact ties | decimal half-way points | powers of ten, half a last digit below | fast range's edge, beyond the field | total |
|---|---|---|---|---|---|---|---|---|---|
| `F12.6`, `T_s` | 2 | 4 | 22 | 112 | 1 122 | 1 496 | 884 | 130 | 3 772 |
| `F10.6` | 2 | 4 | 22 | 112 | 748 | 1 190 | 748 | 130 | 2 956 |
| `E13.6E2` | 2 | 4 | 22 | — | 1 088 | 2 550 | 1 904 | 376 | 5 946 |
| `ES16.8E2` | 2 | 4 | 22 | — | 782 | 2 550 | 1 904 | 376 | 5 640 |

Exact ties are odd multiples of 2^−(d+1) for `F` and halves of d-digit integers times powers of ten for `E`/`ES`;
the half-way points are those of the last digit at every decade of the fast range and one beyond; the edges include
exponents of two digits and beyond. The integer descriptors (`I3`, `I5`, `I8`, `I10`) take 0 and 10^k − 1, 10^k,
10^k + 1 at every decade, both signs, to the int32 extremes (66 values, every width's edges among them); `I0` the same
to 2^53 − 1 (107). The constructed values are computed in floating point: 27 092 of the 499 012 values of all the
lists differ by one ulp between ifx and gfortran: without the near-tie window a C part's mismatched records differ by
up to 6 between the two compilers, the mismatched fields per descriptor do not, nor do they with the −0.0 rules
exchanged. The test also prints the records per second of both routes on records of in-range values, one thread.

The full size, by hand: `build/verif/tests/test_fmt 3000000` (3·10^6 random records per layout): 0 of its
39 065 513 records differ under ifx and under gfortran.

## RED

Each mutation changes one line of `src/lib/fmt.f90`; the counts are mismatched fields of the fixed parts (P1 and the
constructed values, the same at any record count), at the registered size, ifx and gfortran alike where both are
given.

- **No near-tie window** (a scaled value rounded by its fraction alone): P1 fails on 355 828 of its 3 000 000
  records, every layout's C part on 1 251 to 1 619 of its 5 946 (gfortran 1 249 to 1 625); `F12.6` 114 929,
  `E13.6E2` 331 672, `ES16.8E2` 768, `F10.6` 784, `T_s` 368 (both compilers); the integer descriptors pass.
- **The two compilers' −0.0 rules exchanged**: P1 fails on its 117 647 records with −0.0, and every −0.0 of the C
  parts: `F12.6` 117 730 (83 of them in C), `E13.6E2` 26, `ES16.8E2` 8, `F10.6` 10, `T_s` 2 (both compilers).
- **The sign written only when the value does not round to zero**: P1 fails on 4 778 records; `F12.6` 12 223,
  `F10.6` 820, `T_s` 164 (ifx); the `E`, `ES` and integer descriptors pass.
- **The integers' sign dropped**: `I8` 18 924, `I3` 1 082, `I5` 4 324, `I10` 4 864, `I0` 214 (ifx); P1 (positive
  IDs) and the real descriptors pass.
