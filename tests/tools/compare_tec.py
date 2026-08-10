#!/usr/bin/env python3
"""Tolerance comparison of two IGLOO ASCII `.tec` field files.

    compare_tec.py <a.tec> <b.tec> [tol]        # CLI, exit 0 == within tolerance
    from compare_tec import compare_tec         # (ok, worst, ndiff, nval) tuple

THE METRIC IS SCALE-RELATIVE, and that choice is load-bearing rather than stylistic:

    err = max|a_i - b_i| / max|a_i|

i.e. each difference is normalised by the FIELD's magnitude, not by the local value. A per-value
relative error `|a-b|/|a|` is ill-conditioned on these fields. Measured 2026-08-10 on `vie-plait`:
one cell holds 6.03e-12 against a field scale of 6.0, and an absolute difference of 7.3e-24 there
reads as a 1.2e-12 "relative" error — apparently at tolerance, actually 24 orders below the field.
Field-wide the worst absolute difference on that comparison was 1.06e-22, i.e. 1.8e-23 of scale.

Both wrong moves are available with the per-value metric: calling 1.2e-12 a regression (a false
positive), or loosening the tolerance to swallow it — which would then hide a REAL defect of the
same magnitude somewhere the local value is O(1). Normalising by field scale removes the choice.

Why any tolerance at all: these fields are accumulated with `!$OMP ATOMIC UPDATE`, and floating
point addition is not associative, so the value depends on thread interleaving — and, once MPI is
live, on how the rank count partitions the sum. Both are the same FP-reassociation class. Measured
floors on this host: run-to-run at fixed thread count <= 4e-15; across rank counts, bit-identical on
coupled-body / khrt-e2e / db-2daxi and 1.8e-23 on vie-plait. DEFAULT_TOL sits decades above those and
decades below the 5.2e-11 that the F7 accumulator-shape bug produced.

`tests/tools/check_twosweep.py` carries its own equivalent comparator. It is deliberately NOT
refactored to import this one: its PASS stdout is byte-stable on purpose (see the note at its
`cmp_tec`, and BUGS.md O2), and seven green gates depend on that. Duplicating ~20 lines is the
cheaper side of that trade.
"""
import sys

DEFAULT_TOL = 1.0e-12


def tec_values(path):
    """One float per line; non-numeric lines (headers) become None placeholders so the two files
    stay index-aligned even where they carry text."""
    vals = []
    with open(path, errors='replace') as fh:
        for line in fh:
            try:
                vals.append(float(line.strip()))
            except ValueError:
                vals.append(None)
    return vals


def compare_tec(a, b, tol=DEFAULT_TOL):
    """(ok, worst_scale_relative, n_values_differing, n_values). Length mismatch => worst = inf."""
    va, vb = tec_values(a), tec_values(b)
    if len(va) != len(vb):
        return False, float('inf'), abs(len(va) - len(vb)), len(va)
    scale = max((abs(x) for x in va if x is not None), default=1.0) or 1.0
    worst, ndiff = 0.0, 0
    for x, y in zip(va, vb):
        if x is None or y is None or x == y:
            continue
        ndiff += 1
        worst = max(worst, abs(x - y) / scale)
    return worst <= tol, worst, ndiff, len(va)


def main():
    if len(sys.argv) < 3:
        print("[FAIL] usage: compare_tec.py <a.tec> <b.tec> [tol]")
        return 2
    tol = float(sys.argv[3]) if len(sys.argv) > 3 else DEFAULT_TOL
    ok, worst, ndiff, nval = compare_tec(sys.argv[1], sys.argv[2], tol)
    if ok:
        print(f"[PASS] within tolerance (scale-relative <= {tol:.0e}), {nval} values")
        return 0
    print(f"[FAIL] scale-relative error {worst:.3e} exceeds {tol:.0e} "
          f"({ndiff} of {nval} values differ)")
    return 1


if __name__ == "__main__":
    sys.exit(main())
