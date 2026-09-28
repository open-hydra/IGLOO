#!/usr/bin/env python3
"""A23/S4-E robustness gate: heavy KH shedding must stay bounded, mass-exact and crash-free.

This is a STRESS case, not a validation case — there is no oracle and no paper. It runs the
breakup/khrt-e2e fixture with `mShedLim` at 0.01 instead of the 0.03 default, so a parent
sheds after a 1 % accumulated mass loss rather than 3 %. Measured effect: ~700 children
instead of ~240, with the deepest single parent shedding ~66 times, which drives the growable
per-parent shed list through ~7 doublings.

What it is actually protecting:

  * the unbounded shed path cannot be exercised to any depth by khrt-e2e alone (~22 sheds from
    the deepest parent there), so list growth under repeated `move_alloc` is otherwise untested;
  * mass must stay exact no matter how many times a parcel splits — every shed both strips the
    parent and creates the child, so an error in that pairing compounds with shed count and
    shows up here long before it would in the production case;
  * the runaway guard (`maxShed` in Lib_Integration) must NOT fire on a merely-aggressive case.
    If it does, either the guard is set too low or the case has become genuinely divergent —
    both are findings, so the warning is gated as a failure rather than ignored.

The ctest command supplies the rc>=128 signal gate, so a crash or OOM fails independently.
"""
import sys

OUTLOC  = "OUTPUT/outloc-A.dat"
RUN_OUT = "run_out.txt"

MDOT_INJECTED = 3.0909090909e-01   # identical fixture to khrt-e2e (25 cells x 1.23636364e-02)
MDOT_TOL      = 1.0e-4

MIN_CHILDREN  = 500   # observed 701; a floor, so a collapse of the shed path fails loudly
MIN_MAX_SHED  = 20    # observed 66; proves the per-parent list really did grow deep


def main():
    try:
        rows = [l.split() for l in open(OUTLOC)]
    except FileNotFoundError:
        print(f"[FAIL] {OUTLOC} not found -- did the solver run?")
        return 1

    tot, n = 0.0, 0
    for c in rows:
        if len(c) == 9 and c[0][0] in "0123456789-":
            try:
                tot += float(c[6]); n += 1
            except ValueError:
                continue

    nkids, maxshed, capped = 0, 0, 0
    for line in open(RUN_OUT):
        if "number of children" in line:
            nkids += int(line.split("=")[-1])
        elif "sheds from one parcel" in line:
            maxshed = max(maxshed, int(line.split("max")[-1].split()[0]))
        elif "shed cap" in line:
            capped += 1

    rel = abs(tot - MDOT_INJECTED) / MDOT_INJECTED
    ok_mass = rel <= MDOT_TOL
    ok_kids = nkids >= MIN_CHILDREN
    ok_deep = maxshed >= MIN_MAX_SHED
    ok_cap  = capped == 0

    print("KHRT heavy-shedding robustness gate (mShedLim=0.01):")
    print(f"  parcels exited     : {n} ({nkids} shed children, need >= {MIN_CHILDREN}) "
          f"{'PASS' if ok_kids else 'FAIL'}")
    print(f"  deepest single parcel: {maxshed} sheds (need >= {MIN_MAX_SHED}) "
          f"{'PASS' if ok_deep else 'FAIL'}")
    print(f"  mass conservation  : mdot_out={tot:.8e} vs mdot_in={MDOT_INJECTED:.8e}, "
          f"rel={rel:.2e} (tol {MDOT_TOL:.0e}) {'PASS' if ok_mass else 'FAIL'}")
    print(f"  runaway guard      : fired {capped} time(s), expected 0 "
          f"{'PASS' if ok_cap else 'FAIL'}")

    if not ok_kids:
        print("[FAIL] the shed path collapsed -- far fewer children than this case produces.")
    if not ok_deep:
        print("[FAIL] no parcel shed deeply; the growable list is not being exercised, so "
              "this case is no longer a stress test (a one-shed-per-parent cap may be back).")
    if not ok_mass:
        print("[FAIL] mass is not conserved under heavy shedding -- the per-shed "
              "strip/create pairing is wrong; the error compounds with shed count.")
    if not ok_cap:
        print("[FAIL] the runaway guard fired on a merely-aggressive case: either maxShed is "
              "set too low, or shedding has genuinely diverged.")

    if ok_mass and ok_kids and ok_deep and ok_cap:
        print("\n[PASS] heavy shedding stays bounded, mass-exact and crash-free.")
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
