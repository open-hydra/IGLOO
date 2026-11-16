#!/usr/bin/env python3
"""burnout-cold: burn-box's stream with T-ign 2000 K and d-burnout 60 um, twice the injected 30 um. In the 600 K gas
the parcels never ignite: a burning model removes a parcel only while it burns, so none is removed and the mass
source of the material is zero in every cell (each segment's mass in equals its mass out). Removed cold at the
first step, the stream would hand the gas all its mass in the inlet cells."""
import sys

SRC = "OUTPUT/source.tec"


def main():
    lines = open(SRC).read().split("\n")
    names = lines[0].split("=", 1)[1].split()
    zone = lines[1]
    dims = [int(zone.split(k + "=")[1].split(",")[0]) for k in ("I", "J", "K")]
    nnode = dims[0] * dims[1] * dims[2]
    ncell = (dims[0] - 1) * (dims[1] - 1) * (dims[2] - 1)
    nums = []
    for line in lines[2:]:
        if not line.strip() or line.lstrip()[0].isalpha():
            continue
        nums += [float(x) for x in line.split()]
    k = names.index("wdot(A)")
    w = nums[3 * nnode + (k - 3) * ncell: 3 * nnode + (k - 2) * ncell]
    total, peak = sum(w), max(abs(x) for x in w)
    print("burnout-cold: cells %d, sum of wdot(A) %.6e kg/s, largest |wdot(A)| %.6e kg/s" % (len(w), total, peak))
    if peak != 0.0:
        print("[FAIL] a cold parcel was removed: the stream lost mass")
        return 1
    print("[PASS] no cold parcel removed: the mass source is zero in every cell")
    return 0


if __name__ == "__main__":
    sys.exit(main())
