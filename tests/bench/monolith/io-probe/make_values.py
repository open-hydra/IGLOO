#!/usr/bin/env python3
"""The random sets of the I/O probe (igloo-io.slurm, step read): N doubles written raw (values.bin, little-endian)
and as text in the box's layout, one value per line "  {v:.<p>E}" for p = 15, 16 and 17 (v15E.txt, v16E.txt,
v17E.txt). Of every 1000 values: 900 with |v| log-uniform over [1e-25, 1e25], 97 with random bits over every finite
double (subnormals included), one zero of each sign and one power of ten 1e-22 .. 1e22; signs random, the seed
fixed.

    python3 make_values.py <dir> <N>
"""
import os
import random
import struct
import sys


def values(n, seed=12345):
    rnd = random.Random(seed)
    pack, unpack = struct.Struct("<Q").pack, struct.Struct("<d").unpack
    out = []
    for i in range(n):
        r = i % 1000
        if r == 997:
            v = 0.0
        elif r == 998:
            v = -0.0
        elif r == 999:
            v = float("1e%d" % rnd.randint(-22, 22))
        elif r >= 900:
            while True:
                v = unpack(pack(rnd.getrandbits(64)))[0]
                if v == v and abs(v) != float("inf"):
                    break
        else:
            v = 10.0 ** rnd.uniform(-25.0, 25.0)
            if rnd.getrandbits(1):
                v = -v
        out.append(v)
    return out


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    d, n = sys.argv[1], int(sys.argv[2])
    os.makedirs(d, exist_ok=True)
    v = values(n)
    with open(os.path.join(d, "values.bin"), "wb") as f:
        f.write(struct.pack("<%dd" % n, *v))
    for p in (15, 16, 17):
        fmt = "  {:.%dE}\n" % p
        with open(os.path.join(d, "v%dE.txt" % p), "w") as f:
            for k in range(0, n, 100000):
                f.write("".join(fmt.format(x) for x in v[k:k + 100000]))
    print("make_values: %d values in %s" % (n, d))


if __name__ == "__main__":
    main()
