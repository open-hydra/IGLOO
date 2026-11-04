#!/usr/bin/env python3
"""Known answers of the per-parcel stream for test_rng_stream's RS4, computed independently of IGLOO: splitmix64 on
Python's unbounded integers reduced modulo 2**64 (the published constants), the top 53 bits as the uniform, the
Marsaglia polar normal and the chi-square sum as IGLOO_Lib_Statistics draws them, the pin-time seed
mix64(rng_seed + famID*GOLDEN + ID*IDMIX). Prints the constants as signed 64-bit integers and as reals."""
import math

M64 = (1 << 64) - 1
GOLDEN, IDMIX = 0x9E3779B97F4A7C15, 0xC2B2AE3D27D4EB4F
MIX_A, MIX_B = 0xBF58476D1CE4E5B9, 0x94D049BB133111EB
RNG_SEED = 42
NDRAW, NNORM, NCHI, KCHI = 1000, 1000, 300, 3


def mix64(z):
    z = ((z ^ (z >> 30)) * MIX_A) & M64
    z = ((z ^ (z >> 27)) * MIX_B) & M64
    return z ^ (z >> 31)


def seed_for(fam, pid):
    return mix64((RNG_SEED + fam * GOLDEN + pid * IDMIX) & M64)


class Stream:
    def __init__(self, state):
        self.s = state & M64

    def uniform(self):
        while True:
            self.s = (self.s + GOLDEN) & M64
            u = (mix64(self.s) >> 11) * 2.0 ** -53
            if u > 0.0:
                return u

    def normal(self):
        while True:
            u1 = 2.0 * self.uniform() - 1.0
            u2 = 2.0 * self.uniform() - 1.0
            s = u1 * u1 + u2 * u2
            if 0.0 < s < 1.0:
                return u1 * math.sqrt(-2.0 * math.log(s) / s)

    def chisq(self, k):
        x = 0.0
        for _ in range(k):
            z = self.normal()
            x = x + z * z
        return x


def signed(x):
    return x - (1 << 64) if x >> 63 else x


def main():
    s1, s2 = seed_for(1, 1), seed_for(4, 200000)
    print(f"rngSeedFor(1, 1)      = {signed(s1)}")
    print(f"rngSeedFor(4, 200000) = {signed(s2)}")
    st = Stream(s1)
    for _ in range(NDRAW):
        u = st.uniform()
    print(f"from rngSeedFor(1, 1): state after {NDRAW} draws {signed(st.s)}, last draw {u!r}")
    for name, start in (("rngSeedFor(1, 1)", s1), ("-1", M64)):
        st = Stream(start)
        tot = sum(st.normal() for _ in range(NNORM))
        print(f"from {name}: state after {NNORM} normals {signed(st.s)}, their sum {tot!r}")
        st = Stream(start)
        tot = sum(st.chisq(KCHI) for _ in range(NCHI))
        print(f"from {name}: state after {NCHI} chi-squares({KCHI}) {signed(st.s)}, their sum {tot!r}")


if __name__ == "__main__":
    main()
