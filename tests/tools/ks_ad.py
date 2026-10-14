#!/usr/bin/env python3
"""Goodness of fit of a sample to a fully specified CDF: Kolmogorov-Smirnov and Anderson-Darling.

    ks_d(sample, cdf)                 D_n = max_i max(i/n - F(x_(i)), F(x_(i)) - (i-1)/n)
    ks_d_weighted(sample, w, cdf)     the same on the weighted empirical CDF (cumulative w / sum w)
    ad_a2(sample, cdf)                A^2 = -n - (1/n) sum_i (2i-1) [ln F(x_(i)) + ln(1 - F(x_(n+1-i)))]
    ks_pvalue(d, n)                   asymptotic Kolmogorov P(sqrt(n) D_n > sqrt(n) d)
    ad_pvalue(a2)                     P(A^2 > a2), case 0, from Marsaglia & Marsaglia's ADinf

`cdf` is the prescribed law F (no parameter estimated from the sample), so both statistics follow
their case-0 distributions. F is clamped to [1e-15, 1 - 1e-15] before the logarithms: an inverse-CDF
sampler maps u = 0 onto the lower end of the support, where F = 0.

Upper points used by the gates (alpha = 0.05 / 0.01 / 0.001):
    Kolmogorov K_alpha (D_crit = K_alpha / sqrt(n), n >= 35):  1.358 / 1.628 / 1.949
    Anderson-Darling, case 0 (asymptotic, any n >= 10):         2.492 / 3.878 / 5.969
The Kolmogorov points solve 1 - 2 sum_k (-1)^(k-1) exp(-2 k^2 K^2) = 1 - alpha; the Anderson-Darling
points are those of the exact limiting law (Anderson & Darling's series; 3.878 is the 1 % point, tables
printing 3.857 are wrong), which ADinf reproduces to 3e-3. Running this file prints both sets.

References: Stephens, J. Am. Stat. Assoc. 69 (1974) 730; Marsaglia & Marsaglia, J. Stat. Softw. 9(2) (2004).
"""
import math

K_CRIT = {0.05: 1.35810, 0.01: 1.62762, 0.001: 1.94947}
A2_CRIT = {0.05: 2.492, 0.01: 3.878, 0.001: 5.969}
F_CLAMP = 1.0e-15


def _sorted_cdf(sample, cdf):
    return [cdf(x) for x in sorted(sample)]


def ks_d(sample, cdf):
    F = _sorted_cdf(sample, cdf)
    n = len(F)
    return max(max((i + 1) / n - f, f - i / n) for i, f in enumerate(F))


def ks_d_weighted(sample, weights, cdf):
    pairs = sorted(zip(sample, weights))
    tot = sum(w for _, w in pairs)
    d, cum = 0.0, 0.0
    for x, w in pairs:
        f = cdf(x)
        lo = cum / tot
        cum += w
        d = max(d, cum / tot - f, f - lo)
    return d


def ad_a2(sample, cdf):
    F = [min(max(f, F_CLAMP), 1.0 - F_CLAMP) for f in _sorted_cdf(sample, cdf)]
    n = len(F)
    s = 0.0
    for i in range(n):
        s += (2 * i + 1) * (math.log(F[i]) + math.log(1.0 - F[n - 1 - i]))
    return -n - s / n


def kolmogorov_cdf(x, kmax=200):
    if x <= 0.0:
        return 0.0
    s = sum((-1) ** (k - 1) * math.exp(-2.0 * k * k * x * x) for k in range(1, kmax + 1))
    return 1.0 - 2.0 * s


def ks_pvalue(d, n):
    return 1.0 - kolmogorov_cdf(math.sqrt(n) * d)


def adinf(z):
    """Marsaglia & Marsaglia (2004): the limiting CDF P(A^2 <= z) of the case-0 statistic."""
    if z <= 0.0:
        return 0.0
    if z < 2.0:
        return (math.exp(-1.2337141 / z) / math.sqrt(z)
                * (2.00012 + (.247105 - (.0649821 - (.0347962 - (.011672 - .00168691 * z) * z) * z) * z) * z))
    return math.exp(-math.exp(1.0776 - (2.30695 - (.43424 - (.082433 - (.008056 - .0003146 * z) * z) * z) * z) * z))


def ad_pvalue(a2):
    return 1.0 - adinf(a2)


def _bisect(f, lo, hi):
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if (f(mid) > 0) == (f(lo) > 0):
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


if __name__ == "__main__":
    for a in sorted(K_CRIT, reverse=True):
        k = _bisect(lambda x: 1.0 - kolmogorov_cdf(x) - a, 0.5, 3.0)
        z = _bisect(lambda x: 1.0 - adinf(x) - a, 0.5, 15.0)
        print(f"alpha = {a:<6}  Kolmogorov K = {k:.5f} (used {K_CRIT[a]})   "
              f"Anderson-Darling z = {z:.4f} (used {A2_CRIT[a]})")
