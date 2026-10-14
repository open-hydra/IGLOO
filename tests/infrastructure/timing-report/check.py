#!/usr/bin/env python3
"""[IGLOO-timing] report gate (timing = on): the report's structure and its counts against oracles
that do not come from the timers.

  1. one header, then exactly one line per phase, in order, each a finite time >= 0;
  2. one line per thread of rank 0 and one line per rank (the rank count from the solver's own
     `MPI ranks = n` banner, 1 when absent);
  3. parcel counts: rank r integrated exactly the generation-0 parcels the strided ownership gives
     it (ip = r+1, r+1+n, ...), the ranks sum to the solver's `number of particles` lines, and
     rank 0's threads sum to rank 0's count;
  4. the rank imbalance max/mean is >= 1 (exactly 1 at one rank);
  5. the allocated sizes equal the closed form of the fixture's mesh (gas-order 2, 3D, source
     output only): gas (nsc+9)(Ni+2)(Nj+2)(Nk+2); accumulators (nm+4)(Ni+1)(Nj+1)(Nk+1);
     geometry = geo block 4(N+1)^3 + 16 N^3, gas block 3(N+2)^3 + (N+1)^3, source block copy
     3(N+1)^3 + 7 N^3 (N^3 = Ni Nj Nk, (N+k)^3 = (Ni+k)(Nj+k)(Nk+k)).

PASS stdout carries no times (they vary run to run).
"""
import math
import re
import sys

PHASES = ['setup:input', 'setup:gas-read', 'setup:blocks', 'setup:gas-import', 'setup:pinning',
          'accumulators', 'ghost-fill', 'solve-other', 'parcel-loop', 'census',
          'reduce', 'finalize', 'smoothing', 'writeout+merge']

fails = []


def fail(msg):
    fails.append(msg)


def main():
    text = open('run_out.txt', errors='replace').read()
    tl = [l.strip() for l in text.splitlines() if l.strip().startswith('[IGLOO-timing]')]
    m = re.search(r'MPI ranks = (\d+)', text)
    nranks = int(m.group(1)) if m else 1
    nparc = sum(int(x) for x in re.findall(r'number of particles =\s+(\d+)', text))
    nchild = sum(int(x) for x in re.findall(r'number of children =\s+(\d+)', text))
    if nparc <= 0:
        print('[FAIL] the solver reported no particles: nothing to account for')
        return 1
    if nchild:
        fail(f'fixture sheds {nchild} children: the strided-ownership oracle covers generation 0 only')
    if not tl:
        print('[FAIL] no [IGLOO-timing] line in run_out.txt (timing = on is in input.ini)')
        return 1

    # 1. header + phases
    hdr = [l for l in tl if l.startswith('[IGLOO-timing] sweep ')]
    if len(hdr) != 1:
        fail(f'{len(hdr)} header lines, expected 1')
        nthr = 0
    else:
        h = re.search(r'sweep (\d+): (\d+) solve\(s\), (\d+) rank\(s\) x (\d+) thread\(s\)', hdr[0])
        if not h:
            fail(f'unreadable header: {hdr[0]}')
            nthr = 0
        else:
            nthr = int(h.group(4))
            if int(h.group(1)) != 0 or int(h.group(2)) != 1:
                fail(f'header says sweep {h.group(1)}, {h.group(2)} solve(s); expected sweep 0, 1 solve')
            if int(h.group(3)) != nranks:
                fail(f'header says {h.group(3)} rank(s), the solver banner {nranks}')
    ph = [re.match(r'\[IGLOO-timing\] phase (\S+)\s+(\S+) s\s+ranks min\s+(\S+) max\s+(\S+)', l)
          for l in tl if l.startswith('[IGLOO-timing] phase ')]
    names = [p.group(1) if p else None for p in ph]
    if names != PHASES:
        fail(f'phase lines {names} differ from the expected {PHASES}')
    for p in ph:
        if not p:
            continue
        t, lo, hi = (float(p.group(k)) for k in (2, 3, 4))
        if not all(math.isfinite(v) and v >= 0 for v in (t, lo, hi)) or not (lo <= t <= hi):
            fail(f'phase {p.group(1)}: time {t} min {lo} max {hi} is not a finite ordered triple >= 0')
    loop = [float(p.group(2)) for p in ph if p and p.group(1) == 'parcel-loop']
    if loop and loop[0] <= 0:
        fail('parcel-loop time is 0 although parcels were integrated')

    # 2-3. threads and ranks
    th = [re.match(r'\[IGLOO-timing\] thread\s+(\d+)\s+parcels\s+(\d+)\s+busy\s+(\S+) s', l)
          for l in tl if l.startswith('[IGLOO-timing] thread ')]
    if [t.group(1) if t else None for t in th] != [str(t) for t in range(nthr)]:
        fail(f'{len(th)} thread lines, expected {nthr} (threads 0..{nthr - 1})')
    rk = [re.match(r'\[IGLOO-timing\] rank\s+(\d+)\s+parcel-loop\s+(\S+) s\s+parcels\s+(\d+)', l)
          for l in tl if l.startswith('[IGLOO-timing] rank ')]
    if [r.group(1) if r else None for r in rk] != [str(r) for r in range(nranks)]:
        fail(f'{len(rk)} rank lines, expected {nranks} (ranks 0..{nranks - 1})')
    else:
        counts = [int(r.group(3)) for r in rk]
        expect = [len(range(r, nparc, nranks)) for r in range(nranks)]
        if counts != expect:
            fail(f'per-rank parcel counts {counts}, strided ownership of {nparc} parcels gives {expect}')
        if th and all(th) and sum(int(t.group(2)) for t in th) != counts[0]:
            fail(f'rank 0 threads integrated {sum(int(t.group(2)) for t in th)} parcels, rank 0 {counts[0]}')

    # 4. imbalance
    im = [re.search(r'max/mean\s+(\S+)', l) for l in tl if l.startswith('[IGLOO-timing] imbalance')]
    if len(im) != 1 or not im[0]:
        fail(f'{len(im)} imbalance lines, expected 1')
    else:
        v = float(im[0].group(1))
        if not (v >= 0.9995 and (nranks > 1 or abs(v - 1) < 1e-3)):
            fail(f'imbalance max/mean {v} (must be >= 1, and 1 at one rank)')

    # 5. sizes against the closed form of the mesh
    #> every rank prints its block sizes: one entry per block index
    blocks = list({int(b[0]): tuple(int(x) for x in b[1:]) for b in
                   re.findall(r'- Block\s+(\d+) size =\s+(\d+)\s+(\d+)\s+(\d+)', text)}.values())
    gasfile = re.search(r'gas-file\s*=\s*(\S+)', open('input.ini').read()).group(1)
    with open(gasfile, errors='replace') as fh:
        head = ''.join(fh.readline() for _ in range(5))
    nsc = len(re.findall(r'(?:Roi|rho)\(', head))
    mats = [l.split() for l in open('INPUT/phase.txt').read().splitlines()[1:] if l.strip()]
    nm = len(mats)
    sz = [re.search(r'gas (\d+)\s+geometry (\d+)\s+accumulators (\d+)', l) for l in tl if l.startswith('[IGLOO-timing] sizes')]
    if len(sz) != 1 or not sz[0]:
        fail(f'{len(sz)} sizes lines, expected 1')
    elif not blocks or any(nk < 2 for _, _, nk in blocks):
        fail(f'closed form written for 3D blocks, the solver reported {blocks}')
    else:
        gas = geo = acc = 0
        for ni, nj, nk in blocks:
            n3 = ni * nj * nk
            p1 = (ni + 1) * (nj + 1) * (nk + 1)
            p2 = (ni + 2) * (nj + 2) * (nk + 2)
            gas += (nsc + 9) * p2
            acc += (nm + 4) * p1
            geo += (4 * p1 + 16 * n3) + (3 * p2 + p1) + (3 * p1 + 7 * n3)
        got = tuple(int(sz[0].group(k)) for k in (1, 2, 3))
        if got != (gas, geo, acc):
            fail(f'sizes gas/geometry/accumulators {got}, closed form {(gas, geo, acc)}')

    if fails:
        for f in fails:
            print(f'  [FAIL] {f}')
        return 1
    print(f'  [PASS] {len(PHASES)} phase lines, {nthr} thread line(s), {nranks} rank line(s); '
          f'{nparc} parcels split as the strided ownership gives; sizes match the closed form')
    return 0


if __name__ == '__main__':
    sys.exit(main())
