#!/usr/bin/env python3
"""Rewrite a dispersed-phase BC file (IGLOO's INPUT/<prefix>bc.txt) to the 401 record of ATLAS's mixture rule over
the dispersed phases: field 1 krho/(1 - krho_tot) and G as a 10th field (the writers of the fixtures call it after
writing the former record, field 1 krho; a one-off conversion of a committed fixture is the same call).

usage: python3 -B bc401.py [--check] [--families N] <bc.txt> [<bc.txt> ...]
       from bc401 import convert; convert(path, families=N, quiet=True)   (a generator, after writing)

Records as IGLOO reads them: a header of 6 integers (block i j k face code), then one payload line for 101, 201 and
401-403, a counts line and n1 + n2 donor lines for 102, nothing for the other codes. The table may be repeated once per
family (copies); per cell, over its copies in file order: krhoTot = the sum of field 1 of its 401 records (a one-copy
file: its value added once per family, IGLOO's rule; the families from --families or the sibling <prefix>phase.txt,
the sum of the group counts of its material lines, or 1 without that file), G = the sum of field 1 of its 402/403 records (likewise). Each 401
payload becomes k' = krho/(1 - krhoTot) in place of krho (17 significant digits: the double read back is the one
computed), the law name and ds when the line lacks them (Dirac, 0: IGLOO's defaults), and G last. Every other line is
kept as it is; a file whose 401 lines have their 10th field already is left alone (a second call changes nothing).
A cell whose copies carry different codes is reported (only there G differs from 0). With --check
nothing is written. Exit 1 on a file that cannot be read as records, or a krhoTot of 1 or more.
"""
import os
import re
import sys


def is_header(tokens):
    return len(tokens) == 6 and all(re.fullmatch(r'-?\d+', t) for t in tokens)


def walk(lines):
    """[(header line index, cell key, code, payload line indices)], or raise."""
    recs, i = [], 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        t = lines[i].split()
        if not is_header(t):
            raise ValueError(f'line {i + 1}: not a header of 6 integers: {lines[i]!r}')
        code = int(t[5])
        if code in (101, 201, 401, 402, 403):
            pay = [i + 1]
        elif code == 102:
            n1, n2 = (int(x) for x in lines[i + 1].replace(',', ' ').split()[:2])
            pay = list(range(i + 1, i + 2 + n1 + n2))
        else:
            pay = []
        recs.append((i, tuple(int(x) for x in t[:5]), code, pay))
        i = (pay[-1] if pay else i) + 1
    return recs


def families_of(path):
    d, b = os.path.split(path)
    phase = os.path.join(d, b[:-len('bc.txt')] + 'phase.txt')
    if not os.path.exists(phase):
        return 1                                  # a generator writing bc.txt before phase.txt: one family
    lines = [l.split() for l in open(phase).read().splitlines()[1:] if l.strip() and not l.lstrip().startswith(('#', '!', ';'))]
    return sum(int(l[1]) for l in lines if len(l) > 1 and re.fullmatch(r'\d+', l[1]))


def set_field(line, k, text):
    """line with its k-th blank-separated field (1-based) replaced by text, the spacing kept (a generator editing
    a payload written under the rule)."""
    parts = re.split(r'(\s+)', line)
    n = 0
    for i, part in enumerate(parts):
        if part and not part.isspace():
            n += 1
            if n == k:
                parts[i] = text
                return ''.join(parts)
    raise ValueError(f'a line of {n} fields has no field {k}: {line!r}')


def convert(path, check=False, families=None, quiet=False):
    lines = open(path).read().split('\n')
    recs = walk(lines)
    keys = {}
    for _, key, code, _ in recs:
        keys.setdefault(key, []).append(code)
    n401 = [len(lines[pay[0]].split()) for _, _, code, pay in recs if code == 401]
    if n401 and min(n401) >= 10:
        return False                                  # written under the rule already (a 10th field, G): as it is
    if n401 and max(n401) >= 10:
        raise ValueError('401 lines with and without the 10th field G')
    ncopies = len(recs) // len(keys) if keys else 1
    if any(len(v) != ncopies for v in keys.values()):
        raise ValueError('cells with different copy counts')
    mixed = [k for k, v in keys.items() if len(set(v)) > 1]
    reps = 1
    if ncopies == 1 and any(c in (401, 402, 403) for v in keys.values() for c in v):
        reps = families if families else families_of(path)
    ktot, gtot = {}, {}
    for _, key, code, pay in recs:                    # file order: the order IGLOO sums in
        if code in (401, 402, 403):
            v = float(lines[pay[0]].replace(',', ' ').split()[0])
            acc = ktot if code == 401 else gtot
            for _ in range(reps):
                acc[key] = acc.get(key, 0.0) + v
    bad = [k for k, s in ktot.items() if not s < 1.0]
    n401 = 0
    for _, key, code, pay in recs:
        if code != 401:
            continue
        n401 += 1
        line = lines[pay[0]]
        m = re.match(r'(\s*)(\S+)(.*)$', line)
        lead, first, rest = m.groups()
        k = float(first.rstrip(','))
        kk = k / (1.0 - ktot[key])
        ntok = len(line.split())
        if ntok < 7:
            raise ValueError(f'line {pay[0] + 1}: a 401 payload of {ntok} tokens (7 at least)')
        extra = ''
        if ntok < 8:
            extra += '   Dirac'
        if ntok < 9:
            extra += f'   {0.0:.16E}'
        lines[pay[0]] = f'{lead}{kk:.16E}{rest.rstrip()}{extra}   {gtot.get(key, 0.0):.16E}'
    if not check and n401:
        with open(path, 'w') as fh:
            fh.write('\n'.join(lines))
    if not quiet:
        print(f'{path}: {len(recs)} records, {ncopies} cop(y/ies), {n401} of 401 converted'
              f'{" (" + str(reps) + " families on one copy)" if reps > 1 else ""}; '
              f'max krhoTot {max(ktot.values(), default=0.0):.6g}, cells with G > 0: {sum(1 for v in gtot.values() if v)}'
              f'{", MIXED codes in " + str(len(mixed)) + " cell(s)" if mixed else ""}'
              f'{", krhoTot >= 1 in " + str(len(bad)) + " cell(s)" if bad else ""}')
    return bool(bad)


def main(argv):
    check = '--check' in argv
    fam = None
    if '--families' in argv:
        fam = int(argv[argv.index('--families') + 1])
    files = [a for i, a in enumerate(argv) if not a.startswith('--') and not (i and argv[i - 1] == '--families')]
    rc = 0
    for f in files:
        try:
            rc |= convert(f, check, fam)
        except (ValueError, OSError, IndexError) as e:
            print(f'{f}: CANNOT CONVERT: {e}')
            rc = 1
    return rc


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
