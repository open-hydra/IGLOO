#!/usr/bin/env python3
"""Writes the three boundary-file refusal fixtures from drag-stokes's bc.txt:

  bc-legacy/INPUT/bc.txt         every record's code in the legacy single-digit schema (401 -> 4, 400 -> 1);
  bc-empty-payload/INPUT/bc.txt  the first inlet record's payload line left empty;
  krho-negative/INPUT/bc.txt     every inlet record's krho (the payload's first token) negated.

Run from tests/infrastructure/refusals/; the output is committed."""
import os

SRC = os.path.join("..", "..", "standard", "drag-stokes", "INPUT", "bc.txt")
LEGACY = {401: 4, 400: 1}


def legacy(lines):
    out, n = [], 0
    for ln in lines:
        t = ln.split()
        if len(t) == 6 and all(x.isdigit() for x in t) and int(t[5]) in LEGACY:
            ln = ln[:ln.rindex(t[5])] + f"{LEGACY[int(t[5])]:>{len(t[5])}}" + "\n"
            n += 1
        out.append(ln)
    return out, n


def empty_payload(lines):
    out, done = [], False
    for i, ln in enumerate(lines):
        if not done and i > 0 and lines[i - 1].split()[-1:] == ["401"] and len(lines[i - 1].split()) == 6:
            out.append("\n")
            done = True
            continue
        out.append(ln)
    return out, done


def negative_krho(lines):
    out, n = [], 0
    for i, ln in enumerate(lines):
        if i > 0 and lines[i - 1].split()[-1:] == ["401"] and len(lines[i - 1].split()) == 6:
            k = ln.index(ln.split()[0])
            assert k > 0 and float(ln.split()[0]) > 0.0
            ln = ln[:k - 1] + "-" + ln[k:]
            n += 1
        out.append(ln)
    return out, n


def main():
    src = open(SRC).readlines()
    lg, n = legacy(src)
    os.makedirs(os.path.join("bc-legacy", "INPUT"), exist_ok=True)
    open(os.path.join("bc-legacy", "INPUT", "bc.txt"), "w").writelines(lg)
    print(f"bc-legacy/INPUT/bc.txt: {n} record codes rewritten to the legacy schema")
    ep, done = empty_payload(src)
    assert done
    os.makedirs(os.path.join("bc-empty-payload", "INPUT"), exist_ok=True)
    open(os.path.join("bc-empty-payload", "INPUT", "bc.txt"), "w").writelines(ep)
    print("bc-empty-payload/INPUT/bc.txt: the first inlet payload line emptied")
    nk, n = negative_krho(src)
    assert n > 0
    os.makedirs(os.path.join("krho-negative", "INPUT"), exist_ok=True)
    open(os.path.join("krho-negative", "INPUT", "bc.txt"), "w").writelines(nk)
    print(f"krho-negative/INPUT/bc.txt: {n} inlet krho values negated")


if __name__ == "__main__":
    main()
