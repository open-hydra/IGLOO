#!/usr/bin/env python3
"""History census of the test pages: tests/**/INFO.md, tests/VERIFICATION_MATRIX.md and the comment
lines of tests/CMakeLists.txt.

The pages describe the tests, their oracles and their measured results; ledger ids, plan tags, dates,
commit hashes and chronology belong to git and to the working notes. The census lists every line that
carries one of those markers (the five drop classes below) and exits 1 while one is left.

The census is necessary, not sufficient: it sees the listed patterns only. A page can carry history in
other words ("since", "until", "earlier", a past tense), so a page at 0 still has to be read.

Plan tags include the pointers to the working notes (FINDINGS, its D-<topic>-<n> ids) and the plan's
open-question ids (OQ<n>, Q<n>, Q-<letter>). Chronology includes the RED-on-the-base forms, pre-fix,
"was fixed", "used to", "first validation", "had shipped", "this case found" and "promoted". "used to"
also matches the present passive ("is used to"): reword such a hit or justify it.

The allowlist keeps by rule the tokens the drop patterns also match. What each rule admits:
  numeral            any number; it also admits a hash made of digits only or of digits around one 'e'
  citation locator   App./Appendix/Eq./Fig./Table/§ followed by an id
  input key B0, B1   the Reitz KH constants of [IGLOO-Models] (Register_IGLOO.f90), written as a key
                     (in backticks or followed by '=') or on a line that names the KH model
  species O2         O2 inside a list of two or more formulas joined by ',', '/' or '+'
  check-id families  the A/B/O-letter ids of one checker, on the pages and matrix rows of the cases that
                     checker defines, on any line there: a ledger id of the same letters on such a line
                     is admitted too unless a never-exempt rule below catches it
A drop match is exempt only when an allow rule of its class covers its whole span on that line. A
ledger-like token is never exempt right after "bug", "bugs" or "ledger", right before "fix", "fixes" or
"fixed", or anywhere on a line that names the ledger ("ledger", "BUGS").

Usage: history_census.py [--ref REV] [--raw] [-v] [--plants FILE]
  --ref REV      read the pages at a git revision (default: the working tree, tracked pages only)
  --raw          no allowlist (the plain pattern census)
  -v             print every line that keeps a drop marker
  --plants FILE  refute the census instead: FILE holds one plant per line, "<class>TAB<page>TAB<text>";
                 each text is checked as a line of the page it names and must be caught by <class>
                 ("none": must not be caught at all); exits 1 if any plant fails
"""
import argparse
import collections
import fnmatch
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MATRIX = "tests/VERIFICATION_MATRIX.md"
CMAKE = "tests/CMakeLists.txt"

# drop classes
DROP = collections.OrderedDict([
    ("ledger", re.compile(r"(?<![A-Za-z0-9_-])(?:A\d{1,2}[a-h]?|O\d{1,2}|B\d{1,2})(?![A-Za-z0-9_])")),
    ("plantag", re.compile(r"\b(?:B-VAL-?\d*|E-VAL-?\d*[a-z]?|W-plan|W plan|Tier-[12]|tier-[12]|phase M\d|S4-D\d|"
                           r"plan-bucket|BUGS\.md|FINDINGS(?:\.md)?|D-[A-Z]+-\d+|OQ\d+|Q-[A-Z]|Q\d+|"
                           r"psat-scalar-cleanup|SESSION|skill)\b")),
    ("date", re.compile(r"\b20\d\d-\d\d-\d\d\b")),
    ("hash", re.compile(r"(?<![A-Za-z0-9])[0-9a-f]{7,10}(?![A-Za-z0-9])")),
    ("chrono", re.compile(r"\b(?:[Ff]ixed in|[Ff]ound by|[Pp]reviously|used to|[Bb]efore the fix|"
                          r"RED on the (?:old|previous|pre|base)|was removed|no longer|first cut|"
                          r"(?:on|with|from) the (?:previous|old|earlier|original) "
                          r"(?:binary|build|tree|fork|commit|code|version)|"
                          r"[Pp]re-fix|[Pp]ost-fix|(?:was|were|been) fixed|(?:was|were) found|"
                          r"[Ff]irst validation|[Ff]irst coverage|for the first time|"
                          r"(?:had|was|were) shipped|shipped (?:unexercised|untested)|"
                          r"(?:[Tt]his|[Tt]he) (?:case|test|gate|check) (?:found|FOUND)|"
                          r"[Pp]romot(?:ed|ion)|[Hh]istorically)\b")),
])

# never exempt: a ledger-like token after bug/bugs/ledger, before fix/fixes/fixed, or on a ledger line
_ID = r"\**[ABO]\d{1,2}[a-h]?\**"
FORCED = [
    re.compile(r"\b(?:[Bb]ugs?|[Ll]edger)\s+(" + _ID + r"(?:\s*(?:/|,|\+|and)\s*" + _ID + r")*)"),
    re.compile(r"(?<![A-Za-z0-9_-])([ABO]\d{1,2}[a-h]?)(?=[ -](?:fix|fixes|fixed)\b)"),
]
FORCED_LINE = re.compile(r"\b(?:[Ll]edger|BUGS)\b")

# allow rules: (name, drop classes it exempts, regex, page globs, matrix-row globs, line regex, source);
# page globs "*" = every page, matrix-row globs "*" = every matrix line; line regex None = any line, else the
# rule applies only on a line it matches; an empty class set exempts nothing (a kept-list family that no
# drop pattern matches, listed so the printout shows the whole kept list)
ALL = ("*",)
KH_LINE = r"KH|Reitz|[Ss]tripping|breakupOde|dStable|\[IGLOO-Models\]"
_SPECIES = r"(?:N2|H2O|CO2|CO|H2|Ar|O2)"
ALLOW = [
    ("numeral", {"hash"}, r"[-+−]?\d*\.?\d+(?:[eE][-+−]?\d+)?", ALL, ALL, None,
     "kept: constants, tolerances, measured results"),
    ("citation locator", {"ledger"}, r"(?:App\.|Appendix|[Ee]q\.|Fig\.|Table|§)\s*[A-Z]?\d+[a-z]?", ALL, ALL, None,
     "kept: citation tags and their locators"),
    ("input key B0, B1 (as a key)", {"ledger"}, r"`B[01]`|(?<![A-Za-z0-9_-])B[01](?=\s*=)", ALL, ALL, None,
     "src/lib/config/Register_IGLOO.f90: [IGLOO-Models] B0, B1 (the Reitz KH constants)"),
    ("input key B0, B1 (on a KH line)", {"ledger"}, r"(?<![A-Za-z0-9_-])B[01](?![A-Za-z0-9_])", ALL, ALL, KH_LINE,
     "src/lib/config/Register_IGLOO.f90: [IGLOO-Models] B0, B1 (the Reitz KH constants)"),
    ("species O2 (in a species list)", {"ledger"}, _SPECIES + r"(?:\s*[,/+]\s*" + _SPECIES + r")+", ALL, ALL, None,
     "kept: a chemical formula"),
    ("drag unit checks A1-A4", {"ledger"}, r"A[1-4]′?", ("tests/standard/drag/INFO.md",), ("test_drag",), None,
     "tests/standard/drag/test_drag.f90"),
    ("body-source checks B0-B5", {"ledger"}, r"B[0-5]'?", ("tests/infrastructure/body-*/INFO.md",),
     ("body-*", "check-order"), None, "tests/tools/body_source_check.py"),
    ("sym-plane-bounce checks B1-B4", {"ledger"}, r"B[1-4]", ("tests/infrastructure/sym-plane-bounce*/INFO.md",),
     ("sym-plane-bounce*",), None, "tests/infrastructure/sym-plane-bounce/check.py"),
    ("evap-breakup-box checks B1-B3", {"ledger"}, r"B[1-3]", ("tests/evaporation/evap-breakup-box/INFO.md",),
     ("evap-breakup-box",), None, "tests/evaporation/evap-breakup-box/check.py"),
    ("temperature unit checks B1-B3c", {"ledger"}, r"B[1-3][a-c]?", ("tests/standard/temperature/INFO.md",),
     ("test_temperature",), None, "tests/standard/temperature/test_temperature.f90"),
    ("solid-box-2mat checks B1-B2", {"ledger"}, r"B[12]", ("tests/solidification/solid-box-2mat/INFO.md",),
     ("solid-box-2mat",), None, "tests/solidification/solid-box-2mat/check.py"),
    ("body-force-z checks B0-B2", {"ledger"}, r"B[0-2]", ("tests/standard/body-force-z/INFO.md",),
     ("body-force-z",), None, "tests/standard/body-force-z/check.py"),
    ("output checks O1-O2", {"ledger"}, r"O[12]",
     ("tests/standard/out-off/INFO.md", "tests/standard/vie-plait-dense/INFO.md"), ("out-off", "vie-plait-dense"),
     None, "tests/standard/out-off/check.py, tests/standard/vie-plait-dense/check.py"),
    ("kept check-id families", set(),
     r"\b(?:XD[1-5]|XE[1-3]|XH1|PE[1-4]|LK[1-6]|TC[1-7]|CB[1-4]|E[1-9]|IP[1-5]|DC[1-3]|ET[1-4][a-d]?|KH[1-3]|"
     r"PR\d{1,2}|SG\d{1,2})\b", ALL, ALL, None, "kept: check ids inside a family"),
    ("paper model names", set(), r"\bMHB98\b|\bM[1-8]\b", ALL, ALL, None, "kept: paper model names (MHB98 M1-M8)"),
]
ALLOW = [(n, c, re.compile(r), p, w, re.compile(l) if l else None, s) for n, c, r, p, w, l, s in ALLOW]
ROW = re.compile(r"^\|\s*\*\*([^*|]+)\*\*\s*\|")


def git(*args):
    return subprocess.run(["git", "-C", ROOT] + list(args), stdout=subprocess.PIPE,
                          universal_newlines=True, check=True).stdout


def pages(ref):
    names = git("ls-tree", "-r", "--name-only", ref, "tests") if ref else git("ls-files", "tests")
    out = [f for f in names.split("\n")
           if f and (os.path.basename(f) == "INFO.md" or f in (MATRIX, CMAKE))]
    return sorted(set(out))   # an unmerged path is listed once per stage


def text(path, ref):
    if ref:
        return git("show", "%s:%s" % (ref, path))
    with open(os.path.join(ROOT, path), encoding="utf-8") as fh:
        return fh.read()


def applies(globs, key):
    return key is not None and any(g == "*" or fnmatch.fnmatchcase(key, g) for g in globs)


def kind(path):
    return "INFO.md" if path.endswith("INFO.md") else ("MATRIX" if path == MATRIX else "CMake")


def line_census(path, line, raw, exempt):
    """The drop matches of one line of `path` that no allow rule exempts: [(class, token)]."""
    if path == CMAKE and not line.lstrip().startswith("#"):
        return []
    m = ROW.match(line) if path == MATRIX else None
    row = m.group(1).strip() if m else None
    if FORCED_LINE.search(line):
        forced = [(0, len(line))]
    else:
        forced = [f.span(1) for rx in FORCED for f in rx.finditer(line)]
    spans = []
    if not raw:
        for name, classes, rx, pglobs, rglobs, lrx, _ in ALLOW:
            if path == MATRIX:
                ok = (pglobs == ALL and rglobs == ALL) or applies(rglobs, row)
            else:
                ok = applies(pglobs, path)
            if not ok or (lrx is not None and not lrx.search(line)):
                continue
            spans += [(name, classes, a.start(), a.end()) for a in rx.finditer(line)]
    out = []
    for cls, rx in DROP.items():
        for d in rx.finditer(line):
            s, e = d.span()
            if cls == "ledger" and any(a <= s and e <= b for a, b in forced):
                cover = None
            else:
                cover = next((n for n, c, a, b in spans if cls in c and a <= s and e <= b), None)
            if cover:
                exempt[cover] += 1
            else:
                out.append((cls, d.group(0)))
    return out


def census(ref, raw):
    hits = []                       # (path, lineno, class, token, line)
    exempt = collections.Counter()  # allow rule -> spans exempted
    for path in pages(ref):
        for no, line in enumerate(text(path, ref).split("\n"), 1):
            for cls, tok in line_census(path, line, raw, exempt):
                hits.append((path, no, cls, tok, line.strip()))
    return hits, exempt


def plants(fname, ref):
    """Each plant is a line of the page it names; it must be caught by its class ('none': by no class)."""
    bad = n = 0
    known = set(pages(ref))
    for raw in open(fname, encoding="utf-8").read().split("\n"):
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        want, path, planted = raw.split("\t", 2)
        if path not in known:
            print("ERROR   %s: no such page in the census" % path)
            bad += 1
            continue
        n += 1
        got = sorted(set(c for c, _ in line_census(path, planted, False, collections.Counter())))
        ok = (not got) if want == "none" else (want in got)
        bad += not ok
        print("%-7s %-8s %s | %s | caught by: %s" % ("ok" if ok else "FAILED", want, path, planted[:110],
                                                    ",".join(got) or "nothing"))
    print("plants: %d, failed: %d" % (n, bad))
    return 1 if bad else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--ref", default=None)
    ap.add_argument("--raw", action="store_true")
    ap.add_argument("-v", action="store_true")
    ap.add_argument("--plants", default=None)
    a = ap.parse_args()
    if a.plants:
        return plants(a.plants, a.ref)
    hits, exempt = census(a.ref, a.raw)
    allp = pages(a.ref)
    ninfo = sum(1 for p in allp if p.endswith("INFO.md"))
    print("history census of %s (%s): %d INFO.md pages, %s, the comment lines of %s"
          % ("the working tree" if not a.ref else a.ref, "raw" if a.raw else "allowlist on", ninfo, MATRIX, CMAKE))
    per = collections.Counter((kind(p), c) for p, _, c, _, _ in hits)
    lines = collections.Counter(kind(p) for p, _ in set((p, n) for p, n, _, _, _ in hits))
    nfiles = len(set(p for p, _, _, _, _ in hits if p.endswith("INFO.md")))
    print("%-8s %8s %8s %8s" % ("class", "INFO.md", "MATRIX", "CMake"))
    for c in DROP:
        print("%-8s %8d %8d %8d" % (c, per[("INFO.md", c)], per[("MATRIX", c)], per[("CMake", c)]))
    print("lines with a drop marker: INFO.md %d (%d pages) | VERIFICATION_MATRIX %d | CMakeLists comments %d"
          % (lines["INFO.md"], nfiles, lines["MATRIX"], lines["CMake"]))
    if not a.raw:
        print("allowlist (rule [classes exempted] pages / matrix rows / line <- source: spans exempted):")
        for name, classes, rx, pglobs, rglobs, lrx, src in ALLOW:
            print("  %s [%s] /%s/ pages %s rows %s line %s <- %s: %d"
                  % (name, ",".join(sorted(classes)) or "-", rx.pattern, ",".join(pglobs), ",".join(rglobs),
                     "/%s/" % lrx.pattern if lrx else "any", src, exempt[name]))
        for rx in FORCED:
            print("  never exempt: /%s/" % rx.pattern)
        print("  never exempt: every ledger-like token on a line matching /%s/" % FORCED_LINE.pattern)
    if a.v:
        for p, n, c, t, l in hits:
            print("%s:%d: [%s] %s | %s" % (p, n, c, t, l[:160]))
    print("drop markers left: %d" % len(hits))
    return 1 if hits else 0


if __name__ == "__main__":
    sys.exit(main())
