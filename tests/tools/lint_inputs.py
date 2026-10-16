#!/usr/bin/env python3
"""lint-test-inputs: the INI files of the test cases hold value lists FiNeR counts right, in ASCII only.

Rules:
  blanks  a value (the text after the first `=` of an option line, an inline `;`, `#` or `!` comment
          removed) separates two tokens by anything but one blank. FiNeR counts the values of a list as
          its blanks plus one, so `vp = -1.0   -3.3` reads as four values (count_values = 4; a get into
          four reals returns no error and leaves the last two unset; FiNeR v2.2.0, measured), and a tab
          is no separator at all.
  ascii   a byte above 127 anywhere in the file (FiNeR v2.2.0 still finds the keys after such a line;
          ASCII keeps the inputs portable).
Every *.ini under the tests directory is read once (a symlink through its target), skipping the
directories runs and checks write: OUTPUT/, ref-*/, tec_n*/, _run-*/.

Usage: lint_inputs.py <tests dir>
"""
import os
import re
import sys

SKIP_DIRS = ("OUTPUT", "__pycache__")
SKIP_PREFIXES = ("ref-", "tec_n", "_run-")
COMMENT_LINE = ";#!"
INLINE = re.compile(r"[;#!]")
GAP = re.compile(r"(?<=\S)[ \t]+(?=\S)")


def ini_files(root):
    seen, out = set(), []
    for dp, dn, fn in os.walk(root):
        dn[:] = sorted(d for d in dn if d not in SKIP_DIRS and not d.startswith(SKIP_PREFIXES))
        for f in sorted(fn):
            if f.endswith(".ini"):
                p = os.path.join(dp, f)
                real = os.path.realpath(p)
                if real not in seen and os.path.isfile(real):
                    seen.add(real)
                    out.append(p)
    return out


def lint(path):
    bad = []
    with open(path, "rb") as fh:
        raw = fh.read()
    for n, line in enumerate(raw.split(b"\n"), 1):
        if any(b > 127 for b in line):
            bad.append((n, "ascii", line.decode("utf-8", "replace").strip()))
        text = line.decode("ascii", "replace").strip()
        if not text or text[0] in COMMENT_LINE or text.startswith("[") or "=" not in text:
            continue
        value = INLINE.split(text.split("=", 1)[1], 1)[0].strip()
        if any(g.group(0) != " " for g in GAP.finditer(value)):
            bad.append((n, "blanks", text))
    return bad


def main():
    if len(sys.argv) != 2:
        print("[FAIL] usage: lint_inputs.py <tests dir>")
        return 1
    root = os.path.abspath(sys.argv[1])
    files = ini_files(root)
    if not files:
        print("[FAIL] no .ini file under %s" % root)
        return 1
    hits = []
    for p in files:
        real, rroot = os.path.realpath(p), os.path.realpath(root)
        shown = os.path.relpath(real, rroot) if real.startswith(rroot + os.sep) else os.path.relpath(p, root)
        for n, rule, text in lint(p):
            hits.append("%s:%d: %s: %s" % (shown, n, rule, text[:100]))
    print("lint-test-inputs: %d INI files, %d findings" % (len(files), len(hits)))
    for h in hits:
        print("  " + h)
    if hits:
        print("[FAIL] fix the lines above (one blank between list values, ASCII only)")
        return 1
    print("[PASS] every value list single-blank separated, every file ASCII")
    return 0


if __name__ == "__main__":
    sys.exit(main())
