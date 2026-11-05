#!/usr/bin/env python3
"""History census of the tests' text: the test pages and the test sources, scripts and inputs.

The pages are tests/**/INFO.md, tests/VERIFICATION_MATRIX.md and the comment lines of tests/CMakeLists.txt. The
sources are every other tracked .py, .f90, .F90, .sh, .cmake, .ini and .md file under tests/, except symlinks, the
data files under a directory whose name starts with INPUT (their .md files are read) and this script, whose patterns
name the markers. A source is read for the text a reader meets, not for its code:
  .py        comments, docstrings and string literals (tokenize); an f-string's replacement fields are code and are
             blanked by a brace-depth scan that keeps {{ and }} as text; from Python 3.12 an f-string arrives as
             FSTRING_START ... FSTRING_END with its fields as ordinary tokens, so its source text is sliced from the
             first to the matching last and blanked the same way: every interpreter reads the same text
  .f90 .F90  comments and string literals, quote-aware (a doubled quote stays inside its string), a string continued
             by & read across its lines; !$ sentinels, !DIR$/!DEC$ directives and cpp lines are code
  .sh        comments (not the #! line) and quoted strings
  .cmake     comments
  .ini       comment lines and inline comments, as tools/lint_inputs.py reads them (;, # or !)
  .md        every line

The pages and the sources describe the tests, their oracles and their measured results; ledger ids, plan tags,
dates, commit hashes and chronology belong to git and to the working notes. The census lists every line that
carries one of those markers (the five drop classes below) and exits 1 while one is left.

The census is necessary, not sufficient: it sees the listed patterns only. A text can carry history in other
words ("since", "until", "earlier", a past tense), so a file at 0 still has to be read.

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
A source is matched against these rules as itself and as its directory's INFO.md, so the check ids a case's page
admits are admitted in that case's checker too; the printout counts the spans only this mapping admits. Rules of
the sources only:
  body-source checks B0-B5   in tests/tools/body_source_check.py, which defines them
  KH constants B0, B1        the [IGLOO-Models] constants in formulas and lists of the KH tests (reitz-khrt,
                             khrt-*, breakup-per-material)
  DOI                        a DOI in any source: its hash-like or ledger-like pieces are no commit and no id
  RT-path quadrants Q1-Q4    the check labels of tests/breakup/reitz-khrt/test_khrt_interaction.f90
  print-format placeholder   d.dddddddd, an ES16.8 field written in digit placeholders: its d-run is no commit
  working-notes exclude      a quoted string that is exactly /plan-bucket in tests/bench/monolith/sync-src-monolith.sh:
                             the rsync filter that keeps the untracked working notes off the cluster
A drop match is exempt only when an allow rule of its class covers its whole span on that line. A
ledger-like token is never exempt right after "bug", "bugs" or "ledger", right before "fix", "fixes" or
"fixed", or anywhere on a line that names the ledger ("ledger", "BUGS").

Usage: history_census.py [--scope pages|sources|all] [--ref REV] [--raw] [-v] [--plants FILE]
  --scope SCOPE  the pages, the sources or all of them (default: all)
  --ref REV      read the files at a git revision (default: the working tree, tracked files only)
  --raw          no allowlist (the plain pattern census)
  -v             print every line that keeps a drop marker
  --plants FILE  refute the census instead: FILE holds one plant per line, "<class>TAB<file>TAB<text>";
                 each text is checked as a line of the page it names, or as text of the source it names passed
                 through that file's extraction ("\\n" in the text is a line break, "\\\\" a backslash), and must
                 be caught by <class> ("none": must not be caught at all); exits 1 if any plant fails
The files are listed by git: outside a git checkout (no .git at the top of the tree, as in an export) the census
prints "[SKIP] not a git checkout" and exits 77, CTest's skip code of the history-census test.
"""
import argparse
import collections
import fnmatch
import io
import os
import re
import subprocess
import sys
import tokenize

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MATRIX = "tests/VERIFICATION_MATRIX.md"
CMAKE = "tests/CMakeLists.txt"
SELF = "tests/tools/history_census.py"
KINDS = ("py", "f90", "F90", "sh", "cmake", "ini", "md")

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
    ("plane-back checks B0-B4", {"ledger"}, r"B[0-4]", ("tests/infrastructure/plane-back/INFO.md",),
     ("plane-back",), None, "tests/infrastructure/plane-back/check.py"),
    ("beckstead-air checks B0-B4", {"ledger"}, r"B[0-4]", ("tests/combustion/beckstead-air/INFO.md",),
     ("beckstead-air",), None, "tests/combustion/beckstead-air/check.py"),
    ("burnout-threshold checks B0-B2", {"ledger"}, r"B[0-2]", ("tests/evaporation/burnout-threshold/INFO.md",),
     ("burnout-threshold",), None, "tests/evaporation/burnout-threshold/check.py"),
    ("output checks O1-O2", {"ledger"}, r"O[12]",
     ("tests/standard/out-off/INFO.md", "tests/standard/vie-plait-dense/INFO.md"), ("out-off", "vie-plait-dense"),
     None, "tests/standard/out-off/check.py, tests/standard/vie-plait-dense/check.py"),
    ("erosion-law checks O1-O4", {"ledger"}, r"O[1-4]", ("tests/erosion/laws/INFO.md",), ("test_erosion_laws",),
     None, "tests/erosion/laws/test_erosion_laws.f90"),
    ("kept check-id families", set(),
     r"\b(?:XD[1-5]|XE[1-3]|XH1|PE[1-4]|LK[1-6]|TC[1-7]|CB[1-4]|E[1-9]|IP[1-5]|DC[1-3]|ET[1-4][a-d]?|KH[1-3]|"
     r"PR\d{1,2}|SG\d{1,2})\b", ALL, ALL, None, "kept: check ids inside a family"),
    ("paper model names", set(), r"\bMHB98\b|\bM[1-8]\b", ALL, ALL, None, "kept: paper model names (MHB98 M1-M8)"),
]
ALLOW = [(n, c, re.compile(r), p, w, re.compile(l) if l else None, s) for n, c, r, p, w, l, s in ALLOW]
ROW = re.compile(r"^\|\s*\*\*([^*|]+)\*\*\s*\|")

# allow rules of the sources only: (name, drop classes it exempts, regex, file globs, source)
SOURCE_ALLOW = [
    ("body-source checks B0-B5 (the checker)", {"ledger"}, r"B[0-5]'?", ("tests/tools/body_source_check.py",),
     "tests/tools/body_source_check.py"),
    ("KH constants B0, B1 (the KH tests)", {"ledger"}, r"(?<![A-Za-z0-9_-])B[01](?![A-Za-z0-9_])",
     ("tests/breakup/reitz-khrt/*", "tests/breakup/khrt-*/*", "tests/breakup/breakup-per-material/*"),
     "src/lib/config/Register_IGLOO.f90: [IGLOO-Models] B0, B1 (the Reitz KH constants)"),
    ("DOI", {"hash", "ledger"}, r"(?:DOI:?|doi:?|doi\.org/)\s*10\.\d{4,9}/[^\s|)]+", ALL,
     "kept: a citation identifier"),
    ("RT-path quadrants Q1-Q4", {"plantag"}, r"\bQ[1-4]\b", ("tests/breakup/reitz-khrt/test_khrt_interaction.f90",),
     "tests/breakup/reitz-khrt/test_khrt_interaction.f90"),
    ("print-format placeholder d.dddddddd", {"hash"}, r"(?<![A-Za-z0-9.])d\.d{7,10}(?![A-Za-z0-9])", ALL,
     "kept: a print format in digit placeholders (an ES16.8 field)"),
    ("working-notes exclude", {"plantag"}, r"^/plan-bucket$", ("tests/bench/monolith/sync-src-monolith.sh",),
     "tests/bench/monolith/sync-src-monolith.sh: the rsync filter keeping the untracked working notes off the cluster"),
]
SOURCE_ALLOW = [(n, c, re.compile(r), g, s) for n, c, r, g, s in SOURCE_ALLOW]
MAPPING = "directory mapping"


def git(*args):
    return subprocess.run(["git", "-C", ROOT] + list(args), stdout=subprocess.PIPE,
                          universal_newlines=True, check=True).stdout


def tracked(ref):
    """{path: mode} of every tracked file under tests/, at a revision or in the index."""
    rows = git("ls-tree", "-r", "-z", ref, "tests") if ref else git("ls-files", "-s", "-z", "tests")
    out = {}
    for row in rows.split("\0"):
        if row:
            meta, path = row.split("\t", 1)
            out[path] = meta.split()[0]
    return out


def pages(ref):
    return sorted(f for f in tracked(ref) if os.path.basename(f) == "INFO.md" or f in (MATRIX, CMAKE))


def is_page(path):
    return os.path.basename(path) == "INFO.md" or path in (MATRIX, CMAKE)


def ext_of(path):
    base = os.path.basename(path)
    return base.rsplit(".", 1)[1] if "." in base else ""


def sources(ref):
    """The sources in scope, sorted, and the count of each skipped kind."""
    out, skipped = [], collections.Counter()
    for path, mode in tracked(ref).items():
        ext = ext_of(path)
        if ext not in KINDS:
            continue
        if is_page(path):
            skipped["pages"] += 1
        elif mode == "120000":
            skipped["symlinks"] += 1
        elif mode == "160000":
            skipped["submodules"] += 1
        elif ext != "md" and any(d.startswith("INPUT") for d in path.split("/")[:-1]):
            skipped["INPUT data"] += 1
        elif path == SELF:
            skipped["the census itself"] += 1
        else:
            out.append(path)
    return sorted(out), skipped


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


# ---- the sources: what a reader meets in each kind of file -------------------------------------------------------

_FSTART = tuple(t for t in (getattr(tokenize, "FSTRING_START", None), getattr(tokenize, "TSTRING_START", None))
                if t is not None)
_FEND = tuple(t for t in (getattr(tokenize, "FSTRING_END", None), getattr(tokenize, "TSTRING_END", None))
              if t is not None)


def blank_fields(s):
    """An f-string's text with each replacement field {...} replaced by {} ({{ and }} are text)."""
    out, i, n = [], 0, len(s)
    while i < n:
        if s[i] == "{":
            if i + 1 < n and s[i + 1] == "{":
                out.append("{{")
                i += 2
                continue
            depth, j = 1, i + 1
            while j < n and depth:
                depth += {"{": 1, "}": -1}.get(s[j], 0)
                j += 1
            out.append("{}")
            i = j
            continue
        out.append(s[i])
        i += 1
    return "".join(out)


def _split(start_line, s):
    return [(start_line + k, part) for k, part in enumerate(s.split("\n"))]


def py_pieces(src):
    """(line, text) of every comment, docstring and string literal, and whether tokenize failed (every line read)."""
    lines = io.StringIO(src).readlines()
    out = []
    try:
        toks = list(tokenize.generate_tokens(io.StringIO(src).readline))
    except (tokenize.TokenError, IndentationError, SyntaxError):
        return list(enumerate(src.split("\n"), 1)), True
    i = 0
    while i < len(toks):
        t = toks[i]
        if t.type == tokenize.COMMENT:
            out.append((t.start[0], t.string))
        elif t.type == tokenize.STRING:
            s = t.string
            if "f" in re.match(r"[A-Za-z]*", s).group(0).lower():
                s = blank_fields(s)
            out += _split(t.start[0], s)
        elif t.type in _FSTART:
            depth, j = 1, i + 1
            while j < len(toks) and depth:
                depth += 1 if toks[j].type in _FSTART else (-1 if toks[j].type in _FEND else 0)
                j += 1
            (r1, c1), (r2, c2) = t.start, toks[j - 1].end
            s = lines[r1 - 1][c1:] + "".join(lines[r1:r2 - 1]) + lines[r2 - 1][:c2] if r2 > r1 else lines[r1 - 1][c1:c2]
            out += _split(r1, blank_fields(s))
            i = j
            continue
        i += 1
    return out, False


def f90_pieces(src):
    out, q, buf = [], None, ""
    for n, line in enumerate(src.split("\n"), 1):
        i, end = 0, len(line)
        if q is None:
            s = line.lstrip()
            if s.lower().startswith(("!$", "!dir$", "!dec$")) or s.startswith("#"):
                continue                             # OpenMP sentinels, compiler directives, cpp lines
        else:                                        # a string continued from the line above
            j = end - len(line.lstrip())
            i = j + 1 if j < end and line[j] == "&" else 0
            buf = ""
        while i < end:
            c = line[i]
            if q:
                if c == q:
                    if i + 1 < end and line[i + 1] == q:  # a doubled quote inside the string
                        buf += c
                        i += 2
                        continue
                    out.append((n, buf))
                    q, buf = None, ""
                else:
                    buf += c
            elif c in "'\"":
                q, buf = c, ""
            elif c == "!":
                out.append((n, line[i:]))
                break
            i += 1
        if q:
            if buf.rstrip().endswith("&"):           # continued on the next line
                out.append((n, buf.rstrip()[:-1]))
                buf = ""
            else:
                out.append((n, buf))
                q, buf = None, ""
    return out


def sh_pieces(src):
    out = []
    for n, line in enumerate(src.split("\n"), 1):
        if n == 1 and line.startswith("#!"):
            continue
        i, end, q, buf = 0, len(line), None, ""
        while i < end:
            c = line[i]
            if q == "'":
                if c == "'":
                    out.append((n, buf))
                    q, buf = None, ""
                else:
                    buf += c
            elif q == '"':
                if c == "\\" and i + 1 < end:
                    buf += line[i:i + 2]
                    i += 2
                    continue
                if c == '"':
                    out.append((n, buf))
                    q, buf = None, ""
                else:
                    buf += c
            else:
                if c == "\\":
                    i += 2
                    continue
                if c in "'\"":
                    q, buf = c, ""
                elif c == "#" and (i == 0 or line[i - 1] in " \t;|&("):
                    out.append((n, line[i:]))
                    break
            i += 1
        if q and buf:
            out.append((n, buf))
    return out


def cmake_pieces(src):
    out = []
    for n, line in enumerate(src.split("\n"), 1):
        q = False
        for i, c in enumerate(line):
            if c == '"' and (i == 0 or line[i - 1] != "\\"):
                q = not q
            elif c == "#" and not q:
                out.append((n, line[i:]))
                break
    return out


_INI_INLINE = re.compile(r"[;#!]")


def ini_pieces(src):
    out = []
    for n, line in enumerate(src.split("\n"), 1):
        s = line.strip()
        if not s:
            continue
        if s[0] in ";#!":
            out.append((n, line))
            continue
        head = s.index("]") + 1 if s.startswith("[") and "]" in s else (s.index("=") + 1 if "=" in s else None)
        if head is None:
            out.append((n, line))                    # neither an option nor a section: read whole
            continue
        m = _INI_INLINE.search(s, head)
        if m:
            out.append((n, s[m.start():]))
    return out


def pieces(path, src):
    """(line, text) pieces of a source and whether its .py failed to tokenize."""
    ext = ext_of(path)
    if ext == "py":
        return py_pieces(src)
    reader = {"f90": f90_pieces, "F90": f90_pieces, "sh": sh_pieces, "cmake": cmake_pieces, "ini": ini_pieces}
    if ext in reader:
        return reader[ext](src), False
    return list(enumerate(src.split("\n"), 1)), False


def source_census(path, piece, raw, exempt, rules=None, mapping=True):
    """The drop matches of one piece of the source `path` that no allow rule exempts: [(class, token)]. A page rule
    applies to the source as itself or as its directory's INFO.md; exempt[MAPPING] counts the matches only the
    second admits."""
    rules = SOURCE_ALLOW if rules is None else rules
    page = os.path.join(os.path.dirname(path), "INFO.md")
    if FORCED_LINE.search(piece):
        forced = [(0, len(piece))]
    else:
        forced = [f.span(1) for rx in FORCED for f in rx.finditer(piece)]
    own, mapped = [], []
    if not raw:
        for name, classes, rx, pglobs, _, lrx, _ in ALLOW:
            if lrx is not None and not lrx.search(piece):
                continue
            if applies(pglobs, path):
                own += [(name, classes, a.start(), a.end()) for a in rx.finditer(piece)]
            elif mapping and applies(pglobs, page):
                mapped += [(name, classes, a.start(), a.end()) for a in rx.finditer(piece)]
        for name, classes, rx, globs, _ in rules:
            if applies(globs, path):
                own += [(name, classes, a.start(), a.end()) for a in rx.finditer(piece)]
    out = []
    for cls, rx in DROP.items():
        for d in rx.finditer(piece):
            s, e = d.span()
            cover = via = None
            if not (cls == "ledger" and any(a <= s and e <= b for a, b in forced)):
                cover = next((n for n, c, a, b in own if cls in c and a <= s and e <= b), None)
                if cover is None:
                    cover = via = next((n for n, c, a, b in mapped if cls in c and a <= s and e <= b), None)
            if cover:
                exempt[cover] += 1
                if via:
                    exempt[MAPPING] += 1
            else:
                out.append((cls, d.group(0)))
    return out


def census_sources(ref, raw, rules=None, mapping=True):
    """hits [(path, line, class, token, piece)], exempted spans per rule, the files in scope, the skipped counts and
    the .py files that failed to tokenize (read whole)."""
    files, skipped = sources(ref)
    hits, exempt, unparsable = [], collections.Counter(), []
    for path in files:
        pcs, bad = pieces(path, text(path, ref))
        if bad:
            unparsable.append(path)
        for no, piece in pcs:
            for cls, tok in source_census(path, piece, raw, exempt, rules, mapping):
                hits.append((path, no, cls, tok, piece.strip()))
    return hits, exempt, files, skipped, unparsable


def source_kind(path):
    return ("tools/" if path.startswith("tests/tools/") else "") + "*." + ext_of(path)


# ---- plants and the printout ----------------------------------------------------------------------------------

def unescape(s):
    out, i = [], 0
    while i < len(s):
        if s[i] == "\\" and i + 1 < len(s) and s[i + 1] in "n\\":
            out.append("\n" if s[i + 1] == "n" else "\\")
            i += 2
        else:
            out.append(s[i])
            i += 1
    return "".join(out)


def plants(fname, ref):
    """Each plant is a line of the page it names or text of the source it names; it must be caught by its class
    ('none': by no class)."""
    bad = n = 0
    known = set(pages(ref))
    srcs = set(sources(ref)[0])
    for raw in open(fname, encoding="utf-8").read().split("\n"):
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        want, path, planted = raw.split("\t", 2)
        planted = unescape(planted)
        if path in known:
            got = set(c for line in planted.split("\n") for c, _ in line_census(path, line, False,
                                                                                 collections.Counter()))
        elif path in srcs:
            got = set(c for _, piece in pieces(path, planted)[0]
                      for c, _ in source_census(path, piece, False, collections.Counter()))
        else:
            print("ERROR   %s: no such page or source in the census" % path)
            bad += 1
            continue
        got = sorted(got)
        n += 1
        ok = (not got) if want == "none" else (want in got)
        bad += not ok
        print("%-7s %-8s %s | %s | caught by: %s" % ("ok" if ok else "FAILED", want, path,
                                                    planted.replace("\n", "\\n")[:110], ",".join(got) or "nothing"))
    print("plants: %d, failed: %d" % (n, bad))
    return 1 if bad else 0


def report_pages(a):
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
    return len(hits)


def report_sources(a):
    hits, exempt, files, skipped, unparsable = census_sources(a.ref, a.raw)
    nk = collections.Counter(source_kind(p) for p in files)
    print("history census of the sources of %s (%s): %d files (%s); not read: %s"
          % ("the working tree" if not a.ref else a.ref, "raw" if a.raw else "allowlist on", len(files),
             ", ".join("%s %d" % kv for kv in sorted(nk.items())),
             ", ".join("%s %d" % kv for kv in sorted(skipped.items())) or "none"))
    per = collections.Counter((source_kind(p), c) for p, _, c, _, _ in hits)
    hit_lines = set((p, n) for p, n, _, _, _ in hits)
    print("%-12s %7s %7s %7s %7s %7s %7s %7s" % ("kind", *DROP, "lines", "files"))
    for k in sorted(nk):
        kl = [x for x in hit_lines if source_kind(x[0]) == k]
        print("%-12s %7d %7d %7d %7d %7d %7d %7d" % ((k,) + tuple(per[(k, c)] for c in DROP)
                                                       + (len(kl), len(set(p for p, _ in kl)))))
    print("%-12s %7d %7d %7d %7d %7d %7d %7d" % (("total",) + tuple(sum(per[(k, c)] for k in nk) for c in DROP)
                                                   + (len(hit_lines), len(set(p for p, _ in hit_lines)))))
    top = collections.Counter(p.split("/")[1] for p, _ in hit_lines)
    print("lines with a drop marker by directory: %s" % (", ".join("%s %d" % kv for kv in
                                                                   sorted(top.items(), key=lambda kv: (-kv[1], kv[0])))
                                                        or "none"))
    if not a.raw:
        print("allowlist in the sources (rule [classes exempted] files <- source: spans exempted):")
        print("  the pages' rules above, each source matched as itself and as its directory's INFO.md: %s"
              % (", ".join("%s %d" % (n, exempt[n]) for n, _, _, _, _, _, _ in ALLOW if exempt[n]) or "none"))
        print("  %s (admitted only as the directory's INFO.md): %d" % (MAPPING, exempt[MAPPING]))
        for name, classes, rx, globs, src in SOURCE_ALLOW:
            print("  %s [%s] /%s/ files %s <- %s: %d" % (name, ",".join(sorted(classes)), rx.pattern,
                                                         ",".join(globs), src, exempt[name]))
    if unparsable:
        print("not tokenized, every line read: %s" % ", ".join(unparsable))
    if a.v:
        for p, n, c, t, l in hits:
            print("%s:%d: [%s] %s | %s" % (p, n, c, t, l[:160]))
    return len(hits)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--scope", choices=("pages", "sources", "all"), default="all")
    ap.add_argument("--ref", default=None)
    ap.add_argument("--raw", action="store_true")
    ap.add_argument("-v", action="store_true")
    ap.add_argument("--plants", default=None)
    a = ap.parse_args()
    if not os.path.exists(os.path.join(ROOT, ".git")):  # an export: no index to list the files
        print("[SKIP] not a git checkout")
        return 77
    if a.plants:
        return plants(a.plants, a.ref)
    left = {}
    if a.scope in ("pages", "all"):
        left["pages"] = report_pages(a)
    if a.scope in ("sources", "all"):
        left["sources"] = report_sources(a)
    total = sum(left.values())
    if a.scope == "all":
        print("drop markers left: %d (pages %d, sources %d)" % (total, left["pages"], left["sources"]))
    else:
        print("drop markers left: %d" % total)
    return 1 if total else 0


if __name__ == "__main__":
    sys.exit(main())
