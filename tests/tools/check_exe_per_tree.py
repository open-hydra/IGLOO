#!/usr/bin/env python3
"""exe-per-tree: the tests of a build tree run that tree's own executables, never the shared copy.

<source>/bin/IGLOO and <source>/bin/DocGen are copies each build publishes for users; a test that ran
them could meet a relink of another tree in progress (a launch then finds a missing or half-written
file). Fails on
  (0) a function(igloo_<name> ...) of tests/CMakeLists.txt whose body runs $<TARGET_FILE:IGLOO> and is
      not one of the helpers rule (2) checks;
  (1) a registered test whose command or properties name <source>/bin/IGLOO or <source>/bin/DocGen;
  (2) a test registered through igloo_e2e_case, igloo_e2e_case_mpi, igloo_refusal_case or
      igloo_threads_case whose ENVIRONMENT lacks IGLOO_EXE (the executable its check.py runs for a
      sibling case);
  (3) a tests/**/*.py that builds the path bin/IGLOO anywhere but in the default of an IGLOO_EXE read
      (os.environ.get("IGLOO_EXE", <path>) or os.getenv("IGLOO_EXE", <path>)).
The MPI gates that take the executable as an argument (mpi-consistency, mpi-consistency-two-mat,
mpi-consistency-psd, mpi-determinism) are not helper registrations and need no IGLOO_EXE.

Usage: check_exe_per_tree.py <build>/tests/CTestTestfile.cmake <source dir>
"""
import ast
import os
import re
import subprocess
import sys

HELPERS = ("igloo_e2e_case", "igloo_e2e_case_mpi", "igloo_refusal_case", "igloo_threads_case")
SHARED = ("IGLOO", "DocGen")
SKIP_DIRS = ("OUTPUT", "__pycache__")
SKIP_PREFIXES = ("ref-", "tec_n", "_run-")


def parse_ctest(path):
    """{name: (command text, properties text)} from a generated CTestTestfile.cmake."""
    tests = {}
    with open(path) as fh:
        for line in fh:
            m = re.match(r"add_test\(\[=\[(.+?)\]=\]\s*(.*)\)\s*$", line)
            if m:
                tests.setdefault(m.group(1), ["", ""])[0] = m.group(2)
                continue
            m = re.match(r"set_tests_properties\(\[=\[(.+?)\]=\]\s+PROPERTIES\s*(.*)\)\s*$", line)
            if m:
                tests.setdefault(m.group(1), ["", ""])[1] = m.group(2)
    return tests


def env_of(props):
    m = re.search(r'(?:^|\s)ENVIRONMENT\s+"((?:[^"\\]|\\.)*)"', props)
    return m.group(1) if m else ""


def helper_names(tests, srcdir):
    """Names registered through the helpers: from the generated backtraces, else from tests/CMakeLists.txt."""
    names, have_bt = set(), False
    for name, (_, props) in tests.items():
        if "_BACKTRACE_TRIPLES" in props:
            have_bt = True
            if any(";%s;" % h in props for h in HELPERS):
                names.add(name)
    if have_bt:
        return names, "backtraces"
    text = open(os.path.join(srcdir, "tests", "CMakeLists.txt")).read()
    for h, n in re.findall(r"^\s*(%s)\(\s*([A-Za-z0-9_.+-]+)" % "|".join(HELPERS), text, re.M):
        names.add(n)
    return names & set(tests), "tests/CMakeLists.txt"


def solver_helpers(srcdir):
    """The function(igloo_<name> ...) definitions of tests/CMakeLists.txt whose body runs $<TARGET_FILE:IGLOO>."""
    runners, cur = [], None
    for line in open(os.path.join(srcdir, "tests", "CMakeLists.txt")):
        s = line.strip()
        m = re.match(r"function\(\s*(igloo_[A-Za-z0-9_]+)", s)
        if m:
            cur = m.group(1)
            continue
        if re.match(r"endfunction\s*\(", s):
            cur = None
            continue
        if cur and not s.startswith("#") and "$<TARGET_FILE:IGLOO>" in s and cur not in runners:
            runners.append(cur)
    return runners


def strval(node):
    if isinstance(node, getattr(ast, "Constant", ())) and isinstance(getattr(node, "value", None), str):
        return node.value
    if hasattr(ast, "Str") and isinstance(node, ast.Str):
        return node.s
    return None


def reads_igloo_exe(call):
    f = call.func
    getter = (isinstance(f, ast.Attribute) and f.attr == "get" and isinstance(f.value, ast.Attribute)
              and f.value.attr == "environ") or (isinstance(f, ast.Attribute) and f.attr == "getenv")
    return getter and call.args and strval(call.args[0]) == "IGLOO_EXE"


def builds_shared_path(node):
    """A path literal ending in bin/IGLOO, a join of "bin", "IGLOO", or a pathlib "bin" / "IGLOO"."""
    tail = "bin" + "/" + "IGLOO"
    s = strval(node)
    if s is not None:
        return s == tail or s.endswith("/" + tail)
    if isinstance(node, ast.JoinedStr) and node.values:
        last = strval(node.values[-1])
        return last is not None and (last == "/" + tail or last.endswith("/" + tail))
    if isinstance(node, ast.Call):
        vals = [strval(a) for a in node.args]
        return any(a == "bin" and b == "IGLOO" for a, b in zip(vals, vals[1:]))
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div) and strval(node.right) == "IGLOO":
        return isinstance(node.left, ast.BinOp) and strval(node.left.right) == "bin"
    return False


def unguarded(path):
    try:
        tree = ast.parse(open(path, encoding="utf-8", errors="replace").read(), path)
    except SyntaxError as exc:
        return ["line %s: not parsable (%s)" % (exc.lineno, exc.msg)]
    guarded = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and reads_igloo_exe(node):
            for sub in list(node.args[1:]) + [k.value for k in node.keywords]:
                guarded.update(id(n) for n in ast.walk(sub))
    hits = set()
    for node in ast.walk(tree):
        if id(node) not in guarded and builds_shared_path(node):
            hits.add(node.lineno)
    return ["line %d" % n for n in sorted(hits)]


def python_files(srcdir):
    """Tracked tests/**/*.py (git), else every .py under tests/ outside run directories."""
    try:
        out = subprocess.check_output(["git", "-C", srcdir, "ls-files", "-z", "--", "tests/*.py"],
                                      stderr=subprocess.DEVNULL)
        files = [f for f in out.decode().split("\0") if f]
        if files:
            return sorted(files), "git ls-files"
    except (OSError, subprocess.CalledProcessError):
        pass
    files = []
    for dp, dn, fn in os.walk(os.path.join(srcdir, "tests")):
        dn[:] = [d for d in dn if d not in SKIP_DIRS and not d.startswith(SKIP_PREFIXES)]
        files += [os.path.relpath(os.path.join(dp, f), srcdir) for f in fn if f.endswith(".py")]
    return sorted(files), "walk"


def main():
    if len(sys.argv) != 3:
        print("[FAIL] usage: check_exe_per_tree.py <CTestTestfile.cmake> <source dir>")
        return 1
    ctest, srcdir = sys.argv[1], os.path.abspath(sys.argv[2])
    tests = parse_ctest(ctest)
    if not tests:
        print("[FAIL] no registered test found in %s" % ctest)
        return 1
    me = os.path.basename(__file__)
    shared = re.compile(r"%s/bin/(%s)(?![A-Za-z0-9_.+-])" % (re.escape(srcdir), "|".join(SHARED)))

    runners = solver_helpers(srcdir)
    bad0 = sorted(f for f in runners if f not in HELPERS)
    bad1 = sorted(n for n, (cmd, props) in tests.items()
                  if me not in cmd and (shared.search(cmd) or shared.search(props)))
    helpers, source = helper_names(tests, srcdir)
    bad2 = sorted(n for n in helpers if not re.search(r"(^|;)IGLOO_EXE=", env_of(tests[n][1])))
    pyfiles, how = python_files(srcdir)
    bad3 = []
    for rel in pyfiles:
        if os.path.basename(rel) == me:
            continue
        for where in unguarded(os.path.join(srcdir, rel)):
            bad3.append("%s %s" % (rel, where))

    print("exe-per-tree: %d registered tests, %d helper registrations (%s), %d test scripts (%s)"
          % (len(tests), len(helpers), source, len(pyfiles), how))
    print("  (0) helpers that run the IGLOO executable and are not checked: %d" % len(bad0))
    for n in bad0:
        print("      " + n)
    print("  (1) tests naming %s/bin/{%s}: %d" % ("<source>", ",".join(SHARED), len(bad1)))
    for n in bad1:
        print("      " + n)
    print("  (2) helper registrations without IGLOO_EXE: %d" % len(bad2))
    for n in bad2:
        print("      " + n)
    print("  (3) scripts building bin/IGLOO outside an IGLOO_EXE default: %d" % len(bad3))
    for b in bad3:
        print("      " + b)
    if not helpers:
        print("[FAIL] no helper registration found: the IGLOO_EXE check would be vacuous")
        return 1
    if not runners:
        print("[FAIL] no helper of tests/CMakeLists.txt runs $<TARGET_FILE:IGLOO>: rule (0) would be vacuous")
        return 1
    if bad0 or bad1 or bad2 or bad3:
        print("[FAIL] a test can run the shared executable another build tree publishes")
        return 1
    print("[PASS] every test runs its build tree's executable")
    return 0


if __name__ == "__main__":
    sys.exit(main())
