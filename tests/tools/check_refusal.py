#!/usr/bin/env python3
"""
check_refusal.py -- gate for a case the solver must REFUSE at setup.

    check_refusal.py "<expected message substring>" [more substrings...]

Run from the case directory after `IGLOO 1>run_out.txt 2>run_err.txt; echo $? > rc.txt`.
PASS requires ALL of:
  1. exit code exactly IGLOO_ERROR_STOP_RC -- what this compiler's `error stop 'msg'` returns,
     probed at configure time by tools/stop_code_probe.f90 and set by ctest: 128 under
     ifx/ifort (and why the e2e harness, which reads rc >= 128 as "died on signal", cannot host
     a refusal), 1 under gfortran. Unset (a hand-run outside ctest) means 128. 0 = did not
     refuse; 129-255 = a signal death (a SIGSEGV after the message is not a refusal); a plain
     `stop` exits 0 and is therefore also caught (ledger O26);
  2. every expected substring in run_err.txt -- the `error stop` PAYLOAD, i.e. the string
     of the stop that actually ended the run. A message merely printed to stdout before some
     other death does not count, so key the expected text on the stop string, not on the
     [ERROR] line that precedes it;
  3. nothing was set up or integrated past the choke point: no "Placing particles" (injection),
     "Compute particles dynamics" or "Stop condition" line, and no OUTPUT/*.dat.
A green solver (rc 0, output written) FAILS this gate on every count: proven on drag-stokes
(see tests/infrastructure/refusals/INFO.md).
"""
import glob
import os
import sys


def main(argv):
    if len(argv) < 2:
        print("[FAIL] usage: check_refusal.py <expected substring> [...]")
        return 2
    expected = argv[1:]
    stop_rc = int(os.environ.get("IGLOO_ERROR_STOP_RC", "128"))
    rc = 0
    try:
        code = int(open("rc.txt").read().split()[0])
    except (FileNotFoundError, ValueError, IndexError):
        print("[FAIL] rc.txt missing or unreadable -- the harness did not record the exit code")
        return 1
    log, err = "", ""
    try:
        log = open("run_out.txt", errors="replace").read()
    except FileNotFoundError:
        pass
    try:
        err = open("run_err.txt", errors="replace").read()
    except FileNotFoundError:
        pass
    if code == 0:
        print("[FAIL] solver exited 0 -- it did not refuse (a plain `stop` exits 0 too)")
        rc = 1
    elif code != stop_rc:
        print(f"[FAIL] solver exited {code}, not {stop_rc} (this compiler's `error stop`) -- "
              "a signal death or some other exit, not a refusal")
        rc = 1
    for s in expected:
        if s in err:
            print(f"refused with: {s!r} PASS")
        else:
            where = " (it is in run_out.txt, but the run did not END on it)" if s in log else ""
            print(f"[FAIL] expected error-stop payload not found in run_err.txt: {s!r}{where}")
            rc = 1
    for marker in ("Placing particles", "Compute particles dynamics", "Stop condition"):
        if marker in log:
            print(f"[FAIL] '{marker}' in the log -- the solver integrated before/without refusing")
            rc = 1
    dats = glob.glob(os.path.join("OUTPUT", "*.dat"))
    if dats:
        print(f"[FAIL] {len(dats)} OUTPUT/*.dat written -- the run produced output instead of refusing")
        rc = 1
    if rc == 0:
        print(f"[PASS] refused at setup (exit {code}), nothing integrated")
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv))
