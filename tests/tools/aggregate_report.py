#!/usr/bin/env python3
"""T8/T9 closure for the IGLOO verification suite.

T8 — aggregate every per-test verif_*.csv (written into the ctest working dir)
     into a single markdown report with a per-row PASS/FAIL table and totals.
T9 — consistency gate, nonzero exit if:
     (a) any CSV row is FAIL;
     (b) any test executable's family directory lacks an INFO.md;
     (c) any CSV exists with zero parsed rows (report corruption).

Usage: aggregate_report.py <build_tests_dir> <source_tests_dir> [-o report.md]
"""
import argparse
import glob
import os
import sys

# family dir of each test exe (kept in sync with CMakeLists.txt)
FAMILY_DIRS = [
    "support",
    "infrastructure/gas_reconstruction", "infrastructure/ini_pipeline",
    "standard/drag", "standard/temperature",
    "evaporation", "evaporation/interface-neq", "evaporation/tc-analytic",
    "combustion",
    "breakup/tab", "breakup/etab", "breakup/pilch-erdman",
    "breakup/reitz-diwakar",
]


def parse_csv(path):
    rows = []
    with open(path) as f:
        header = f.readline()
        if not header.lower().startswith("case"):
            return rows
        for line in f:
            c = [t.strip() for t in line.split(",")]
            if len(c) < 8:
                continue
            # variable names may themselves contain commas (e.g. "y,yDot"):
            # csv schema is case,variable...,Linf,L2,p_obs,p_expected,tol,result
            rows.append({
                "case": c[0],
                "variable": ",".join(c[1:len(c) - 6]),
                "Linf": c[-6], "tol": c[-2], "result": c[-1],
            })
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("build_dir")
    ap.add_argument("source_dir")
    ap.add_argument("-o", "--out", default=None)
    a = ap.parse_args()
    out = a.out or os.path.join(a.build_dir, "verification_report.md")

    bad = []

    # ---- T9(b): every family dir carries an INFO.md (support exempt) ----
    for d in FAMILY_DIRS:
        if d == "support":
            continue
        info = os.path.join(a.source_dir, d, "INFO.md")
        if not os.path.isfile(info):
            bad.append(f"missing INFO.md: {d}/")

    # ---- T8: aggregate ----
    csvs = sorted(glob.glob(os.path.join(a.build_dir, "verif_*.csv")))
    lines = ["# IGLOO verification suite — aggregated report", ""]
    total = passed = 0
    lines += ["| csv | case | variable | Linf | tol | result |",
              "|---|---|---|---|---|---|"]
    for path in csvs:
        rows = parse_csv(path)
        if not rows:
            bad.append(f"no parsable rows: {os.path.basename(path)}")
            continue
        for r in rows:
            total += 1
            ok = r["result"].upper() == "PASS"
            passed += ok
            if not ok:
                bad.append(f"FAIL row: {os.path.basename(path)}:{r['case']}")
            lines.append(f"| {os.path.basename(path)} | {r['case']} | "
                         f"{r['variable']} | {r['Linf']} | {r['tol']} | {r['result']} |")
    lines += ["", f"**{passed}/{total} rows PASS** across {len(csvs)} report files.",
              "", "Bug-transcription pins emit no CSV: test_drag_probes, "
              "test_heat_probes, test_evap_probes. They were written as xfail probes, "
              "but every bug they cover is fixed, so they exit 0 and run as ORDINARY "
              "gates — no WILL_FAIL property is set on any test.", ""]
    if bad:
        lines += ["## CONSISTENCY FAILURES", ""] + [f"- {b}" for b in bad]

    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"[report] {out}: {passed}/{total} rows PASS, {len(csvs)} csv files")
    for b in bad:
        print(f"[consistency-FAIL] {b}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
