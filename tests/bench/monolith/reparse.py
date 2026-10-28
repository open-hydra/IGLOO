#!/usr/bin/env python3
"""Parse again the repetitions of a timing job whose report fields are missing, and write the corrected .out.

    python3 reparse.py <job.out> <logs dir> > <corrected .out>

A repetition whose IREP line carries no report fields (no `reports=`: the parser failed on that run's output) is
parsed again from the log igloo_run kept (<logs dir>/<tag>-r<rep>.log: the run's stdout, a `---- stderr ----`
line, its stderr) with igloo_parse of _igloo_common.sh beside this script, under the job's wait policy (its
WAITPOLICY line). The fields are appended to its IREP line, its IAFF line is inserted before it, and the IOK line
of its point is computed again from the point's IREP lines with igloo_summary. Every other line is copied
unchanged; an IFAIL line is decided from the log by igloo_run itself, never by the parser, and stays as it is.
A summary goes to stderr.
"""
import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))


def snippet(text, name):
    """The python of the shell function `name` of _igloo_common.sh (its <<'PY' heredoc)."""
    m = re.search(r"\n%s \(\) \{\n[^\n]*<<'PY'\n(.*?)\nPY\n\}" % re.escape(name), text, re.S)
    if not m:
        sys.exit("reparse: no %s in _igloo_common.sh" % name)
    return m.group(1) + "\n"


def run(code, args, env):
    p = subprocess.run([sys.executable, "-"] + args, input=code, stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, universal_newlines=True, env=env)
    if p.returncode != 0:
        raise RuntimeError(p.stderr.strip().splitlines()[-1] if p.stderr.strip() else "rc %d" % p.returncode)
    return p.stdout


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    out_path, logs = sys.argv[1:3]
    common = open(os.path.join(HERE, "_igloo_common.sh")).read()
    parse, summary = snippet(common, "igloo_parse"), snippet(common, "igloo_summary")
    lines = open(out_path, errors="surrogateescape").read().splitlines(True)

    env = dict(os.environ)
    for line in lines:
        if line.startswith("WAITPOLICY "):
            for k, v in re.findall(r"(\w+)=(\S+)", line):
                if v == "unset":
                    env.pop(k, None)
                else:
                    env[k] = v
            break

    fixed, iaff, tags = {}, {}, []
    tmp = tempfile.mkdtemp(prefix="reparse.")
    try:
        for i, line in enumerate(lines):
            if not line.startswith("IREP ") or " reports=" in line:
                continue
            d = dict(re.findall(r"(\w+)=(\S*)", line))
            tag, rep = d.get("tag", ""), d.get("rep", "")
            log = os.path.join(logs, "%s-r%s.log" % (tag, rep))
            try:
                text = open(log, errors="surrogateescape").read()
            except OSError as e:
                sys.stderr.write("reparse: %s r%s: %s\n" % (tag, rep, e))
                continue
            sout, sep, serr = text.partition("---- stderr ----\n")
            fo, fe = os.path.join(tmp, "out"), os.path.join(tmp, "err")
            with open(fo, "w", errors="surrogateescape") as f:
                f.write(sout)
            with open(fe, "w", errors="surrogateescape") as f:
                f.write(serr if sep else "")
            try:
                res = run(parse, [tag, rep, fo, fe], env).splitlines()
            except RuntimeError as e:
                sys.stderr.write("reparse: %s r%s: the parser failed again: %s\n" % (tag, rep, e))
                continue
            if len(res) < 2 or " reports=" not in " " + res[1]:
                sys.stderr.write("reparse: %s r%s: the log holds no report\n" % (tag, rep))
                continue
            body = line.rstrip("\n")
            fixed[i] = body + ("" if body.endswith(" ") else " ") + res[1] + "\n"
            iaff[i] = res[0] + "\n"
            if tag not in tags:
                tags.append(tag)

        new = list(lines)
        for i in fixed:
            new[i] = fixed[i]
        ioks = {}
        for tag in tags:
            reps = [l for l in new if l.startswith("IREP tag=%s " % tag)]
            iok = [j for j, l in enumerate(new) if l.startswith("IOK tag=%s " % tag)]
            if len(iok) != 1:
                sys.stderr.write("reparse: %s: %d IOK lines, left as they are\n" % (tag, len(iok)))
                continue
            asked = re.search(r" reps=(\S+)", new[iok[0]])
            fr = os.path.join(tmp, "reps")
            with open(fr, "w", errors="surrogateescape") as f:
                f.writelines(reps)
            ioks[iok[0]] = run(summary, [fr, asked.group(1) if asked else str(len(reps))], env)
        for j in ioks:
            new[j] = ioks[j]
    finally:
        for f in ("out", "err", "reps"):
            try:
                os.remove(os.path.join(tmp, f))
            except OSError:
                pass
        os.rmdir(tmp)

    w = sys.stdout.buffer
    for i, line in enumerate(new):
        if i in iaff:
            w.write(iaff[i].encode("utf-8", "surrogateescape"))
        w.write(line.encode("utf-8", "surrogateescape"))
    w.flush()
    sys.stderr.write("reparse: %d IREP line(s) completed, %d IAFF line(s) inserted, %d IOK line(s) recomputed (%s)\n"
                     % (len(fixed), len(iaff), len(ioks), " ".join(tags) or "-"))


if __name__ == "__main__":
    main()
