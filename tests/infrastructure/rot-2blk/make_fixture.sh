#!/usr/bin/env bash
# rot-2blk fixture. Gas of the rotated two-block box (INPUT/solfile.tec) and of its one-block sibling (INPUT-1blk/),
# the ATLAS BCB deck's mesh (INPUT-atlas/mesh.tec), all from tools/make_box_case.py; INPUT/drop-bc.txt is ATLAS BCB's
# output on INPUT-atlas/ (ATLAS a3689c3, open-hydra/ATLAS main, gfortran build), regenerated only when ATLASDIR
# names an ATLAS tree:  ATLASDIR=<atlas> ./make_fixture.sh
set -euo pipefail
cd "$(dirname "$0")"
GEN=../../tools/make_box_case.py
G=2000   # gas temperature gradient along z [K/m]

python3 "$GEN" INPUT-atlas --mesh-only --blocks 2 --rotate2 y
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
python3 "$GEN" "$tmp/gas" --blocks 2 --rotate2 y --kv 0.1 --kt 0.5 --zgradient $G
cp "$tmp/gas/solfile.tec" INPUT/solfile.tec
python3 "$GEN" INPUT-1blk --kv 0.1 --kt 0.5 --zgradient $G

if [ -n "${ATLASDIR:-}" ]; then
    cp INPUT-atlas/input.ini INPUT-atlas/drop-phase.txt INPUT-atlas/mesh.tec "$tmp/"
    (cd "$tmp" && "$ATLASDIR/ATLAS.sh" BCB)
    cp "$tmp/fromATLAStoSolver/drop-bc.txt" INPUT/drop-bc.txt
    echo "INPUT/drop-bc.txt written by ATLAS BCB ($ATLASDIR)"
fi
