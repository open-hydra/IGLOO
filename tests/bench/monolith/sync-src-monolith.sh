#!/usr/bin/env bash
# Ship THIS IGLOO tree (the checkout this script lives in, libraries included) to the cluster as files and submit
# its build there as a short job; the job installs $IGLOO_PERF_ROOT/bin/IGLOO-<cfg> and its .commit file.
#   IGLOO_PERF_REMOTE=<ssh host> IGLOO_PERF_ROOT=<harness root there> sync-src-monolith.sh [cfg ...]   (default: hyb)
#   sync-src-monolith.sh --local <dir>     the same checks and transfer into <dir>/src and <dir>/scripts, no build
# Optional: IGLOO_PERF_PARTITION (default free) for the build job.
#
# The tree travels without .git, except FiNeR's five dependencies: IGLOO's configure refuses a dependency directory
# that is not a git checkout and fetches from GitHub when its HEAD is off the pin (ensure_finer_dependencies in
# CMakeLists.txt), so they travel as checkouts at their pins. The three library submodules travel as plain files
# (their CMakeLists.txt is all the configure asks), each checked first to be a checkout of its own at its gitlink.
# OSLO's tracked lib/Intel-ODE/lib/intel64/libiode_intel64.a travels too (OSLO links it; *.a is excluded otherwise).
# IGLOO-commit.txt (ignored by git) names the commits and dirty-file counts of what was shipped.
set -eu
HERE=$(cd "$(dirname "$0")" && pwd -P)
SRC=$(cd "$HERE/../../.." && pwd -P)
DEPS="BeFoR64 FACE FLAP PENF StringiFor"
LIBS="lib/ORION lib/OSLO lib/third_party/FiNeR"
ODE=lib/OSLO/lib/Intel-ODE/lib/intel64/libiode_intel64.a
die () { echo "sync: $*" >&2; exit 1; }

cd "$SRC"
pin_of () { sed -n "s/^ *set(_FINER_DEP_PIN_$1  *\([0-9a-f]\{40\}\))$/\1/p" "$SRC/CMakeLists.txt"; }
for L in $LIBS; do
  link=$(git -C "$SRC" ls-tree HEAD "$L" | awk '{print $3}')
  [ -n "$link" ] || die "no gitlink for $L"
  [ -e "$SRC/$L/.git" ] || die "$L is not a git checkout of its own: clone it over https at $link"
  [ "$(cd "$SRC/$L" && git rev-parse --show-toplevel)" = "$SRC/$L" ] || die "$L is not the top of its own checkout"
  head=$(git -C "$SRC/$L" rev-parse HEAD)
  [ "$head" = "$link" ] || die "$L is at $head, its gitlink is $link"
done
for d in $DEPS; do
  pin=$(pin_of "$d"); dir=$SRC/lib/third_party/FiNeR/src/third_party/$d
  [ -n "$pin" ] || die "no _FINER_DEP_PIN_$d in CMakeLists.txt"
  [ -e "$dir/.git" ] || die "FiNeR dependency $d is not a git checkout: configure the tree once locally"
  head=$(git -C "$dir" rev-parse HEAD)
  [ "$head" = "$pin" ] || die "FiNeR dependency $d is at $head, the pin is $pin: configure the tree once locally"
done
[ -f "$SRC/$ODE" ] || die "$ODE is missing"

{
  echo "IGLOO $(git -C "$SRC" rev-parse HEAD) dirty=$(git -C "$SRC" status --porcelain | wc -l)"
  for L in $LIBS; do
    echo "$L gitlink=$(git -C "$SRC" ls-tree HEAD "$L" | awk '{print $3}') head=$(git -C "$SRC/$L" rev-parse HEAD)" \
         "dirty=$(git -C "$SRC/$L" status --porcelain | wc -l)"
  done
  for d in $DEPS; do echo "FiNeR-dependency $d head=$(git -C "$SRC/lib/third_party/FiNeR/src/third_party/$d" rev-parse HEAD)"; done
  echo "synced $(date -u +%Y-%m-%dT%H:%M:%SZ) from $(hostname):$SRC"
} > "$SRC/IGLOO-commit.txt"
cat "$SRC/IGLOO-commit.txt"

# rsync takes the first matching rule: the two includes come before the excludes they override.
FILTERS=(
  --include="/$ODE"
  --include='/lib/third_party/FiNeR/src/third_party/*/.git/***'
  --exclude='.git' --exclude='.claude' --exclude='.claude-context' --exclude='/plan-bucket' --exclude='/site'
  --exclude='__pycache__' --exclude='/build*' --exclude='build/' --exclude='bin' --exclude='bin-*' --exclude='modules'
  --exclude='*.o' --exclude='*.mod' --exclude='*.a'
  --exclude='/lib/ORION/lib/TecIO/tecio*-build-*' --exclude='/lib/ORION/lib/TecIO/tecio*-install-*'
  --exclude='/tests/*/*/OUTPUT' --exclude='/tests/*/*/ref-*'
)

if [ "${1:-}" = --local ]; then
  D=${2:?--local <dir>}
  mkdir -p "$D/src" "$D/scripts" "$D/bin" "$D/logs" "$D/case"
  rsync -a --delete "${FILTERS[@]}" "$SRC/" "$D/src/"
  rsync -a --delete "$HERE/" "$D/scripts/"
  cmp -s "$SRC/$ODE" "$D/src/$ODE" || die "$ODE did not arrive intact"
  for d in $DEPS; do
    [ "$(git -C "$D/src/lib/third_party/FiNeR/src/third_party/$d" rev-parse HEAD)" = "$(pin_of "$d")" ] \
      || die "FiNeR dependency $d did not arrive at its pin"
  done
  echo "rehearsal copy in $D: nothing submitted"
  exit 0
fi

REMOTE=${IGLOO_PERF_REMOTE:?IGLOO_PERF_REMOTE: the ssh host of the cluster}
ROOT=${IGLOO_PERF_ROOT:?IGLOO_PERF_ROOT: the harness root on the cluster}
PART=${IGLOO_PERF_PARTITION:-free}

ssh "$REMOTE" "mkdir -p $ROOT/src $ROOT/scripts $ROOT/bin $ROOT/logs $ROOT/case"
rsync -a --delete "${FILTERS[@]}" "$SRC/" "$REMOTE:$ROOT/src/"
rsync -a --delete "$HERE/" "$REMOTE:$ROOT/scripts/"

# What the configure there needs arrived: OSLO's archive byte for byte, the dependencies at their pins.
lmd5=$(md5sum < "$SRC/$ODE" | cut -d' ' -f1)
rmd5=$(ssh "$REMOTE" "md5sum < $ROOT/src/$ODE" | cut -d' ' -f1)
[ "$lmd5" = "$rmd5" ] || die "$ODE: md5 $rmd5 on $REMOTE, $lmd5 here"
echo "ODE archive md5 $rmd5 (both sides)"
pins=""; for d in $DEPS; do pins="$pins $d=$(pin_of "$d")"; done
ssh "$REMOTE" "bash -lc 'for p in $pins; do d=\${p%%=*}; h=\$(git -C $ROOT/src/lib/third_party/FiNeR/src/third_party/\$d rev-parse HEAD); [ \"\$h\" = \"\${p#*=}\" ] && echo \"dependency \$d at its pin\" || { echo \"dependency \$d at \$h, not at its pin\"; exit 1; }; done'" \
  || die "FiNeR's dependencies did not arrive at their pins"

# The build runs in a job on the timing partition, never on the login node; its log is logs/igloo-build-<cfg>-<id>.out.
# A login profile may print to stdout: the job id is the line sbatch --parsable writes.
for cfg in "${@:-hyb}"; do
  id=$(ssh "$REMOTE" "bash -lc 'cd $ROOT && sbatch --parsable --job-name=igloo-build-$cfg --partition=$PART --nodes=1 --ntasks=1 --cpus-per-task=16 --time=00:30:00 --output=$ROOT/logs/%x-%j.out --wrap=\"JOBS=16 bash -l $ROOT/scripts/build-cfg.sh $cfg $ROOT/src $ROOT/bin\"'" \
       | grep -E '^[0-9]+(;[^ ]*)?$' | tail -1)
  [ -n "$id" ] || die "the build job of $cfg was not submitted (no job id from sbatch)"
  echo "BUILD cfg=$cfg job=${id%%;*} submitted $(date -u +%Y-%m-%dT%H:%M:%SZ) log=$ROOT/logs/igloo-build-$cfg-${id%%;*}.out"
done
