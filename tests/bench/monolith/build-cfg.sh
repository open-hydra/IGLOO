#!/usr/bin/env bash
# Build one IGLOO configuration into <IGLOO root>/build-<cfg> and, given an install directory, copy the executable
# there as IGLOO-<cfg> beside the tree's IGLOO-commit.txt (IGLOO-<cfg>.commit); the probe configuration builds and
# copies the I/O probe's executables instead (io-probe, io-probe-wr[-main|-mod|-both], io-probe.commit).
#   build-cfg.sh <cfg> [<IGLOO root> [<install dir>]]     cfg: hyb (MPI + OpenMP, the campaign binary) | omp |
#                                                              probe (the hyb tree with IGLOO_IO_PROBE=ON)
#
# Runs inside a job on the partition the timing jobs use: -xHost (cmake/SetFortranFlags.cmake) ties the binary to
# the build node's instruction set, so the build node's CPU model and flag digest are printed beside its md5.
# TecIO is off (the MPI tree would link ORION's teciompi otherwise) and the tree publishes nothing into
# <IGLOO root>/bin. Intel MPI's mpiifort wraps ifx through I_MPI_F90.
set -eu
CFG=${1:?configuration: hyb or omp}
SRC=${2:-$(cd "$(dirname "$0")/../../.." && pwd)}
INST=${3:-}
FC=${FC:-$(command -v ifx || true)}
MPIFC=${MPIFC:-$(command -v mpiifort || true)}
CC=${CC:-$(command -v icx || command -v gcc || true)}
CXX=${CXX:-$(command -v icpx || command -v g++ || true)}
EXES=IGLOO
case $CFG in
  hyb) OPTS="-DUSE_MPI=ON -DUSE_OPENMP=ON" ;;
  omp) OPTS="-DUSE_MPI=OFF -DUSE_OPENMP=ON" ;;
  probe) OPTS="-DUSE_MPI=ON -DUSE_OPENMP=ON -DIGLOO_IO_PROBE=ON"
         EXES="io-probe io-probe-wr io-probe-wr-main io-probe-wr-mod io-probe-wr-both" ;;
  *) echo "unknown configuration $CFG (hyb, omp or probe)" >&2; exit 2 ;;
esac
[ -n "$FC" ] && [ -n "$CC" ] && [ -n "$CXX" ] || { echo "build-cfg: FC, CC and CXX needed (FC=$FC CC=$CC CXX=$CXX)" >&2; exit 2; }
[ "$CFG" = omp ] || [ -n "$MPIFC" ] || { echo "build-cfg: no mpiifort for $CFG" >&2; exit 2; }
export I_MPI_F90=$FC

echo "BUILDHOST $(hostname) | CPU $(lscpu | sed -n 's/^Model name: *//p') | flags-md5 $(lscpu | grep '^Flags' | md5sum | cut -c1-12)"
echo "COMPILERS FC=$FC ($("$FC" --version 2>&1 | head -1)) MPIFC=${MPIFC:-none} CC=$CC CXX=$CXX"
echo "CMAKE $(cmake --version | head -1)"

cmake -S "$SRC" -B "$SRC/build-$CFG" \
  -DCMAKE_BUILD_TYPE=RELEASE \
  -DCMAKE_Fortran_COMPILER="$FC" \
  -DCMAKE_C_COMPILER="$CC" \
  -DCMAKE_CXX_COMPILER="$CXX" \
  ${MPIFC:+-DMPI_Fortran_COMPILER="$MPIFC"} \
  -DUSE_TECIO=OFF -DUSE_SUNDIALS=OFF -DIGLOO_PUBLISH_BIN=OFF \
  $OPTS

if [ "$CFG" = probe ]; then cmake --build "$SRC/build-$CFG" -j "${JOBS:-16}" --target $EXES
else cmake --build "$SRC/build-$CFG" -j "${JOBS:-16}"; fi

for e in $EXES; do
  B=$SRC/build-$CFG/bin/$e
  echo "== $CFG: $B"
  md5sum "$B"
done
B=$SRC/build-$CFG/bin/${EXES%% *}
echo "libmpi $(ldd "$B" | grep -c libmpi || true)"
ldd "$B" | grep -E 'libmpi|libiomp5' || echo "   (no libmpi/libiomp5 linked)"

if [ -n "$INST" ]; then
  mkdir -p "$INST"
  if [ "$CFG" = probe ]; then
    for e in $EXES; do
      cp -f "$SRC/build-$CFG/bin/$e" "$INST/$e.tmp" && mv -f "$INST/$e.tmp" "$INST/$e"
      echo "INSTALLED $(md5sum "$INST/$e")"
    done
    C=$INST/io-probe.commit
  else
    cp -f "$B" "$INST/IGLOO-$CFG.tmp" && mv -f "$INST/IGLOO-$CFG.tmp" "$INST/IGLOO-$CFG"
    echo "INSTALLED $(md5sum "$INST/IGLOO-$CFG")"
    C=$INST/IGLOO-$CFG.commit
  fi
  if [ -f "$SRC/IGLOO-commit.txt" ]; then cp -f "$SRC/IGLOO-commit.txt" "$C"
  else echo "unknown (no IGLOO-commit.txt in $SRC)" > "$C"; fi
fi
