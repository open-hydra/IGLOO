#!/usr/bin/env bash
#===============================================================================
#         FILE: test.sh
#        USAGE: ./test.sh [all|<category>|<test-name>|clean] [-- <cmake args>]
#  DESCRIPTION: MOSE-style runner for the IGLOO test suite. Configures + builds
#               build/verif/ (build/mpi/ when USE_MPI=ON), whose tests run that
#               tree's own IGLOO, publishes its IGLOO and DocGen to bin/, then
#               dispatches to ctest by label or name.
#               Categories: standard evaporation combustion solidification breakup infrastructure unit e2e
#===============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IGLOO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CATEGORIES="standard evaporation combustion solidification breakup infrastructure repeatability mpi unit e2e"

usage() {
    echo "Usage: ./test.sh [all|<category>|<test-name>|clean] [-- <extra cmake args>]"
    echo "  all            build + run the full suite + aggregated report"
    echo "  <category>     one of: ${CATEGORIES}"
    echo "  mpi            USE_MPI=ON ./test.sh mpi -- -DUSE_TECIO=OFF (builds build/mpi, see mpi/INFO.md)"
    echo "  <test-name>    single ctest name (e.g. conv-nu, test_breakup_tab)"
    echo "  clean          remove what the tests write in their directories and build/verif/ (build/mpi/ with USE_MPI=ON)"
    exit 1
}

[ $# -ge 1 ] || usage
TARGET="$1"; shift
EXTRA_CMAKE_ARGS=()
if [[ $# -gt 0 && "$1" == "--" ]]; then shift; EXTRA_CMAKE_ARGS=("$@"); fi

# USE_MPI defaults OFF. USE_MPI=ON (or -- -DUSE_MPI=ON) builds build/mpi, so the serial tree is
# never reconfigured as an MPI build.
MPI="${USE_MPI:-OFF}"
for a in "${EXTRA_CMAKE_ARGS[@]}"; do
    case "${a}" in -DUSE_MPI=*|-DUSE_MPI:BOOL=*) MPI="${a#*=}" ;; esac
done
case "${MPI^^}" in ON|TRUE|YES|Y|1) MPI=ON ;; *) MPI=OFF ;; esac
BUILD_DIR="${IGLOO_ROOT}/build/verif"
if [[ "${MPI}" == ON ]]; then BUILD_DIR="${IGLOO_ROOT}/build/mpi"; fi

if [[ "${TARGET}" == clean ]]; then
    # every directory a test runs in holds an input.ini; remove what the runs and their checks write there
    while IFS= read -r d; do
        rm -rf "${d}/OUTPUT" "${d}"/run_out*.txt "${d}"/run_err*.txt "${d}/rc.txt" "${d}"/run_n*.txt \
               "${d}"/err_n*.txt "${d}"/tec_n*/ "${d}"/ref-*/ "${d}"/threads-omp*/ "${d}/diag_threads.txt"
        echo "cleaned ${d#"${SCRIPT_DIR}"/}"
    done < <(find "${SCRIPT_DIR}" -name input.ini -not -path '*/OUTPUT/*' -not -path '*/ref-*/*' \
                  -not -path '*/_run-*/*' -not -path '*/threads-omp*/*' -printf '%h\n' | sort -u)
    rm -rf "${BUILD_DIR}"
    exit 0
fi

cmake -B "${BUILD_DIR}" -S "${IGLOO_ROOT}" \
    -DUSE_OPENMP=ON \
    -DUSE_MPI=${MPI} \
    -DUSE_SUNDIALS=OFF \
    -DBUILD_VERIFICATION=ON \
    -DCMAKE_BUILD_TYPE=RELEASE \
    "${EXTRA_CMAKE_ARGS[@]}"

cmake --build "${BUILD_DIR}" -j

# bin/ holds the executables of the tree tested last: copy, then rename (a launch never sees a partial file).
mkdir -p "${IGLOO_ROOT}/bin"
for exe in IGLOO DocGen; do
    part="${IGLOO_ROOT}/bin/.${exe}.test-sh.$$.part"
    cp -f "${BUILD_DIR}/bin/${exe}" "${part}"
    mv -f "${part}" "${IGLOO_ROOT}/bin/${exe}"
done

if [[ "${TARGET}" == all ]]; then
    ctest --test-dir "${BUILD_DIR}" --output-on-failure
    # aggregate per-test CSVs into one report + consistency gate
    python3 -B "${SCRIPT_DIR}/tools/aggregate_report.py" \
        "${BUILD_DIR}/tests" "${SCRIPT_DIR}"
    # MOSE-style: refresh each case's OUTPUT/<case>.svg overlay (non-gating;
    # verify.py is a no-op where matplotlib is absent)
    for v in "${SCRIPT_DIR}"/*/*/verify.py; do
        [ -f "${v}" ] || continue
        ( cd "$(dirname "${v}")" && python3 -B verify.py ) || true
    done
    # unit-family overlays from the verif_dump curves written at ctest time
    python3 -B "${SCRIPT_DIR}/tools/plot_curves.py" \
        --build-dir "${BUILD_DIR}/tests" || true
    # MOSE-style site pipeline: sync the overlays into docs/vv/images/ where the
    # per-case V&V pages {% include %} them (ignore-missing until first sync).
    mkdir -p "${IGLOO_ROOT}/docs/vv/images"
    for s in "${SCRIPT_DIR}"/*/*/OUTPUT/*.svg "${BUILD_DIR}"/tests/curves-svg/*.svg; do
        [ -f "${s}" ] || continue
        cp -f "${s}" "${IGLOO_ROOT}/docs/vv/images/"
        echo "svg -> docs/vv/images/$(basename "${s}")"
    done
elif grep -qw "${TARGET}" <<< "${CATEGORIES}"; then
    ctest --test-dir "${BUILD_DIR}" -L "^${TARGET}\$" --output-on-failure
else
    ctest --test-dir "${BUILD_DIR}" -R "^${TARGET}\$" --output-on-failure
fi
