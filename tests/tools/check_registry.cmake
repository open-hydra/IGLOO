# Gate: docs/user/registry.md is GENERATED output that is also TRACKED, so it can
# silently disagree with the registry it documents. Regenerate into the build tree
# and diff against the committed file.
#
# Invoked by ctest with -DDOCGEN=<exe> -DTRACKED=<path> -DGENERATED=<path>.
# Twice in one session (2026-07-30) the tracked file drifted: once when a new
# `blowing` key was added without re-running DocGen, and the same regeneration
# exposed `TC` missing from the `evaporation` allowed list since the F2 merge.

execute_process(COMMAND "${DOCGEN}" "${GENERATED}" RESULT_VARIABLE rc
                OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(NOT rc EQUAL 0)
    message(FATAL_ERROR "DocGen failed (exit ${rc}):\n${out}${err}")
endif()

if(NOT EXISTS "${GENERATED}")
    message(FATAL_ERROR "DocGen wrote nothing to ${GENERATED}")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files "${TRACKED}" "${GENERATED}"
                RESULT_VARIABLE differ OUTPUT_QUIET ERROR_QUIET)

if(differ EQUAL 0)
    message(STATUS "registry-docs: docs/user/registry.md matches the registry")
else()
    message(FATAL_ERROR
        "docs/user/registry.md is STALE — it does not match what DocGen produces "
        "from src/lib/config/Register_IGLOO.f90.\n"
        "Regenerate and commit:\n"
        "    cmake --build <build-dir> --target DocGen\n"
        "    ./bin/DocGen\n"
        "(run DocGen from the repo root; it writes docs/user/registry.md)\n"
        "Tracked:   ${TRACKED}\n"
        "Generated: ${GENERATED}")
endif()
