# rc-gate — the e2e harness fails a non-zero exit

**Purpose.** Every `igloo_e2e_case` runs one command (`igloo_e2e_command` in `tests/CMakeLists.txt`): the
executable in the case directory, the exit-code gate (`igloo_rc_gate`), `tools/check_finite.py`, then the case's
`check.py`. The gate fails the case on any non-zero exit before the oracle can read a partial `OUTPUT/`: a signal
death (exit above 128) prints `[FAIL] solver died on signal (exit N)`, every other non-zero exit `[FAIL] solver exit
N`. A runtime's crash codes cannot be listed by value: a severe run-time error exits with its error number (ifx), a
run-time library error with 2 (gfortran), a caught signal with 1, an integer `error stop` with its code, ifx's
character `error stop` with 128. IGLOO itself never exits non-zero on a run that completes.
This case pins that rule on the very string the e2e cases run.

**Probe.** `tools/crash_code_probe.f90`, built as `crash_code_probe`: a list-directed integer read from a blank
internal record, which ends the program on a run-time error ("end-of-file during read") with an exit code below
128. It writes no `OUTPUT/` file, so `check_finite.py` passes and only the exit-code gate can stop the command.
No stdin is read.

**Gate** (ctest's regular expressions, independent of the gate's code). `e2e-rc-gate` runs
`igloo_e2e_command` around the probe in this directory; it passes only when the output holds `[FAIL] solver exit
N` (`PASS_REGULAR_EXPRESSION`) and not this directory's `check.py` line `[PASS] rc-gate oracle reached`
(`FAIL_REGULAR_EXPRESSION`, which wins). Measured: `[FAIL] solver exit 24` under ifx, `[FAIL] solver exit 2`
under gfortran; the number is reported, not asserted.

**RED.** The gate failing only signal deaths (`rc >= 128`): the probe's 24 passes it, `check_finite.py` finds no
file, `check.py` prints `[PASS] rc-gate oracle reached` and the test fails on the FAIL expression.
