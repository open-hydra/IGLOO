! Configure-time probe (tests/CMakeLists.txt try_run): the exit code this compiler's runtime gives
! `error stop` with a character stop-code -- the form of every IGLOO refusal. tools/check_refusal.py
! gates on exactly this code (IGLOO_ERROR_STOP_RC).
program stop_code_probe
  error stop 'stop_code_probe'
end program stop_code_probe
