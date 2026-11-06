! e2e-rc-gate's probe: a list-directed read past a blank internal record (a run-time error, exit code below 128).
program crash_code_probe
  implicit none
  character(len=8) :: rec = ' '
  integer :: v
  read(rec, *) v
  print *, v
end program crash_code_probe
