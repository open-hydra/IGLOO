!>
!> Two-sweep repeatability harness (work package C).
!>
!> The 49-entry suite performs exactly ONE sweep, so it can never detect a state leak
!> across a second solve(). This driver does what an embedding does -- integrate the same
!> pinned population twice -- and the gate is:
!>
!>     for a STEADY gas, sweep 2's output must be byte-identical to sweep 1's.
!>
!> That one statement catches every leak by construction: a surviving child changes the
!> record count, a stale `d` changes the mass column, an unreset `brkupVar` changes the
!> exit position, a corrupted accumulator changes source.tec.
!>
!> OUTPUT/ is moved aside between sweeps rather than suffixed, so the harness needs no
!> production filename change (that is P7's sweep tag, for the embedding's benefit).
!>
!> Two modes:
!>   twosweep                 steady gas, 2 sweeps. Sweep 2 must reproduce sweep 1.
!>   twosweep <gasfile.tec>   gas-refresh cycle, 3 sweeps against the external_gas hook:
!>                              sweep 1  original field
!>                              sweep 2  field with U doubled   -> MUST differ from sweep 1
!>                              sweep 3  original field again   -> MUST match sweep 1
!>                            Sweep 2 proves the refresh reaches the solver at all; sweep 3
!>                            proves it fully overwrites, leaving no residue of sweep 2.
program twosweep
  use IGLOO_module,     only: obj_IGLOO
  use IGLOO_IO,         only: read_TECsolfile
  use Lib_ORION_data
  implicit none
  type(obj_IGLOO)      :: IGLOOsolver
  type(orion_data)     :: gas0, gasFast
  character(len=512)   :: gasfile
  integer              :: sweep, nargs

  nargs = command_argument_count()

  if (nargs == 0) then
    call IGLOOsolver%setup_static()
    do sweep = 1, 2
      write(*,'(/A,I0,A/)') ' ======== twosweep: SWEEP ', sweep, ' ========'
      call IGLOOsolver%reset_state()
      call IGLOOsolver%solve()
      call IGLOOsolver%writeout()
      if (sweep == 1) call park('OUTPUT.sweep1')
    enddo

  else
    call get_command_argument(1, gasfile)
    call read_TECsolfile(trim(gasfile), gas0)
    call read_TECsolfile(trim(gasfile), gasFast)
    call scaleU(gasFast, 2._8)

    call IGLOOsolver%setup_static(gas0)

    write(*,'(/A/)') ' ======== gas-cycle: SWEEP 1 (original field) ========'
    call IGLOOsolver%reset_state()
    call IGLOOsolver%solve();  call IGLOOsolver%writeout();  call park('OUTPUT.sweep1')

    write(*,'(/A/)') ' ======== gas-cycle: SWEEP 2 (U doubled) ========'
    call IGLOOsolver%reset_state(gasFast)
    call IGLOOsolver%solve();  call IGLOOsolver%writeout();  call park('OUTPUT.sweep2')

    write(*,'(/A/)') ' ======== gas-cycle: SWEEP 3 (original field again) ========'
    call IGLOOsolver%reset_state(gas0)
    call IGLOOsolver%solve();  call IGLOOsolver%writeout()
  endif

contains

  !> Park a finished OUTPUT/ under `dest` so the next sweep writes a fresh tree.
  subroutine park(dest)
    character(len=*), intent(in) :: dest
    integer :: st
    call execute_command_line('rm -rf '//dest//' && mv OUTPUT '//dest//' && mkdir -p OUTPUT', &
                              exitstat=st)
    if (st /= 0) error stop 'twosweep: could not park OUTPUT'
  end subroutine park

  !> Scale the streamwise gas velocity. varnames carries X,Y,Z first, so the variable at
  !  varnames(v) lives in vars(v-3) -- the same offset import_gas uses.
  subroutine scaleU(gas, factor)
    type(orion_data), intent(inout) :: gas
    real(8),          intent(in)    :: factor
    integer :: ib, v, nhit
    nhit = 0
    do ib = 1, size(gas%block)
      do v = 1, size(gas%varnames)
        if (trim(adjustl(gas%varnames(v))) == 'U') then
          gas%block(ib)%vars(v-3,:,:,:) = factor * gas%block(ib)%vars(v-3,:,:,:)
          nhit = nhit + 1
        endif
      enddo
    enddo
    if (nhit == 0) error stop 'twosweep: no U variable found to perturb'
    write(*,'(A,F5.2,A,I0,A)') ' [gas-cycle] scaled U by ', factor, ' in ', nhit, ' block(s)'
  end subroutine scaleU

end program twosweep
