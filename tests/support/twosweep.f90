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
program twosweep
  use IGLOO_module, only: obj_IGLOO
  implicit none
  type(obj_IGLOO) :: IGLOOsolver
  integer :: sweep

  call IGLOOsolver%setup_static()

  do sweep = 1, 2
    write(*,'(/A,I0,A/)') ' ======== twosweep: SWEEP ', sweep, ' ========'
    call IGLOOsolver%reset_state()
    call IGLOOsolver%solve()
    call IGLOOsolver%writeout()
    if (sweep == 1) call parkSweep1()
  enddo

contains

  !> Park sweep 1's OUTPUT/ so sweep 2 writes a fresh tree to compare against.
  subroutine parkSweep1()
    integer :: st
    call execute_command_line( &
      'rm -rf OUTPUT.sweep1 && mv OUTPUT OUTPUT.sweep1 && mkdir -p OUTPUT', exitstat=st)
    if (st /= 0) error stop 'twosweep: could not park sweep-1 OUTPUT'
  end subroutine parkSweep1

end program twosweep
