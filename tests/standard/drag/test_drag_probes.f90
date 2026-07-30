program test_drag_probes
    !
    ! Bug-transcription pins. Written as xfail probes (exit 1 while a bug is live),
    ! but A1/A5/A6 + the Wen-Yu C-flag are all FIXED, so this exits 0 and ctest runs
    ! it as an ORDINARY gate: no WILL_FAIL property is set, and a regression turns it
    ! red the normal way.
    !
    !   XD1  Crowe finite value  (bug A1): exp coded in the denominator ->
    !        Cd ~ hfun/toll ~ 1e20 at Re=1000, Ma=2 (literature: Cd = O(1))
    !   XD2  Hermsen finite value (A1, same coded term)
    !   XD3  Putnam plateau = 0.4392 SOURCE-FAITHFUL (Shimada2006 eq. 17; user
    !        decision 2026-07-16: originals unavailable, Shimada authoritative).
    !        The 3.6% C0 jump at Re=1000 is a property of the published form.
    !   XD4  Henderson C0 continuity at Ma=1.75 (A6): blend slope 4/3 kept --
    !        Shimada eq. 22's 3/4 cannot reach the M>=1.75 branch (typo)
    !   XD5  Wen-Yu high-Re plateau = 0.43 SOURCE-FAITHFUL (Shimada eq. 16; same
    !        decision; the ~2% step at Re=1000 is inherent to the published form)
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64
    use IGLOO_Lib_Drag, only: drag
    implicit none

    integer, parameter :: WENYU=6, PUTNAM=7, HENDERSON=11, CROWE=12, HERMSEN=13
    real(R8), parameter :: G = 1.4_R8, Tr = 1._R8, eps = 1.0e-9_R8
    real(R8) :: cd, clo, chi, jump
    integer  :: nbad

    nbad = 0

    cd = drag(1000._R8, 2._R8, G, Tr, CROWE)
    call probe('XD1 Crowe Cd(Re=1e3,Ma=2)', cd, cd < 10._R8, nbad)

    cd = drag(1000._R8, 2._R8, G, Tr, HERMSEN)
    call probe('XD2 Hermsen Cd(Re=1e3,Ma=2)', cd, cd < 10._R8, nbad)

    chi = drag(1000._R8*(1._R8+eps), 0._R8, G, Tr, PUTNAM)
    call probe('XD3 Putnam plateau=0.4392 @Re>=1000', chi, &
               abs(chi-0.4392_R8) < 1.0e-12_R8, nbad)

    clo = drag(100._R8, 1.75_R8*(1._R8-eps), G, Tr, HENDERSON)
    chi = drag(100._R8, 1.75_R8*(1._R8+eps), G, Tr, HENDERSON)
    jump = abs(chi-clo)/cli_safe(clo)
    call probe('XD4 Henderson jump @Ma=1.75', jump, jump < 1.0e-6_R8, nbad)

    clo = drag(1000._R8*(1._R8-eps), 0._R8, G, Tr, WENYU)
    chi = drag(1000._R8*(1._R8+eps), 0._R8, G, Tr, WENYU)
    jump = abs(chi-clo)/clo
    call probe('XD5 Wen-Yu plateau=0.43 @Re>1000', chi, &
               abs(chi-0.43_R8) < 1.0e-12_R8 .and. jump < 2.0e-2_R8, nbad)

    if (nbad > 0) then
        write(*,'(a,i0,a)') 'test_drag_probes: ', nbad, ' known bug(s) still present (expected)'
        stop 1
    else
        write(*,'(a)') 'test_drag_probes: all known bugs FIXED — promote probes to the gate'
        stop 0
    end if

contains

    pure function cli_safe(x) result(y)
        real(R8), intent(in) :: x
        real(R8) :: y
        y = max(abs(x), 1.0e-30_R8)
    end function cli_safe

    subroutine probe(name, val, fixed, nbad)
        character(len=*), intent(in)    :: name
        real(R8),         intent(in)    :: val
        logical,          intent(in)    :: fixed
        integer,          intent(inout) :: nbad
        if (fixed) then
            write(*,'(a,a,es12.4)') '  [FIXED]        ', name//' =', val
        else
            nbad = nbad + 1
            write(*,'(a,a,es12.4)') '  [STILL-BROKEN] ', name//' =', val
        end if
    end subroutine probe

end program test_drag_probes
