program test_heat_probes
    !
    ! A7 RESOLVED 2026-07-15 as SOURCE-FAITHFUL (was an xfail bug probe):
    ! Shimada2006 eq. 45 (quoting NASA-SP-8039) IS  Nu = 2.5 Re^0.15 + 0.04 Re,
    ! with no +2 conduction floor — the code transcribes the source exactly.
    ! The missing stagnant-sphere limit is a property of the PUBLISHED
    ! correlation (documented limitation: prefer other models for Re -> 0).
    !
    !   XH1  JAXA1 transcription pin: Nu(Re) = 2.5 Re^0.15 + 0.04 Re exactly,
    !        including the as-published Nu(Re->0) -> 0 limit.
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64
    use IGLOO_Lib_Heat, only: heat
    implicit none

    real(R8), parameter :: Res(4) = [1.0e-6_R8, 1._R8, 100._R8, 1.0e4_R8]
    real(R8) :: nu, ref, err
    integer  :: i, nbad

    nbad = 0
    err  = 0._R8
    do i = 1, 4
        nu  = heat(Res(i), 0.7_R8, 0._R8, 1)     ! JAXA1
        ref = 2.5_R8*Res(i)**0.15_R8 + 0.04_R8*Res(i)
        err = max(err, abs(nu-ref)/ref)
    end do
    if (err > 1.0e-14_R8) then
        nbad = nbad + 1
        write(*,'(a,es12.4)') '  [XH1 FAIL] JAXA1 deviates from Shimada2006 eq.45, rel err =', err
    else
        write(*,'(a,es12.4)') '  [XH1 PASS] JAXA1 == 2.5 Re^0.15 + 0.04 Re, rel err =', err
    end if

    if (nbad > 0) then
        write(*,'(a)') 'test_heat_probes: FAIL'
        stop 1
    else
        write(*,'(a)') 'test_heat_probes: PASS (JAXA1 source-faithful; no +2 floor by design of the source)'
        stop 0
    end if

end program test_heat_probes
