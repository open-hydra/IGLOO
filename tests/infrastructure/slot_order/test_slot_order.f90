program test_slot_order
    !
    ! The packed property slots of the right-hand sides. setupRHS lays out the auxiliary array
    ! (Lib_RHS ind_*), packAuxVars writes each per-material property at its slot, and the readers
    ! take two contiguous slices by their OWN named indices: the evaporation parameters
    ! aux(ind_Mv : ind_Mv+nep-1) as ep(iMv..imuLiq) of IGLOO_Lib_Evaporation (rhsEvaporation,
    ! rhsEvapBreakup), and the metal block aux(ind_mb : ind_mb+nmp-1) as mp(imKb..imQc) of
    ! IGLOO_Lib_Combustion (rhsAlCombustion; imHf and imCps in rhsSolidification). Every property gets
    ! a distinct value, and every named read must return the property of that name.
    !
    !   S1  model 2 (CEM): the ten ep slots; nep = 10; the slice inside aux; ind_Lv.. agree with ep
    !   S2  model 4 (CEM + Reitz-Diwakar): the same ep slots
    !   S3  model 5 (Beckstead): the eleven metal slots; nmp = nmetal = 11; the slice inside aux
    !   S4  model 6 (solidification): the same metal block, h-fus and cp-solid where its RHS reads them
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64
    use IGLOO_particles,       only: obj_particle
    use Lib_RHS,               only: setupRHS, packAuxVars, nauxvar, nmetal, ind_mb, &
                                     ind_Mv, ind_Lv, ind_Cv, ind_Le, ind_Yi, ind_e1, ind_e2, ind_e3, ind_e4, ind_e5
    use IGLOO_Lib_Evaporation, only: nep, iMv, iLv, icpv, iLe, iYinf, iLvMvOverRu, iinvTboil, &
                                     ialphaE, ikLiq, imuLiq
    use IGLOO_Lib_Combustion,  only: nmp, imKb, imN, imXe, imBeta, imXi, imTig, imTm, imHf, imTn, imCps, imQc
    implicit none

    type(obj_particle)    :: p
    real(R8), allocatable :: aux(:)
    logical  :: propFlags(5)
    integer  :: nAux, nAuxSt, nEvV, nbad

    nbad = 0
    propFlags = .false.
    !> distinct values, one per property
    p%cp = 1.0_R8;  p%rho = 2.0_R8;  p%m = 3.0_R8;  p%d = 4.0_R8;  p%mdot = 5.0_R8
    p%sigma = 6.0_R8;  p%mup = 7.0_R8
    p%Mv = 101._R8;  p%Lv = 102._R8;  p%cpv = 103._R8;  p%Le = 104._R8;  p%Yinf = 105._R8
    p%LvMvOverRu = 106._R8;  p%invTboil = 107._R8;  p%alphaE = 108._R8;  p%kLiq = 109._R8
    p%muLiq = 110._R8
    p%Kburn = 201._R8;  p%nBurn = 202._R8;  p%Xeff = 203._R8;  p%betaPart = 204._R8
    p%xiCap = 205._R8;  p%Tign = 206._R8;  p%Tmelt = 207._R8;  p%hFus = 208._R8
    p%Tnuc = 209._R8;  p%cpSol = 210._R8;  p%qComb = 211._R8

    ! ---- S1, S2: the evaporation slice --------------------------------------------------------
    call setupRHS(2, 0, 2, 0, 0, 0, 0, 0, propFlags, [real(R8)::], 0, 0._R8, nAux, nAuxSt, nEvV)
    call pack()
    call check(nep == 10, 'S1 nep = 10 evaporation slots', nbad)
    call checkEvap('S1 model 2')
    call setupRHS(4, 2, 2, 0, 0, 0, 0, 0, propFlags, [0._R8], 1, 1._R8, nAux, nAuxSt, nEvV)
    call pack()
    call checkEvap('S2 model 4')

    ! ---- S3, S4: the metal block ----------------------------------------------------------------
    call setupRHS(5, 0, 0, 0, 0, 0, 1, 0, propFlags, [real(R8)::], 0, 0._R8, nAux, nAuxSt, nEvV)
    call pack()
    call check(nmp == nmetal .and. nmp == 11, 'S3 nmp = nmetal = 11 metal slots', nbad)
    call checkMetal('S3 model 5')
    call setupRHS(6, 0, 0, 0, 0, 0, 0, 1, propFlags, [real(R8)::], 0, 0._R8, nAux, nAuxSt, nEvV)
    call pack()
    call checkMetal('S4 model 6')
    call check(aux(ind_mb+imHf-1) == p%hFus .and. aux(ind_mb+imCps-1) == p%cpSol, &
               'S4 model 6: h-fus and cp-solid where rhsSolidification reads them', nbad)

    write(*,*)
    if (nbad == 0) then
        write(*,'(A)') '[PASS] test_slot_order: every property is read by its own name'
    else
        write(*,'(A,I0,A)') '[FAIL] test_slot_order: ', nbad, ' check(s) failed'
        error stop 1
    endif

contains

    subroutine pack()
        if (allocated(aux)) deallocate(aux)
        allocate(aux(nauxvar))
        call packAuxVars(p, nauxvar, aux)
    end subroutine pack

    !> ep(name) = aux(ind_Mv + name - 1) for the ten names; the slice lies inside aux; the
    !  separately named ind_Lv .. ind_e5 point at the same slots.
    subroutine checkEvap(tag)
        character(len=*), intent(in) :: tag
        real(R8) :: ep(nep)
        call check(ind_Mv > 0 .and. ind_Mv + nep - 1 <= nauxvar, tag//': the ep slice lies inside aux', nbad)
        if (.not. (ind_Mv > 0 .and. ind_Mv + nep - 1 <= nauxvar)) return
        ep = aux(ind_Mv:ind_Mv+nep-1)
        call check(ep(iMv) == p%Mv .and. ep(iLv) == p%Lv .and. ep(icpv) == p%cpv .and. ep(iLe) == p%Le &
                   .and. ep(iYinf) == p%Yinf, tag//': ep Mv, Lv, cpv, Le, Yinf', nbad)
        call check(ep(iLvMvOverRu) == p%LvMvOverRu .and. ep(iinvTboil) == p%invTboil .and. &
                   ep(ialphaE) == p%alphaE .and. ep(ikLiq) == p%kLiq .and. ep(imuLiq) == p%muLiq, &
                   tag//': ep LvMv/Ru, 1/Tboil, alpha-e, k-liq, mu-liq', nbad)
        call check(all([ind_Lv, ind_Cv, ind_Le, ind_Yi, ind_e1, ind_e2, ind_e3, ind_e4, ind_e5] - ind_Mv + 1 == &
                       [iLv, icpv, iLe, iYinf, iLvMvOverRu, iinvTboil, ialphaE, ikLiq, imuLiq]), &
                   tag//': ind_Lv..ind_e5 name the same slots as ep', nbad)
        if (nbad > 0) write(*,'(A,10F7.1)') '        ep = ', ep
    end subroutine checkEvap

    !> mp(name) = aux(ind_mb + name - 1) for the eleven names; the slice lies inside aux.
    subroutine checkMetal(tag)
        character(len=*), intent(in) :: tag
        real(R8) :: mp(nmp)
        call check(ind_mb > 0 .and. ind_mb + nmp - 1 <= nauxvar, tag//': the metal slice lies inside aux', nbad)
        if (.not. (ind_mb > 0 .and. ind_mb + nmp - 1 <= nauxvar)) return
        mp = aux(ind_mb:ind_mb+nmp-1)
        call check(mp(imKb) == p%Kburn .and. mp(imN) == p%nBurn .and. mp(imXe) == p%Xeff .and. &
                   mp(imBeta) == p%betaPart .and. mp(imXi) == p%xiCap .and. mp(imQc) == p%qComb, &
                   tag//': mp K-burn, n-burn, X-eff, beta-part, xi-cap, q-comb', nbad)
        call check(mp(imTig) == p%Tign .and. mp(imTm) == p%Tmelt .and. mp(imHf) == p%hFus .and. &
                   mp(imTn) == p%Tnuc .and. mp(imCps) == p%cpSol, &
                   tag//': mp T-ign, T-melt, h-fus, T-nuc, cp-solid', nbad)
        if (nbad > 0) write(*,'(A,11F7.1)') '        mp = ', mp
    end subroutine checkMetal

    subroutine check(cond, name, nbad)
        logical,          intent(in)    :: cond
        character(len=*), intent(in)    :: name
        integer,          intent(inout) :: nbad
        if (cond) then
            write(*,'(a,a)') '  [PASS] ', name
        else
            write(*,'(a,a)') '  [FAIL] ', name
            nbad = nbad + 1
        end if
    end subroutine check

end program test_slot_order
