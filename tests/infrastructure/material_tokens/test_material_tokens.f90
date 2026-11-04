program test_material_tokens
    !
    ! The material-line tokens of the phase file (key=value after "<name> <groups>", read by
    ! IGLOO_IO_INI::apply_material_key): six model words, the breakup word, fourteen reals and the
    ! twenty-three breakup constants; their defaults (set_material_defaults), the finalize rules
    ! (finalize_material_models) and the breakup constants' resolution (apply_breakup_tokens), driven on
    ! an obj_material built here: no input file.
    !
    !   MT1  the 14 reals after set_material_defaults + finalize, on a material whose 14 reals were
    !        poisoned first (a dropped default cannot hide behind the type's initializer)
    !   MT2  each real token alone: its field equals the literal, the other 13 keep their defaults
    !   MT3  each model word, every allowed value, from the opposite preset: its selector takes the
    !        value, the other five selectors, evapWord (but for evaporation) and the 14 reals stay
    !   MT4  finalize: T-nuc = 0.8 T-melt from T-melt alone; an explicit T-nuc kept; combustion with
    !        evaporation disables the evaporation, in either token order
    !   MT5  a repeated key: the last value wins
    !   MT6  each breakup constant token alone, on a material of each of the five models (its constants
    !        poisoned, bpMethod 2): on its own model exactly its slot takes the value (the two WeCrit
    !        doubled, TAB-method sets bpMethod), on the four others nothing moves; the 14 reals and the
    !        selectors keep theirs. The token table is this test's own, never Lib_Breakup's
    !   MT7  TAB's product-size scale: TAB-n=4 gives Gamma(1.5)/Gamma(1.75), with TAB-method=1 it is 1
    !
    ! Equality is exact: the reader's list-directed read and the compiler's literal both round the same
    ! decimal string correctly.
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64
    use IGLOO_data_phases, only: obj_material
    use IGLOO_IO_INI,      only: set_material_defaults, apply_material_key, finalize_material_models, &
                                 apply_breakup_tokens
    implicit none

    integer,  parameter :: NR = 14
    character(len=9),  parameter :: RKEY(NR) = [character(len=9) :: 'alpha-e', 'k-liq', 'mu-liq', 'K-burn', &
        'n-burn', 'X-eff', 'beta-part', 'xi-cap', 'T-ign', 'q-comb', 'T-melt', 'h-fus', 'T-nuc', 'cp-solid']
    character(len=12), parameter :: RVAL(NR) = [character(len=12) :: '0.75', '0.125', '0.0009765625', &
        '4.5e-07', '1.5', '0.5', '0.25', '0.375', '2100', '31000000', '2300', '1070000', '1900', '1300']
    real(R8), parameter :: RLIT(NR) = [0.75_R8, 0.125_R8, 0.0009765625_R8, 4.5e-07_R8, 1.5_R8, 0.5_R8, &
        0.25_R8, 0.375_R8, 2100._R8, 31000000._R8, 2300._R8, 1070000._R8, 1900._R8, 1300._R8]
    !> the defaults before finalize (T-nuc 0: finalize resolves it)
    real(R8), parameter :: DEF0(NR) = [1._R8, 0._R8, 0._R8, 0._R8, 1.8_R8, 1._R8, 0._R8, 0._R8, &
        2350._R8, 0._R8, 2327._R8, 0._R8, 0._R8, 0._R8]
    real(R8), parameter :: POISON = -999._R8

    !> model words: key, value, selector (1 evaporation, 2 liquid-conduction, 3 interface, 4 boiling,
    !  5 combustion, 6 solidification), expected value
    integer, parameter :: NW = 15
    character(len=17), parameter :: WKEY(NW) = [character(len=17) :: &
        'evaporation', 'evaporation', 'evaporation', 'evaporation', 'evaporation', 'evaporation', &
        'liquid-conduction', 'liquid-conduction', 'interface', 'interface', 'boiling', 'boiling', &
        'combustion', 'solidification', 'solidification']
    character(len=13), parameter :: WVAL(NW) = [character(len=13) :: 'NoEvaporation', 'd2-law', 'CEM', &
        'CEM-B', 'ASM', 'TC', 'ITC', 'P2T', 'VLE', 'LK', 'clamp', 'ZGR', 'Beckstead', 'on', 'off']
    integer, parameter :: WSEL(NW) = [1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 6, 6]
    integer, parameter :: WEXP(NW) = [0, 1, 2, 3, 4, 5, 0, 1, 0, 1, 0, 1, 1, 1, 0]
    integer, parameter :: SENT(6)  = [-11, -12, -13, -14, -15, -16]

    !> breakup constant tokens: key, model (1 Pilch-Erdman, 2 Reitz-Diawakar, 3 Reitz-KHRT, 4 TAB, 5 ETAB),
    !  slot of the model's constants (0: bpMethod), value, stored value (the two WeCrit doubled)
    integer, parameter :: NT = 23
    character(len=13), parameter :: TKEY(NT) = [character(len=13) :: 'PE-Cd', 'PE-B', 'RD-WeBag', 'RD-Cb', &
        'RD-Cstrip', 'RD-Cs', 'KHRT-B0', 'KHRT-B1', 'KHRT-Ctau', 'KHRT-CRT', 'KHRT-mShedLim', 'KHRT-WeLimit', &
        'TAB-Comega', 'TAB-Cmu', 'TAB-WeCrit', 'TAB-method', 'TAB-n', 'ETAB-k1', 'ETAB-k2', 'ETAB-WeCrit', &
        'ETAB-WeTrans', 'ETAB-Comega', 'ETAB-Cmu']
    integer, parameter :: TMOD(NT)  = [1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5]
    integer, parameter :: TSLOT(NT) = [1, 2, 1, 2, 3, 4, 1, 2, 3, 4, 5, 6, 1, 2, 3, 0, 4, 1, 2, 3, 4, 5, 6]
    character(len=6), parameter :: TVAL(NT) = [character(len=6) :: '1.25', '0.125', '7.5', '3.25', '0.625', &
        '8', '0.5', '40', '1.5', '0.25', '0.0625', '12', '10', '4.5', '7', '1', '2.75', '0.3', '0.15', '9', &
        '100', '6.5', '3.5']
    real(R8), parameter :: TEXP(NT) = [1.25_R8, 0.125_R8, 7.5_R8, 3.25_R8, 0.625_R8, 8._R8, 0.5_R8, 40._R8, &
        1.5_R8, 0.25_R8, 0.0625_R8, 12._R8, 10._R8, 4.5_R8, 14._R8, 1._R8, 2.75_R8, 0.3_R8, 0.15_R8, 18._R8, &
        100._R8, 6.5_R8, 3.5_R8]
    character(len=14), parameter :: MWORD(5) = [character(len=14) :: 'Pilch-Erdman', 'Reitz-Diawakar', &
        'Reitz-KHRT', 'TAB', 'ETAB']
    integer, parameter :: NBP(5) = [2, 4, 6, 4, 6]

    type(obj_material) :: mat
    real(R8) :: v(NR), ex(NR), tmelt, tnuc
    real(R8), allocatable :: bpex(:)
    integer  :: nbad, i, k, s, mex
    logical  :: ok, okAll

    nbad = 0

    ! ---- MT1: the defaults ---------------------------------------------------------------------
    call fresh(mat)
    call finalize_material_models(1, mat)
    tmelt = DEF0(11)
    tnuc  = 0.8_R8*tmelt
    ex = DEF0;  ex(13) = tnuc
    v  = reals(mat)
    call check(all(v == ex), 'MT1 the 14 reals after set_material_defaults + finalize', nbad)
    if (any(v /= ex)) call show(v, ex)

    ! ---- MT2: each real token alone -------------------------------------------------------------
    do i = 1, NR
        call fresh(mat)
        call apply_material_key(1, mat, trim(RKEY(i)), trim(RVAL(i)))
        v  = reals(mat)
        ex = DEF0;  ex(i) = RLIT(i)
        call check(all(v == ex), 'MT2 '//trim(RKEY(i))//'='//trim(RVAL(i))// &
                   ': the field takes the value, the other 13 keep their defaults', nbad)
        if (any(v /= ex)) call show(v, ex)
    end do

    ! ---- MT3: each model word from the opposite preset --------------------------------------------
    do i = 1, NW
        call fresh(mat)
        do k = 1, 6
            call setSel(mat, k, SENT(k))
        end do
        call setSel(mat, WSEL(i), opposite(WSEL(i), WEXP(i)))
        mat%evapWord = 'unset'
        call apply_material_key(1, mat, trim(WKEY(i)), trim(WVAL(i)))
        ok = getSel(mat, WSEL(i)) == WEXP(i)
        do k = 1, 6
            if (k /= WSEL(i)) ok = ok .and. getSel(mat, k) == SENT(k)
        end do
        if (WSEL(i) == 1) then
            ok = ok .and. trim(mat%evapWord) == trim(WVAL(i))
        else
            ok = ok .and. trim(mat%evapWord) == 'unset'
        end if
        ok = ok .and. all(reals(mat) == DEF0)
        call check(ok, 'MT3 '//trim(WKEY(i))//'='//trim(WVAL(i))//' -> '//itoa(WEXP(i)), nbad)
        if (.not. ok) write(*,'(a,6(1x,i0),2a)') '        selectors:', (getSel(mat, k), k = 1, 6), &
                                                  '  evapWord: ', trim(mat%evapWord)
    end do

    ! ---- MT4: finalize ----------------------------------------------------------------------------
    call fresh(mat)
    call apply_material_key(1, mat, 'T-melt', '2300')
    call finalize_material_models(1, mat)
    tmelt = RLIT(11)
    tnuc  = 0.8_R8*tmelt
    call check(mat%Tmelt == tmelt .and. mat%Tnuc == tnuc, 'MT4 T-melt=2300 alone: T-nuc = 0.8 T-melt', nbad)
    call fresh(mat)
    call apply_material_key(1, mat, 'T-nuc', '1900')
    call finalize_material_models(1, mat)
    call check(mat%Tnuc == RLIT(13) .and. mat%Tmelt == DEF0(11), 'MT4 T-nuc=1900 kept by finalize', nbad)
    call fresh(mat)
    call apply_material_key(1, mat, 'evaporation', 'CEM')
    call finalize_material_models(1, mat)
    call check(mat%evapSelect == 2 .and. mat%combSelect == 0, 'MT4 evaporation=CEM alone kept by finalize', nbad)
    call fresh(mat)
    call apply_material_key(1, mat, 'evaporation', 'CEM')
    call apply_material_key(1, mat, 'combustion', 'Beckstead')
    call finalize_material_models(1, mat)
    call check(mat%evapSelect == 0 .and. mat%combSelect == 1, &
               'MT4 evaporation=CEM then combustion=Beckstead: evaporation disabled', nbad)
    call fresh(mat)
    call apply_material_key(1, mat, 'combustion', 'Beckstead')
    call apply_material_key(1, mat, 'evaporation', 'CEM')
    call finalize_material_models(1, mat)
    call check(mat%evapSelect == 0 .and. mat%combSelect == 1, &
               'MT4 combustion=Beckstead then evaporation=CEM: evaporation disabled', nbad)

    ! ---- MT5: a repeated key ----------------------------------------------------------------------
    call fresh(mat)
    call apply_material_key(1, mat, 'alpha-e', '0.75')
    call apply_material_key(1, mat, 'alpha-e', '0.5')
    call check(mat%alphaE == 0.5_R8, 'MT5 alpha-e=0.75 then alpha-e=0.5: the last wins', nbad)
    call fresh(mat)
    call apply_material_key(1, mat, 'evaporation', 'CEM')
    call apply_material_key(1, mat, 'evaporation', 'ASM')
    call check(mat%evapSelect == 4 .and. trim(mat%evapWord) == 'ASM', &
               'MT5 evaporation=CEM then evaporation=ASM: the last wins', nbad)

    ! ---- MT6: each breakup constant token alone, on a material of each model -----------------------
    do i = 1, NT
        okAll = .true.
        do s = 1, 5
            call breaking(mat, s)
            call apply_material_key(1, mat, trim(TKEY(i)), trim(TVAL(i)))
            call apply_breakup_tokens(1, mat)
            bpex = spread(POISON, 1, NBP(s))
            mex  = 2
            if (TMOD(i) == s) then
                if (TSLOT(i) == 0) then
                    mex = nint(TEXP(i))
                else
                    bpex(TSLOT(i)) = TEXP(i)
                end if
            end if
            ok = size(mat%bp) == NBP(s)
            if (ok) ok = all(mat%bp == bpex)
            ok = ok .and. mat%bpMethod == mex .and. mat%brkupSelect == s .and. all(reals(mat) == DEF0)
            do k = 1, 6
                ok = ok .and. getSel(mat, k) == 0
            end do
            if (.not. ok) write(*,'(3a,*(1x,es12.4))') '        on ', trim(MWORD(s)), ': bp =', mat%bp
            okAll = okAll .and. ok
        end do
        call check(okAll, 'MT6 '//trim(TKEY(i))//'='//trim(TVAL(i))//': its slot of '//trim(MWORD(TMOD(i)))// &
                   ' only, the four other models unmoved', nbad)
    end do

    ! ---- MT7: TAB's product-size scale ------------------------------------------------------------
    call breaking(mat, 4)
    mat%bpScale = POISON
    call apply_material_key(1, mat, 'TAB-n', '4')
    call apply_breakup_tokens(1, mat)
    call check(abs(mat%bpScale - gamma(1.5_R8)/gamma(1.75_R8)) <= 1e-15_R8, &
               'MT7 TAB-n=4 under method 2: bpScale = Gamma(1.5)/Gamma(1.75)', nbad)
    call breaking(mat, 4)
    mat%bpScale = POISON
    call apply_material_key(1, mat, 'TAB-n', '4')
    call apply_material_key(1, mat, 'TAB-method', '1')
    call apply_breakup_tokens(1, mat)
    call check(mat%bpScale == 1._R8 .and. mat%bpMethod == 1, 'MT7 TAB-n=4 with TAB-method=1: bpScale = 1', nbad)

    write(*,*)
    if (nbad == 0) then
        write(*,'(a)') '[PASS] test_material_tokens: every token sets its own field, defaults and finalize as pinned'
    else
        write(*,'(a,i0,a)') '[FAIL] test_material_tokens: ', nbad, ' check(s) failed'
        error stop 1
    end if

contains

    !> A material as the reader starts one: the 14 reals poisoned, the selectors at their type defaults,
    !  then set_material_defaults.
    subroutine fresh(m)
        type(obj_material), intent(inout) :: m
        call setReals(m, spread(POISON, 1, NR))
        m%evapSelect = 0;  m%liqSelect = 0;  m%intfSelect = 0
        m%boilSelect = 0;  m%combSelect = 0; m%solidSelect = 0
        m%evapWord = ' '
        m%brkupSelect = 0;  m%brkupWord = ' '
        m%bpGiven = .false.;  m%bpTok = 0._R8
        if (allocated(m%bp)) deallocate(m%bp)
        m%bpMethod = 1;  m%bpScale = 1._R8
        call set_material_defaults(m)
    end subroutine fresh

    !> A fresh material of breakup model s (by its token), its constants poisoned and bpMethod 2.
    subroutine breaking(m, s)
        type(obj_material), intent(inout) :: m
        integer,            intent(in)    :: s
        call fresh(m)
        call apply_material_key(1, m, 'breakup', trim(MWORD(s)))
        allocate(m%bp(NBP(s)))
        m%bp = POISON
        m%bpMethod = 2
    end subroutine breaking

    !> The 14 reals in the reader's token order.
    function reals(m) result(r)
        type(obj_material), intent(in) :: m
        real(R8) :: r(NR)
        r = [m%alphaE, m%kLiq, m%muLiq, m%Kburn, m%nBurn, m%Xeff, m%betaPart, m%xiCap, &
             m%Tign, m%qComb, m%Tmelt, m%hFus, m%Tnuc, m%cpSol]
    end function reals

    subroutine setReals(m, r)
        type(obj_material), intent(inout) :: m
        real(R8),           intent(in)    :: r(NR)
        m%alphaE = r(1);  m%kLiq  = r(2);  m%muLiq    = r(3);  m%Kburn = r(4);  m%nBurn = r(5)
        m%Xeff   = r(6);  m%betaPart = r(7);  m%xiCap = r(8);  m%Tign  = r(9);  m%qComb = r(10)
        m%Tmelt  = r(11); m%hFus  = r(12); m%Tnuc     = r(13); m%cpSol = r(14)
    end subroutine setReals

    integer function getSel(m, k)
        type(obj_material), intent(in) :: m
        integer,            intent(in) :: k
        select case (k)
        case (1); getSel = m%evapSelect
        case (2); getSel = m%liqSelect
        case (3); getSel = m%intfSelect
        case (4); getSel = m%boilSelect
        case (5); getSel = m%combSelect
        case default; getSel = m%solidSelect
        end select
    end function getSel

    subroutine setSel(m, k, val)
        type(obj_material), intent(inout) :: m
        integer,            intent(in)    :: k, val
        select case (k)
        case (1); m%evapSelect  = val
        case (2); m%liqSelect   = val
        case (3); m%intfSelect  = val
        case (4); m%boilSelect  = val
        case (5); m%combSelect  = val
        case default; m%solidSelect = val
        end select
    end subroutine setSel

    !> A valid selector value other than the expected one (evaporation 0..5; combustion 0 = off).
    integer function opposite(k, expected)
        integer, intent(in) :: k, expected
        select case (k)
        case (1); opposite = mod(expected + 3, 6)
        case (5); opposite = 0
        case default; opposite = 1 - expected
        end select
    end function opposite

    function itoa(n) result(s)
        integer, intent(in) :: n
        character(len=:), allocatable :: s
        character(len=12) :: buf
        write(buf,'(i0)') n
        s = trim(buf)
    end function itoa

    subroutine show(got, want)
        real(R8), intent(in) :: got(NR), want(NR)
        integer :: j
        do j = 1, NR
            if (got(j) /= want(j)) write(*,'(a,a9,a,es24.16,a,es24.16)') '        ', RKEY(j), ': ', got(j), &
                                                                      '  expected ', want(j)
        end do
    end subroutine show

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

end program test_material_tokens
