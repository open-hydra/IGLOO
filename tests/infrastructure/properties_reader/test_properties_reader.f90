program test_properties_reader
    !
    ! The production reader of INPUT/<prefix>properties.dat (read_cdp_properties) on the fixtures
    ! of the working directory; the first argument names the directory's leg set:
    !
    !   hexadecane    PR1 legacy- (T 280..620) against legacy1- (T 1..620): same cp, rho(T), h(T)
    !                 PR2 psat-: the Psat column tabulated on [Tmin, Tmax], psat(560 K) = 1 atm
    !                 PR3 perm-: columns found by name, in any order
    !                 PR6 legacy-: relative enthalpy datum, hOff = 0
    !                 PR7 abs-: "Enthalpy_abs" is the enthalpy column, with its datum
    !                 PR4 validate_psat_column, PR5 check_table_nodes, PR8 classify_table_tokens,
    !                 PR11 check_table_columns, PR12 tableValue, PR13 scan_rows (scan-*.dat)
    !                 PR14 comp_TfromTab: the inverse of lookupTab on tables from 250 K, 300 K and 1 K
    !                 PR15 lookupTab at a non-finite temperature: NaN, and no invalid flag
    !                 PR-V2 validate_vapour_columns: each code once
    !   two-material  PR9 nm-: zone A tabulated, zone B all zero (Clausius-Clapeyron)
    !   atlas-water   PR10 atlas-: the table ATLAS GPB writes for water with psat-vapour = H2O
    !                 PR-V1 atlas-: its Lv, Mv, Cpv columns by name, Lv(T) and cpv(T) through tableValue
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64
    use IGLOO_IO_INI,         only: read_IGLOO_input
    use IGLOO_IO,             only: read_cdp_properties
    use IGLOO_data_phases,    only: obj_material
    use IGLOO_Lib_Properties, only: Tmin, Tmax, tableValue, lookupTab, comp_TfromTab, &
                                    validate_psat_column, check_table_nodes, &
                                    classify_table_tokens, check_table_columns, scan_rows, nzone_max, &
                                    TAB_OK, TAB_NO_TEMPERATURE, TAB_NO_CP, TAB_NO_DENSITY, TAB_NO_ENTHALPY, &
                                    TAB_TWO_ENTHALPY, TAB_DUPLICATE, TAB_FEW_ROWS, TAB_OFF_NODE, TAB_NONFINITE, &
                                    TAB_RHO_NONPOSITIVE, TAB_CP_NONPOSITIVE, TAB_H_NONMONOTONE, &
                                    TAB_H_CP_MISMATCH, TAB_DATUM_MISMATCH, TAB_NEGATIVE_T, &
                                    PSAT_OK, PSAT_ABSENT, PSAT_NONFINITE, PSAT_NEGATIVE, PSAT_NONMONOTONE, &
                                    PSAT_CONSTANT, PSAT_TBOIL_RANGE, PSAT_TBOIL_MISMATCH, &
                                    validate_vapour_columns, VAP_OK, VAP_ABSENT, VAP_NONFINITE, VAP_PARTIAL, &
                                    VAP_NONPOSITIVE, VAP_MV_VARIES, VAP_WITHOUT_PSAT, VAP_CC_MISMATCH
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan, ieee_positive_inf, &
                                             ieee_negative_inf
    use, intrinsic :: ieee_exceptions, only: ieee_set_flag, ieee_get_flag, ieee_invalid
    use IGLOO_variables,      only: llen, nm
    implicit none

    real(R8), parameter :: Patm = 101325._R8
    character(len=2)      :: method
    logical               :: srcSwitch, eulSwitch
    character(len=llen)   :: gasfile
    real(R8), allocatable :: pos0(:,:), vel0(:,:), mdot(:), diam(:), temp0(:)
    character(len=32)     :: mode
    integer               :: nbad

    call get_command_argument(1, mode)
    call read_IGLOO_input(method, srcSwitch, eulSwitch, gasfile, pos0, vel0, temp0, mdot, diam)
    nbad = 0

    select case (trim(mode))
    case ('hexadecane')
        call hexadecane_legs(nbad)
        call function_legs(nbad)
        call inversion_legs(nbad)
        call nonfinite_legs(nbad)
    case ('two-material')
        call two_material_legs(nbad)
    case ('atlas-water')
        call atlas_water_legs(nbad)
    case default
        write(*,'(a)') 'test_properties_reader: unknown mode "'//trim(mode)//'"'
        nbad = 1
    end select

    if (nbad == 0) then
        write(*,'(a)') 'test_properties_reader '//trim(mode)//': OVERALL PASS'; stop 0
    else
        write(*,'(a,i0,a)') 'test_properties_reader '//trim(mode)//': OVERALL FAIL (', nbad, ' assertion(s))'
        stop 1
    end if

contains

    subroutine hexadecane_legs(nbad)
        integer, intent(inout) :: nbad
        type(obj_material), allocatable :: mat(:)
        real(R8), allocatable :: rho1(:), h1(:), rhoL(:), hL(:), psatP(:), rhoP(:), hP(:)
        character(len=32) :: prefix

        ! regression half: the table from 1 K, as every e2e fixture
        prefix = 'legacy1-'
        call read_cdp_properties(prefix, mat)
        call check(Tmin == 1 .and. Tmax == 620,                  'PR1 legacy1-: Tmin = 1, Tmax = 620', nbad)
        call check(.not. mat(1)%cpVariable .and. mat(1)%cp == 2800._R8, 'PR1 legacy1-: constant cp = 2800', nbad)
        call check(mat(1)%rhoVariable,                           'PR1 legacy1-: rho varies', nbad)
        call check(mat(1)%hDatum == 'relative' .and. mat(1)%hOff == 0._R8, 'PR6 legacy1-: relative datum, hOff = 0', nbad)
        call check(.not. mat(1)%psatVariable .and. .not. allocated(mat(1)%psatTab), &
                   'PR1 legacy1-: no Psat column, no psat table', nbad)
        if (allocated(mat(1)%rhoTab)) rho1 = mat(1)%rhoTab
        h1 = mat(1)%hTab

        ! PR1 / PR6: the same rows from 280 K
        prefix = 'legacy-'
        call read_cdp_properties(prefix, mat)
        call check(Tmin == 280 .and. Tmax == 620,                'PR1 legacy-: Tmin = 280, Tmax = 620', nbad)
        call check(.not. mat(1)%cpVariable .and. mat(1)%cp == 2800._R8, 'PR1 legacy-: constant cp = 2800', nbad)
        call check(mat(1)%rhoVariable,                           'PR1 legacy-: rho varies', nbad)
        if (mat(1)%rhoVariable .and. allocated(rho1)) then
            call check(lbound(mat(1)%rhoTab,1) == 280 .and. ubound(mat(1)%rhoTab,1) == 620, &
                       'PR1 legacy-: rho table on [280, 620]', nbad)
            call check(all(mat(1)%rhoTab == rho1(280:620)),      'PR1 legacy-: rho(T) = legacy1- rho(T), bitwise', nbad)
        endif
        call check(all(mat(1)%hTab == h1(280:620)),              'PR6 legacy-: h(T) = legacy1- h(T), bitwise', nbad)
        call check(mat(1)%hDatum == 'relative' .and. mat(1)%hOff == 0._R8, 'PR6 legacy-: relative datum, hOff = 0', nbad)
        call check(.not. mat(1)%psatVariable .and. .not. allocated(mat(1)%psatTab), &
                   'PR1 legacy-: no Psat column, no psat table', nbad)
        if (allocated(mat(1)%rhoTab)) rhoL = mat(1)%rhoTab
        hL = mat(1)%hTab

        ! PR7: the absolute datum
        prefix = 'abs-'
        call read_cdp_properties(prefix, mat)
        call check(.not. mat(1)%cpVariable .and. mat(1)%cp == 2800._R8, 'PR7 abs-: constant cp = 2800', nbad)
        call check(mat(1)%hDatum == 'absolute',                  'PR7 abs-: "Enthalpy_abs" read as the absolute datum', nbad)
        call check(mat(1)%hOff == -2.0e6_R8,                     'PR7 abs-: hOff = -2.0e6 J/kg', nbad)
        call check(mat(1)%hTab(280) == -1216000._R8,             'PR7 abs-: h(280 K) = cp*T + hOff', nbad)

        ! PR2: the Psat column
        prefix = 'psat-'
        call read_cdp_properties(prefix, mat)
        call check(mat(1)%psatVariable,                          'PR2 psat-: psat tabulated', nbad)
        if (allocated(mat(1)%psatTab)) then
            call check(lbound(mat(1)%psatTab,1) == 280 .and. ubound(mat(1)%psatTab,1) == 620, &
                       'PR2 psat-: psat table on [280, 620]', nbad)
            call check(abs(mat(1)%psatTab(560) - Patm) <= 1.e-9_R8*Patm, 'PR2 psat-: psat(560 K) = 1 atm', nbad)
            call check(mat(1)%psatTab(474) == 1.017887e4_R8,     'PR2 psat-: psat(474 K) = the file''s value', nbad)
            psatP = mat(1)%psatTab
        else
            call check(.false.,                                  'PR2 psat-: psat table allocated', nbad)
        endif
        call check(.not. mat(1)%cpVariable .and. mat(1)%cp == 2800._R8, 'PR2 psat-: cp unchanged', nbad)
        if (mat(1)%rhoVariable .and. allocated(rhoL)) &
            call check(all(mat(1)%rhoTab == rhoL),               'PR2 psat-: rho(T) unchanged, bitwise', nbad)
        call check(all(mat(1)%hTab == hL),                       'PR2 psat-: h(T) unchanged, bitwise', nbad)
        if (allocated(mat(1)%rhoTab)) rhoP = mat(1)%rhoTab
        hP = mat(1)%hTab

        ! PR3: the same columns in another order
        prefix = 'perm-'
        call read_cdp_properties(prefix, mat)
        call check(mat(1)%psatVariable .and. allocated(psatP),   'PR3 perm-: psat tabulated', nbad)
        if (mat(1)%psatVariable .and. allocated(psatP)) &
            call check(all(mat(1)%psatTab == psatP),             'PR3 perm-: psat(T) = psat-, bitwise', nbad)
        call check(.not. mat(1)%cpVariable .and. mat(1)%cp == 2800._R8, 'PR3 perm-: cp found by name', nbad)
        call check(mat(1)%rhoVariable,                           'PR3 perm-: rho found by name', nbad)
        if (mat(1)%rhoVariable .and. allocated(rhoP)) &
            call check(all(mat(1)%rhoTab == rhoP),               'PR3 perm-: rho(T) = psat-, bitwise', nbad)
        call check(all(mat(1)%hTab == hP),                       'PR3 perm-: h(T) = psat-, bitwise', nbad)
        call check(mat(1)%hDatum == 'relative' .and. mat(1)%hOff == 0._R8, 'PR3 perm-: relative datum, hOff = 0', nbad)
    end subroutine hexadecane_legs


    subroutine two_material_legs(nbad)
        integer, intent(inout) :: nbad
        type(obj_material), allocatable :: mat(:)
        character(len=32) :: prefix

        prefix = 'nm-'
        call read_cdp_properties(prefix, mat)
        call check(nm == 2 .and. size(mat) == 2,                 'PR9 nm-: two materials', nbad)
        call check(Tmin == 1 .and. Tmax == 620,                  'PR9 nm-: Tmin = 1, Tmax = 620', nbad)
        call check(mat(1)%psatVariable,                          'PR9 nm-: zone A psat tabulated', nbad)
        if (allocated(mat(1)%psatTab)) then
            call check(abs(mat(1)%psatTab(560) - Patm) <= 1.e-9_R8*Patm, 'PR9 nm-: zone A psat(560 K) = 1 atm', nbad)
            call check(mat(1)%psatTab(1) == 0._R8,               'PR9 nm-: zone A psat(1 K) underflows to 0', nbad)
        endif
        call check(.not. mat(2)%psatVariable .and. .not. allocated(mat(2)%psatTab), &
                   'PR9 nm-: zone B all zero, Clausius-Clapeyron', nbad)
        call check(mat(1)%cp == 2800._R8 .and. mat(2)%cp == 2800._R8, 'PR9 nm-: cp of both zones', nbad)
        if (mat(1)%rhoVariable .and. mat(2)%rhoVariable) then
            call check(all(mat(1)%rhoTab == mat(2)%rhoTab),      'PR9 nm-: rho(T) of both zones, bitwise', nbad)
        else
            call check(.false.,                                  'PR9 nm-: rho varies in both zones', nbad)
        endif
    end subroutine two_material_legs


    subroutine atlas_water_legs(nbad)
        integer, intent(inout) :: nbad
        type(obj_material), allocatable :: mat(:)
        character(len=32) :: prefix

        prefix = 'atlas-'
        call read_cdp_properties(prefix, mat)
        call check(trim(mat(1)%matName) == 'H2O(L)',             'PR10 atlas-: material H2O(L)', nbad)
        call check(Tmin == 280 .and. Tmax == 380,                'PR10 atlas-: Tmin = 280, Tmax = 380', nbad)
        call check(.not. mat(1)%cpVariable .and. mat(1)%cp == 4184._R8, 'PR10 atlas-: constant cp = 4184', nbad)
        call check(.not. mat(1)%rhoVariable .and. mat(1)%rho == 997._R8, 'PR10 atlas-: constant rho = 997', nbad)
        call check(mat(1)%hDatum == 'absolute',                  'PR10 atlas-: "Enthalpy_abs" read as the absolute datum', nbad)
        call check(abs(mat(1)%hOff + 17112459.6_R8) <= 1.e-9_R8*17112459.6_R8, &
                   'PR10 atlas-: hOff = h(280 K) - cp*280 = -17112459.6 J/kg', nbad)
        call check(mat(1)%psatVariable,                          'PR10 atlas-: psat tabulated', nbad)
        if (allocated(mat(1)%psatTab)) then
            call check(lbound(mat(1)%psatTab,1) == 280 .and. ubound(mat(1)%psatTab,1) == 380, &
                       'PR10 atlas-: psat table on [280, 380]', nbad)
            call check(mat(1)%psatTab(300) == 3.533623e3_R8,     'PR10 atlas-: psat(300 K) = the file''s 3.533623e+03', nbad)
            call check(mat(1)%psatTab(373) == 9.930803e4_R8 .and. mat(1)%psatTab(374) == 1.028750e5_R8, &
                       'PR10 atlas-: psat(373 K), psat(374 K) = the file''s values', nbad)
            call check(tableValue(mat(1)%psatTab, Tmin, Tmax, 275._R8) == mat(1)%psatTab(280) .and. &
                       tableValue(mat(1)%psatTab, Tmin, Tmax, 385._R8) == mat(1)%psatTab(380), &
                       'PR10 atlas-: psat 5 K outside the table = the end values', nbad)
        endif
        call atlas_vapour_legs(mat, nbad)
    end subroutine atlas_water_legs


    !> PR-V1: Lv, Mv, Cpv found by name; Lv and cpv through tableValue at six temperatures (below the table,
    !  the datum temperature, a node, between nodes, the boiling temperature, above) equal to the file's
    !  own numbers by the same formula (exact on a node and outside, 4 eps between nodes); Lv(298.15 K)
    !  the material's scalar latent heat.
    subroutine atlas_vapour_legs(mat, nbad)
        type(obj_material), intent(in)    :: mat(:)
        integer,            intent(inout) :: nbad
        real(R8), parameter :: Tq(6) = [250._R8, 298.15_R8, 300._R8, 337.5_R8, 373.15_R8, 400._R8]
        logical,  parameter :: onNode(6) = [.true., .false., .true., .false., .false., .true.]
        real(R8) :: T(200), Lvf(200), Mvf(200), Cpvf(200), row(8), eLv, eCpv, tol
        integer  :: u, ios, n, k
        character(len=512) :: line

        call check(mat(1)%LvVariable .and. mat(1)%cpvVariable,   'PR-V1 atlas-: Lv and Cpv tabulated', nbad)
        call check(mat(1)%Mv == 18.015_R8,                       'PR-V1 atlas-: Mv = 18.015 from its column', nbad)
        if (.not. (allocated(mat(1)%LvTab) .and. allocated(mat(1)%cpvTab))) then
            call check(.false.,                                  'PR-V1 atlas-: Lv and Cpv tables allocated', nbad)
            return
        endif
        open(newunit=u, file='INPUT/atlas-properties.dat', status='old', action='read')
        n = 0
        do
            read(u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            read(line, *, iostat=ios) row
            if (ios /= 0) cycle
            n = n + 1
            T(n) = row(1); Lvf(n) = row(6); Mvf(n) = row(7); Cpvf(n) = row(8)
        enddo
        close(u)
        call check(n == Tmax - Tmin + 1 .and. nint(T(1)) == Tmin,  'PR-V1 atlas-: the file read here spans the table', nbad)
        call check(all(Mvf(1:n) == 18.015_R8),                   'PR-V1 atlas-: the Mv column is 18.015 on every row', nbad)
        do k = 1, 6
            eLv  = rel(tableValue(mat(1)%LvTab, Tmin, Tmax, Tq(k)), file_value(Lvf, Tq(k)))
            eCpv = rel(tableValue(mat(1)%cpvTab, Tmin, Tmax, Tq(k)), file_value(Cpvf, Tq(k)))
            tol  = merge(0._R8, 4._R8*epsilon(1._R8), onNode(k))
            write(*,'(a,f7.2,a,es9.2,a,es9.2)') '    T = ', Tq(k), ' K: rel Lv ', eLv, ', rel Cpv ', eCpv
            call check(eLv <= tol .and. eCpv <= tol, 'PR-V1 atlas-: Lv and cpv at T = the file''s by the same formula', nbad)
        enddo
        call check(rel(mat(1)%Lv, file_value(Lvf, 298.15_R8)) <= 4._R8*epsilon(1._R8), &
                   'PR-V1 atlas-: the scalar Lv = Lv(298.15 K) of the column', nbad)
    end subroutine atlas_vapour_legs


    !> A column read from the file at T: linear between the nodes, the end value outside (tableValue's rule).
    pure function file_value(col, Tx) result(v)
        real(R8), intent(in) :: col(:), Tx
        real(R8) :: v, Tc
        integer  :: i
        Tc = min(max(Tx, real(Tmin, R8)), real(Tmax, R8))
        if (Tc >= real(Tmax, R8)) then
            v = col(Tmax - Tmin + 1)
        else
            i = int(Tc)
            v = col(i - Tmin + 1) + (col(i - Tmin + 2) - col(i - Tmin + 1))*(Tc - real(i, R8))
        endif
    end function file_value


    pure function rel(a, b) result(e)
        real(R8), intent(in) :: a, b
        real(R8) :: e
        e = abs(a - b)/abs(b)
    end function rel


    !> The table checks and the psat lookup on synthetic input, one outcome per code.
    subroutine function_legs(nbad)
        integer, intent(inout) :: nbad
        real(R8) :: tab(100:110), bad(100:110), T(11), cp(11), rho(11), h(11), nan
        character(len=64) :: tok4(4), tok5(5)
        integer :: k, icp, irho, ih, ips, code, ilv, imv, icpv
        integer :: nzone, nsize, ndata, ntrail, badline, nannounced(nzone_max), nrows(nzone_max)
        logical :: relative

        nan = ieee_value(1._R8, ieee_quiet_nan)
        do k = 100, 110
            tab(k) = Patm*exp(-5000._R8*(1._R8/real(k, R8) - 1._R8/105._R8))
        enddo

        ! PR4: every Psat code, in the validator's order of precedence
        call check(validate_psat_column(tab, 100, 110, 105._R8, Patm) == PSAT_OK,      'PR4 rising column, psat(Tboil) = 1 atm: OK', nbad)
        bad = 0._R8
        call check(validate_psat_column(bad, 100, 110, 105._R8, Patm) == PSAT_ABSENT,  'PR4 all zero: ABSENT', nbad)
        bad = tab; bad(103) = nan
        call check(validate_psat_column(bad, 100, 110, 105._R8, Patm) == PSAT_NONFINITE, 'PR4 a NaN: NONFINITE', nbad)
        bad = tab; bad(100) = -1._R8
        call check(validate_psat_column(bad, 100, 110, 105._R8, Patm) == PSAT_NEGATIVE, 'PR4 a negative value: NEGATIVE', nbad)
        bad = tab; bad(103) = tab(104); bad(104) = tab(103)
        call check(validate_psat_column(bad, 100, 110, 105._R8, Patm) == PSAT_NONMONOTONE, 'PR4 a decrease: NONMONOTONE', nbad)
        bad = Patm
        call check(validate_psat_column(bad, 100, 110, 105._R8, Patm) == PSAT_CONSTANT, 'PR4 constant: CONSTANT', nbad)
        call check(validate_psat_column(tab, 100, 110, 110._R8, Patm) == PSAT_TBOIL_RANGE .and. &
                   validate_psat_column(tab, 100, 110, 99.5_R8, Patm) == PSAT_TBOIL_RANGE, &
                   'PR4 Tboil at Tmax or below Tmin: TBOIL_RANGE', nbad)
        call check(validate_psat_column(3._R8*tab, 100, 110, 105._R8, Patm) == PSAT_TBOIL_MISMATCH .and. &
                   validate_psat_column(0.4_R8*tab, 100, 110, 105._R8, Patm) == PSAT_TBOIL_MISMATCH, &
                   'PR4 psat(Tboil) = 3 atm or 0.4 atm: TBOIL_MISMATCH', nbad)
        call check(validate_psat_column(1.9_R8*tab, 100, 110, 105._R8, Patm) == PSAT_OK .and. &
                   validate_psat_column(0.6_R8*tab, 100, 110, 105._R8, Patm) == PSAT_OK, &
                   'PR4 psat(Tboil) = 1.9 atm or 0.6 atm: OK (a database curve)', nbad)
        bad = tab; bad(100:102) = 0._R8
        call check(validate_psat_column(bad, 100, 110, 105._R8, Patm) == PSAT_OK, &
                   'PR4 equal neighbours (low-T underflow to 0): OK', nbad)

        ! PR5: the nodes
        T = [(real(279 + k, R8), k = 1, 11)]
        call check(check_table_nodes(T) == TAB_OK,                                'PR5 280..290 K: OK', nbad)
        call check(check_table_nodes(T + 1.e-7_R8) == TAB_OK,                     'PR5 1e-7 K off the nodes: OK', nbad)
        call check(check_table_nodes(T + 0.4_R8) == TAB_OFF_NODE,                 'PR5 all 0.4 K off the nodes: OFF_NODE', nbad)
        call check(check_table_nodes([T(1:5), T(7:11)]) == TAB_OFF_NODE,          'PR5 a missing node: OFF_NODE', nbad)
        call check(check_table_nodes(T(1:1)) == TAB_FEW_ROWS,                     'PR5 one row: FEW_ROWS', nbad)
        call check(check_table_nodes([T(1:5), nan, T(7:11)]) == TAB_NONFINITE,    'PR5 a NaN node: NONFINITE', nbad)
        call check(check_table_nodes(T - 281._R8) == TAB_NEGATIVE_T,              'PR5 from -1 K: NEGATIVE_T', nbad)

        ! PR8: the header names
        tok4 = [character(len=64) :: 'Temperature', 'Cp', 'Density', 'Enthalpy']
        call classify_table_tokens(tok4, icp, irho, ih, ips, relative, code, ilv, imv, icpv)
        call check(code == TAB_OK .and. icp == 1 .and. irho == 2 .and. ih == 3 .and. ips == 0 .and. relative, &
                   'PR8 the 4-column header: OK, relative, no Psat', nbad)
        tok5 = [character(len=64) :: 'temperature', 'PSAT', 'enthalpy_abs', 'cp', 'density']
        call classify_table_tokens(tok5, icp, irho, ih, ips, relative, code, ilv, imv, icpv)
        call check(code == TAB_OK .and. icp == 3 .and. irho == 4 .and. ih == 2 .and. ips == 1 .and. .not. relative, &
                   'PR8 aliases in another order: OK, absolute, Psat first', nbad)
        tok5 = [character(len=64) :: 'Temperature', 'Cp', 'Density', 'Enthalpy', 'Enthalpy_abs']
        call classify_table_tokens(tok5, icp, irho, ih, ips, relative, code, ilv, imv, icpv)
        call check(code == TAB_TWO_ENTHALPY,                                     'PR8 both enthalpy names: TWO_ENTHALPY', nbad)
        tok5 = [character(len=64) :: 'Temperature', 'Cp', 'Density', 'Enthalpy', 'Cp']
        call classify_table_tokens(tok5, icp, irho, ih, ips, relative, code, ilv, imv, icpv)
        call check(code == TAB_DUPLICATE,                                        'PR8 Cp twice: DUPLICATE', nbad)
        tok4 = [character(len=64) :: 'Cp', 'Temperature', 'Density', 'Enthalpy']
        call classify_table_tokens(tok4, icp, irho, ih, ips, relative, code, ilv, imv, icpv)
        call check(code == TAB_NO_TEMPERATURE,                                   'PR8 Temperature not first: NO_TEMPERATURE', nbad)
        tok4 = [character(len=64) :: 'Temperature', 'Density', 'Enthalpy', 'Psat']
        call classify_table_tokens(tok4, icp, irho, ih, ips, relative, code, ilv, imv, icpv)
        call check(code == TAB_NO_CP,                                            'PR8 no Cp: NO_CP', nbad)
        tok4 = [character(len=64) :: 'Temperature', 'Cp', 'Enthalpy', 'Psat']
        call classify_table_tokens(tok4, icp, irho, ih, ips, relative, code, ilv, imv, icpv)
        call check(code == TAB_NO_DENSITY,                                       'PR8 no Density: NO_DENSITY', nbad)
        tok4 = [character(len=64) :: 'Temperature', 'Cp', 'Density', 'Psat']
        call classify_table_tokens(tok4, icp, irho, ih, ips, relative, code, ilv, imv, icpv)
        call check(code == TAB_NO_ENTHALPY,                                      'PR8 no enthalpy: NO_ENTHALPY', nbad)

        ! PR11: the columns
        cp = 2800._R8; rho = 700._R8; h = cp*T
        call check(check_table_columns(T, cp, rho, h, .true.) == TAB_OK,          'PR11 constant cp, h = cp*T, relative: OK', nbad)
        call check(check_table_columns(T, cp, rho, h - 2.e6_R8, .false.) == TAB_OK, 'PR11 absolute with an offset: OK', nbad)
        call check(check_table_columns(T, cp, rho, h - 2.e6_R8, .true.) == TAB_DATUM_MISMATCH, &
                   'PR11 relative with an offset: DATUM_MISMATCH', nbad)
        h(6) = h(6) + 10._R8
        call check(check_table_columns(T, cp, rho, h, .true.) == TAB_H_CP_MISMATCH, 'PR11 a row off cp*T: H_CP_MISMATCH', nbad)
        h = cp*T; h(6) = h(5)
        call check(check_table_columns(T, cp, rho, h, .true.) == TAB_H_NONMONOTONE, 'PR11 h not increasing: H_NONMONOTONE', nbad)
        h = cp*T; rho(3) = 0._R8
        call check(check_table_columns(T, cp, rho, h, .true.) == TAB_RHO_NONPOSITIVE, 'PR11 rho = 0: RHO_NONPOSITIVE', nbad)
        rho = 700._R8; cp(3) = -1._R8
        call check(check_table_columns(T, cp, rho, h, .true.) == TAB_CP_NONPOSITIVE, 'PR11 cp < 0: CP_NONPOSITIVE', nbad)
        cp = 2800._R8; rho(4) = nan
        call check(check_table_columns(T, cp, rho, h, .true.) == TAB_NONFINITE,   'PR11 a NaN density: NONFINITE', nbad)
        rho = 700._R8
        cp = [(1250._R8 + 10._R8*real(k, R8), k = 1, 11)]
        h(1) = cp(1)*T(1)
        do k = 2, 11
            h(k) = h(k-1) + 0.5_R8*(cp(k-1) + cp(k))
        enddo
        call check(check_table_columns(T, cp, rho, h, .true.) == TAB_OK,          'PR11 varying cp, h its trapezoid sum: OK', nbad)
        do k = 2, 11
            h(k) = h(k-1) + cp(k)
        enddo
        call check(check_table_columns(T, cp, rho, h, .true.) == TAB_H_CP_MISMATCH, &
                   'PR11 varying cp, h its right Riemann sum: H_CP_MISMATCH', nbad)

        ! PR12: psat through tableValue
        call check(abs(tableValue(tab, 100, 110, 104.25_R8) - (0.75_R8*tab(104) + 0.25_R8*tab(105))) <= 1.e-14_R8*tab(105), &
                   'PR12 inside: linear between the nodes', nbad)
        call check(tableValue(tab, 100, 110, 95._R8) == tab(100) .and. tableValue(tab, 100, 110, 115._R8) == tab(110), &
                   'PR12 outside: the end values', nbad)
        call check(tableValue(tab, 100, 110, 110._R8) == tab(110),                'PR12 at Tmax: tab(Tmax)', nbad)
        call check(ieee_is_nan(tableValue(tab, 100, 110, nan)),                   'PR12 a NaN temperature: NaN', nbad)

        ! PR13: the row scan
        call scan_rows('scan-short.dat', 4, nzone, nsize, ndata, nannounced, nrows, ntrail, badline)
        call check(badline == 7 .and. nrows(1) == 5 .and. ntrail == 0, 'PR13 a short row: badline = its line 7', nbad)
        call scan_rows('scan-trail.dat', 4, nzone, nsize, ndata, nannounced, nrows, ntrail, badline)
        call check(ntrail == 1 .and. badline == 0 .and. nrows(1) == 3, 'PR13 text after the rows: ntrail = 1', nbad)
        call scan_rows('scan-count.dat', 4, nzone, nsize, ndata, nannounced, nrows, ntrail, badline)
        call check(nannounced(1) == 4 .and. nrows(1) == 3,             'PR13 I=4 over 3 rows: announced 4, held 3', nbad)
        call scan_rows('../two-material/INPUT/nm-properties.dat', 5, nzone, nsize, ndata, nannounced, nrows, ntrail, badline)
        call check(nzone == 2 .and. nsize == 2 .and. ndata == 2 .and. all(nannounced(1:2) == 620) .and. &
                   all(nrows(1:2) == 620) .and. ntrail == 0 .and. badline == 0, 'PR13 two zones of 620 rows: clean', nbad)

        call vapour_code_legs(tab, nan, nbad)
    end subroutine function_legs


    !> PR-V2: validate_vapour_columns returns each code once, on columns beside the PR4 Psat column
    !  (ln p = -5000/T + const, so Lv*Mv/Ru = 5000 K is its Clausius-Clapeyron slope).
    subroutine vapour_code_legs(ps, nan, nbad)
        real(R8), intent(in)    :: ps(100:110), nan
        integer,  intent(inout) :: nbad
        real(R8), parameter :: RuT = 8314.46_R8
        real(R8) :: Lv(100:110), Mv(100:110), cv(100:110), zero(100:110), bad(100:110)

        Lv = 5000._R8*RuT/18._R8; Mv = 18._R8; cv = 1900._R8; zero = 0._R8
        call check(validate_vapour_columns(Lv, Mv, cv, ps, .true., 100, 110, 105._R8, RuT, 0._R8) == VAP_OK, &
                   'PR-V2 the three columns on the Psat slope: OK', nbad)
        call check(validate_vapour_columns(zero, zero, zero, ps, .true., 100, 110, 105._R8, RuT, 18._R8) == VAP_ABSENT, &
                   'PR-V2 all three zero: ABSENT', nbad)
        bad = Lv; bad(103) = nan
        call check(validate_vapour_columns(bad, Mv, cv, ps, .true., 100, 110, 105._R8, RuT, 0._R8) == VAP_NONFINITE, &
                   'PR-V2 a NaN in Lv: NONFINITE', nbad)
        bad = cv; bad(100:102) = 0._R8
        call check(validate_vapour_columns(Lv, Mv, bad, ps, .true., 100, 110, 105._R8, RuT, 0._R8) == VAP_PARTIAL, &
                   'PR-V2 Cpv zero on three rows: PARTIAL', nbad)
        call check(validate_vapour_columns(Lv, Mv, -cv, ps, .true., 100, 110, 105._R8, RuT, 0._R8) == VAP_NONPOSITIVE, &
                   'PR-V2 Cpv negative: NONPOSITIVE', nbad)
        bad = Mv; bad(104) = 18.5_R8
        call check(validate_vapour_columns(Lv, bad, cv, ps, .true., 100, 110, 105._R8, RuT, 0._R8) == VAP_MV_VARIES, &
                   'PR-V2 Mv not constant: MV_VARIES', nbad)
        call check(validate_vapour_columns(Lv, Mv, cv, zero, .false., 100, 110, 105._R8, RuT, 0._R8) == VAP_WITHOUT_PSAT, &
                   'PR-V2 no Psat column: WITHOUT_PSAT', nbad)
        call check(validate_vapour_columns(Lv/1000._R8, Mv, cv, ps, .true., 100, 110, 105._R8, RuT, 0._R8) == &
                   VAP_CC_MISMATCH, 'PR-V2 Lv in kJ/kg: CC_MISMATCH', nbad)
        ! the slope with Mv from the key: the column absent, MvKey 18 (OK) or 18/2.5 (ratio 2.5)
        call check(validate_vapour_columns(Lv, zero, cv, ps, .true., 100, 110, 105._R8, RuT, 18._R8) == VAP_OK .and. &
                   validate_vapour_columns(Lv, zero, cv, ps, .true., 100, 110, 105._R8, RuT, 18._R8/2.5_R8) == &
                   VAP_CC_MISMATCH, 'PR-V2 Mv from the key: 18 OK, 18/2.5 CC_MISMATCH', nbad)
    end subroutine vapour_code_legs


    !> comp_TfromTab against lookupTab on tables that start above 1 K, and ICE's energy rule on the same table.
    subroutine inversion_legs(nbad)
        integer, intent(inout) :: nbad
        real(R8) :: hA(250:260), hB(300:301), hC(1:5), Ts(2), nan, hOff, eLo, eHi, e, Tice, h, T
        integer  :: k, lo0, hi0
        logical  :: ok

        lo0 = Tmin; hi0 = Tmax
        nan = ieee_value(1._R8, ieee_quiet_nan)

        ! 250..260 K: cp = 1000 + 10 (T - 250), h its trapezoid sum from h(250) = 2.5e5
        hA(250) = 2.5e5_R8
        do k = 251, 260
            hA(k) = hA(k-1) + 0.5_R8*((1000._R8 + 10._R8*real(k-251, R8)) + (1000._R8 + 10._R8*real(k-250, R8)))
        enddo
        Tmin = 250; Tmax = 260
        ok = .true.
        do k = 250, 260
            ok = ok .and. comp_TfromTab(hA, hA(k)) == real(k, R8)
        enddo
        call check(ok,                                                  'PR14 250..260 K: every node, exactly', nbad)
        call check(roundtrip(hA, [250.25_R8, 254.5_R8, 259.999_R8]),    'PR14 250..260 K: between the nodes', nbad)
        call check(roundtrip(hA, [245.5_R8, 200._R8]),                  'PR14 250..260 K: below Tmin, the first segment extended', nbad)
        call check(roundtrip(hA, [263.25_R8, 300._R8]),                 'PR14 250..260 K: above Tmax, the last segment extended', nbad)
        call check(ieee_is_nan(comp_TfromTab(hA, nan)),                 'PR14 a NaN enthalpy: NaN', nbad)

        ! ICE's e = h - hOff, hOff = h(Tmin) - Tmin (h(Tmin+1) - h(Tmin)): the line from e(0) = 0 below Tmin,
        ! the last segment above Tmax, forward and inverse
        hOff = hA(250) - 250._R8*(hA(251) - hA(250))
        eLo = hA(250) - hOff
        eHi = hA(260) - hOff
        ok = .true.
        Ts = [245.5_R8, 200._R8]
        do k = 1, 2
            e    = eLo*(Ts(k)/250._R8)
            h    = lookupTab(hA, Ts(k))
            Tice = 250._R8*((h - hOff)/eLo)
            T    = comp_TfromTab(hA, h)
            ok = ok .and. abs(e + hOff - h) <= 4._R8*epsilon(1._R8)*abs(h) .and. &
                 abs(T - Tice) <= 4._R8*epsilon(1._R8)*(abs(h)/(hA(251) - hA(250)) + abs(T))
        enddo
        Ts = [263.25_R8, 300._R8]
        do k = 1, 2
            e    = eHi + (eHi - (hA(259) - hOff))*(Ts(k) - 260._R8)
            h    = lookupTab(hA, Ts(k))
            Tice = 260._R8 + ((h - hOff) - eHi)/(eHi - (hA(259) - hOff))
            T    = comp_TfromTab(hA, h)
            ok = ok .and. abs(e + hOff - h) <= 4._R8*epsilon(1._R8)*abs(h) .and. &
                 abs(T - Tice) <= 4._R8*epsilon(1._R8)*(abs(h)/(hA(260) - hA(259)) + abs(T))
        enddo
        call check(ok,                                                  'PR14 250..260 K: ICE''s energy rule outside the table', nbad)

        ! two rows, 300..301 K
        hB = [3.0e5_R8, 3.0e5_R8 + 1005._R8]
        Tmin = 300; Tmax = 301
        call check(comp_TfromTab(hB, hB(300)) == 300._R8 .and. comp_TfromTab(hB, hB(301)) == 301._R8, &
                   'PR14 300..301 K: both nodes, exactly', nbad)
        call check(roundtrip(hB, [300.5_R8, 299._R8, 290._R8, 305._R8]), 'PR14 300..301 K: between, below and above', nbad)

        ! 1..5 K, an enthalpy below h(1 K)
        hC(1) = 1000._R8
        do k = 2, 5
            hC(k) = hC(k-1) + 0.5_R8*((1000._R8 + 10._R8*real(k-2, R8)) + (1000._R8 + 10._R8*real(k-1, R8)))
        enddo
        Tmin = 1; Tmax = 5
        call check(roundtrip(hC, [0.5_R8, 2.75_R8]),                    'PR14 1..5 K: below 1 K and between the nodes', nbad)

        Tmin = lo0; Tmax = hi0
    end subroutine inversion_legs


    !> lookupTab on the hexadecane enthalpy at T = +Inf, -Inf and NaN: NaN, and the NaN lookup raises no invalid flag.
    subroutine nonfinite_legs(nbad)
        integer, intent(inout) :: nbad
        type(obj_material), allocatable :: mat(:)
        character(len=32) :: prefix
        real(R8), volatile :: Tin, v
        logical :: flag

        prefix = 'legacy-'
        call read_cdp_properties(prefix, mat)
        call check(ieee_is_nan(lookupTab(mat(1)%hTab, ieee_value(1._R8, ieee_positive_inf))), &
                   'PR15 legacy-: h at T = +Inf is NaN', nbad)
        call check(ieee_is_nan(lookupTab(mat(1)%hTab, ieee_value(1._R8, ieee_negative_inf))), &
                   'PR15 legacy-: h at T = -Inf is NaN', nbad)
        Tin = ieee_value(1._R8, ieee_quiet_nan)
        call ieee_set_flag(ieee_invalid, .false.)
        v = lookupTab(mat(1)%hTab, Tin)
        call check(ieee_is_nan(v),                                      'PR15 legacy-: h at T = NaN is NaN', nbad)
        call ieee_get_flag(ieee_invalid, flag)
        call check(.not. flag,                                          'PR15 legacy-: the NaN lookup raises no invalid flag', nbad)
    end subroutine nonfinite_legs


    !> lookupTab then comp_TfromTab returns every T within 4 eps (|h|/s + |T|), s the table's smallest slope.
    logical function roundtrip(tab, Ts)
        real(R8), intent(in) :: tab(Tmin:Tmax), Ts(:)
        real(R8) :: h, s, T
        integer  :: k
        s = minval(tab(Tmin+1:Tmax) - tab(Tmin:Tmax-1))
        roundtrip = .true.
        do k = 1, size(Ts)
            h = lookupTab(tab, Ts(k))
            T = comp_TfromTab(tab, h)
            if (.not. abs(T - Ts(k)) <= 4._R8*epsilon(1._R8)*(abs(h)/s + abs(Ts(k)))) then
                write(*,'(a,es24.16,a,es24.16)') '         T = ', Ts(k), ' returned ', T
                roundtrip = .false.
            endif
        enddo
    end function roundtrip


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

end program test_properties_reader
