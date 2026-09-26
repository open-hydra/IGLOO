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
    !   two-material  PR9 nm-: zone A tabulated, zone B all zero (Clausius-Clapeyron)
    !   atlas-water   PR10 atlas-: the table ATLAS GPB writes for water with psat-vapour = H2O
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64
    use IGLOO_IO_INI,         only: read_IGLOO_input
    use IGLOO_IO,             only: read_cdp_properties
    use IGLOO_data_phases,    only: obj_material
    use IGLOO_Lib_Properties, only: Tmin, Tmax, tableValue, validate_psat_column, check_table_nodes, &
                                    classify_table_tokens, check_table_columns, scan_rows, nzone_max, &
                                    TAB_OK, TAB_NO_TEMPERATURE, TAB_NO_CP, TAB_NO_DENSITY, TAB_NO_ENTHALPY, &
                                    TAB_TWO_ENTHALPY, TAB_DUPLICATE, TAB_FEW_ROWS, TAB_OFF_NODE, TAB_NONFINITE, &
                                    TAB_RHO_NONPOSITIVE, TAB_CP_NONPOSITIVE, TAB_H_NONMONOTONE, &
                                    TAB_H_CP_MISMATCH, TAB_DATUM_MISMATCH, TAB_NEGATIVE_T, &
                                    PSAT_OK, PSAT_ABSENT, PSAT_NONFINITE, PSAT_NEGATIVE, PSAT_NONMONOTONE, &
                                    PSAT_CONSTANT, PSAT_TBOIL_RANGE, PSAT_TBOIL_MISMATCH
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
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
    end subroutine atlas_water_legs


    !> The table checks and the psat lookup on synthetic input, one outcome per code.
    subroutine function_legs(nbad)
        integer, intent(inout) :: nbad
        real(R8) :: tab(100:110), bad(100:110), T(11), cp(11), rho(11), h(11), nan
        character(len=64) :: tok4(4), tok5(5)
        integer :: k, icp, irho, ih, ips, code
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
        call classify_table_tokens(tok4, icp, irho, ih, ips, relative, code)
        call check(code == TAB_OK .and. icp == 1 .and. irho == 2 .and. ih == 3 .and. ips == 0 .and. relative, &
                   'PR8 the 4-column header: OK, relative, no Psat', nbad)
        tok5 = [character(len=64) :: 'temperature', 'PSAT', 'enthalpy_abs', 'cp', 'density']
        call classify_table_tokens(tok5, icp, irho, ih, ips, relative, code)
        call check(code == TAB_OK .and. icp == 3 .and. irho == 4 .and. ih == 2 .and. ips == 1 .and. .not. relative, &
                   'PR8 aliases in another order: OK, absolute, Psat first', nbad)
        tok5 = [character(len=64) :: 'Temperature', 'Cp', 'Density', 'Enthalpy', 'Enthalpy_abs']
        call classify_table_tokens(tok5, icp, irho, ih, ips, relative, code)
        call check(code == TAB_TWO_ENTHALPY,                                     'PR8 both enthalpy names: TWO_ENTHALPY', nbad)
        tok5 = [character(len=64) :: 'Temperature', 'Cp', 'Density', 'Enthalpy', 'Cp']
        call classify_table_tokens(tok5, icp, irho, ih, ips, relative, code)
        call check(code == TAB_DUPLICATE,                                        'PR8 Cp twice: DUPLICATE', nbad)
        tok4 = [character(len=64) :: 'Cp', 'Temperature', 'Density', 'Enthalpy']
        call classify_table_tokens(tok4, icp, irho, ih, ips, relative, code)
        call check(code == TAB_NO_TEMPERATURE,                                   'PR8 Temperature not first: NO_TEMPERATURE', nbad)
        tok4 = [character(len=64) :: 'Temperature', 'Density', 'Enthalpy', 'Psat']
        call classify_table_tokens(tok4, icp, irho, ih, ips, relative, code)
        call check(code == TAB_NO_CP,                                            'PR8 no Cp: NO_CP', nbad)
        tok4 = [character(len=64) :: 'Temperature', 'Cp', 'Enthalpy', 'Psat']
        call classify_table_tokens(tok4, icp, irho, ih, ips, relative, code)
        call check(code == TAB_NO_DENSITY,                                       'PR8 no Density: NO_DENSITY', nbad)
        tok4 = [character(len=64) :: 'Temperature', 'Cp', 'Density', 'Psat']
        call classify_table_tokens(tok4, icp, irho, ih, ips, relative, code)
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
    end subroutine function_legs


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
