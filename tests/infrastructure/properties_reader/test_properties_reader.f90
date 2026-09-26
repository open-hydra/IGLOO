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
    !   two-material  PR9 nm-: zone A tabulated, zone B all zero (Clausius-Clapeyron)
    !   atlas-water   PR10 atlas-: the table ATLAS GPB writes for water with psat-vapour = H2O
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64
    use IGLOO_IO_INI,         only: read_IGLOO_input
    use IGLOO_IO,             only: read_cdp_properties
    use IGLOO_data_phases,    only: obj_material
    use IGLOO_Lib_Properties, only: Tmin, Tmax
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
        endif
    end subroutine atlas_water_legs


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
