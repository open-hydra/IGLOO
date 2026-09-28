program test_bc_families
    !
    ! A particle boundary file carrying several families.
    !
    ! ATLAS writes <name>-bc.txt with one copy of every mesh block's boundary table per
    ! (material, population): block-outermost, material-major, population-minor. That is IGLOO's
    ! family order (famID counts the groups material by material), so copy c of a block is
    ! family c. A file with one copy -- what ATLAS writes for a single-population phase, and every
    ! box fixture -- feeds every family; any other record count is refused at setup.
    !
    ! Each case writes its own INPUT/famtest-<case>-bc.txt in the ATLAS layout (header
    ! "b i j k f code", f = 1..6, n outer, m inner; face 1 an inlet with a payload line, every
    ! other face 100) and reads it through read_cdp_bc_file. Blocks 4 x 3 x 1 (38 boundary
    ! faces) and 5 x 3 x 1 (46), H = 0.25. mdotGas is seeded on every face cell, so initMdotGas
    ! returns at once and the case stays a parser test. Two families throughout.
    !
    !   F1  one copy (control)       fam 1 = fam 2 = (0.1, 1e-6), krhoTot 0.2          RED the same
    !   F2  two copies, one block    fam 1 (0.1, 1e-6), fam 2 (0.2, 2e-6), krhoTot 0.3
    !                                RED fam 2 (0.1, 1e-6), krhoTot 0.2 (copy 2 never read)
    !       run twice: two materials of one group (F2a), one material of two groups (F2b)
    !   F3  two copies, two blocks   b1 as F2; b2 (0.3, 3e-6)/(0.4, 4e-6), krhoTot 0.7, three 401
    !                                cells per block, all on face 1
    !                                RED block 2 reads block 1's copy 2: b2 face 1 (0.2, 2e-6) for
    !                                both families, three more 401 cells on b2 face 6
    !   F4  two copies, 402          mdotPart = (0.1 + 0.2) H^2                         RED 0.2 H^2
    !   F5  two copies, ds on copy 2 properties(2,9) = 0.01, dsSwitch on               RED 0, off
    !
    ! A pair (a, b) is a family's payload columns 1 and 6: krho (401) or gp (402), and rp.
    !
    use, intrinsic :: iso_fortran_env, only: real64
    use IGLOO_IO,          only: read_cdp_bc_file
    use IGLOO_data_block,  only: obj_block, obj_flowblock, obj_sourceblock, obj_eulerblock
    use IGLOO_data_phases, only: obj_material
    use IGLOO_variables,   only: nb, nm, mesh2D, dsSwitch
    implicit none

    real(real64), parameter :: H = 0.25_real64
    real(real64), parameter :: TOL = 1.0e-15_real64
    real(real64), parameter :: MDOTGAS = 1.2e-3_real64
    integer :: nfail

    nfail = 0
    call execute_command_line('mkdir -p INPUT', wait=.true.)

    call f1_one_copy
    call f2_two_copies('F2a', [1, 1])
    call f2_two_copies('F2b', [2])
    call f3_two_blocks
    call f4_mass_flux
    call f5_ds

    if (nfail == 0) then
        print '(a)', ''
        print '(a)', '[PASS] test_bc_families: copy c of every block is family c; one copy feeds every family'
    else
        print '(a)', ''
        print '(a,i0,a)', '[FAIL] test_bc_families: ', nfail, ' assertion(s) failed'
        error stop 1
    endif

contains

    !> F1 -- the control: one copy is every family's, and krhoTot sums it once per family.
    subroutine f1_one_copy
        type(obj_block)                 :: geo(1)
        type(obj_flowblock)             :: gas(1)
        type(obj_sourceblock)           :: src(1)
        type(obj_eulerblock)            :: eul(1,1)
        type(obj_material), allocatable :: mat(:)
        integer :: u
        call setup(geo, reshape([4, 3, 1], [3, 1]), mat, [1, 1])
        open(newunit=u, file='INPUT/famtest-F1-bc.txt', status='replace', action='write')
        call writeCopy(u, 1, 4, 3, 1, 401, 0.1_real64, 1.0e-6_real64, 0.0_real64)
        close(u)
        call read_cdp_bc_file('famtest-F1-', mat, geo, gas, src, eul, .false., .false.)
        call expectInlet('F1', geo(1), [0.1_real64, 0.1_real64], [1.0e-6_real64, 1.0e-6_real64], 0.2_real64)
        call expectTags('F1', geo(1), 3)
    end subroutine f1_one_copy

    !> F2 -- two copies on one block: copy 2 is family 2, whether the two families are two
    !  materials or two groups of one material.
    subroutine f2_two_copies(lab, groups)
        character(*), intent(in) :: lab
        integer,      intent(in) :: groups(:)
        type(obj_block)                 :: geo(1)
        type(obj_flowblock)             :: gas(1)
        type(obj_sourceblock)           :: src(1)
        type(obj_eulerblock)            :: eul(1,1)
        type(obj_material), allocatable :: mat(:)
        integer :: u
        call setup(geo, reshape([4, 3, 1], [3, 1]), mat, groups)
        open(newunit=u, file='INPUT/famtest-'//lab//'-bc.txt', status='replace', action='write')
        call writeCopy(u, 1, 4, 3, 1, 401, 0.1_real64, 1.0e-6_real64, 0.0_real64)
        call writeCopy(u, 1, 4, 3, 1, 401, 0.2_real64, 2.0e-6_real64, 0.0_real64)
        close(u)
        call read_cdp_bc_file('famtest-'//lab//'-', mat, geo, gas, src, eul, .false., .false.)
        call expectInlet(lab, geo(1), [0.1_real64, 0.2_real64], [1.0e-6_real64, 2.0e-6_real64], 0.3_real64)
        call expectTags(lab, geo(1), 3)
    end subroutine f2_two_copies

    !> F3 -- two blocks: the copies of block 1 come first, then those of block 2.
    subroutine f3_two_blocks
        type(obj_block)                 :: geo(2)
        type(obj_flowblock)             :: gas(2)
        type(obj_sourceblock)           :: src(2)
        type(obj_eulerblock)            :: eul(2,2)
        type(obj_material), allocatable :: mat(:)
        integer :: u
        call setup(geo, reshape([4, 3, 1, 5, 3, 1], [3, 2]), mat, [1, 1])
        open(newunit=u, file='INPUT/famtest-F3-bc.txt', status='replace', action='write')
        call writeCopy(u, 1, 4, 3, 1, 401, 0.1_real64, 1.0e-6_real64, 0.0_real64)
        call writeCopy(u, 1, 4, 3, 1, 401, 0.2_real64, 2.0e-6_real64, 0.0_real64)
        call writeCopy(u, 2, 5, 3, 1, 401, 0.3_real64, 3.0e-6_real64, 0.0_real64)
        call writeCopy(u, 2, 5, 3, 1, 401, 0.4_real64, 4.0e-6_real64, 0.0_real64)
        close(u)
        call read_cdp_bc_file('famtest-F3-', mat, geo, gas, src, eul, .false., .false.)
        call expectInlet('F3 block 1', geo(1), [0.1_real64, 0.2_real64], [1.0e-6_real64, 2.0e-6_real64], 0.3_real64)
        call expectTags('F3 block 1', geo(1), 3)
        call expectInlet('F3 block 2', geo(2), [0.3_real64, 0.4_real64], [3.0e-6_real64, 4.0e-6_real64], 0.7_real64)
        call expectTags('F3 block 2', geo(2), 3)
    end subroutine f3_two_blocks

    !> F4 -- a 402 inlet: mdotPart sums gp * area over the families' own gp.
    subroutine f4_mass_flux
        type(obj_block)                 :: geo(1)
        type(obj_flowblock)             :: gas(1)
        type(obj_sourceblock)           :: src(1)
        type(obj_eulerblock)            :: eul(1,1)
        type(obj_material), allocatable :: mat(:)
        integer :: u, m
        call setup(geo, reshape([4, 3, 1], [3, 1]), mat, [1, 1])
        open(newunit=u, file='INPUT/famtest-F4-bc.txt', status='replace', action='write')
        call writeCopy(u, 1, 4, 3, 1, 402, 0.1_real64, 1.0e-6_real64, 0.0_real64)
        call writeCopy(u, 1, 4, 3, 1, 402, 0.2_real64, 2.0e-6_real64, 0.0_real64)
        close(u)
        call read_cdp_bc_file('famtest-F4-', mat, geo, gas, src, eul, .false., .false.)
        do m = 1, 3
            associate(cell => geo(1)%face(1)%cell(m, 1))
            call expectEq('F4 cell '//itoa(m)//' mdotPart', cell%mdotPart, 0.3_real64*H*H)
            call expectEq('F4 cell '//itoa(m)//' krhoTot', cell%krhoTot, 0.0_real64)
            end associate
        enddo
    end subroutine f4_mass_flux

    !> F5 -- the per-cell spacing (column 9) of family 2 reaches family 2 and switches ds on.
    subroutine f5_ds
        type(obj_block)                 :: geo(1)
        type(obj_flowblock)             :: gas(1)
        type(obj_sourceblock)           :: src(1)
        type(obj_eulerblock)            :: eul(1,1)
        type(obj_material), allocatable :: mat(:)
        integer :: u
        call setup(geo, reshape([4, 3, 1], [3, 1]), mat, [1, 1])
        open(newunit=u, file='INPUT/famtest-F5-bc.txt', status='replace', action='write')
        call writeCopy(u, 1, 4, 3, 1, 401, 0.1_real64, 1.0e-6_real64, 0.0_real64)
        call writeCopy(u, 1, 4, 3, 1, 401, 0.2_real64, 2.0e-6_real64, 0.01_real64)
        close(u)
        call read_cdp_bc_file('famtest-F5-', mat, geo, gas, src, eul, .false., .false.)
        call expectEq('F5 fam 1 ds', geo(1)%face(1)%cell(1, 1)%properties(1, 9), 0.0_real64)
        call expectEq('F5 fam 2 ds', geo(1)%face(1)%cell(1, 1)%properties(2, 9), 0.01_real64)
        if (dsSwitch) then
            print '(a)', '  [ok]   famtest F5 dsSwitch: got T want T'
        else
            print '(a)', '  [FAIL] famtest F5 dsSwitch: got F want T'
            nfail = nfail + 1
        endif
    end subroutine f5_ds


    ! ------------------------------------------------------------------ fixture

    !> Blocks of the given sizes through the production geometry path, mdotGas seeded on every
    !  face cell; materials with the given group counts.
    subroutine setup(geo, dims, mat, groups)
        type(obj_block),                 intent(inout) :: geo(:)
        integer,                         intent(in)    :: dims(:,:)
        type(obj_material), allocatable, intent(out)   :: mat(:)
        integer,                         intent(in)    :: groups(:)
        integer :: b, f
        nb = size(dims, 2); nm = size(groups)
        mesh2D = all(dims(3,:) == 1)
        dsSwitch = .false.
        do b = 1, nb
            call buildGeo(geo(b), dims(1,b), dims(2,b), dims(3,b))
            do f = 1, 6
                geo(b)%face(f)%cell(:,:)%mdotGas = MDOTGAS
            enddo
        enddo
        allocate(mat(nm))
        mat(:)%ngroups = groups
    end subroutine setup

    subroutine buildGeo(blk, nx, ny, nz)
        type(obj_block), intent(inout) :: blk
        integer,         intent(in)    :: nx, ny, nz
        integer :: i, j, k
        blk%Nx = nx; blk%Ny = ny; blk%Nz = nz
        allocate(blk%node(3, 0:nx, 0:ny, 0:nz))
        do k = 0, nz; do j = 0, ny; do i = 0, nx
            blk%node(:, i, j, k) = [i*H, j*H, k*H]
        enddo; enddo; enddo
        allocate(blk%center(3, 1:nx, 1:ny, 1:nz))
        call blk%compute_geometry                 ! allocates face(f)%cell, fills centre + normal
    end subroutine buildGeo

    !> One copy of block b's boundary table in the ATLAS layout: header "b i j k f code" per face
    !  cell (fmn2ijk's face rule), f = 1..6, n outer, m inner; face 1 carries code1 and a payload
    !  line in drag-stokes's form, every other face 100.
    subroutine writeCopy(u, b, nx, ny, nz, code1, col1, rp, ds)
        integer,      intent(in) :: u, b, nx, ny, nz, code1
        real(real64), intent(in) :: col1, rp, ds
        integer :: f, m, n, i, j, k, mend(6), nend(6)
        mend(1:2) = ny; nend(1:2) = nz
        mend(3:4) = nx; nend(3:4) = nz
        mend(5:6) = nx; nend(5:6) = ny
        do f = 1, 6; do n = 1, nend(f); do m = 1, mend(f)
            select case (f)
            case (1); i = 1;  j = m;  k = n
            case (2); i = nx; j = m;  k = n
            case (3); i = m;  j = 1;  k = n
            case (4); i = m;  j = ny; k = n
            case (5); i = m;  j = n;  k = 1
            case (6); i = m;  j = n;  k = nz
            end select
            if (f == 1) then
                write(u, '(6i8)') b, i, j, k, f, code1
                write(u, '(2(2x,es12.5),2(3x,a),3(2x,es12.5),3x,a,2x,es12.5)') col1, 0.1_real64, &
                    'normal,', 'normal,', 1.0_real64, rp, 0.0_real64, 'Dirac', ds
            else
                write(u, '(6i8)') b, i, j, k, f, 100
            endif
        enddo; enddo; enddo
    end subroutine writeCopy


    ! ------------------------------------------------------------------ report

    !> Face-1 inlet cells of blk: code 401, (column 1, column 6) per family, krhoTot.
    subroutine expectInlet(lab, blk, col1, rp, krhoTot)
        character(*),    intent(in) :: lab
        type(obj_block), intent(in) :: blk
        real(real64),    intent(in) :: col1(:), rp(:), krhoTot
        integer :: m, fam
        do m = 1, blk%Ny
            associate(cell => blk%face(1)%cell(m, 1))
            if (cell%bcdef /= 401) then
                print '(a,i0,a)', '  [FAIL] famtest '//lab//' face 1 cell '//itoa(m)//' bcdef: got ', &
                    cell%bcdef, ' want 401'
                nfail = nfail + 1
            else
                do fam = 1, size(col1)
                    call expectEq(lab//' cell '//itoa(m)//' fam '//itoa(fam)//' col1', cell%properties(fam, 1), col1(fam))
                    call expectEq(lab//' cell '//itoa(m)//' fam '//itoa(fam)//' rp',   cell%properties(fam, 6), rp(fam))
                enddo
                call expectEq(lab//' cell '//itoa(m)//' krhoTot', cell%krhoTot, krhoTot)
            endif
            end associate
        enddo
    end subroutine expectInlet

    !> Inlet cells of blk: exactly n401 cells tagged 401, all on face 1; every other cell 100.
    subroutine expectTags(lab, blk, n401)
        character(*),    intent(in) :: lab
        type(obj_block), intent(in) :: blk
        integer,         intent(in) :: n401
        integer :: f, per(6), got100, want100
        got100 = 0; want100 = 0
        do f = 1, 6
            per(f)  = count(blk%face(f)%cell(:,:)%bcdef == 401)
            got100  = got100 + count(blk%face(f)%cell(:,:)%bcdef == 100)
            want100 = want100 + size(blk%face(f)%cell)
        enddo
        want100 = want100 - n401
        if (per(1) == n401 .and. sum(per(2:6)) == 0 .and. got100 == want100) then
            print '(a,2(i0,a))', '  [ok]   famtest '//lab//' tags: ', per(1), ' x 401 on face 1, ', got100, ' x 100'
        else
            print '(a,6(i0,1x),a,i0,a,2(i0,a))', '  [FAIL] famtest '//lab//' tags: 401 per face ', per, &
                '; ', got100, ' x 100; want ', n401, ' x 401 on face 1, ', want100, ' x 100'
            nfail = nfail + 1
        endif
    end subroutine expectTags

    subroutine expectEq(label, got, want)
        character(*), intent(in) :: label
        real(real64), intent(in) :: got, want
        if (abs(got - want) <= TOL*max(1.0_real64, abs(want))) then
            print '(a,a,es13.6,a,es13.6)', '  [ok]   famtest ', label//': got ', got, ' want ', want
        else
            print '(a,a,es13.6,a,es13.6)', '  [FAIL] famtest ', label//': got ', got, ' want ', want
            nfail = nfail + 1
        endif
    end subroutine expectEq

    pure function itoa(i) result(s)
        integer, intent(in) :: i
        character(len=:), allocatable :: s
        character(len=12) :: buf
        write(buf, '(i0)') i
        s = trim(buf)
    end function itoa

end program test_bc_families
