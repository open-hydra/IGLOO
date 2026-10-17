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
    ! other face 400) and reads it through read_cdp_bc_file. Blocks 4 x 3 x 1 (38 boundary
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
    !   F6  chimera (102) faces, two blocks, two copies: face 2 of block 1 and face 1 of block 2
    !       carry ATLAS's payload (the donor counts of the two ghost layers, then one line per
    !       donor: block, i, j, k, weight), repeated in copy 2. Each facelet keeps its counts and
    !       its layer-1 donors and weights in file order (one, two, or none with a layer-2 donor);
    !       layer 2 is dropped; every other face 400. RED: a reader without 102 stops at pass 0
    !       on the first chimera record; with the layer-1 lines stored in reverse order the
    !       two-donor facelet fails (0.75, 0.25)
    !   F6b the field mollifier's link table (buildInterfaceLinks) on F6's blocks: one link per 102
    !       facelet with a layer-1 donor (5 of 6: the facelet without one keeps zero flux), each
    !       normal to index direction 1, its donors and weights in file order in the flat donor
    !       list (6 entries). RED: a table of the 101/201 records only holds no link
    !   F7  a one-layer (2D) block with faces 3/4 periodic (201): face 3's connection lines carry
    !       ATLAS's nine integers (block, i, j, k, face, d11..d22), face 4's the eight-integer form
    !       (block, i, j, face, d11..d22); each cell keeps its partner cell, partner face and the
    !       four orientation integers. RED: a reader taking four index integers on a one-layer mesh
    !       stores k (1) as face 3's partner face and (4, 1, 0, 0) as its orientation
    !   F8  the field mollifier's link table (buildInterfaceLinks) on F7's block: 8 links, each face-3
    !       cell joined to the face-4 cell of its column and back, both normal to index direction 2.
    !       RED: a table of the 101 records only holds no link
    !
    ! A pair (a, b) is a family's payload columns 1 and 6: krho (401) or gp (402), and rp.
    !
    use, intrinsic :: iso_fortran_env, only: real64
    use IGLOO_IO,          only: read_cdp_bc_file
    use IGLOO_data_block,  only: obj_block, obj_flowblock, obj_sourceblock, obj_eulerblock, buildInterfaceLinks
    use IGLOO_Lib_Mollify, only: cellLink, linkDonors
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
    call f6_chimera
    call f7_connection_2d

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
        type(obj_material), allocatable :: mat(:)
        integer :: u
        call setup(geo, reshape([4, 3, 1], [3, 1]), mat, [1, 1])
        open(newunit=u, file='INPUT/famtest-F1-bc.txt', status='replace', action='write')
        call writeCopy(u, 1, 4, 3, 1, 401, 0.1_real64, 1.0e-6_real64, 0.0_real64)
        close(u)
        call read_cdp_bc_file('famtest-F1-', mat, geo, gas)
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
        type(obj_material), allocatable :: mat(:)
        integer :: u
        call setup(geo, reshape([4, 3, 1], [3, 1]), mat, groups)
        open(newunit=u, file='INPUT/famtest-'//lab//'-bc.txt', status='replace', action='write')
        call writeCopy(u, 1, 4, 3, 1, 401, 0.1_real64, 1.0e-6_real64, 0.0_real64)
        call writeCopy(u, 1, 4, 3, 1, 401, 0.2_real64, 2.0e-6_real64, 0.0_real64)
        close(u)
        call read_cdp_bc_file('famtest-'//lab//'-', mat, geo, gas)
        call expectInlet(lab, geo(1), [0.1_real64, 0.2_real64], [1.0e-6_real64, 2.0e-6_real64], 0.3_real64)
        call expectTags(lab, geo(1), 3)
    end subroutine f2_two_copies

    !> F3 -- two blocks: the copies of block 1 come first, then those of block 2.
    subroutine f3_two_blocks
        type(obj_block)                 :: geo(2)
        type(obj_flowblock)             :: gas(2)
        type(obj_material), allocatable :: mat(:)
        integer :: u
        call setup(geo, reshape([4, 3, 1, 5, 3, 1], [3, 2]), mat, [1, 1])
        open(newunit=u, file='INPUT/famtest-F3-bc.txt', status='replace', action='write')
        call writeCopy(u, 1, 4, 3, 1, 401, 0.1_real64, 1.0e-6_real64, 0.0_real64)
        call writeCopy(u, 1, 4, 3, 1, 401, 0.2_real64, 2.0e-6_real64, 0.0_real64)
        call writeCopy(u, 2, 5, 3, 1, 401, 0.3_real64, 3.0e-6_real64, 0.0_real64)
        call writeCopy(u, 2, 5, 3, 1, 401, 0.4_real64, 4.0e-6_real64, 0.0_real64)
        close(u)
        call read_cdp_bc_file('famtest-F3-', mat, geo, gas)
        call expectInlet('F3 block 1', geo(1), [0.1_real64, 0.2_real64], [1.0e-6_real64, 2.0e-6_real64], 0.3_real64)
        call expectTags('F3 block 1', geo(1), 3)
        call expectInlet('F3 block 2', geo(2), [0.3_real64, 0.4_real64], [3.0e-6_real64, 4.0e-6_real64], 0.7_real64)
        call expectTags('F3 block 2', geo(2), 3)
    end subroutine f3_two_blocks

    !> F4 -- a 402 inlet: mdotPart sums gp * area over the families' own gp.
    subroutine f4_mass_flux
        type(obj_block)                 :: geo(1)
        type(obj_flowblock)             :: gas(1)
        type(obj_material), allocatable :: mat(:)
        integer :: u, m
        call setup(geo, reshape([4, 3, 1], [3, 1]), mat, [1, 1])
        open(newunit=u, file='INPUT/famtest-F4-bc.txt', status='replace', action='write')
        call writeCopy(u, 1, 4, 3, 1, 402, 0.1_real64, 1.0e-6_real64, 0.0_real64)
        call writeCopy(u, 1, 4, 3, 1, 402, 0.2_real64, 2.0e-6_real64, 0.0_real64)
        close(u)
        call read_cdp_bc_file('famtest-F4-', mat, geo, gas)
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
        type(obj_material), allocatable :: mat(:)
        integer :: u
        call setup(geo, reshape([4, 3, 1], [3, 1]), mat, [1, 1])
        open(newunit=u, file='INPUT/famtest-F5-bc.txt', status='replace', action='write')
        call writeCopy(u, 1, 4, 3, 1, 401, 0.1_real64, 1.0e-6_real64, 0.0_real64)
        call writeCopy(u, 1, 4, 3, 1, 401, 0.2_real64, 2.0e-6_real64, 0.01_real64)
        close(u)
        call read_cdp_bc_file('famtest-F5-', mat, geo, gas)
        call expectEq('F5 fam 1 ds', geo(1)%face(1)%cell(1, 1)%properties(1, 9), 0.0_real64)
        call expectEq('F5 fam 2 ds', geo(1)%face(1)%cell(1, 1)%properties(2, 9), 0.01_real64)
        if (dsSwitch) then
            print '(a)', '  [ok]   famtest F5 dsSwitch: got T want T'
        else
            print '(a)', '  [FAIL] famtest F5 dsSwitch: got F want T'
            nfail = nfail + 1
        endif
    end subroutine f5_ds


    !> F6 -- chimera (102) facelets: face 2 of block 1 (4 x 3 x 1) names cells of block 2, face 1
    !  of block 2 (5 x 3 x 1) cells of block 1; two copies, the second repeating the first.
    subroutine f6_chimera
        type(obj_block)                 :: geo(2)
        type(obj_flowblock)             :: gas(2)
        type(obj_material), allocatable :: mat(:)
        integer :: u, c, m, b, f, n102(2), n400(2)
        call setup(geo, reshape([4, 3, 1, 5, 3, 1], [3, 2]), mat, [1, 1])
        open(newunit=u, file='INPUT/famtest-F6-bc.txt', status='replace', action='write')
        do c = 1, 2
            call writeChimCopy(u, 1, 4, 3, 1, 2)
        enddo
        do c = 1, 2
            call writeChimCopy(u, 2, 5, 3, 1, 1)
        enddo
        close(u)
        call read_cdp_bc_file('famtest-F6-', mat, geo, gas)
        call expectChim('F6 block 1 face 2 m=1', geo(1)%face(2)%cell(1,1), [1, 1], &
                        reshape([2, 1, 1, 1], [1, 4]), [1.0_real64])
        call expectChim('F6 block 1 face 2 m=2', geo(1)%face(2)%cell(2,1), [2, 1], &
                        reshape([2, 2, 1, 2, 2, 2, 1, 1], [2, 4]), [0.25_real64, 0.75_real64])
        call expectChim('F6 block 1 face 2 m=3', geo(1)%face(2)%cell(3,1), [0, 1], &
                        reshape([integer ::], [0, 4]), [real(real64) ::])
        do m = 1, 3
            call expectChim('F6 block 2 face 1 m='//itoa(m), geo(2)%face(1)%cell(m,1), [1, 0], &
                            reshape([1, 4, m, 1], [1, 4]), [1.0_real64])
        enddo
        do b = 1, 2
            n102(b) = 0; n400(b) = 0
            do f = 1, 6
                n102(b) = n102(b) + count(geo(b)%face(f)%cell(:,:)%bcdef == 102)
                n400(b) = n400(b) + count(geo(b)%face(f)%cell(:,:)%bcdef == 400)
            enddo
        enddo
        if (all(n102 == 3) .and. n400(1) == 35 .and. n400(2) == 43) then
            print '(a)', '  [ok]   famtest F6 tags: 3 x 102 and the rest 400 on both blocks'
        else
            print '(a,4(i0,1x))', '  [FAIL] famtest F6 tags (102 per block, 400 per block): ', n102, n400
            nfail = nfail + 1
        endif
        call f6b_chimera_links(geo)
    end subroutine f6_chimera

    !> F6b -- the mollifier's link table on F6's blocks: one link per 102 facelet with layer-1 donors,
    !  receiver (block, face, cell), dir 1, nd and d0 into the flat donor list, in face-walk order.
    subroutine f6b_chimera_links(geo)
        type(obj_block), intent(in) :: geo(:)
        type(cellLink), allocatable :: links(:)
        type(linkDonors) :: donors
        !> expected links: block, face, i, j, nd, d0; expected donors: block, i, j, k, weight
        integer, parameter :: EL(6, 5) = reshape([1, 2, 4, 1, 1, 0,   1, 2, 4, 2, 2, 1, &
                                                  2, 1, 1, 1, 1, 3,   2, 1, 1, 2, 1, 4,   2, 1, 1, 3, 1, 5], [6, 5])
        integer, parameter :: ED(4, 6) = reshape([2, 1, 1, 1,   2, 1, 2, 1,   2, 2, 2, 1, &
                                                  1, 4, 1, 1,   1, 4, 2, 1,   1, 4, 3, 1], [4, 6])
        real(real64), parameter :: EW(6) = [1.0_real64, 0.25_real64, 0.75_real64, 1.0_real64, 1.0_real64, 1.0_real64]
        integer :: l, nbad
        call buildInterfaceLinks(geo, links, donors)
        nbad = 0
        if (size(links) /= 5 .or. size(donors%w) /= 6) then
            nbad = 1
        else
            do l = 1, 5
                associate(c => links(l))
                if (any([c%b, c%f, c%i, c%j, c%nd, c%d0] /= EL(:, l)) .or. c%k /= 1 .or. c%dir /= 1) nbad = nbad + 1
                end associate
            enddo
            do l = 1, 6
                if (any([donors%b(l), donors%i(l), donors%j(l), donors%k(l)] /= ED(:, l)) .or. donors%w(l) /= EW(l)) &
                    nbad = nbad + 1
            enddo
        endif
        if (nbad == 0) then
            print '(a)', '  [ok]   famtest F6b link table: 5 links of the 102 facelets with a layer-1 donor, dir 1, '// &
                '6 donors and weights in file order'
        else
            print '(a,i0,a,i0,a,i0,a)', '  [FAIL] famtest F6b link table: ', size(links), ' links (expected 5), ', &
                size(donors%w), ' donors (expected 6), ', nbad, ' entries off'
            nfail = nfail + 1
        endif
    end subroutine f6b_chimera_links


    !> F7 -- 201 connection lines on a one-layer block (4 x 3 x 1): nine integers on face 3, eight on
    !  face 4; partner cell, face and orientation as written.
    subroutine f7_connection_2d
        type(obj_block)                 :: geo(1)
        type(obj_flowblock)             :: gas(1)
        type(obj_material), allocatable :: mat(:)
        integer :: u, f, m, n, i, j, mend(6), nend(6)
        call setup(geo, reshape([4, 3, 1], [3, 1]), mat, [1, 1])
        mend(1:2) = 3; nend(1:2) = 1; mend(3:4) = 4; nend(3:4) = 1; mend(5:6) = 4; nend(5:6) = 3
        open(newunit=u, file='INPUT/famtest-F7-bc.txt', status='replace', action='write')
        do f = 1, 6; do n = 1, nend(f); do m = 1, mend(f)
            select case (f)
            case (1); i = 1; j = m
            case (2); i = 4; j = m
            case (3); i = m; j = 1
            case (4); i = m; j = 3
            case (5, 6); i = m; j = n
            end select
            if (f == 1) then
                write(u, '(6i8)') 1, i, j, 1, f, 401
                write(u, '(2(2x,es12.5),2(3x,a),3(2x,es12.5),3x,a,2x,es12.5)') 0.1_real64, 0.1_real64, &
                    'normal,', 'normal,', 1.0_real64, 1.0e-6_real64, 0.0_real64, 'Dirac', 0.0_real64
            elseif (f == 3) then
                write(u, '(6i8)') 1, i, j, 1, f, 201
                write(u, '(9i8)') 1, m, 3, 1, 4, 1, 0, 0, 1
            elseif (f == 4) then
                write(u, '(6i8)') 1, i, j, 1, f, 201
                write(u, '(8i8)') 1, m, 1, 3, 1, 0, 0, 1
            else
                write(u, '(6i8)') 1, i, j, 1, f, 400
            endif
        enddo; enddo; enddo
        close(u)
        call read_cdp_bc_file('famtest-F7-', mat, geo, gas)
        do m = 1, 4
            call expectConn('F7 face 3 m='//itoa(m)//' (nine integers)', geo(1)%face(3)%cell(m,1), &
                            [1, m, 3, 1], 4, [1, 0, 0, 1])
            call expectConn('F7 face 4 m='//itoa(m)//' (eight integers)', geo(1)%face(4)%cell(m,1), &
                            [1, m, 1, 1], 3, [1, 0, 0, 1])
        enddo
        call f8_periodic_links(geo)
    end subroutine f7_connection_2d

    !> F8 -- the mollifier's link table on F7's block: one link per 201 face cell, joined to the cell
    !  of its column on the other face, both faces normal to index direction 2.
    subroutine f8_periodic_links(geo)
        type(obj_block), intent(in) :: geo(:)
        type(cellLink), allocatable :: links(:)
        type(linkDonors) :: donors
        integer :: l, nbad, n3, n4
        call buildInterfaceLinks(geo, links, donors)
        nbad = 0; n3 = 0; n4 = 0
        do l = 1, size(links)
            associate(c => links(l))
            if (c%b /= 1 .or. c%pb /= 1 .or. c%dir /= 2 .or. c%k /= 1 .or. c%pk /= 1 .or. c%pi /= c%i) then
                nbad = nbad + 1
            elseif (c%f == 3 .and. c%j == 1 .and. c%pf == 4 .and. c%pj == 3) then
                n3 = n3 + 1
            elseif (c%f == 4 .and. c%j == 3 .and. c%pf == 3 .and. c%pj == 1) then
                n4 = n4 + 1
            else
                nbad = nbad + 1
            endif
            end associate
        enddo
        if (size(links) == 8 .and. n3 == 4 .and. n4 == 4 .and. nbad == 0) then
            print '(a)', '  [ok]   famtest F8 link table: 8 links, face 3 (i,1) <-> face 4 (i,3) of each column, dir 2'
        else
            print '(a,4(i0,a))', '  [FAIL] famtest F8 link table: ', size(links), ' links (expected 8), ', n3, &
                ' face 3 -> 4, ', n4, ' face 4 -> 3, ', nbad, ' other'
            nfail = nfail + 1
        endif
    end subroutine f8_periodic_links


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
    !  line in drag-stokes's form, every other face 400.
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
                write(u, '(6i8)') b, i, j, k, f, 400
            endif
        enddo; enddo; enddo
    end subroutine writeCopy


    !> One copy of block b's table with face cf a chimera face: its cells carry code 102 and the
    !  payload of F6 (block 1: one, two and no layer-1 donor in block 2; block 2: one donor in
    !  block 1, layer 2 empty); every other face 400.
    subroutine writeChimCopy(u, b, nx, ny, nz, cf)
        integer, intent(in) :: u, b, nx, ny, nz, cf
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
            if (f /= cf) then
                write(u, '(6i8)') b, i, j, k, f, 400
                cycle
            endif
            write(u, '(6i8)') b, i, j, k, f, 102
            if (b == 2) then
                write(u, '(2i8)') 1, 0
                write(u, '(4i8,es24.16)') 1, 4, m, 1, 1.0_real64
            elseif (m == 1) then
                write(u, '(2i8)') 1, 1
                write(u, '(4i8,es24.16)') 2, 1, 1, 1, 1.0_real64
                write(u, '(4i8,es24.16)') 2, 2, 1, 1, 1.0_real64
            elseif (m == 2) then
                write(u, '(2i8)') 2, 1
                write(u, '(4i8,es24.16)') 2, 1, 2, 1, 0.25_real64
                write(u, '(4i8,es24.16)') 2, 2, 2, 1, 0.75_real64
                write(u, '(4i8,es24.16)') 2, 3, 2, 1, 1.0_real64
            else
                write(u, '(2i8)') 0, 1
                write(u, '(4i8,es24.16)') 2, 1, 3, 1, 1.0_real64
            endif
        enddo; enddo; enddo
    end subroutine writeChimCopy


    ! ------------------------------------------------------------------ report

    !> A 101/201 cell: partner cell, partner face and orientation integers.
    subroutine expectConn(lab, cell, conn, face, orient)
        use IGLOO_data_block, only: obj_bc_cell
        character(*),      intent(in) :: lab
        type(obj_bc_cell), intent(in) :: cell
        integer,           intent(in) :: conn(4), face, orient(4)
        if (all(cell%connection == conn) .and. cell%connectionFace == face .and. all(cell%orient == orient)) then
            print '(a,4(i0,1x),a,i0,a,4(i0,1x))', '  [ok]   famtest '//lab//': partner ', conn, 'face ', face, &
                ', orientation ', orient
        else
            print '(a,4(i0,1x),a,i0,a,4(i0,1x))', '  [FAIL] famtest '//lab//': partner ', cell%connection, 'face ', &
                cell%connectionFace, ', orientation ', cell%orient
            nfail = nfail + 1
        endif
    end subroutine expectConn

    !> A chimera facelet: code 102, the two layer counts, and the layer-1 donors and weights in order.
    subroutine expectChim(lab, cell, nchi, donor, w)
        use IGLOO_data_block, only: obj_bc_cell
        character(*),      intent(in) :: lab
        type(obj_bc_cell), intent(in) :: cell
        integer,           intent(in) :: nchi(2), donor(:,:)
        real(real64),      intent(in) :: w(:)
        logical :: ok
        ok = cell%bcdef == 102 .and. all(cell%chimN == nchi)
        if (ok) then
            if (nchi(1) == 0) then
                ok = .not. allocated(cell%chimDonor) .and. .not. allocated(cell%chimW)
            else
                ok = allocated(cell%chimDonor) .and. allocated(cell%chimW)
                if (ok) ok = all(shape(cell%chimDonor) == shape(donor)) .and. size(cell%chimW) == size(w)
                if (ok) ok = all(cell%chimDonor == donor) .and. all(abs(cell%chimW - w) <= TOL)
            endif
        endif
        if (ok) then
            print '(a,2(i0,a))', '  [ok]   famtest '//lab//': 102, counts ', nchi(1), ' and ', nchi(2), &
                ', layer-1 donors and weights as written'
        else
            print '(a,i0,a,2(i0,1x))', '  [FAIL] famtest '//lab//': bcdef ', cell%bcdef, ', counts ', cell%chimN
            if (allocated(cell%chimDonor)) print '(a,*(i0,1x))', '         donors ', cell%chimDonor
            if (allocated(cell%chimW))     print '(a,*(es12.5,1x))', '         weights ', cell%chimW
            nfail = nfail + 1
        endif
    end subroutine expectChim


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

    !> Inlet cells of blk: exactly n401 cells tagged 401, all on face 1; every other cell 400.
    subroutine expectTags(lab, blk, n401)
        character(*),    intent(in) :: lab
        type(obj_block), intent(in) :: blk
        integer,         intent(in) :: n401
        integer :: f, per(6), got400, want400
        got400 = 0; want400 = 0
        do f = 1, 6
            per(f)  = count(blk%face(f)%cell(:,:)%bcdef == 401)
            got400  = got400 + count(blk%face(f)%cell(:,:)%bcdef == 400)
            want400 = want400 + size(blk%face(f)%cell)
        enddo
        want400 = want400 - n401
        if (per(1) == n401 .and. sum(per(2:6)) == 0 .and. got400 == want400) then
            print '(a,2(i0,a))', '  [ok]   famtest '//lab//' tags: ', per(1), ' x 401 on face 1, ', got400, ' x 400'
        else
            print '(a,6(i0,1x),a,i0,a,2(i0,a))', '  [FAIL] famtest '//lab//' tags: 401 per face ', per, &
                '; ', got400, ' x 400; want ', n401, ' x 401 on face 1, ', want400, ' x 400'
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
