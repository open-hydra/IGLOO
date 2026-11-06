program test_mollify_blocks
    !
    ! The field mollifier over several blocks joined by 101, 201 and 102 links (mollifyBlocks).
    !
    ! A one-block grid with a random field q in [-1, 1] and random cell volumes V in [0.5, 1.5] is
    ! cut into blocks along whole index planes; every face cell on a cut is linked to the cell
    ! across it (both ways, as buildInterfaceLinks lists a conforming 101 interface). The blocks
    ! go through mollifyBlocks, the uncut grid through binomial_smooth, 8 passes each.
    !
    !   U1a  every cut of the grid: each block's cells equal the one-block cells to 1e-13 of the
    !        field scale; 3D cuts along x, y and z, three blocks along x, 2 x 2 blocks in x-y, a
    !        grid above the OpenMP threshold, and 2D (Nz = 1) cuts along x and y. The arithmetic
    !        is the one-block arithmetic: built with -fp-model strict every cut is bit-identical;
    !        the default build evaluates the update differently in the two loops (a few ULP)
    !   U1b  Sum(q V) over the blocks after the passes equals the sum before to 1e-13 of
    !        Sum(|q| V)
    !   U1c  no aligned link: mollifyBlocks equals binomial_smooth block by block, bit for bit,
    !        with no link at all and with links between faces normal to different directions
    !        (dir = 0, zero flux)
    !   U1d  a periodic pair (201): every face-3 cell of a block linked to the face-4 cell of its
    !        column in the same block and back; the oracle is the test's own smoother with j
    !        wrapped (zero flux in i and k), to 1e-13 of the field scale: one block in 3D and in 2D,
    !        and the 3D grid cut along x into two blocks, each with its own pair (101 and 201 links
    !        together); Sum(q V) to 1e-13 of Sum(|q| V)
    !   U-X1 the exact intersection of two hexahedra (24 tetrahedra each, clipped pairwise): two cubes of edge
    !        H offset by s H give (1 - s_x)(1 - s_y)(1 - s_z) H^3 to 1e-14 of H^3, touching and apart 0
    !   U-X2 a unit cube and itself turned by 45 degrees about z: 2(sqrt 2 - 1) to 1e-14
    !   U-X3 a warped 3 x 3 x 3 block and a warped cell inside it: the overlaps with the 27 cells, and with
    !        their 216 octants, sum to the cell's V24 to 1e-13 (V24 - cellVol printed, not asserted)
    !   U-X4 a block shifted by 0.3 cell: an inside cell's weights 0.3 and 0.7 on the right targets (1e-14)
    !   U-X5 coincident cells: W = 1 exactly, one target; one node moved by an ULP: W = 1 by the snap
    !   U-X6 chordInside: the fraction of a chord inside the other block's cells (closed forms, 1e-14)
    !   U-X7 overlapStatus: 0 valid, 1 a non-positive volume (an inverted cell), 2 a block's cells overlapping
    !   U-X8 timing only: every candidate pair of two NG^3 boxes offset by (0.3, 0.2, 0.1) cells (NG = 12,
    !        IGLOO_UX8_N sets it), the pair count and the wall time
    !   U-T  the overlap map's tie rule at a chimera face: two blocks of a regular grid overlapping
    !        along x, block 2 shifted by half a cell, so that a column of each block has its centroids
    !        on the other block's chimera (102) face: neither is covered by the other block (coverage
    !        count 1), while centroids on the other block's interior faces are (inclusive); moved inside
    !        by a quarter of the tie tolerance tau = max(1e-12 h, 64 eps L) of the cell holding them
    !        they are still on the face, by 2 tau they are covered
    !   U-M  (i) the chimera overlap mirror (IGLOO_overlap): one 3D grid of warped cells cut into two
    !        blocks that overlap by four columns (the overlap map built by buildOverlapMap from the
    !        geometry), random owner-only source and eulerian fields (each overlap cell's deposit in
    !        the block that owns it, 0 in the other); after mirrorSource each block's source equals
    !        the grid's on its range bit for bit, after mirrorEuler its eulerian densities to 1e-14 of
    !        the field scale (the two copies of a warped cell can differ in volume by an ULP, and the
    !        amount moves as q V_source/V_target), the coverage count is 2 in the overlap and 1
    !        elsewhere, the mirrored totals equal the other copy's owner deposits and their sum equals
    !        Sum (n - 1) S over the owner-only fields, to 1e-14 of Sum |S|
    !        (ii) the grid cut into two blocks that overlap by whole columns, each holding the grid's
    !        values on its range (the mirror's result on coincident cells); block 1's face-2 cells and
    !        block 2's face-1 cells are chimera (102) links with one donor of weight 1, the cell across
    !        the face in the other block: each block equals the uncut grid smoothed by binomial_smooth
    !        to 1e-13 of the field scale (3D, 2D, a one-column overlap, above the OpenMP threshold)
    !        (iii) the guard's identity: each block's Sum(q V) changes by exactly what its 102 cells
    !        took from the donor halo (mollifyBlocks' exchanged), to 1e-13 of Sum(|q| V)
    !
    use, intrinsic :: iso_fortran_env, only: real64
    use IGLOO_Lib_Mollify, only: binomial_smooth, mollifyBlocks, cellLink, blockField, linkDonors
    implicit none

    integer, parameter :: NPASS = 8
    real(real64), parameter :: TOL_CUT = 1.0e-13_real64, TOL_SUM = 1.0e-13_real64
    integer(8) :: seed
    integer :: nfail

    nfail = 0
    seed = 20261014_8

    !               label           Nx  Ny  Nz  cuts in i   cuts in j   cuts in k
    call cutCase('3D x',           12,  9,  8, [6],        [integer::], [integer::])
    call cutCase('3D y',           12,  9,  8, [integer::], [4],        [integer::])
    call cutCase('3D z',           12,  9,  8, [integer::], [integer::], [5])
    call cutCase('3D x x',         12,  9,  8, [4, 8],     [integer::], [integer::])
    call cutCase('3D x-y',         12,  9,  8, [5],        [3],        [integer::])
    call cutCase('3D x OpenMP',    80, 12, 12, [40],       [integer::], [integer::])
    call cutCase('2D x',           14, 11,  1, [7],        [integer::], [integer::])
    call cutCase('2D y',           14, 11,  1, [integer::], [6],        [integer::])
    call noLinkCases()
    !                    label            Nx  Ny  Nz  cuts in i
    call periodicCase('3D periodic y',     10,  7,  6, [integer::])
    call periodicCase('2D periodic y',     12,  9,  1, [integer::])
    call periodicCase('3D periodic y, x',  12,  7,  6, [6])
    call clipCases()
    call tieCase()
    call mirrorCase()
    !                   label                Nx  Ny  Nz  overlap columns
    call chimeraCase('3D',                   12,  9,  8,  5,  8)
    call chimeraCase('2D',                   14, 11,  1,  6,  9)
    call chimeraCase('3D one-column overlap', 10,  7,  6,  5,  5)
    call chimeraCase('3D OpenMP',            80, 12, 12, 37, 44)

    if (nfail == 0) then
        print '(a)', ''
        print '(a)', '[PASS] test_mollify_blocks: blocks joined by 101, 201 and 102 links smooth as one block'
    else
        print '(a)', ''
        print '(a,i0,a)', '[FAIL] test_mollify_blocks: ', nfail, ' assertion(s) failed'
        error stop 1
    endif

contains

    !> U1a + U1b on one cut of an nx x ny x nz grid (block-local index ranges from the cut planes).
    subroutine cutCase(lab, nx, ny, nz, ci, cj, ck)
        character(*), intent(in) :: lab
        integer,      intent(in) :: nx, ny, nz, ci(:), cj(:), ck(:)
        real(real64), allocatable :: qg(:,:,:), vg(:,:,:), q1(:,:,:)
        type(blockField), allocatable :: bq(:), bv(:)
        type(cellLink),   allocatable :: links(:)
        integer, allocatable :: lo(:,:), hi(:,:)
        integer :: ndim, b
        real(real64) :: scale, worst, s0, s1, l1
        integer :: nbit

        call randomGrid(nx, ny, nz, qg, vg)
        ndim = merge(2, 3, nz == 1)
        call cutBlocks(nx, ny, nz, ci, cj, ck, lo, hi)
        call linkCuts(lo, hi, links)
        allocate(bq(size(lo,2)), bv(size(lo,2)))
        s0 = 0._real64
        do b = 1, size(lo,2)
            bq(b)%q   = qg(lo(1,b):hi(1,b), lo(2,b):hi(2,b), lo(3,b):hi(3,b))
            bv(b)%q   = vg(lo(1,b):hi(1,b), lo(2,b):hi(2,b), lo(3,b):hi(3,b))
            s0 = s0 + sum(bq(b)%q*bv(b)%q)
        enddo
        q1 = qg
        call binomial_smooth(q1, vg, NPASS, ndim)
        call mollifyBlocks(bq, bv, links, NPASS, ndim)

        scale = maxval(abs(q1)); worst = 0._real64; nbit = 0; s1 = 0._real64; l1 = 0._real64
        do b = 1, size(lo,2)
            associate(r => q1(lo(1,b):hi(1,b), lo(2,b):hi(2,b), lo(3,b):hi(3,b)))
            worst = max(worst, maxval(abs(bq(b)%q - r))/scale)
            nbit  = nbit + count(bq(b)%q /= r)
            end associate
            s1 = s1 + sum(bq(b)%q*bv(b)%q)
            l1 = l1 + sum(abs(bq(b)%q)*bv(b)%q)
        enddo
        if (worst <= TOL_CUT) then
            print '(a,i0,a,i0,a,es9.2,a,i0,a)', '  [ok]   U1a '//lab//': ', size(lo,2), ' blocks, ', size(links), &
                ' links; worst ', worst, ' of the field scale (', nbit, ' cells not bit-identical)'
        else
            print '(a,i0,a,i0,a,es9.2,a,i0,a)', '  [FAIL] U1a '//lab//': ', size(lo,2), ' blocks, ', size(links), &
                ' links; worst ', worst, ' of the field scale (', nbit, ' cells differ)'
            nfail = nfail + 1
        endif
        if (abs(s1 - s0) <= TOL_SUM*l1) then
            print '(a,es9.2,a)', '  [ok]   U1b '//lab//': Sum(q V) drift ', abs(s1 - s0)/l1, ' of Sum(|q| V)'
        else
            print '(a,es9.2,a)', '  [FAIL] U1b '//lab//': Sum(q V) drift ', abs(s1 - s0)/l1, ' of Sum(|q| V)'
            nfail = nfail + 1
        endif
    end subroutine cutCase

    !> U1c: two blocks along x, no aligned link -- each block smoothed on its own.
    subroutine noLinkCases()
        real(real64), allocatable :: qg(:,:,:), vg(:,:,:)
        type(blockField), allocatable :: bq(:), bv(:), br(:)
        type(cellLink),   allocatable :: links(:)
        integer, allocatable :: lo(:,:), hi(:,:)
        integer :: b, pass, l, nbit

        do pass = 1, 2
            call randomGrid(12, 9, 8, qg, vg)
            call cutBlocks(12, 9, 8, [6], [integer::], [integer::], lo, hi)
            if (pass == 1) then
                allocate(links(0))
            else
                call linkCuts(lo, hi, links)
                do l = 1, size(links)
                    links(l)%dir = 0          ! faces normal to different directions: zero flux
                enddo
            endif
            allocate(bq(2), bv(2), br(2))
            do b = 1, 2
                bq(b)%q   = qg(lo(1,b):hi(1,b), lo(2,b):hi(2,b), lo(3,b):hi(3,b))
                bv(b)%q = vg(lo(1,b):hi(1,b), lo(2,b):hi(2,b), lo(3,b):hi(3,b))
                br(b)%q   = bq(b)%q
                call binomial_smooth(br(b)%q, bv(b)%q, NPASS, 3)
            enddo
            call mollifyBlocks(bq, bv, links, NPASS, 3)
            nbit = count(bq(1)%q /= br(1)%q) + count(bq(2)%q /= br(2)%q)
            if (nbit == 0) then
                print '(a,i0,a)', '  [ok]   U1c '//trim(merge('no link         ', 'non-aligned link', pass == 1))// &
                    ': bit-identical to binomial_smooth per block (', size(links), ' links)'
            else
                print '(a,i0,a)', '  [FAIL] U1c '//trim(merge('no link         ', 'non-aligned link', pass == 1))// &
                    ': ', nbit, ' cells differ from binomial_smooth per block'
                nfail = nfail + 1
            endif
            deallocate(links, bq, bv, br)
        enddo
    end subroutine noLinkCases

    !> U1d: the grid cut along x after the listed planes, each block's face-3 and face-4 cells linked
    !  to each other (a periodic pair in j), against the uncut grid smoothed with j wrapped.
    subroutine periodicCase(lab, nx, ny, nz, ci)
        character(*), intent(in) :: lab
        integer,      intent(in) :: nx, ny, nz, ci(:)
        real(real64), allocatable :: qg(:,:,:), vg(:,:,:), q1(:,:,:)
        type(blockField), allocatable :: bq(:), bv(:)
        type(cellLink),   allocatable :: links(:), cut(:)
        integer, allocatable :: lo(:,:), hi(:,:)
        integer :: ndim, b, i, k, n, nbit
        real(real64) :: scale, worst, s0, s1, l1

        call randomGrid(nx, ny, nz, qg, vg)
        ndim = merge(2, 3, nz == 1)
        call cutBlocks(nx, ny, nz, ci, [integer::], [integer::], lo, hi)
        call linkCuts(lo, hi, cut)
        allocate(links(size(cut) + 2*sum((hi(1,:)-lo(1,:)+1)*(hi(3,:)-lo(3,:)+1))))
        links(1:size(cut)) = cut
        n = size(cut)
        do b = 1, size(lo,2)
            do k = 1, hi(3,b)-lo(3,b)+1; do i = 1, hi(1,b)-lo(1,b)+1
                n = n + 1
                links(n) = cellLink(b=b, f=3, i=i, j=1,  k=k, pb=b, pf=4, pi=i, pj=ny, pk=k, dir=2)
                n = n + 1
                links(n) = cellLink(b=b, f=4, i=i, j=ny, k=k, pb=b, pf=3, pi=i, pj=1,  pk=k, dir=2)
            enddo; enddo
        enddo
        allocate(bq(size(lo,2)), bv(size(lo,2)))
        s0 = 0._real64
        do b = 1, size(lo,2)
            bq(b)%q = qg(lo(1,b):hi(1,b), lo(2,b):hi(2,b), lo(3,b):hi(3,b))
            bv(b)%q = vg(lo(1,b):hi(1,b), lo(2,b):hi(2,b), lo(3,b):hi(3,b))
            s0 = s0 + sum(bq(b)%q*bv(b)%q)
        enddo
        q1 = qg
        call periodicSmooth(q1, vg, NPASS, ndim)
        call mollifyBlocks(bq, bv, links, NPASS, ndim)

        scale = maxval(abs(q1)); worst = 0._real64; nbit = 0; s1 = 0._real64; l1 = 0._real64
        do b = 1, size(lo,2)
            associate(r => q1(lo(1,b):hi(1,b), lo(2,b):hi(2,b), lo(3,b):hi(3,b)))
            worst = max(worst, maxval(abs(bq(b)%q - r))/scale)
            nbit  = nbit + count(bq(b)%q /= r)
            end associate
            s1 = s1 + sum(bq(b)%q*bv(b)%q)
            l1 = l1 + sum(abs(bq(b)%q)*bv(b)%q)
        enddo
        if (worst <= TOL_CUT .and. abs(s1 - s0) <= TOL_SUM*l1) then
            print '(a,i0,a,i0,a,es9.2,a,i0,a,es9.2,a)', '  [ok]   U1d '//lab//': ', size(lo,2), ' block(s), ', &
                size(links), ' links; worst ', worst, ' of the field scale (', nbit, &
                ' cells not bit-identical); Sum(q V) drift ', abs(s1 - s0)/l1, ' of Sum(|q| V)'
        else
            print '(a,i0,a,i0,a,es9.2,a,i0,a,es9.2,a)', '  [FAIL] U1d '//lab//': ', size(lo,2), ' block(s), ', &
                size(links), ' links; worst ', worst, ' of the field scale (', nbit, &
                ' cells differ); Sum(q V) drift ', abs(s1 - s0)/l1, ' of Sum(|q| V)'
            nfail = nfail + 1
        endif
    end subroutine periodicCase

    !> The oracle of U1d: binomial_smooth's passes with j wrapped (cell 1 and cell Ny neighbours), zero
    !  flux in i and k; the flux theta*min(V, V')*(q' - q) of every face.
    subroutine periodicSmooth(q, V, npass, ndim)
        real(real64), intent(inout) :: q(:,:,:)
        real(real64), intent(in)    :: V(:,:,:)
        integer,      intent(in)    :: npass, ndim
        real(real64), allocatable :: s(:,:,:)
        real(real64), parameter :: th = 0.25_real64
        real(real64) :: phiL, phiR
        integer :: nx, ny, nz, i, j, k, p, d, jp, jm
        nx = size(q,1); ny = size(q,2); nz = size(q,3)
        do p = 1, npass
            do d = 1, ndim
                s = q
                do k = 1, nz; do j = 1, ny; do i = 1, nx
                    phiR = 0._real64; phiL = 0._real64
                    select case (d)
                    case (1)
                        if (i < nx) phiR = th*min(V(i,j,k), V(i+1,j,k))*(s(i+1,j,k) - s(i,j,k))
                        if (i > 1)  phiL = th*min(V(i-1,j,k), V(i,j,k))*(s(i,j,k) - s(i-1,j,k))
                    case (2)
                        jp = merge(1, j+1, j == ny); jm = merge(ny, j-1, j == 1)
                        phiR = th*min(V(i,j,k), V(i,jp,k))*(s(i,jp,k) - s(i,j,k))
                        phiL = th*min(V(i,jm,k), V(i,j,k))*(s(i,j,k) - s(i,jm,k))
                    case (3)
                        if (k < nz) phiR = th*min(V(i,j,k), V(i,j,k+1))*(s(i,j,k+1) - s(i,j,k))
                        if (k > 1)  phiL = th*min(V(i,j,k-1), V(i,j,k))*(s(i,j,k) - s(i,j,k-1))
                    end select
                    q(i,j,k) = s(i,j,k) + phiR/V(i,j,k) - phiL/V(i,j,k)
                enddo; enddo; enddo
            enddo
        enddo
    end subroutine periodicSmooth

    !> U-X1..U-X8: the exact intersection of two hexahedra (IGLOO_overlap), its weights, guards and chord split.
    subroutine clipCases()
        use IGLOO_data_block, only: obj_block
        use IGLOO_overlap,    only: cellOverlapVolume, cellVolume24, cellTargets, overlapStatus, chordInside, &
                                    lexVertices, octantVertices, cutRecord, bucketGrid, bucketBuild, bucketCandidates
        use IGLOO_variables,  only: mesh2D
        real(real64), parameter :: H = 2.5e-3_real64, X0(3) = [0.0625_real64, 0.0125_real64, 0.005_real64]
        real(real64), parameter :: PI = 3.141592653589793_real64
        real(real64) :: va(3,8), vb(3,8), vc(3,8,2), w(2), s(3), ov, ex, err, v24, sumOv, rot(2,2)
        real(real64) :: ex6(7), got6(7), tsec
        real(real64), allocatable :: blo(:,:), bhi(:,:)
        type(obj_block) :: blk, one
        type(cutRecord) :: rec
        type(bucketGrid) :: grid
        logical, allocatable :: mark(:)
        integer, allocatable :: cand(:)
        integer :: n, i, j, k, a, b, e, nbad, nc, npair, ng, c0, c1, crate
        character(len=40) :: env

        mesh2D = .false.
        !> U-X1: two cubes of edge H offset by s H
        nbad = 0
        do n = 1, 5
            select case (n)
            case (1); s = [0.3_real64, 0._real64, 0._real64]
            case (2); s = [0.5_real64, 0.5_real64, 0._real64]
            case (3); s = [0.3_real64, 0.2_real64, 0.1_real64]
            case (4); s = [1._real64, 0._real64, 0._real64]
            case (5); s = [1.2_real64, 0._real64, 0._real64]
            end select
            va = boxNodes(X0, [H, H, H]); vb = boxNodes(X0 + s*H, [H, H, H])
            ov = cellOverlapVolume(va, vb)
            ex = product(max(0._real64, 1._real64 - s))*H**3
            err = abs(ov - ex)/H**3
            print '(a,3f5.2,a,es23.15,a,es9.2)', '  [info] U-X1 s =', s, ': |a ∩ b|/H^3 = ', ov/H**3, ', off ', err
            if (err > 1.0e-14_real64) nbad = nbad + 1
        enddo
        call verdict('U-X1 two cubes offset by (0.3,0,0), (0.5,0.5,0), (0.3,0.2,0.1), touching, apart: the closed form '// &
                     'to 1e-14 of H^3', nbad == 0, nbad)

        !> U-X2: a unit cube and itself rotated by 45 degrees about its z axis
        va = boxNodes([-0.5_real64, -0.5_real64, -0.5_real64], [1._real64, 1._real64, 1._real64])
        rot = reshape([cos(PI/4), sin(PI/4), -sin(PI/4), cos(PI/4)], [2, 2])
        vb = va
        do e = 1, 8
            vb(1:2, e) = matmul(rot, va(1:2, e))
        enddo
        ov = cellOverlapVolume(va, vb)
        call verdict('U-X2 unit cube and itself turned by 45 degrees: '//es15(ov)//' against 2(sqrt 2 - 1) = '// &
                     es15(2._real64*(sqrt(2._real64) - 1._real64))//' (1e-14)', &
                     abs(ov - 2._real64*(sqrt(2._real64) - 1._real64)) <= 1.0e-14_real64, 1)

        !> U-X3: a warped 3 x 3 x 3 block and a warped cell inside it; the tiling, by the cells and by their octants
        call warpedBlock(blk, 3, H, 0.1_real64*H, 0)
        va = boxNodes([0.8_real64*H, 0.8_real64*H, 0.8_real64*H], [1.4_real64*H, 1.4_real64*H, 1.4_real64*H])
        do e = 1, 8
            va(:, e) = va(:, e) + 0.1_real64*H*warp(va(:, e)/H, 1)
        enddo
        v24 = cellVolume24(va)
        sumOv = 0._real64; ov = 0._real64
        do k = 1, 3; do j = 1, 3; do i = 1, 3
            call lexVertices(blk, i, j, k, vb)
            sumOv = sumOv + cellOverlapVolume(va, vb)
            do c0 = 0, 7
                call octantVertices(blk, i, j, k, iand(c0, 1), iand(c0/2, 1), iand(c0/4, 1), vb)
                ov = ov + cellOverlapVolume(va, vb)
            enddo
        enddo; enddo; enddo
        one%Nx = 1; one%Ny = 1; one%Nz = 1
        allocate(one%node(3, 0:1, 0:1, 0:1))
        do k = 0, 1; do j = 0, 1; do i = 0, 1
            one%node(:, i, j, k) = va(:, 1 + i + 2*j + 4*k)
        enddo; enddo; enddo
        call one%precomputeMetric
        print '(a,es9.2,a)', '  [info] U-X3 V24 - cellVol (hexVolume''s five tetrahedra) = ', (v24 - one%cellVol(1,1,1))/v24, &
            ' of V24 (not asserted)'
        call verdict('U-X3 the 27 warped cells tile the warped cell: Sum |d ∩ c| = V24(d) to '//es(abs(sumOv - v24)/v24)// &
                     ' (1e-13)', abs(sumOv - v24) <= 1.0e-13_real64*v24, 1)
        call verdict('U-X3 their 216 octants tile it too: '//es(abs(ov - v24)/v24)//' (1e-13)', &
                     abs(ov - v24) <= 1.0e-13_real64*v24, 1)

        !> U-X4: block 2 shifted by 0.3 cell along x: an inside cell's two weights
        va = boxNodes(X0 + [5._real64, 0._real64, 0._real64]*H, [H, H, H])
        vc(:,:,1) = boxNodes(X0 + [4.3_real64, 0._real64, 0._real64]*H, [H, H, H])
        vc(:,:,2) = boxNodes(X0 + [5.3_real64, 0._real64, 0._real64]*H, [H, H, H])
        call cellTargets(va, vc, w, sumOv)
        nbad = count(abs(w - [0.3_real64, 0.7_real64]) > 1.0e-14_real64)
        vb = vc(:,:,2)
        vc(:,:,1) = boxNodes(X0 + [5._real64, 0._real64, 0._real64]*H, [H, H, H])
        vc(:,:,2) = boxNodes(X0 + [6._real64, 0._real64, 0._real64]*H, [H, H, H])
        call cellTargets(vb, vc, w, ov)
        nbad = nbad + count(abs(w - [0.7_real64, 0.3_real64]) > 1.0e-14_real64)
        call verdict('U-X4 a 0.3-cell shift: a block-1 cell''s weights 0.3, 0.7 into block 2''s cells k - 1, k, a '// &
                     'block-2 cell''s 0.7, 0.3 into block 1''s (1e-14)', nbad == 0, nbad)

        !> U-X5: coincident cells, then one node moved by an ULP
        va = boxNodes(X0, [H, H, H])
        vc(:,:,1) = va
        vc(:,:,2) = boxNodes(X0 + [H, 0._real64, 0._real64], [H, H, H])
        call cellTargets(va, vc, w, sumOv)
        nbad = merge(0, 1, w(1) == 1._real64 .and. w(2) == 0._real64)
        va(1, 8) = nearest(va(1, 8), 1._real64)
        call cellTargets(va, vc, w, sumOv)
        nbad = nbad + merge(0, 1, w(1) == 1._real64 .and. w(2) == 0._real64)
        print '(a,2es24.16)', '  [info] U-X5 a node moved by an ULP: W = ', w
        call verdict('U-X5 coincident cells: W = 1 exactly, one target; a node moved by an ULP: W = 1 exactly by the '// &
                     'snap, one target', nbad == 0, nbad)

        !> U-X6: the chord split against one cell, then against two cells with a gap
        call makeRecord(rec, [X0], [1])
        ex6 = [1._real64, 0._real64, 0.7_real64, 0._real64, 1._real64, 0._real64, 0.5_real64]
        got6(1) = chordInside(X0 + [0.2_real64, 0.5_real64, 0.5_real64]*H, X0 + [0.8_real64, 0.5_real64, 0.5_real64]*H, rec)
        got6(2) = chordInside(X0 + [1.2_real64, 0.5_real64, 0.5_real64]*H, X0 + [1.8_real64, 0.5_real64, 0.5_real64]*H, rec)
        got6(3) = chordInside(X0 + [-0.3_real64, 0.5_real64, 0.5_real64]*H, X0 + [0.7_real64, 0.5_real64, 0.5_real64]*H, rec)
        got6(4) = chordInside(X0 + [0._real64, 0.5_real64, 0.5_real64]*H + [0.5_real64*rec%tol, 0._real64, 0._real64], &
                              X0 + [-0.5_real64, 0.5_real64, 0.5_real64]*H, rec)
        got6(5) = chordInside(X0 + [0.5_real64, 0.5_real64, 0.5_real64]*H, X0 + [0.5_real64, 0.5_real64, 0.5_real64]*H, rec)
        got6(6) = chordInside(X0 + [1.5_real64, 0.5_real64, 0.5_real64]*H, X0 + [1.5_real64, 0.5_real64, 0.5_real64]*H, rec)
        call makeRecord(rec, [X0, X0 + [2._real64*H, 0._real64, 0._real64]], [1, 2])
        got6(7) = chordInside(X0 + [0.5_real64, 0.5_real64, 0.5_real64]*H, X0 + [2.5_real64, 0.5_real64, 0.5_real64]*H, rec)
        nbad = count(abs(got6 - ex6) > 1.0e-14_real64)
        print '(a,7f8.4)', '  [info] U-X6 lambda: ', got6
        call verdict('U-X6 chordInside: inside 1, outside 0, crossing at 0.3 of the chord 0.7, from the plane outward 0, '// &
                     'zero length inside 1 and outside 0, out and back in 0.5 (1e-14)', nbad == 0, nbad)

        !> U-X7: the guards
        va = boxNodes(X0, [2._real64*H, H, H])
        vc(:,:,1) = boxNodes(X0, [H, H, H])
        vc(:,:,2) = boxNodes(X0 + [H, 0._real64, 0._real64], [H, H, H])
        nbad = merge(0, 1, overlapStatus(cellVolume24(va), cellOverlapVolume(va, vc(:,:,1)) &
                                                          + cellOverlapVolume(va, vc(:,:,2))) == 0)
        vb = va
        do e = 1, 8, 2
            vb(:, e) = va(:, e+1); vb(:, e+1) = va(:, e)
        enddo
        nbad = nbad + merge(0, 1, overlapStatus(cellVolume24(vb), cellOverlapVolume(vb, vc(:,:,1)) &
                                                                 + cellOverlapVolume(vb, vc(:,:,2))) == 1)
        vc(:,:,2) = boxNodes(X0 + [0.7_real64*H, 0._real64, 0._real64], [1.3_real64*H, H, H])
        nbad = nbad + merge(0, 1, overlapStatus(cellVolume24(va), cellOverlapVolume(va, vc(:,:,1)) &
                                                                 + cellOverlapVolume(va, vc(:,:,2))) == 2)
        call verdict('U-X7 overlapStatus: a valid pair 0, an inverted cell 1, two overlapping cells of one block 2', &
                     nbad == 0, nbad)

        !> U-X8 (timing, not asserted): every candidate pair of two NG^3 boxes offset by (0.3, 0.2, 0.1) cells
        ng = 12
        call get_environment_variable('IGLOO_UX8_N', env)
        if (len_trim(env) > 0) read(env, *) ng
        allocate(blo(3, ng**3), bhi(3, ng**3), mark(ng**3))
        do k = 1, ng; do j = 1, ng; do i = 1, ng
            n = i + ng*(j - 1 + ng*(k - 1))
            blo(:, n) = ([i, j, k] - 1 + [0.3_real64, 0.2_real64, 0.1_real64])*H
            bhi(:, n) = blo(:, n) + H
        enddo; enddo; enddo
        call bucketBuild(blo, bhi, grid)
        mark = .false.
        npair = 0; ov = 0._real64
        call system_clock(c0, crate)
        do k = 1, ng; do j = 1, ng; do i = 1, ng
            va = boxNodes(([i, j, k] - 1)*H*[1._real64, 1._real64, 1._real64], [H, H, H])
            call bucketCandidates(grid, va(:,1), va(:,8), mark, cand, nc)
            do a = 1, nc
                b = cand(a)
                ov = ov + cellOverlapVolume(va, boxNodes(blo(:, b), [H, H, H]))
            enddo
            npair = npair + nc
        enddo; enddo; enddo
        call system_clock(c1)
        tsec = real(c1 - c0, real64)/real(crate, real64)
        ex = (real(ng, real64) - 0.3_real64)*(real(ng, real64) - 0.2_real64)*(real(ng, real64) - 0.1_real64)*H**3
        print '(a,i0,a,i0,a,f8.2,a,es9.2,a,es9.2,a)', '  [info] U-X8 ', ng, '^3 cells: ', npair, ' candidate pairs, ', tsec, &
            ' s on one core (', tsec/max(npair, 1), ' s per pair); Sum |a ∩ b| off the closed form by ', abs(ov - ex)/ex, &
            ' (not asserted; IGLOO_UX8_N sets the size)'

    end subroutine clipCases

    !> The lexicographic nodes of the box [lo, lo + ext].
    pure function boxNodes(lo, ext) result(v)
        real(real64), intent(in) :: lo(3), ext(3)
        real(real64) :: v(3,8)
        integer :: a, b, c
        do c = 0, 1; do b = 0, 1; do a = 0, 1
            v(:, 1 + a + 2*b + 4*c) = lo + [a, b, c]*ext
        enddo; enddo; enddo
    end function boxNodes

    !> A smooth displacement field of unit amplitude (two variants).
    pure function warp(p, m) result(d)
        real(real64), intent(in) :: p(3)
        integer,      intent(in) :: m
        real(real64) :: d(3)
        if (m == 0) then
            d = [sin(1.3_real64*p(2) + 0.4_real64*p(3)), cos(0.9_real64*p(3) + 1.1_real64*p(1)), &
                 sin(0.7_real64*p(1) - 1.2_real64*p(2))]
        else
            d = [cos(0.6_real64*p(3) - 1.4_real64*p(2)), sin(1.2_real64*p(1) + 0.5_real64*p(3)), &
                 cos(1.1_real64*p(2) + 0.8_real64*p(1))]
        endif
    end function warp

    !> An n^3 block of cells of edge h, every node moved by amp times warp(x/h, m).
    subroutine warpedBlock(g, n, h, amp, m)
        use IGLOO_data_block, only: obj_block
        type(obj_block), intent(inout) :: g
        integer,         intent(in)    :: n, m
        real(real64),    intent(in)    :: h, amp
        integer :: i, j, k
        g%Nx = n; g%Ny = n; g%Nz = n
        allocate(g%node(3, 0:n, 0:n, 0:n), g%center(3, n, n, n))
        do k = 0, n; do j = 0, n; do i = 0, n
            g%node(:, i, j, k) = [i, j, k]*h
            g%node(:, i, j, k) = g%node(:, i, j, k) + amp*warp(g%node(:, i, j, k)/h, m)
        enddo; enddo; enddo
        call g%compute_geometry
        call g%precomputeMetric
    end subroutine warpedBlock

    !> A cut record whose other block is the cells of edge 2.5e-3 m at lo(:, c): their six face planes (normals
    !  into the cell) and 12-plane data; the tolerance of the first.
    subroutine makeRecord(r, lo, ids)
        use IGLOO_overlap,               only: cutRecord
        use IGLOO_RayFaceIntersection3D, only: guide, computeFacePlane, precomputeHex12
        type(cutRecord), intent(out) :: r
        real(real64),    intent(in)  :: lo(:)
        integer,         intent(in)  :: ids(:)
        real(real64) :: v(3,8), g(3,8), pl(3,2)
        integer :: c, f, nc
        nc = size(ids)
        allocate(r%pn(3, 6*nc), r%pp(3, 6*nc), r%hn(3,2,6,nc), r%hc(3,2,6,nc), r%hd(2,6,nc), r%ht(nc))
        do c = 1, nc
            v = boxNodes(lo(3*c-2:3*c), [2.5e-3_real64, 2.5e-3_real64, 2.5e-3_real64])
            !> getVertices' order from the lexicographic one
            g(:,1) = v(:,1); g(:,2) = v(:,5); g(:,3) = v(:,7); g(:,4) = v(:,3)
            g(:,5) = v(:,2); g(:,6) = v(:,6); g(:,7) = v(:,8); g(:,8) = v(:,4)
            call precomputeHex12(g, r%hn(:,:,:,c), r%hc(:,:,:,c), r%hd(:,:,c), r%ht(c))
            do f = 1, 6
                pl = computeFacePlane(g(:, guide(f,:)))
                r%pn(:, 6*(c-1) + f) = -pl(:,1)
                r%pp(:, 6*(c-1) + f) = pl(:,2)
            enddo
        enddo
        r%t = 2
        r%tol = r%ht(1)
    end subroutine makeRecord

    function es15(x) result(r)
        real(real64), intent(in) :: x
        character(:), allocatable :: r
        character(24) :: buf
        write(buf, '(es23.16)') x
        r = trim(adjustl(buf))
    end function es15

    !> U-T: block 1 = columns 1..I1 of a regular grid of edge H, block 2 = NX2 columns from x = (I0 - 1/2) H
    !  - delta; face 2 of block 1 and face 1 of block 2 are chimera faces. delta = 0: block 1's column I0 and
    !  block 2's column I1 - I0 + 1 have their centroids on the other block's chimera face.
    subroutine tieCase()
        use IGLOO_data_block, only: obj_block
        use IGLOO_overlap,    only: overlapMap, buildOverlapMap
        use IGLOO_variables,  only: mesh2D
        integer, parameter :: NY = 3, NZ = 2, I0 = 5, I1 = 8, NX2 = 6
        real(real64), parameter :: H = 1.0e-2_real64
        type(obj_block), allocatable :: geo(:)
        type(overlapMap) :: ovl
        integer, allocatable :: npair(:,:)
        real(real64) :: tau, delta, x0(2), v(3,8)
        integer :: icase, b, i, j, k, f, nfr, nexp, nbad
        character(len=60) :: lab

        mesh2D = .false.
        !> the tie tolerance of the two cells that hold the fringe centroids (block 2's first, block 1's last)
        tau = 0._real64
        do b = 1, 2
            v(1, :) = merge(I1 - 1, I0 - 1, b == 1)*H + H*[0, 0, 0, 0, 1, 1, 1, 1]
            v(2, :) = H*[0, 0, 1, 1, 0, 0, 1, 1]
            v(3, :) = H*[0, 1, 1, 0, 0, 1, 1, 0]
            if (b == 2) v(1, :) = v(1, :) + 0.5_real64*H
            tau = max(tau, max(1.e-12_real64*norm2(v(:,7) - v(:,1)), 64._real64*epsilon(1._real64)*maxval(abs(v))))
        enddo
        do icase = 1, 3
            select case (icase)
            case (1); delta = 0._real64;          lab = 'on the face'
            case (2); delta = 0.25_real64*tau;    lab = 'inside by tau/4 (the tie band)'
            case (3); delta = 2._real64*tau;      lab = 'inside by 2 tau'
            end select
            x0 = [0._real64, (I0 - 1 + 0.5_real64)*H - delta]
            if (allocated(geo)) deallocate(geo)
            allocate(geo(2))
            do b = 1, 2
                geo(b)%Nx = merge(I1, NX2, b == 1); geo(b)%Ny = NY; geo(b)%Nz = NZ
                allocate(geo(b)%node(3, 0:geo(b)%Nx, 0:NY, 0:NZ), geo(b)%center(3, geo(b)%Nx, NY, NZ))
                do k = 0, NZ; do j = 0, NY; do i = 0, geo(b)%Nx
                    geo(b)%node(:, i, j, k) = [x0(b) + i*H, j*H, k*H]
                enddo; enddo; enddo
                call geo(b)%compute_geometry
                call geo(b)%precomputeMetric
                call geo(b)%computeBBox
                do f = 1, 6
                    geo(b)%face(f)%cell(:,:)%bcdef = 400
                enddo
            enddo
            geo(1)%face(2)%cell(:,:)%bcdef = 102
            geo(2)%face(1)%cell(:,:)%bcdef = 102
            call buildOverlapMap(geo, ovl, npair)
            !> fringe columns: block 1's column I0, block 2's column I1 - I0 + 1
            nfr = 0
            do k = 1, NZ; do j = 1, NY
                nfr = nfr + ovl%blk(1)%n(I0 + (j-1)*I1 + (k-1)*I1*NY) - 1
                nfr = nfr + ovl%blk(2)%n(I1 - I0 + 1 + (j-1)*NX2 + (k-1)*NX2*NY) - 1
            enddo; enddo
            nexp = merge(2*NY*NZ, 0, icase == 3)
            nbad = abs(nfr - nexp) + count(npair /= reshape([0, (I1 - I0 + merge(1, 0, icase == 3))*NY*NZ, &
                                                              (I1 - I0 + merge(1, 0, icase == 3))*NY*NZ, 0], [2, 2]))
            call verdict('U-T fringe centroids '//trim(lab)//': '//trim(merge('covered    ', 'not covered', icase == 3))// &
                         ' (map pairs '//itoaT(npair(2, 1))//' and '//itoaT(npair(1, 2))//', fringe coverage '// &
                         itoaT(nfr)//' of '//itoaT(2*NY*NZ)//')', nbad == 0, nbad)
        enddo
    end subroutine tieCase

    function itoaT(i) result(r)
        integer, intent(in) :: i
        character(:), allocatable :: r
        character(16) :: buf
        write(buf, '(i0)') i
        r = trim(buf)
    end function itoaT

    !> U-M (i): one grid of NX x NY x NZ warped cells, block 1 = columns 1..I1, block 2 = I0..NX; in the
    !  overlap block 1 owns the columns up to ISPLIT and block 2 the rest.
    subroutine mirrorCase()
        use IGLOO_data_block, only: obj_block, obj_sourceblock, obj_eulerblock
        use IGLOO_overlap,    only: overlapMap, mirroredTotals, buildOverlapMap, mirrorSource, mirrorEuler
        use IGLOO_variables,  only: mesh2D
        integer, parameter :: NX = 12, NY = 5, NZ = 4, NMAT = 2, I0 = 5, I1 = 8, ISPLIT = 6
        real(real64), parameter :: H = 1.0e-2_real64
        type(obj_block)       :: geo(2)
        type(obj_sourceblock) :: src(2)
        type(obj_eulerblock)  :: eul(2)
        type(overlapMap)      :: ovl
        type(mirroredTotals), allocatable :: mir(:)
        integer, allocatable  :: npair(:,:)
        real(real64) :: node(3, 0:NX, 0:NY, 0:NZ)
        real(real64) :: sm(NMAT, NX, NY, NZ), so(3, NX, NY, NZ), se(NX, NY, NZ)
        real(real64) :: ed(NX, NY, NZ), en(NX, NY, NZ), ev(3, NX, NY, NZ), et(NX, NY, NZ)
        real(real64) :: expM(NMAT+4, 2), gotM(NMAT+4, 2), sumNS(NMAT+4), l1(NMAT+4), cellS(NMAT+4), werr, vmis
        integer :: lo(2), hi(2), b, i, j, k, g, c, nsrc, neul, ncov, nm4
        logical :: own

        mesh2D = .false.
        nm4 = NMAT + 4
        do k = 0, NZ; do j = 0, NY; do i = 0, NX
            node(:, i, j, k) = [i, j, k]*H
            if (i > 0 .and. i < NX .and. j > 0 .and. j < NY .and. k > 0 .and. k < NZ) &
                node(:, i, j, k) = node(:, i, j, k) + 0.15_real64*H*[2*uniform() - 1, 2*uniform() - 1, 2*uniform() - 1]
        enddo; enddo; enddo
        do k = 1, NZ; do j = 1, NY; do i = 1, NX
            sm(:, i, j, k) = [2*uniform() - 1, 2*uniform() - 1]
            so(:, i, j, k) = [2*uniform() - 1, 2*uniform() - 1, 2*uniform() - 1]
            se(i, j, k)    = 2*uniform() - 1
            ed(i, j, k)    = 0.5_real64 + uniform()
            en(i, j, k)    = 0.5_real64 + uniform()
            ev(:, i, j, k) = [2*uniform() - 1, 2*uniform() - 1, 2*uniform() - 1]
            et(i, j, k)    = 0.5_real64 + uniform()
        enddo; enddo; enddo

        lo = [1, I0]; hi = [I1, NX]
        do b = 1, 2
            geo(b)%Nx = hi(b) - lo(b) + 1; geo(b)%Ny = NY; geo(b)%Nz = NZ
            allocate(geo(b)%node(3, 0:geo(b)%Nx, 0:NY, 0:NZ), geo(b)%center(3, geo(b)%Nx, NY, NZ))
            geo(b)%node = node(:, lo(b)-1:hi(b), :, :)
            call geo(b)%compute_geometry
            call geo(b)%precomputeMetric
            call geo(b)%computeBBox
            allocate(src(b)%sourceMass(NMAT, geo(b)%Nx, NY, NZ), src(b)%sourceMom(3, geo(b)%Nx, NY, NZ), &
                     src(b)%sourceEn(geo(b)%Nx, NY, NZ))
            allocate(eul(b)%density(geo(b)%Nx, NY, NZ), eul(b)%np(geo(b)%Nx, NY, NZ), &
                     eul(b)%velocity(3, geo(b)%Nx, NY, NZ), eul(b)%temperature(geo(b)%Nx, NY, NZ))
            do k = 1, NZ; do j = 1, NY; do i = 1, geo(b)%Nx
                g = lo(b) + i - 1
                own = g < I0 .or. g > I1 .or. (b == 1 .and. g <= ISPLIT) .or. (b == 2 .and. g > ISPLIT)
                src(b)%sourceMass(:, i, j, k) = merge(sm(:, g, j, k), 0._real64, own)
                src(b)%sourceMom(:, i, j, k)  = merge(so(:, g, j, k), 0._real64, own)
                src(b)%sourceEn(i, j, k)      = merge(se(g, j, k), 0._real64, own)
                eul(b)%density(i, j, k)       = merge(ed(g, j, k), 0._real64, own)
                eul(b)%np(i, j, k)            = merge(en(g, j, k), 0._real64, own)
                eul(b)%velocity(:, i, j, k)   = merge(ev(:, g, j, k), 0._real64, own)
                eul(b)%temperature(i, j, k)   = merge(et(g, j, k), 0._real64, own)
            enddo; enddo; enddo
        enddo

        call buildOverlapMap(geo, ovl, npair)
        !> expected mirrored totals (the other block's owner deposits in this block's overlap columns, in
        !  this block's cell order) and Sum (n - 1) S over the owner-only fields
        expM = 0._real64; sumNS = 0._real64; l1 = 0._real64; ncov = 0
        do b = 1, 2
            do k = 1, NZ; do j = 1, NY; do i = 1, geo(b)%Nx
                g = lo(b) + i - 1
                c = i + (j-1)*geo(b)%Nx + (k-1)*geo(b)%Nx*NY
                cellS = [src(b)%sourceMass(:, i, j, k), src(b)%sourceMom(:, i, j, k), src(b)%sourceEn(i, j, k)]
                sumNS = sumNS + (ovl%blk(b)%n(c) - 1)*cellS
                l1 = l1 + abs(cellS)
                if (ovl%blk(b)%n(c) /= merge(2, 1, g >= I0 .and. g <= I1)) ncov = ncov + 1
                if (g >= I0 .and. g <= I1 .and. (b == 1 .eqv. g > ISPLIT)) &
                    expM(:, b) = expM(:, b) + [sm(:, g, j, k), so(:, g, j, k), se(g, j, k)]
            enddo; enddo; enddo
        enddo
        call mirrorSource(src, ovl, mir)
        call mirrorEuler(eul, geo, ovl)

        nsrc = 0; neul = 0; werr = 0._real64; vmis = 0._real64
        do b = 1, 2
            gotM(:, b) = [mir(b)%mass, mir(b)%mom, mir(b)%en]
            do k = 1, NZ; do j = 1, NY; do i = 1, geo(b)%Nx
                g = lo(b) + i - 1
                if (any(src(b)%sourceMass(:, i, j, k) /= sm(:, g, j, k)) .or. any(src(b)%sourceMom(:, i, j, k) /= so(:, g, j, k)) &
                    .or. src(b)%sourceEn(i, j, k) /= se(g, j, k)) nsrc = nsrc + 1
                if (eul(b)%density(i, j, k) /= ed(g, j, k) .or. eul(b)%np(i, j, k) /= en(g, j, k) .or. &
                    any(eul(b)%velocity(:, i, j, k) /= ev(:, g, j, k)) .or. eul(b)%temperature(i, j, k) /= et(g, j, k)) &
                    neul = neul + 1
                werr = max(werr, abs(eul(b)%density(i, j, k) - ed(g, j, k))/maxval(abs(ed)), &
                           abs(eul(b)%np(i, j, k) - en(g, j, k))/maxval(abs(en)), &
                           maxval(abs(eul(b)%velocity(:, i, j, k) - ev(:, g, j, k)))/maxval(abs(ev)), &
                           abs(eul(b)%temperature(i, j, k) - et(g, j, k))/maxval(abs(et)))
                if (b == 2 .and. g <= I1) &
                    vmis = max(vmis, abs(geo(2)%cellVol(i, j, k) - geo(1)%cellVol(g, j, k))/geo(1)%cellVol(g, j, k))
            enddo; enddo; enddo
        enddo
        call verdict('U-M (i) source: every cell of both blocks equals the grid''s, bit for bit', nsrc == 0, nsrc)
        print '(a,es9.2,a,i0,a,es9.2,a)', '  [info] U-M (i) eulerian: worst ', werr, ' of the field scale (', neul, &
            ' cells not bit-identical; the two copies of a cell differ in volume by up to ', vmis, ')'
        call verdict('U-M (i) eulerian: every cell of both blocks equals the grid''s to 1e-14 of the field scale', &
                     werr <= 1.0e-14_real64, neul)
        call verdict('U-M (i) coverage count: 2 in the overlap, 1 elsewhere; the map pairs 4 x 5 x 4 cells each way', &
                     ncov == 0 .and. all(npair == reshape([0, 80, 80, 0], [2, 2])), ncov)
        nsrc = count(abs(gotM - expM) > 1.0e-14_real64*spread(l1, 2, 2))
        call verdict('U-M (i) mirrored totals: the other block''s owner deposits in the overlap', nsrc == 0, nsrc)
        nsrc = count(abs(gotM(:, 1) + gotM(:, 2) - sumNS) > 1.0e-14_real64*l1)
        call verdict('U-M (i) Sum of the mirrored totals = Sum (n - 1) S of the owner-only fields', nsrc == 0, nsrc)
    end subroutine mirrorCase

    !> U-M (ii) + (iii): block 1 = columns 1..i1, block 2 = i0..nx of an nx x ny x nz grid, each with
    !  the grid's values; block 1's face-2 cells take the halo of block 2's cell in column i1 + 1, block
    !  2's face-1 cells that of block 1's cell in column i0 - 1 (one donor of weight 1 each).
    subroutine chimeraCase(lab, nx, ny, nz, i0, i1)
        character(*), intent(in) :: lab
        integer,      intent(in) :: nx, ny, nz, i0, i1
        real(real64), allocatable :: qg(:,:,:), vg(:,:,:), q1(:,:,:)
        type(blockField) :: bq(2), bv(2)
        type(cellLink),   allocatable :: links(:)
        type(linkDonors) :: donors
        real(real64) :: exch(2), s0(2), s1(2), l1, scale, worst, dev
        integer :: lo(2), hi(2), b, j, k, n, nl, ndim, nbit

        call randomGrid(nx, ny, nz, qg, vg)
        ndim = merge(2, 3, nz == 1)
        lo = [1, i0]; hi = [i1, nx]
        nl = 2*ny*nz
        allocate(links(nl), donors%b(nl), donors%i(nl), donors%j(nl), donors%k(nl), donors%w(nl))
        n = 0
        do k = 1, nz; do j = 1, ny
            n = n + 1
            links(n) = cellLink(b=1, f=2, i=i1, j=j, k=k, dir=1, nd=1, d0=n-1)
            donors%b(n) = 2; donors%i(n) = i1 + 1 - i0 + 1; donors%j(n) = j; donors%k(n) = k
            donors%w(n) = 1._real64
            n = n + 1
            links(n) = cellLink(b=2, f=1, i=1, j=j, k=k, dir=1, nd=1, d0=n-1)
            donors%b(n) = 1; donors%i(n) = i0 - 1; donors%j(n) = j; donors%k(n) = k
            donors%w(n) = 1._real64
        enddo; enddo
        do b = 1, 2
            bq(b)%q = qg(lo(b):hi(b), :, :)
            bv(b)%q = vg(lo(b):hi(b), :, :)
            s0(b) = sum(bq(b)%q*bv(b)%q)
        enddo
        q1 = qg
        call binomial_smooth(q1, vg, NPASS, ndim)
        call mollifyBlocks(bq, bv, links, NPASS, ndim, donors, exch)

        scale = maxval(abs(q1)); worst = 0._real64; nbit = 0; l1 = 0._real64; dev = 0._real64
        do b = 1, 2
            associate(r => q1(lo(b):hi(b), :, :))
            worst = max(worst, maxval(abs(bq(b)%q - r))/scale)
            nbit  = nbit + count(bq(b)%q /= r)
            end associate
            s1(b) = sum(bq(b)%q*bv(b)%q)
            l1 = l1 + sum(abs(bq(b)%q)*bv(b)%q)
        enddo
        do b = 1, 2
            dev = max(dev, abs(s1(b) - s0(b) - exch(b))/l1)
        enddo
        call verdict('U-M (ii) '//lab//': each block equals the uncut grid smoothed, worst '//es(worst)// &
                     ' of the field scale ('//itoa(nbit)//' cells not bit-identical, '//itoa(nl)//' 102 links)', &
                     worst <= TOL_CUT, nbit)
        call verdict('U-M (iii) '//lab//': each block''s Sum(q V) changes by its 102 intake ('//es(exch(1))//', '// &
                     es(exch(2))//'), worst '//es(dev)//' of Sum(|q| V)', dev <= TOL_SUM .and. all(exch /= 0._real64), &
                     count(abs(s1 - s0 - exch) > TOL_SUM*l1))
    end subroutine chimeraCase

    function es(x) result(r)
        real(real64), intent(in) :: x
        character(:), allocatable :: r
        character(16) :: buf
        write(buf, '(es9.2)') x
        r = trim(adjustl(buf))
    end function es

    function itoa(i) result(r)
        integer, intent(in) :: i
        character(:), allocatable :: r
        character(16) :: buf
        write(buf, '(i0)') i
        r = trim(buf)
    end function itoa

    subroutine verdict(lab, ok, nbad)
        character(*), intent(in) :: lab
        logical,      intent(in) :: ok
        integer,      intent(in) :: nbad
        if (ok) then
            print '(a)', '  [ok]   '//lab
        else
            print '(a,i0,a)', '  [FAIL] '//lab//' (', nbad, ' off)'
            nfail = nfail + 1
        endif
    end subroutine verdict

    ! ------------------------------------------------------------------ fixture

    !> Random field in [-1, 1] and volumes in [0.5, 1.5] (a fixed linear congruential stream).
    subroutine randomGrid(nx, ny, nz, Q, V)
        integer, intent(in) :: nx, ny, nz
        real(real64), allocatable, intent(out) :: Q(:,:,:), V(:,:,:)
        integer :: i, j, k
        allocate(Q(nx,ny,nz), V(nx,ny,nz))
        do k = 1, nz; do j = 1, ny; do i = 1, nx
            Q(i,j,k) = 2._real64*uniform() - 1._real64
            V(i,j,k) = 0.5_real64 + uniform()
        enddo; enddo; enddo
    end subroutine randomGrid

    real(real64) function uniform()
        seed = modulo(16807_8*seed, 2147483647_8)
        uniform = real(seed, real64)/2147483647._real64
    end function uniform

    !> Block index ranges (lo/hi(3, nblocks)) of the grid cut after the listed planes.
    subroutine cutBlocks(nx, ny, nz, ci, cj, ck, lo, hi)
        integer, intent(in) :: nx, ny, nz, ci(:), cj(:), ck(:)
        integer, allocatable, intent(out) :: lo(:,:), hi(:,:)
        integer, allocatable :: ei(:), ej(:), ek(:)
        integer :: a, c, d, b
        ei = [0, ci, nx]; ej = [0, cj, ny]; ek = [0, ck, nz]
        allocate(lo(3, (size(ei)-1)*(size(ej)-1)*(size(ek)-1)), hi(3, (size(ei)-1)*(size(ej)-1)*(size(ek)-1)))
        b = 0
        do d = 1, size(ek)-1; do c = 1, size(ej)-1; do a = 1, size(ei)-1
            b = b + 1
            lo(:,b) = [ei(a)+1, ej(c)+1, ek(d)+1]
            hi(:,b) = [ei(a+1), ej(c+1), ek(d+1)]
        enddo; enddo; enddo
    end subroutine cutBlocks

    !> One link per face cell on a cut, both ways: block b's face f cell joined to the block holding
    !  the cell across it, on the opposite face.
    subroutine linkCuts(lo, hi, links)
        integer, intent(in) :: lo(:,:), hi(:,:)
        type(cellLink), allocatable, intent(out) :: links(:)
        type(cellLink), allocatable :: tmp(:)
        integer :: n, b, pb, f, d, i, j, k, g(3), p(3)
        allocate(tmp(2*sum((hi(1,:)-lo(1,:)+1)*(hi(2,:)-lo(2,:)+1) + (hi(1,:)-lo(1,:)+1)*(hi(3,:)-lo(3,:)+1) &
                           + (hi(2,:)-lo(2,:)+1)*(hi(3,:)-lo(3,:)+1))))
        n = 0
        do b = 1, size(lo,2)
            do f = 1, 6
                d = (f+1)/2
                do k = lo(3,b), hi(3,b); do j = lo(2,b), hi(2,b); do i = lo(1,b), hi(1,b)
                    g = [i, j, k]
                    if (mod(f,2) == 1 .and. g(d) /= lo(d,b)) cycle
                    if (mod(f,2) == 0 .and. g(d) /= hi(d,b)) cycle
                    p = g
                    p(d) = g(d) + merge(-1, 1, mod(f,2) == 1)
                    do pb = 1, size(lo,2)
                        if (all(p >= lo(:,pb)) .and. all(p <= hi(:,pb))) exit
                    enddo
                    if (pb > size(lo,2)) cycle
                    n = n + 1
                    tmp(n) = cellLink(b=b, f=f, i=i-lo(1,b)+1, j=j-lo(2,b)+1, k=k-lo(3,b)+1, &
                                      pb=pb, pf=merge(f+1, f-1, mod(f,2) == 1), &
                                      pi=p(1)-lo(1,pb)+1, pj=p(2)-lo(2,pb)+1, pk=p(3)-lo(3,pb)+1, dir=d)
                enddo; enddo; enddo
            enddo
        enddo
        links = tmp(1:n)
    end subroutine linkCuts

end program test_mollify_blocks
