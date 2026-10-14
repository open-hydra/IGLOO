program test_mollify_blocks
    !
    ! The field mollifier over several blocks joined by 101 links (mollifyBlocks).
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
    !
    use, intrinsic :: iso_fortran_env, only: real64
    use IGLOO_Lib_Mollify, only: binomial_smooth, mollifyBlocks, cellLink, blockField
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

    if (nfail == 0) then
        print '(a)', ''
        print '(a)', '[PASS] test_mollify_blocks: blocks joined by 101 links smooth as one block'
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
