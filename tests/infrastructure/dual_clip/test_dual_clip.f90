program test_dual_clip
    !
    ! The ord2 dual (gas) mesh must TILE the domain: its cell volumes must sum to the domain
    ! volume, no more.
    !
    ! The dual node ring is the geo cell CENTRES plus a ghost node reflected through each
    ! boundary face centre (allocation.f90:110-125). So the dual has Nx+1 cells across Nx geo
    ! cells, and the first and last cell in every direction STRADDLE the boundary: half of
    ! each lies outside the domain, where no parcel can deposit. computeEulField divides a
    ! parcel's deposited mass by that volume, so an unclipped boundary cell reports a density
    ! diluted by exactly the outside fraction. Measured on the JPL nozzle, the boundary geo
    ! cells carried only 0.74 of the condensed mass flux -- at the inlet, where that flux is
    ! imposed and must be reproduced exactly.
    !
    ! precomputeDualMetric clips each boundary cell by the fraction of its extent that lies
    ! inside, MEASURED from the mesh via insideFrac. On a uniform mesh with a reflected ghost
    ! that fraction is exactly 1/2, so the analytic oracle is exact and needs no tolerance
    ! beyond accumulated round-off:
    !
    !   sum(dual cellVol) == Lx*Ly*Lz            exact tiling
    !   face cell    = h^3/2      one clipped direction
    !   edge cell    = h^3/4      two
    !   corner cell  = h^3/8      three
    !   interior     = h^3        none
    !
    ! Unclipped, the sum overshoots by (Nx+1)(Ny+1)(Nz+1)/(Nx*Ny*Nz). Measured against the
    ! 2D-only clip: D1 2.500x (60/24), D2 3.375x (27/8), and the corner cell is off by exactly
    ! 8. That factor IS the historical bug -- and it is worst on a COARSE mesh, so a 3D case
    ! small enough to iterate on is where it would have bitten hardest.
    !
    ! D1/D2 run the 3D branch, which had no clip at all until this test was written and which
    ! no e2e case GATES (box fixtures are Nk=5 but ord1, no dual; vie-plait runs the 3D dual at
    ! gas-order 2 with nothing reading its field). D3/D4 run the 2D branch, which shipped
    ! clipped -- they are the regression half, and they also pin that adding k did not disturb i/j.
    !
    ! The fixture is a uniform Cartesian box built through the PRODUCTION geometry path
    ! (compute_geometry + precomputeMetric), with the dual assembled exactly as
    ! allocation.f90 does it, so the ghost placement under test is the shipped rule and not a
    ! restatement of it.
    !
    ! U3 and U3n join two blocks by a 101 face and run the setup pass (assembleInterfaceDual): the
    ! interface dual cell is not clipped and the vertex classes pair coincident vertices. U3p joins a
    ! box's faces 3/4 by 201 and runs the same pass: each seam dual cell is the cell of the box
    ! unrolled over two periods, its vertices paired across the period.
    !
    use, intrinsic :: iso_fortran_env, only: real64
    use IGLOO_data_block, only: obj_block, obj_flowblock, obj_bc_cell, vertexClass, buildInterfaceLinks   ! PRODUCTION types
    use IGLOO_allocation, only: fill_dual_nodes, computeSkewFlag, assembleInterfaceDual
    use IGLOO_Lib_Mollify, only: cellLink, linkDonors
    use IGLOO_variables,  only: ord2, mesh2D
    implicit none

    real(real64), parameter :: H = 0.25_real64      ! uniform spacing, all directions
    integer :: nfail

    nfail = 0

    call runCase('D1 3D  4x3x2', 4, 3, 2, .false.)
    call runCase('D2 3D  2x2x2', 2, 2, 2, .false.)
    call runCase('D3 2D  4x3',   4, 3, 1, .true. )
    call runCase('D4 2D  2x2',   2, 2, 1, .true. )

    call s5_u3_two_block('U3 3D       ', .false.)
    call s5_u3_two_block('U3 one layer', .true.)
    call s5_u3_nonaligned
    call p201_u3_periodic('U3p 3D       ', .false.)
    call p201_u3_periodic('U3p one layer', .true.)

    if (nfail == 0) then
        print '(a)', ''
        print '(a)', '[PASS] test_dual_clip: the dual tiles the domain in 2D and 3D'
    else
        print '(a)', ''
        print '(a,i0,a)', '[FAIL] test_dual_clip: ', nfail, ' assertion(s) failed'
        error stop 1
    endif

contains

    !> Build a uniform Nx x Ny x Nz box, its dual, and check the tiling invariant.
    subroutine runCase(label, nx, ny, nz, is2D)
        character(*), intent(in) :: label
        integer,      intent(in) :: nx, ny, nz
        logical,      intent(in) :: is2D
        type(obj_block) :: geo, dual
        real(real64) :: domVol, dualSum, expect
        integer :: i, j, k, nxe, nye, nze, nclip

        call buildGeo(geo, nx, ny, nz, is2D)
        call buildDual(dual, geo, is2D)
        call dual%precomputeDualMetric(is2D, geo)

        nxe = nx + 1; nye = ny + 1; nze = merge(1, nz + 1, is2D)

        !> (a) exact tiling -- the headline invariant
        domVol  = sum(geo%cellVol)
        dualSum = sum(dual%cellVol)
        call check(label//' : dual sum == domain volume', dualSum, domVol)

        !> (b) per-cell factors: h^3 halved once per clipped direction. Redundant with (a)
        !     only if the errors cannot cancel -- they can (a cell too big beside one too
        !     small), so this is not decoration.
        do k = 1, nze; do j = 1, nye; do i = 1, nxe
            nclip = 0
            if (i == 1 .or. i == nxe) nclip = nclip + 1
            if (j == 1 .or. j == nye) nclip = nclip + 1
            if (.not. is2D) then
                if (k == 1 .or. k == nze) nclip = nclip + 1
            endif
            expect = cellFull(is2D) / 2.0_real64**nclip
            if (abs(dual%cellVol(i,j,k) - expect) > 1.0e-12_real64*expect) then
                print '(a,3(i0,1x),a,es13.6,a,es13.6)', &
                    '  [FAIL] '//label//' cell ', i, j, k, ' got ', dual%cellVol(i,j,k), &
                    ' expected ', expect
                nfail = nfail + 1
                return                              ! one report per case is enough
            endif
        enddo; enddo; enddo
        print '(a,a,es13.6)', '  [ok] ', label//' tiles exactly; sum = ', dualSum
    end subroutine runCase

    !> Volume of one interior dual cell: h^3 in 3D, h^2 * slab in 2D (slab thickness = h,
    !  since buildGeo extrudes a single layer of thickness H).
    pure function cellFull(is2D) result(v)
        logical, intent(in) :: is2D
        real(real64) :: v
        v = H**3
        if (is2D) v = H**3                          ! same: the 2D slab is one H-thick layer
    end function cellFull

    !> Uniform Cartesian geo block through the production geometry path.
    subroutine buildGeo(blk, nx, ny, nz, is2D)
        type(obj_block), intent(out) :: blk
        integer,         intent(in)  :: nx, ny, nz
        logical,         intent(in)  :: is2D
        integer :: i, j, k
        blk%Nx = nx; blk%Ny = ny; blk%Nz = nz
        allocate(blk%node(3, 0:nx, 0:ny, 0:nz))
        do k = 0, nz; do j = 0, ny; do i = 0, nx
            blk%node(:, i, j, k) = [i*H, j*H, k*H]
        enddo; enddo; enddo
        allocate(blk%center(3, 1:nx, 1:ny, 1:nz))
        call blk%compute_geometry
        call blk%precomputeMetric
        if (is2D) continue                          ! mesh2D flattening is not needed: the
                                                    ! 2D branch reads only k=1, and a single
                                                    ! H-thick layer is what it assumes
    end subroutine buildGeo

    !> Dual block assembled exactly as allocation.f90:102-125 does it: interior nodes are geo
    !  cell centres, face ghosts are reflections through the boundary face centre.
    subroutine buildDual(dl, geo, is2D)
        type(obj_block), intent(out) :: dl
        type(obj_block), intent(in)  :: geo
        logical,         intent(in)  :: is2D
        integer :: i, j, k, kg
        dl%Nx = geo%Nx; dl%Ny = geo%Ny; dl%Nz = geo%Nz
        if (is2D) then
            allocate(dl%node(3, 0:dl%Nx+1, 0:dl%Ny+1, 1:1))
        else
            allocate(dl%node(3, 0:dl%Nx+1, 0:dl%Ny+1, 0:dl%Nz+1))
        endif
        kg = merge(1, dl%Nz, is2D)
        do k = 1, kg; do j = 1, dl%Ny; do i = 1, dl%Nx
            dl%node(:, i, j, k) = geo%center(:, i, j, k)
        enddo; enddo; enddo
        do k = 1, kg; do j = 1, dl%Ny
            dl%node(:, 0,       j, k) = 2*geo%face(1)%cell(j,k)%center - geo%center(:,1,      j,k)
            dl%node(:, dl%Nx+1, j, k) = 2*geo%face(2)%cell(j,k)%center - geo%center(:,dl%Nx,  j,k)
        enddo; enddo
        do k = 1, kg; do i = 1, dl%Nx
            dl%node(:, i, 0,       k) = 2*geo%face(3)%cell(i,k)%center - geo%center(:,i,1,      k)
            dl%node(:, i, dl%Ny+1, k) = 2*geo%face(4)%cell(i,k)%center - geo%center(:,i,dl%Ny,  k)
        enddo; enddo
        if (.not. is2D) then
            do j = 1, dl%Ny; do i = 1, dl%Nx
                dl%node(:, i, j, 0)       = 2*geo%face(5)%cell(i,j)%center - geo%center(:,i,j,1)
                dl%node(:, i, j, dl%Nz+1) = 2*geo%face(6)%cell(i,j)%center - geo%center(:,i,j,dl%Nz)
            enddo; enddo
        endif
        call fillEdges(dl, is2D)
    end subroutine buildDual

    !> Edge/corner ghosts by the same cascading extrapolation as allocation.f90:127-155.
    !  computeVolume (3D) and getVertices (2D) read them for the boundary dual cells, so they
    !  must exist. ⚠ The i/j corner block runs in BOTH modes -- see `fill_dual_nodes` in
    !  allocation.f90, whose header spells the three groups out. Omitting it in 2D leaves
    !  node(0,0) undefined and the boundary quad area is then garbage; this test caught
    !  exactly that in its own fixture, which is why that header now exists.
    subroutine fillEdges(dl, is2D)
        type(obj_block), intent(inout) :: dl
        logical,         intent(in)    :: is2D
        integer :: i, j, k, nx, ny, nz
        nx = dl%Nx; ny = dl%Ny; nz = merge(1, dl%Nz, is2D)
        do k = 1, nz
            dl%node(:,0,   0,   k) = 2*dl%node(:,0,   1, k) - dl%node(:,0,   2,   k)
            dl%node(:,nx+1,0,   k) = 2*dl%node(:,nx+1,1, k) - dl%node(:,nx+1,2,   k)
            dl%node(:,0,   ny+1,k) = 2*dl%node(:,0,   ny,k) - dl%node(:,0,   ny-1,k)
            dl%node(:,nx+1,ny+1,k) = 2*dl%node(:,nx+1,ny,k) - dl%node(:,nx+1,ny-1,k)
        enddo
        if (is2D) return
        nz = dl%Nz
        do j = 1, ny
            dl%node(:,0,   j,0)    = 2*dl%node(:,0,   j,1)  - dl%node(:,0,   j,2)
            dl%node(:,nx+1,j,0)    = 2*dl%node(:,nx+1,j,1)  - dl%node(:,nx+1,j,2)
            dl%node(:,0,   j,nz+1) = 2*dl%node(:,0,   j,nz) - dl%node(:,0,   j,nz-1)
            dl%node(:,nx+1,j,nz+1) = 2*dl%node(:,nx+1,j,nz) - dl%node(:,nx+1,j,nz-1)
        enddo
        do i = 1, nx
            dl%node(:,i,0,   0)    = 2*dl%node(:,i,0,   1)  - dl%node(:,i,0,   2)
            dl%node(:,i,ny+1,0)    = 2*dl%node(:,i,ny+1,1)  - dl%node(:,i,ny+1,2)
            dl%node(:,i,0,   nz+1) = 2*dl%node(:,i,0,   nz) - dl%node(:,i,0,   nz-1)
            dl%node(:,i,ny+1,nz+1) = 2*dl%node(:,i,ny+1,nz) - dl%node(:,i,ny+1,nz-1)
        enddo
        dl%node(:,0,   0,   0)    = 2*dl%node(:,0,   0,   1)  - dl%node(:,0,   0,   2)
        dl%node(:,nx+1,0,   0)    = 2*dl%node(:,nx+1,0,   1)  - dl%node(:,nx+1,0,   2)
        dl%node(:,0,   ny+1,0)    = 2*dl%node(:,0,   ny+1,1)  - dl%node(:,0,   ny+1,2)
        dl%node(:,nx+1,ny+1,0)    = 2*dl%node(:,nx+1,ny+1,1)  - dl%node(:,nx+1,ny+1,2)
        dl%node(:,0,   0,   nz+1) = 2*dl%node(:,0,   0,   nz) - dl%node(:,0,   0,   nz-1)
        dl%node(:,nx+1,0,   nz+1) = 2*dl%node(:,nx+1,0,   nz) - dl%node(:,nx+1,0,   nz-1)
        dl%node(:,0,   ny+1,nz+1) = 2*dl%node(:,0,   ny+1,nz) - dl%node(:,0,   ny+1,nz-1)
        dl%node(:,nx+1,ny+1,nz+1) = 2*dl%node(:,nx+1,ny+1,nz) - dl%node(:,nx+1,ny+1,nz-1)
    end subroutine fillEdges

    subroutine check(what, got, expect)
        character(*), intent(in) :: what
        real(real64), intent(in) :: got, expect
        if (abs(got - expect) > 1.0e-12_real64*abs(expect)) then
            print '(a,a,es16.9,a,es16.9,a,es9.2)', '  [FAIL] ', what//' : got ', got, &
                ' expected ', expect, ' rel ', abs(got-expect)/abs(expect)
            nfail = nfail + 1
        endif
    end subroutine check

    !> U3 -- a box pair joined by 101 along x (h_B = 1.5 h_A), the setup pass, against the one-block
    !  box on the same node lines: each interface dual cell is the one-block cell, the dual tiles once.
    subroutine s5_u3_two_block(label, planar)
        character(*), intent(in) :: label
        logical,      intent(in) :: planar
        integer, parameter :: NXA = 3, NXB = 2, NY = 3
        type(obj_block)     :: geo(2), one
        type(obj_flowblock) :: gas(2), gone
        type(cellLink), allocatable    :: links(:)
        type(linkDonors)               :: donors
        type(vertexClass), allocatable :: classes(:)
        real(real64) :: xa(0:NXA), xb(0:NXB), xo(0:NXA+NXB), eA, eB, sI, vo
        integer :: i, j, k, nz, nk, nmoved
        ord2 = .true.; mesh2D = planar
        nz = merge(1, 2, planar); nk = merge(1, nz + 1, planar)
        xa = [(i*H, i = 0, NXA)]
        xb = [(NXA*H + i*1.5_real64*H, i = 0, NXB)]
        xo(0:NXA) = xa; xo(NXA:NXA+NXB) = xb
        call boxNodes(geo(1), xa, NY, nz); call finishBlock(geo(1), gas(1), planar)
        call boxNodes(geo(2), xb, NY, nz); call finishBlock(geo(2), gas(2), planar)
        call boxNodes(one, xo, NY, nz);    call finishBlock(one, gone, planar)
        do k = 1, nz; do j = 1, NY
            call link(geo(1)%face(2)%cell(j,k), [2, 1, j, k], 1, [1, 0, 0, 1])
            call link(geo(2)%face(1)%cell(j,k), [1, NXA, j, k], 2, [1, 0, 0, 1])
        enddo; enddo
        call buildInterfaceLinks(geo, links, donors)
        call assembleInterfaceDual(geo, gas, links, classes, nmoved)
        eA = 0._real64; eB = 0._real64; sI = 0._real64
        do k = 1, nk; do j = 1, NY + 1
            vo = gone%cellVol(NXA+1, j, k)
            eA = max(eA, abs(gas(1)%cellVol(NXA+1, j, k) - vo)/vo)
            eB = max(eB, abs(gas(2)%cellVol(1, j, k) - vo)/vo)
            sI = sI + gas(1)%cellVol(NXA+1, j, k)
        enddo; enddo
        call report(label//' : A interface dual cells vs the one-block cells, max rel', eA, 1.0e-14_real64)
        call report(label//' : B interface dual cells vs the one-block cells, max rel', eB, 1.0e-14_real64)
        call report(label//' : sum_A + sum_B - sum_interface vs sum one-block, rel', &
                    abs(sum(gas(1)%cellVol) + sum(gas(2)%cellVol) - sI - sum(gone%cellVol))/sum(gone%cellVol), &
                    1.0e-14_real64)
    end subroutine s5_u3_two_block

    !> U3n -- block A's face 4 (j-max) joined to block B's face 1 (i-min), B's i along +y and its j
    !  along -x: every class holds two coincident vertices and the one-block box's assembled weight.
    subroutine s5_u3_nonaligned
        integer,      parameter :: NXA = 4, NYA = 2, NXB = 3, NZ = 2
        real(real64), parameter :: HX = 0.25_real64, HYA = 0.3_real64, HYB = 0.45_real64, HZ = 0.25_real64
        type(obj_block)     :: geo(2), one
        type(obj_flowblock) :: gas(2), gone
        type(cellLink), allocatable    :: links(:)
        type(linkDonors)               :: donors
        type(vertexClass), allocatable :: classes(:)
        real(real64) :: dmax, ew, x0(3), lref
        integer :: c, n, i, j, k, nmoved, nbad
        ord2 = .true.; mesh2D = .false.
        geo(1)%Nx = NXA; geo(1)%Ny = NYA; geo(1)%Nz = NZ
        allocate(geo(1)%node(3, 0:NXA, 0:NYA, 0:NZ))
        do k = 0, NZ; do j = 0, NYA; do i = 0, NXA
            geo(1)%node(:, i, j, k) = [i*HX, j*HYA, k*HZ]
        enddo; enddo; enddo
        geo(2)%Nx = NXB; geo(2)%Ny = NXA; geo(2)%Nz = NZ
        allocate(geo(2)%node(3, 0:NXB, 0:NXA, 0:NZ))
        do k = 0, NZ; do j = 0, NXA; do i = 0, NXB
            geo(2)%node(:, i, j, k) = [NXA*HX - j*HX, NYA*HYA + i*HYB, k*HZ]
        enddo; enddo; enddo
        one%Nx = NXA; one%Ny = NYA + NXB; one%Nz = NZ
        allocate(one%node(3, 0:NXA, 0:NYA+NXB, 0:NZ))
        do k = 0, NZ; do j = 0, NYA + NXB; do i = 0, NXA
            one%node(:, i, j, k) = [i*HX, merge(j*HYA, NYA*HYA + (j - NYA)*HYB, j <= NYA), k*HZ]
        enddo; enddo; enddo
        call finishBlock(geo(1), gas(1), .false.)
        call finishBlock(geo(2), gas(2), .false.)
        call finishBlock(one, gone, .false.)
        do k = 1, NZ; do i = 1, NXA
            call link(geo(1)%face(4)%cell(i,k), [2, 1, NXA + 1 - i, k], 1, [-1, 0, 0, 1])
            call link(geo(2)%face(1)%cell(NXA + 1 - i, k), [1, i, NYA, k], 4, [-1, 0, 0, 1])
        enddo; enddo
        call buildInterfaceLinks(geo, links, donors)
        call assembleInterfaceDual(geo, gas, links, classes, nmoved)
        lref = norm2(one%node(:, NXA, NYA+NXB, NZ) - one%node(:, 0, 0, 0))
        nbad = 0; dmax = 0._real64; ew = 0._real64
        do c = 1, size(classes)
            if (size(classes(c)%m, 2) /= 2 .or. count(classes(c)%m(1,:) == 1) /= 1) nbad = nbad + 1
            associate(m1 => classes(c)%m(:, 1))
            x0 = geo(m1(1))%node(:, m1(2) - 1, m1(3) - 1, m1(4) - 1)
            end associate
            do n = 1, size(classes(c)%m, 2)
                associate(mb => classes(c)%m(:, n))
                dmax = max(dmax, norm2(geo(mb(1))%node(:, mb(2) - 1, mb(3) - 1, mb(4) - 1) - x0)/lref)
                ew = max(ew, abs(geo(mb(1))%dualW(mb(2), mb(3), mb(4)) - classes(c)%w)/classes(c)%w)
                if (mb(1) == 1) ew = max(ew, abs(classes(c)%w - one%dualW(mb(2), mb(3), mb(4)))/one%dualW(mb(2), mb(3), mb(4)))
                end associate
            enddo
        enddo
        call report('U3n vertex classes - (NXA+1)(NZ+1)', real(size(classes) - (NXA + 1)*(NZ + 1), real64), 0._real64)
        call report('U3n classes not one vertex of A and one of B', real(nbad, real64), 0._real64)
        call report('U3n class members apart, max distance / diagonal', dmax, 1.0e-14_real64)
        call report('U3n class weights vs the one-block assembled weights, max rel', ew, 1.0e-14_real64)
    end subroutine s5_u3_nonaligned

    !> U3p -- one box whose faces 3/4 are a 201 pair (rows of 0.2, 0.3, 0.3 along y), the setup pass, against the
    !  box unrolled over two periods: each seam dual cell is the unrolled box's interior cell at y = Ly.
    subroutine p201_u3_periodic(label, planar)
        character(*), intent(in) :: label
        logical,      intent(in) :: planar
        integer,      parameter :: NX = 3, NY = 3
        real(real64), parameter :: LY = 0.8_real64
        type(obj_block)     :: geo(1), ref
        type(obj_flowblock) :: gas(1), gref
        type(cellLink), allocatable    :: links(:)
        type(linkDonors)               :: donors
        type(vertexClass), allocatable :: classes(:)
        real(real64) :: xn(0:NX), yn(0:NY), yr(0:2*NY), ty(3), eb, ec, ed, vr, vol, dom
        integer :: i, k, c, n, nz, nk, nd, nc, nbad, nmoved
        ord2 = .true.; mesh2D = planar
        nz = merge(1, 2, planar); nk = merge(1, nz + 1, planar); nd = merge(2, 3, planar)
        xn = [(i*H, i = 0, NX)]
        yn = [0._real64, 0.2_real64, 0.5_real64, LY]
        yr(0:NY) = yn; yr(NY+1:2*NY) = LY + yn(1:NY)
        call boxNodes(geo(1), xn, NY, nz, yn); call finishBlock(geo(1), gas(1), planar)
        call boxNodes(ref, xn, 2*NY, nz, yr);  call finishBlock(ref, gref, planar)
        do k = 1, nz; do i = 1, NX
            call link(geo(1)%face(3)%cell(i,k), [1, i, NY, k], 4, [1, 0, 0, 1], 201)
            call link(geo(1)%face(4)%cell(i,k), [1, i, 1, k], 3, [1, 0, 0, 1], 201)
        enddo; enddo
        call buildInterfaceLinks(geo, links, donors)
        call assembleInterfaceDual(geo, gas, links, classes, nmoved)
        !> (a) one class per seam vertex: dual rows 1 and NY+1 of block 1, the same i and k
        nc = 0
        if (allocated(classes)) nc = size(classes)
        nbad = 0
        do c = 1, nc
            associate(m => classes(c)%m)
            if (size(m, 2) /= 2) then
                nbad = nbad + 1
            else if (any(m(1,:) /= 1) .or. m(2,1) /= m(2,2) .or. m(4,1) /= m(4,2) .or. &
                     minval(m(3,:)) /= 1 .or. maxval(m(3,:)) /= NY + 1) then
                nbad = nbad + 1
            endif
            end associate
        enddo
        call report(label//' : (a) vertex classes - (NX+1)(NZ+1)', real(nc - (NX + 1)*nk, real64), 0._real64)
        call report(label//' : (a) classes not rows 1 and NY+1 of one vertex column', real(nbad, real64), 0._real64)
        !> (b) the face ghosts at the partner centres brought back across the period
        ty = [0._real64, LY, 0._real64]
        eb = 0._real64
        do k = 1, nz; do i = 1, NX
            eb = max(eb, norm2(gas(1)%node(1:nd, i, 0, k) - geo(1)%center(1:nd, i, NY, k) + ty(1:nd))/H, &
                         norm2(gas(1)%node(1:nd, i, NY + 1, k) - geo(1)%center(1:nd, i, 1, k) - ty(1:nd))/H)
        enddo; enddo
        call report(label//' : (b) face ghosts vs the partner centres -+ Ly, max / H', eb, 1.0e-14_real64)
        !> (c) the seam dual cells off the x/z boundaries unclipped: the unrolled box's cells at y = Ly
        ec = 0._real64
        do k = merge(1, 2, planar), merge(1, nz, planar); do i = 2, NX
            vr = gref%cellVol(i, NY + 1, k)
            ec = max(ec, abs(gas(1)%cellVol(i, 1, k) - vr)/vr, abs(gas(1)%cellVol(i, NY + 1, k) - vr)/vr)
        enddo; enddo
        call report(label//' : (c) seam dual cells vs the unrolled cells at y = Ly, max rel', ec, 1.0e-14_real64)
        !> (d) every class's weight and its members' dualW: the unrolled box's dualW at y = Ly
        ed = 0._real64
        do c = 1, nc
            vr = ref%dualW(classes(c)%m(2, 1), NY + 1, classes(c)%m(4, 1))
            ed = max(ed, abs(classes(c)%w - vr)/vr)
            do n = 1, size(classes(c)%m, 2)
                associate(mb => classes(c)%m(:, n))
                ed = max(ed, abs(geo(1)%dualW(mb(2), mb(3), mb(4)) - vr)/vr)
                end associate
            enddo
        enddo
        call report(label//' : (d) class weights vs the unrolled weights at y = Ly, max rel', ed, 1.0e-14_real64)
        !> (e) the dual tiles one period once its seam row NY+1 (a second copy of row 1's cells) is left out
        vol = sum(gas(1)%cellVol)
        do k = 1, nk; do i = 1, NX + 1
            vol = vol - gas(1)%cellVol(i, NY + 1, k)
        enddo; enddo
        dom = NX*H*LY*merge(H, nz*H, planar)
        call report(label//' : (e) sum - seam row NY+1 vs Lx*Ly*Lz, rel', (vol - dom)/dom, 1.0e-14_real64)
    end subroutine p201_u3_periodic

    !> Box nodes on the x node line xn, spacing H in y (or the y node line yn) and z.
    subroutine boxNodes(blk, xn, ny, nz, yn)
        type(obj_block), intent(inout) :: blk
        real(real64),    intent(in)    :: xn(0:)
        integer,         intent(in)    :: ny, nz
        real(real64),    intent(in), optional :: yn(0:)
        integer :: i, j, k
        real(real64) :: y
        blk%Nx = size(xn) - 1; blk%Ny = ny; blk%Nz = nz
        allocate(blk%node(3, 0:blk%Nx, 0:ny, 0:nz))
        do k = 0, nz; do j = 0, ny; do i = 0, blk%Nx
            y = j*H
            if (present(yn)) y = yn(j)
            blk%node(:, i, j, k) = [xn(i), y, k*H]
        enddo; enddo; enddo
    end subroutine boxNodes

    !> Geometry, dual weights and gas dual of a block whose nodes are set, as allocate_blocks builds them.
    subroutine finishBlock(blk, gas, planar)
        type(obj_block),     intent(inout) :: blk
        type(obj_flowblock), intent(inout) :: gas
        logical,             intent(in)    :: planar
        allocate(blk%center(3, 1:blk%Nx, 1:blk%Ny, 1:blk%Nz))
        call blk%compute_geometry
        call blk%precomputeDualWeights
        call blk%precomputeMetric
        gas%Nx = blk%Nx; gas%Ny = blk%Ny; gas%Nz = blk%Nz
        if (planar) then
            allocate(gas%node(3, 0:blk%Nx+1, 0:blk%Ny+1, 1:1))
        else
            allocate(gas%node(3, 0:blk%Nx+1, 0:blk%Ny+1, 0:blk%Nz+1))
        endif
        call fill_dual_nodes(gas, blk, planar)
        call computeSkewFlag(gas, planar)
        call gas%precomputeDualMetric(planar, blk)
        call blk%computeBBox
        call gas%computeBBox
    end subroutine finishBlock

    !> A 101 record (or one of the given code): partner cell [block, i, j, k], partner face, orientation.
    subroutine link(cell, partner, pface, orient, code)
        type(obj_bc_cell), intent(inout) :: cell
        integer,           intent(in)    :: partner(4), pface, orient(4)
        integer,           intent(in), optional :: code
        cell%bcdef = 101
        if (present(code)) cell%bcdef = code
        cell%connection = partner
        cell%connectionFace = pface
        cell%orient = orient
    end subroutine link

    subroutine report(what, got, tol)
        character(*), intent(in) :: what
        real(real64), intent(in) :: got, tol
        if (abs(got) <= tol) then
            print '(a,a,es10.3)', '  [ok] ', what//' = ', got
        else
            print '(a,a,es10.3,a,es9.2)', '  [FAIL] ', what//' = ', got, ' expected |.| <= ', tol
            nfail = nfail + 1
        endif
    end subroutine report

end program test_dual_clip
