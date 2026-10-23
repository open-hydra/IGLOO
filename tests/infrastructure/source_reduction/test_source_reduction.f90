program test_source_reduction
    !
    ! The ord2 source reduction must CONSERVE. Under gas-order = 2 the source accumulators
    ! live on the dual mesh (1..N+1) and finalizeSRC reduces them to geo cells. The shipped
    ! rule was a volume AVERAGE of EXTENSIVE rates,
    !
    !     S_geo(c) = sum_oct subVol(oct,c) * S_dual(c+oct) / sum_oct subVol(oct,c),
    !
    ! which hands each dual cell on with total weight sum_{c touching d} subVol/V_geo(c).
    ! On a uniform Cartesian mesh that is 1 for an interior dual cell but 1/2 on a face row,
    ! 1/4 on an edge and 1/8 on a corner -- so every boundary row leaked most of its
    ! deposit, and the loss scales with how much of the source sits near a boundary. Nothing
    ! printed a conservation total and every telescoping gate ran at gas-order = 1.
    !
    ! The fix divides by the ASSEMBLED dual weight instead:
    !
    !     dualW(d) = sum over every (cell, octant) pair that touches d of subVol,
    !     S_geo(c) = sum_oct subVol(oct,c)/dualW(c+oct) * S_dual(c+oct),
    !
    ! whose octant weights form a partition of unity over the geo cells sharing each dual
    ! cell. Then sum_c S_geo = sum_d S_dual EXACTLY, on any mesh -- the bracket is dualW(d)
    ! by construction. The clipped physical dual volume is deliberately NOT used: enforcing
    ! conservation through an approximate geometry is what bug O15 was.
    !
    ! S1-S6 are uniform Cartesian, where the pre-fix weights are the exact rationals above.
    ! S7 is curvilinear, where the pre-fix error is instead the octant inequality and the
    ! post-fix identity must STILL be exact -- that is the case that separates "conserves on
    ! a nice mesh" from "conserves by construction".
    !
    ! The partition-of-unity invariant is asserted as sum_d dualW == sum_c sum_oct subVol.
    ! It is NOT asserted against sum_c cellVol on S7: cellVol is hexVolume of the parent hex
    ! and hexVolume is a 5-tetrahedron split with abs() per tet, which is NOT additive under
    ! octant subdivision on a hex with non-planar faces. Measured on S7's node formula the
    ! two differ by 5.1e-3 relative (1.1e-2 worst per cell) -- an assertion against cellVol
    ! there would be RED on the CORRECT code.
    !
    ! U2 joins two boxes by a 101 face: a deposit in one block's interface dual cell is shared
    ! between the blocks in proportion to their octant volumes, for the source and the euler field.
    !
    use, intrinsic :: iso_fortran_env, only: real64
    use IGLOO_data_block, only: obj_block, obj_sourceblock, obj_flowblock, obj_eulerblock, &   ! PRODUCTION types
                                vertexClass, buildInterfaceLinks, combineInterfaceDual
    use IGLOO_allocation, only: fill_dual_nodes, computeSkewFlag, assembleInterfaceDual
    use IGLOO_Lib_Mollify, only: cellLink, linkDonors
    use IGLOO_variables,  only: ord2, mesh2D, nm, mollifyPasses
    implicit none

    real(real64), parameter :: H   = 0.25_real64     ! uniform spacing
    real(real64), parameter :: TOL = 1.0e-12_real64  ! post-fix residual measured at <= 2.2e-16
    real(real64), parameter :: TOL14 = 1.0e-14_real64
    integer :: nfail

    nfail = 0
    nm            = 1        ! finalizeSRC sizes accMass(nm)
    mollifyPasses = 0        ! the reduction alone is under test
    ord2          = .true.   ! MUST precede compute_geometry: subVol is allocated only under ord2

    !                label            nx ny nz  2D     curved  deposit (0 = fill every dual cell)
    call runCase('S1 3D 4x3x2 uniform', 4, 3, 2, .false., .false., 0, 0, 0)
    call runCase('S2 3D 2x2x2 uniform', 2, 2, 2, .false., .false., 0, 0, 0)
    call runCase('S3 2D 4x3   uniform', 4, 3, 1, .true. , .false., 0, 0, 0)
    call runCase('S4 2D 2x2   uniform', 2, 2, 1, .true. , .false., 0, 0, 0)

    call runCase('S5 3D interior (2,2,2)', 4, 3, 2, .false., .false., 2, 2, 2)
    call runCase('S5 3D face     (2,2,1)', 4, 3, 2, .false., .false., 2, 2, 1)
    call runCase('S5 3D edge     (2,1,1)', 4, 3, 2, .false., .false., 2, 1, 1)
    call runCase('S5 3D corner   (1,1,1)', 4, 3, 2, .false., .false., 1, 1, 1)

    call runCase('S6 2D interior (2,2,1)', 4, 3, 1, .true. , .false., 2, 2, 1)
    call runCase('S6 2D face     (2,1,1)', 4, 3, 1, .true. , .false., 2, 1, 1)
    call runCase('S6 2D corner   (1,1,1)', 4, 3, 1, .true. , .false., 1, 1, 1)

    call runCase('S7 curvilinear all',      4, 3, 2, .false., .true. , 0, 0, 0)
    call runCase('S7 curvilinear interior', 4, 3, 2, .false., .true. , 2, 2, 2)
    call runCase('S7 curvilinear face',     4, 3, 2, .false., .true. , 2, 2, 1)

    call s5_u2_two_block('U2 101 pair h_B = h_A    ', H)
    call s5_u2_two_block('U2 101 pair h_B = 1.5 h_A', 1.5_real64*H)

    if (nfail == 0) then
        print '(a)', ''
        print '(a)', '[PASS] test_source_reduction: the ord2 reduction conserves on every mesh'
    else
        print '(a)', ''
        print '(a,i0,a)', '[FAIL] test_source_reduction: ', nfail, ' assertion(s) failed'
        error stop 1
    endif

contains

    !> One case: build the geo block, deposit on the dual, reduce, check the totals.
    subroutine runCase(label, nx, ny, nz, is2D, curved, di, dj, dk)
        character(*), intent(in) :: label
        integer,      intent(in) :: nx, ny, nz, di, dj, dk
        logical,      intent(in) :: is2D, curved
        type(obj_block)       :: geo
        type(obj_sourceblock) :: src
        real(real64) :: dualTot, geoTot, ratio, wSum, vSum, cSum
        integer      :: nk, oct, i, j, k, da, db, dc

        mesh2D = is2D                     ! module variable, read by both the weights and finalize
        call buildGeo(geo, nx, ny, nz, curved)
        call geo%precomputeDualWeights

        nk = nz + 1
        if (is2D) nk = 1
        call src%allocate(1, nx+1, ny+1, nk)
        src%Nx = nx; src%Ny = ny; src%Nz = nz
        src%sourceMass = 0._real64; src%sourceMom = 0._real64; src%sourceEn = 0._real64
        if (di == 0) then
            src%sourceMass(1,:,:,:) = 1._real64      ! every dual cell
            src%sourceMom (1,:,:,:) = 2._real64
            src%sourceEn    (:,:,:) = 3._real64
        else
            src%sourceMass(1,di,dj,dk) = 1._real64   ! a single dual cell
            src%sourceMom (1,di,dj,dk) = 2._real64
            src%sourceEn    (di,dj,dk) = 3._real64
        endif
        dualTot = sum(src%sourceMass(1,:,:,:))

        !> partition of unity: the assembled weights are exactly the octant volumes
        wSum = sum(geo%dualW)
        vSum = 0._real64
        cSum = 0._real64
        do k = 1, nz; do j = 1, ny; do i = 1, nx
            do dc = 0, 1; do db = 0, 1; do da = 0, 1
                oct  = 1 + da + 2*db + 4*dc
                vSum = vSum + geo%subVol(oct, i, j, k)
            enddo; enddo; enddo
            cSum = cSum + geo%cellVol(i, j, k)
        enddo; enddo; enddo
        call assertRel(label//' : sum(dualW) == sum(subVol)', wSum, vSum, TOL)
        !> cellVol agrees only where the faces are planar: hexVolume's 5-tet split with
        !> abs() per tet is not additive under octant subdivision on a curved hex.
        if (.not. curved) call assertRel(label//' : sum(dualW) == sum(cellVol)', wSum, cSum, TOL)

        call src%finalize(geo)

        geoTot = sum(src%sourceMass(1,:,:,:))
        ratio  = geoTot / dualTot
        call assertRel(label//' : mass  sum_geo/sum_dual == 1', ratio, 1._real64, TOL)
        call assertRel(label//' : mom   sum_geo/sum_dual == 1', &
                       sum(src%sourceMom(1,:,:,:))/(2._real64*dualTot), 1._real64, TOL)
        call assertRel(label//' : en    sum_geo/sum_dual == 1', &
                       sum(src%sourceEn(:,:,:))/(3._real64*dualTot), 1._real64, TOL)

        call geo%freeBlock
    end subroutine runCase

    !> Uniform or smoothly warped Cartesian box through the PRODUCTION geometry path.
    subroutine buildGeo(blk, nx, ny, nz, curved)
        type(obj_block), intent(out) :: blk
        integer,         intent(in)  :: nx, ny, nz
        logical,         intent(in)  :: curved
        integer :: i, j, k
        blk%Nx = nx; blk%Ny = ny; blk%Nz = nz
        allocate(blk%node(3, 0:nx, 0:ny, 0:nz))
        do k = 0, nz; do j = 0, ny; do i = 0, nx
            if (curved) then
                blk%node(:, i, j, k) = [ i*H + 0.15_real64*H*sin(0.7_real64*j)*cos(0.5_real64*k), &
                                         j*H + 0.12_real64*H*sin(0.9_real64*i),                   &
                                         k*H + 0.10_real64*H*sin(0.6_real64*(i+j)) ]
            else
                blk%node(:, i, j, k) = [i*H, j*H, k*H]
            endif
        enddo; enddo; enddo
        allocate(blk%center(3, 1:nx, 1:ny, 1:nz))
        call blk%compute_geometry
        call blk%precomputeMetric
    end subroutine buildGeo

    !> Two boxes along x joined by 101 (A: 3 cells of H, B: 2 of hb), the setup pass, a unit deposit in
    !  A's interface dual cell of vertex (3,1,1), the combine, then finalize and reduce per block.
    subroutine s5_u2_two_block(label, hb)
        character(*), intent(in) :: label
        real(real64), intent(in) :: hb
        integer, parameter :: NXA = 3, NXB = 2, NY = 3, NZ = 2
        type(obj_block)       :: geo(2)
        type(obj_flowblock)   :: gas(2)
        type(obj_sourceblock) :: src(2)
        type(obj_eulerblock)  :: eul(2,1)
        type(cellLink), allocatable    :: links(:)
        type(linkDonors)               :: donors
        type(vertexClass), allocatable :: classes(:)
        real(real64) :: share, vol, ms(2), qs(2), es(2), eu(4,2)
        integer :: b, i, j, k, nmoved, d(3)

        mesh2D = .false.
        call buildBox(geo(1), gas(1), NXA, NY, NZ, 0.0_real64, H)
        call buildBox(geo(2), gas(2), NXB, NY, NZ, NXA*H, hb)
        call joinX(geo)
        call buildInterfaceLinks(geo, links, donors)
        call assembleInterfaceDual(geo, gas, links, classes, nmoved)

        d = [NXA + 1, 2, 2]
        do b = 1, 2
            call src(b)%allocate(1, geo(b)%Nx + 1, NY + 1, NZ + 1)
            src(b)%sourceMass = 0._real64; src(b)%sourceMom = 0._real64; src(b)%sourceEn = 0._real64
            call eul(b,1)%allocate(geo(b)%Nx + 1, NY + 1, NZ + 1)
            eul(b,1)%density = 0._real64; eul(b,1)%np = 0._real64
            eul(b,1)%velocity = 0._real64; eul(b,1)%temperature = 0._real64
        enddo
        src(1)%sourceMass(1, d(1), d(2), d(3)) = 1._real64
        src(1)%sourceMom (1, d(1), d(2), d(3)) = 2._real64
        src(1)%sourceEn  (   d(1), d(2), d(3)) = 3._real64
        !> euler accumulators hold the deposit over the dual cell's volume
        vol = gas(1)%cellVol(d(1), d(2), d(3))
        eul(1,1)%density    (d(1), d(2), d(3)) = 1._real64/vol
        eul(1,1)%np         (d(1), d(2), d(3)) = 2._real64/vol
        eul(1,1)%velocity(1, d(1), d(2), d(3)) = 3._real64/vol
        eul(1,1)%temperature(d(1), d(2), d(3)) = 4._real64/vol

        call combineInterfaceDual(src, eul, classes, gas)
        do b = 1, 2
            call src(b)%finalize(geo(b))
            call eul(b,1)%reduce(geo(b))
            ms(b) = sum(src(b)%sourceMass(1,:,:,:))
            qs(b) = sum(src(b)%sourceMom(1,:,:,:))
            es(b) = sum(src(b)%sourceEn)
            eu(:,b) = 0._real64
            do k = 1, NZ; do j = 1, NY; do i = 1, geo(b)%Nx
                vol = sum(geo(b)%subVol(:, i, j, k))
                eu(:,b) = eu(:,b) + vol*[eul(b,1)%density(i,j,k), eul(b,1)%np(i,j,k), &
                                         eul(b,1)%velocity(1,i,j,k), eul(b,1)%temperature(i,j,k)]
            enddo; enddo; enddo
        enddo

        share = H/(H + hb)
        call assertRel(label//' : source mass of A == h_A/(h_A+h_B)', ms(1), share, TOL14)
        call assertRel(label//' : source mass of B == h_B/(h_A+h_B)', ms(2), 1._real64 - share, TOL14)
        call assertRel(label//' : mass  total of A + B == 1', ms(1) + ms(2), 1._real64, TOL14)
        call assertRel(label//' : mom   total of A + B == 1', (qs(1) + qs(2))/2._real64, 1._real64, TOL14)
        call assertRel(label//' : en    total of A + B == 1', (es(1) + es(2))/3._real64, 1._real64, TOL14)
        print '(a,a,2es13.6)', '  [..] ', label//' : euler rho total of A, B', eu(1,1), eu(1,2)
        call assertRel(label//' : euler rho   total of A + B == 1', eu(1,1) + eu(1,2), 1._real64, TOL14)
        call assertRel(label//' : euler n     total of A + B == 1', (eu(2,1) + eu(2,2))/2._real64, 1._real64, TOL14)
        call assertRel(label//' : euler rho u total of A + B == 1', (eu(3,1) + eu(3,2))/3._real64, 1._real64, TOL14)
        call assertRel(label//' : euler rho e total of A + B == 1', (eu(4,1) + eu(4,2))/4._real64, 1._real64, TOL14)
        do b = 1, 2
            call geo(b)%freeBlock; call gas(b)%freeBlock
        enddo
    end subroutine s5_u2_two_block

    !> A box of nx x ny x nz cells (x0 + i hx, j H, k H) and its gas dual, built as allocate_blocks does.
    subroutine buildBox(blk, gas, nx, ny, nz, x0, hx)
        type(obj_block),     intent(out) :: blk
        type(obj_flowblock), intent(out) :: gas
        integer,             intent(in)  :: nx, ny, nz
        real(real64),        intent(in)  :: x0, hx
        integer :: i, j, k
        blk%Nx = nx; blk%Ny = ny; blk%Nz = nz
        allocate(blk%node(3, 0:nx, 0:ny, 0:nz))
        do k = 0, nz; do j = 0, ny; do i = 0, nx
            blk%node(:, i, j, k) = [x0 + i*hx, j*H, k*H]
        enddo; enddo; enddo
        allocate(blk%center(3, 1:nx, 1:ny, 1:nz))
        call blk%compute_geometry
        call blk%precomputeDualWeights
        call blk%precomputeMetric
        gas%Nx = nx; gas%Ny = ny; gas%Nz = nz
        allocate(gas%node(3, 0:nx+1, 0:ny+1, 0:nz+1))
        call fill_dual_nodes(gas, blk, .false.)
        call computeSkewFlag(gas, .false.)
        call gas%precomputeDualMetric(.false., blk)
        call blk%computeBBox
        call gas%computeBBox
    end subroutine buildBox

    !> Block 1's face 2 joined to block 2's face 1 by 101, cell for cell.
    subroutine joinX(geo)
        type(obj_block), intent(inout) :: geo(2)
        integer :: j, k
        do k = 1, geo(1)%Nz; do j = 1, geo(1)%Ny
            geo(1)%face(2)%cell(j,k)%bcdef          = 101
            geo(1)%face(2)%cell(j,k)%connection     = [2, 1, j, k]
            geo(1)%face(2)%cell(j,k)%connectionFace = 1
            geo(1)%face(2)%cell(j,k)%orient         = [1, 0, 0, 1]
            geo(2)%face(1)%cell(j,k)%bcdef          = 101
            geo(2)%face(1)%cell(j,k)%connection     = [1, geo(1)%Nx, j, k]
            geo(2)%face(1)%cell(j,k)%connectionFace = 2
            geo(2)%face(1)%cell(j,k)%orient         = [1, 0, 0, 1]
        enddo; enddo
    end subroutine joinX

    !> Relative comparison; prints the measured value so a RED run reports the ratio.
    subroutine assertRel(what, got, want, tol)
        character(*), intent(in) :: what
        real(real64), intent(in) :: got, want, tol
        real(real64) :: err
        err = abs(got - want) / max(abs(want), tiny(1._real64))
        if (err <= tol) then
            print '(a,a,a,es13.6,a,es9.2,a)', '  [PASS] ', what, ' : ', got, '  (err ', err, ')'
        else
            print '(a,a,a,es13.6,a,es9.2,a)', '  [FAIL] ', what, ' : ', got, '  (err ', err, ')'
            nfail = nfail + 1
        endif
    end subroutine assertRel

end program test_source_reduction
