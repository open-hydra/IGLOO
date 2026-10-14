program test_ghost_bc
    !
    ! The ord2 gas ghost ring must know the boundary conditions.
    !
    ! Under gas-order 2 the gas is sampled on a DUAL mesh whose first and last cell in every
    ! direction straddle the domain boundary: half of each lies outside, on a ring of ghost
    ! nodes. fillGhostGradient fills that ring, and until this gate existed it filled it by
    ! blind linear extrapolation, 2q1 - q2, on every face regardless of what the face IS.
    !
    ! At a gas-solid plane that is wrong in the one way that matters to a parcel: the sampled
    ! normal velocity at the plane is (v_ghost + v_1)/2 = (3 v_1 - v_2)/2, not 0, so a parcel
    ! in the boundary dual row is dragged through a plane the gas cannot cross. The bc-aware
    ! fill had been written (ghostState, in the setup-time post-pass of read_cdp_bc_file) but
    ! was overwritten by fillGhostGradient before any parcel read it, and it never ran again
    ! on sweeps >= 1 at all -- import_gas rewrites only the interior.
    !
    ! The fix mirrors the INTERIOR partner's normal component into the ghost,
    !
    !     v_g <- v_g - ((v_g + v_i).n) n        =>   v_g.n = -v_i.n
    !
    ! keeping the tangential part from the linear fill. Mirroring the linear ghost instead
    ! (v_lin - 2(v_lin.n)n) zeroes the plane value only for a uniform normal profile: G2 is
    ! the case that separates them. Keeping the tangential part from the linear fill is what
    ! leaves the swirl-wedge fixtures bit-exact (v.n = 0 there, so the update is a no-op),
    ! where a copy-then-mirror would move their W = omega*y profile.
    !
    ! Cases. Box 4 x 3 (x 2), H = 0.25, so the geo cell centres are y_c = 0.125, 0.375, 0.625
    ! and the dual straddles y = 0 between the ghost row j = 0 and the first interior row.
    ! Every sample is taken AT the plane, through the production path
    ! (gasProperties -> getVertices -> sampleGas2D / interp2ndOrder), never by reading the
    ! ghost array directly -- except G5, where the ghost IS the observable.
    !
    !   G1  2D uniform inflow   f3=300      v=(1,-0.1,0)     RED v_y = -0.1        GREEN 0
    !   G2  2D quadratic        f3=300      v_y=-a y_c^2     RED +3aH^2/4=0.046875 GREEN 0
    !   G3  3D k-plane          f5=300      v=(1,0,-0.1)     RED v_z = -0.1        GREEN 0
    !   G4  2D corner           f1=f3=300   v=(1,-0.1,0)     RED (1,-0.1)          GREEN (0,0)
    !   G5  partner copy        f3/f4=201   T linear         RED T_g = 287.5       GREEN 362.5
    !   G6  positivity guard    f3=300      rho 0.4 / 1.0    RED 0.1               GREEN 0.4
    !   G7  3D corner           f1=f3=f5=300 v=(-.2,-.1,-.1) RED (-.2,-.1,-.1)     GREEN (0,0,0)
    !
    ! Connected faces (two blocks; the ghost IS the observable, compared with the field at the node):
    !
    !   U8  3D, dx 0.3 | 0.6, 101       T, u linear     RED copy: off by b*dx1/2       GREEN exact (1e-12)
    !   U9  3D, dx 0.3 | 0.3, 101       T quadratic     RED snap off: roundoff         GREEN the copy, bit for bit
    !   U10 2D, dx 0.3 | 0.6, 101       T quadratic     RED in-place: order-dependent  GREEN both block orders equal
    !   U11 2D+3D warped, dx 0.3 | 0.6  T, u linear     RED absolute Newton exit       GREEN exact (1e-12)
    !   G6u 2D, overlap, 102            T quadratic     RED own fill: -2 q h^2         GREEN within q h^2/4
    !
    ! G1's T check is a control, not a target: the profile is linear, so both fills reproduce
    ! it exactly and T = 300 at the plane either way. It fails only if the bc pass corrupts a
    ! scalar it must not touch. G6 is its counterpart: there the linear fill DOES go negative
    ! and the guard must fall back to zero-gradient.
    !
    ! ord2 = .true. must be set before any gasProperties call: without it gasProperties reads
    ! the CELL value and every case samples cell (2,1,1), so G1 would read -0.1 both before
    ! and after the fix and the gate could never go green.
    !
    use, intrinsic :: iso_fortran_env, only: real64
    use IGLOO_data_block, only: obj_block, obj_flowblock, obj_bc_cell, fillGhostPartners, buildGhostLinks
    use IGLOO_allocation, only: fill_dual_nodes
    use IGLOO_variables,  only: ord2, mesh2D
    use Lib_Equations,    only: sampleGas2D, interp2ndOrder
    implicit none

    real(real64), parameter :: H = 0.25_real64
    integer,      parameter :: NSP = 9              ! rho u v w T mu gamma R k
    integer :: nfail
    !> field of the two-block cases: T = fa + fb*x + fc*y + fd*z + fq*x^2, u = 1 + (T - fa)/100
    real(real64) :: fa = 300.0_real64, fb = 0.0_real64, fc = 0.0_real64, fd = 0.0_real64, fq = 0.0_real64

    nfail = 0

    call g1_uniform_inflow
    call g2_quadratic
    call g3_kplane
    call g4_corner2D
    call g5_partner
    call g6_guard
    call g7_corner3D
    call u8_stretched
    call u9_snap
    call u10_two_stage
    call u11_warped(.true.)
    call u11_warped(.false.)
    call g6u_chimera

    if (nfail == 0) then
        print '(a)', ''
        print '(a)', '[PASS] test_ghost_bc: the ord2 ghost ring honours 300/301, 201 and the connected faces'
    else
        print '(a)', ''
        print '(a,i0,a)', '[FAIL] test_ghost_bc: ', nfail, ' assertion(s) failed'
        error stop 1
    endif

contains

    !> G1 -- a uniform inflow at a solid plane. The linear fill reproduces the uniform v_y
    !  exactly, so the plane sample is the full -0.1: the crudest form of the defect.
    subroutine g1_uniform_inflow
        type(obj_block)     :: geo
        type(obj_flowblock) :: sol
        real(real64) :: gs(NSP)
        call build2D(geo, sol, 4, 3)
        call tagAll(geo, 400); call tag(geo, 3, 300)
        call setUniform(sol, [1.0_real64, -0.1_real64, 0.0_real64])
        call sol%fillGhostGradient(geo)
        call plane2D(sol, 2, [0.25_real64, 0.0_real64, 0.0_real64], gs)
        call expectLE('G1 v_y at the 300 plane', gs(3), 1.0e-14_real64)
        call expectEq('G1 u   at the 300 plane', gs(2), 1.0_real64,   1.0e-14_real64)
        call expectEq('G1 T   at the 300 plane', gs(5), 300.0_real64, 1.0e-12_real64)
        call freeAll(geo, sol)
    end subroutine g1_uniform_inflow

    !> G2 -- the case that fixes the FORM of the mirror. v_y = -a y^2 is reproduced by neither
    !  fill, and mirroring the linear ghost would leave +a H^2/4 at the plane; mirroring the
    !  interior partner's normal component leaves 0 for any profile.
    subroutine g2_quadratic
        type(obj_block)     :: geo
        type(obj_flowblock) :: sol
        real(real64) :: gs(NSP)
        integer :: j
        call build2D(geo, sol, 4, 3)
        call tagAll(geo, 400); call tag(geo, 3, 300)
        call setUniform(sol, [1.0_real64, 0.0_real64, 0.0_real64])
        do j = 1, sol%Ny
            sol%velocity(2, :, j, 1) = -yc(j)**2                 ! a = 1
        enddo
        call sol%fillGhostGradient(geo)
        call plane2D(sol, 2, [0.25_real64, 0.0_real64, 0.0_real64], gs)
        call expectLE('G2 v_y at the 300 plane (quadratic profile)', gs(3), 1.0e-14_real64)
        call freeAll(geo, sol)
    end subroutine g2_quadratic

    !> G3 -- the same defect on a k-face, which only exists in 3D (both fills skip faces 5/6
    !  when Nz = 1).
    subroutine g3_kplane
        type(obj_block)     :: geo
        type(obj_flowblock) :: sol
        real(real64) :: gs(NSP)
        call build3D(geo, sol, 4, 3, 2)
        call tagAll(geo, 400); call tag(geo, 5, 300)
        call setUniform(sol, [1.0_real64, 0.0_real64, -0.1_real64])
        call sol%fillGhostGradient(geo)
        call plane3D(sol, [2, 2, 1], [0.25_real64, 0.25_real64, 0.0_real64], gs)
        call expectLE('G3 v_z at the 300 k-plane', gs(4), 1.0e-12_real64)
        call expectEq('G3 u   at the 300 k-plane', gs(2), 1.0_real64, 1.0e-12_real64)
        call freeAll(geo, sol)
    end subroutine g3_kplane

    !> G4 -- two solid planes meeting at a corner. The corner ghost is reached by the cascade,
    !  not by the face pass, so it needs the second (edge) mirror pass; without it the corner
    !  dual cell still samples the full interior velocity.
    subroutine g4_corner2D
        type(obj_block)     :: geo
        type(obj_flowblock) :: sol
        real(real64) :: gs(NSP)
        call build2D(geo, sol, 4, 3)
        call tagAll(geo, 400); call tag(geo, 1, 300); call tag(geo, 3, 300)
        call setUniform(sol, [1.0_real64, -0.1_real64, 0.0_real64])
        call sol%fillGhostGradient(geo)
        call plane2D(sol, 1, [0.0_real64, 0.0_real64, 0.0_real64], gs)
        call expectLE('G4 u   at the 2D corner', gs(2), 1.0e-14_real64)
        call expectLE('G4 v_y at the 2D corner', gs(3), 1.0e-14_real64)
        call freeAll(geo, sol)
    end subroutine g4_corner2D

    !> G5 -- 201 partner copy. No ord2 fixture in the tree carries 101 or 201, so this is the
    !  only coverage the partner branch has. The ghost itself is the observable: a partner copy
    !  is not a plane condition and shows up nowhere in a sample at the boundary.
    subroutine g5_partner
        type(obj_block)     :: geo(1)
        type(obj_flowblock) :: sol(1)
        integer :: i
        call build2D(geo(1), sol(1), 4, 3)
        call tagAll(geo(1), 400); call tag(geo(1), 3, 201); call tag(geo(1), 4, 201)
        do i = 1, geo(1)%Nx
            geo(1)%face(3)%cell(i,1)%connection = [1, i, sol(1)%Ny, 1]
            geo(1)%face(4)%cell(i,1)%connection = [1, i, 1,         1]
            geo(1)%face(3)%cell(i,1)%connectionFace = 4
            geo(1)%face(4)%cell(i,1)%connectionFace = 3
        enddo
        call setUniform(sol(1), [1.0_real64, 0.0_real64, 0.0_real64])
        call sol(1)%fillGhostGradient(geo(1))
        call fillGhostPartners(sol, geo)
        call expectEq('G5 T at the 201 ghost row j=0', sol(1)%temperature(2,0,1), &
                      tprof(sol(1)%Ny), 1.0e-12_real64)
        call expectEq('G5 T at the 201 ghost row j=Ny+1', sol(1)%temperature(2,sol(1)%Ny+1,1), &
                      tprof(1), 1.0e-12_real64)
        call freeAll(geo(1), sol(1))
    end subroutine g5_partner

    !> G6 -- the positivity guard. A steep near-wall density makes the linear ghost negative;
    !  a negative density in the sample is worse than a first-order one, so 300/301 falls back
    !  to zero-gradient for the scalars. Only 300/301 does: every other code keeps the linear
    !  fill it has always had.
    subroutine g6_guard
        type(obj_block)     :: geo
        type(obj_flowblock) :: sol
        real(real64) :: gs(NSP)
        call build2D(geo, sol, 4, 3)
        call tagAll(geo, 400); call tag(geo, 3, 300)
        call setUniform(sol, [1.0_real64, 0.0_real64, 0.0_real64])
        sol%density(1,:,1,1) = 0.4_real64
        sol%density(1,:,2,1) = 1.0_real64
        call sol%fillGhostGradient(geo)
        call plane2D(sol, 2, [0.25_real64, 0.0_real64, 0.0_real64], gs)
        call expectEq('G6 rho at the 300 plane (guard -> zero-gradient)', gs(1), &
                      0.4_real64, 1.0e-12_real64)
        call freeAll(geo, sol)
    end subroutine g6_guard

    !> G7 -- three solid planes at a 3D corner. The corner ghost is cascaded from the edge
    !  ghosts, so it is correct only if the edge pass ran before the corner pass; the three
    !  orthogonal mirrors then commute and every component of the 8-node average cancels.
    subroutine g7_corner3D
        type(obj_block)     :: geo
        type(obj_flowblock) :: sol
        real(real64) :: gs(NSP)
        call build3D(geo, sol, 4, 3, 2)
        call tagAll(geo, 400)
        call tag(geo, 1, 300); call tag(geo, 3, 300); call tag(geo, 5, 300)
        call setUniform(sol, [-0.2_real64, -0.1_real64, -0.1_real64])
        call sol%fillGhostGradient(geo)
        call plane3D(sol, [1, 1, 1], [0.0_real64, 0.0_real64, 0.0_real64], gs)
        call expectLE('G7 u   at the 3D corner', gs(2), 1.0e-12_real64)
        call expectLE('G7 v_y at the 3D corner', gs(3), 1.0e-12_real64)
        call expectLE('G7 v_z at the 3D corner', gs(4), 1.0e-12_real64)
        call freeAll(geo, sol)
    end subroutine g7_corner3D


    !> U8 -- a 101 face between cells of 0.3 and 0.6 along x, linear field: every ghost node of
    !  the interface plane (edges and corners too) holds the field at the node -- block 2's from
    !  block 1's interior, block 1's through block 2's own linear fill.
    subroutine u8_stretched
        type(obj_block)     :: geo(2)
        type(obj_flowblock) :: sol(2)
        real(real64) :: eT(2), eU(2)
        fa = 300.0_real64; fb = 100.0_real64; fc = 40.0_real64; fd = 20.0_real64; fq = 0.0_real64
        call buildPair(geo, sol, 4, 0.3_real64, 3, 0.6_real64, 3, 2, 0.0_real64, .false.)
        call fillPair(geo, sol)
        call planeErr(sol(1), sol(1)%Nx+1, eT(1), eU(1))
        call planeErr(sol(2), 0,           eT(2), eU(2))
        call expectLE('U8 T at the fine block''s interface ghosts, |err|/T', eT(1)/fa, 1.0e-12_real64)
        call expectLE('U8 T at the coarse block''s interface ghosts, |err|/T', eT(2)/fa, 1.0e-12_real64)
        call expectLE('U8 u at the interface ghosts of both blocks, |err|', maxval(eU), 1.0e-12_real64)
    end subroutine u8_stretched

    !> U9 -- the uniform pair (0.3 | 0.3): every face ghost lies on the partner's centre to
    !  roundoff; the weights snap to it, so the ghost is the partner's value bit for bit.
    subroutine u9_snap
        type(obj_block)     :: geo(2)
        type(obj_flowblock) :: sol(2)
        real(real64) :: dmax, h
        integer :: j, k, nx, ndiff
        fa = 300.0_real64; fb = 1000.0_real64; fc = 0.0_real64; fd = 0.0_real64; fq = 2000.0_real64
        call buildPair(geo, sol, 4, 0.3_real64, 4, 0.3_real64, 3, 2, 0.0_real64, .false.)
        call fillPair(geo, sol)
        nx = sol(1)%Nx
        h = norm2(sol(2)%node(:,1,1,1) - sol(2)%node(:,0,0,0))
        dmax = 0.0_real64; ndiff = 0
        do k = 1, sol(1)%Nz; do j = 1, sol(1)%Ny
            dmax = max(dmax, norm2(sol(1)%node(:,nx+1,j,k) - sol(2)%node(:,1,j,k))/h, &
                             norm2(sol(2)%node(:,0,j,k)    - sol(1)%node(:,nx,j,k))/h)
            if (sol(1)%temperature(nx+1,j,k) /= sol(2)%temperature(1,j,k)) ndiff = ndiff + 1
            if (sol(2)%temperature(0,j,k)    /= sol(1)%temperature(nx,j,k)) ndiff = ndiff + 1
            if (any(sol(1)%velocity(:,nx+1,j,k) /= sol(2)%velocity(:,1,j,k))) ndiff = ndiff + 1
            if (any(sol(2)%velocity(:,0,j,k)    /= sol(1)%velocity(:,nx,j,k))) ndiff = ndiff + 1
        enddo; enddo
        print '(a,es10.3)', '  [..] U9 face ghost to partner centre, max distance / dual diagonal = ', dmax
        call expectTrue('U9 the snap is exercised: a face ghost off its partner centre by roundoff', dmax > 0.0_real64)
        call expectLE('U9 face ghosts differing from the partner value (T, u)', real(ndiff, real64), 0.0_real64)
    end subroutine u9_snap

    !> U10 -- the two-stage rule on a quadratic field (q h^2 = 1 K, g h = 50 K, h = 0.3): block
    !  1's ghost is evaluated on block 2's stage-1 ring, so both block orders give the same ring.
    subroutine u10_two_stage
        type(obj_block)     :: geoA(2), geoB(2)
        type(obj_flowblock) :: solA(2), solB(2)
        real(real64) :: dfine, dcoarse, x1, x2, want, efine, ecoarse, q
        integer :: j, nf
        q = 1.0_real64/0.09_real64
        fa = 300.0_real64; fb = 50.0_real64/0.3_real64; fc = 0.0_real64; fd = 0.0_real64; fq = q
        call buildPair(geoA, solA, 4, 0.3_real64, 3, 0.6_real64, 3, 1, 0.0_real64, .true.)
        call fillPair(geoA, solA)
        call buildPair(geoB, solB, 4, 0.3_real64, 3, 0.6_real64, 3, 1, 0.0_real64, .true., swap=.true.)
        call fillPair(geoB, solB)
        nf = solA(1)%Nx
        dfine = 0.0_real64; dcoarse = 0.0_real64; efine = 0.0_real64; ecoarse = 0.0_real64
        do j = 1, solA(1)%Ny
            dfine   = max(dfine,   abs(solA(1)%temperature(nf+1,j,1) - solB(2)%temperature(nf+1,j,1)))
            dcoarse = max(dcoarse, abs(solA(2)%temperature(0,j,1)    - solB(1)%temperature(0,j,1)))
            x1 = geoA(2)%center(1,1,j,1); x2 = geoA(2)%center(1,2,j,1)
            want = 1.25_real64*fieldT([x1, 0.0_real64, 0.0_real64]) - 0.25_real64*fieldT([x2, 0.0_real64, 0.0_real64])
            efine = max(efine, abs(solA(1)%temperature(nf+1,j,1) - want))
            ecoarse = max(ecoarse, abs(solA(2)%temperature(0,j,1) - fieldT(solA(2)%node(:,0,j,1))))
        enddo
        print '(a,f9.5,a)', '  [..] U10 fine ghost minus the field at the node: ', &
            solA(1)%temperature(nf+1,1,1) - fieldT(solA(1)%node(:,nf+1,1,1)), ' K'
        call expectLE('U10 fine ghost, block orders 1-2 vs 2-1, |diff|', dfine, 0.0_real64)
        call expectLE('U10 coarse ghost, block orders 1-2 vs 2-1, |diff|', dcoarse, 0.0_real64)
        call expectLE('U10 fine ghost vs 1.25 T(x1) - 0.25 T(x2) of the stage-1 ring, |err|/T', &
                      efine/fa, 1.0e-12_real64)
        call expectLE('U10 coarse ghost vs the field, |err|/(q h^2/4)', ecoarse/(0.25_real64*q*0.09_real64), &
                      1.0_real64 + 1.0e-9_real64)
    end subroutine u10_two_stage

    !> U11 -- the Newton weights: cells of 0.3 | 0.6 along x warped by x(1 + 0.8 y), y(1 + 0.8 x)
    !  (no cell is a parallelogram), linear field: every interface ghost holds the field at the node.
    subroutine u11_warped(planar)
        logical, intent(in) :: planar
        type(obj_block)     :: geo(2)
        type(obj_flowblock) :: sol(2)
        real(real64) :: eT(2), eU(2)
        character(2) :: tag
        fa = 300.0_real64; fb = 100.0_real64; fc = 40.0_real64; fd = 0.0_real64; fq = 0.0_real64
        tag = merge('2D', '3D', planar)
        call buildPair(geo, sol, 4, 0.3_real64, 3, 0.6_real64, 3, merge(1, 2, planar), 0.8_real64, planar)
        call fillPair(geo, sol)
        call planeErr(sol(1), sol(1)%Nx+1, eT(1), eU(1))
        call planeErr(sol(2), 0,           eT(2), eU(2))
        call expectLE('U11 '//tag//' T at the interface ghosts of both blocks, |err|/T', maxval(eT)/fa, 1.0e-12_real64)
        call expectLE('U11 '//tag//' u at the interface ghosts of both blocks, |err|', maxval(eU), 1.0e-12_real64)
    end subroutine u11_warped

    !> G6u -- 102 faces of two overlapping blocks (block 2 shifted by half a cell), quadratic field
    !  (q h^2 = 1 K): each ghost is interpolated between two centres of the other block, within
    !  q h^2/4 of the field; block 1's own linear fill would be off by 2 q h^2. Then the donor
    !  lists pointed at a far cell: the same ring, every 102 node counted as away from its donors.
    subroutine g6u_chimera
        type(obj_block)     :: geo(2)
        type(obj_flowblock) :: sol(2)
        real(real64) :: q, e(2), eU(2), t1(0:4), t2(0:4)
        integer :: nin, nmiss, nfar, j
        q = 1.0_real64/0.09_real64
        fa = 300.0_real64; fb = 50.0_real64/0.3_real64; fc = 0.0_real64; fd = 0.0_real64; fq = q
        call buildOverlap(geo, sol)
        call fillPair(geo, sol)
        call planeErr(sol(1), sol(1)%Nx+1, e(1), eU(1))
        call planeErr(sol(2), 0,           e(2), eU(2))
        call expectLE('G6u 102 ghosts of block 1 vs the field, |err|/(q h^2/4)', e(1)/(0.25_real64*q*0.09_real64), &
                      1.0_real64 + 1.0e-9_real64)
        call expectLE('G6u 102 ghosts of block 2 vs the field, |err|/(q h^2/4)', e(2)/(0.25_real64*q*0.09_real64), &
                      1.0_real64 + 1.0e-9_real64)
        t1 = sol(1)%temperature(sol(1)%Nx+1, 0:4, 1)
        call freePair(geo, sol)
        call buildOverlap(geo, sol)
        do j = 1, geo(1)%Ny
            geo(1)%face(2)%cell(j,1)%chimDonor(:,2) = geo(2)%Nx
        enddo
        call buildGhostLinks(sol, geo, nin, nmiss, nfar)
        call fillPair(geo, sol)
        t2 = sol(1)%temperature(sol(1)%Nx+1, 0:4, 1)
        call expectLE('G6u donors listed far away: ghosts of block 1 changed, max |diff|', maxval(abs(t1 - t2)), &
                      0.0_real64)
        call expectEq('G6u donors listed far away: block 1''s face-2 nodes counted away from their donors', &
                      real(nfar, real64), real(geo(1)%Ny + 2, real64), 0.0_real64)
    end subroutine g6u_chimera


    ! ------------------------------------------------------------------ fixture

    !> y of geo cell row j, and the linear temperature profile on it.
    pure function yc(j) result(y)
        integer, intent(in) :: j
        real(real64) :: y
        y = (real(j,real64) - 0.5_real64)*H
    end function yc

    pure function tprof(j) result(t)
        integer, intent(in) :: j
        real(real64) :: t
        t = 300.0_real64 + 100.0_real64*yc(j)
    end function tprof

    !> 2D fixture: uniform Cartesian geo block through the production geometry path, its dual
    !  node ring through fill_dual_nodes, and a gas block with the ghost ring allocated.
    subroutine build2D(geo, sol, nx, ny)
        type(obj_block),     intent(out) :: geo
        type(obj_flowblock), intent(out) :: sol
        integer,             intent(in)  :: nx, ny
        ord2 = .true.; mesh2D = .true.
        call buildGeo(geo, nx, ny, 1)
        sol%Nx = nx; sol%Ny = ny; sol%Nz = 1
        allocate(sol%node(3, 0:nx+1, 0:ny+1, 1:1))
        call sol%allocate(1, nx, ny, 1, .true.)
        call fill_dual_nodes(sol, geo, .true.)
    end subroutine build2D

    subroutine build3D(geo, sol, nx, ny, nz)
        type(obj_block),     intent(out) :: geo
        type(obj_flowblock), intent(out) :: sol
        integer,             intent(in)  :: nx, ny, nz
        ord2 = .true.; mesh2D = .false.
        call buildGeo(geo, nx, ny, nz)
        sol%Nx = nx; sol%Ny = ny; sol%Nz = nz
        allocate(sol%node(3, 0:nx+1, 0:ny+1, 0:nz+1))
        call sol%allocate(1, nx, ny, nz, .true.)
        call fill_dual_nodes(sol, geo, .false.)
    end subroutine build3D

    subroutine buildGeo(blk, nx, ny, nz)
        type(obj_block), intent(out) :: blk
        integer,         intent(in)  :: nx, ny, nz
        integer :: i, j, k
        blk%Nx = nx; blk%Ny = ny; blk%Nz = nz
        allocate(blk%node(3, 0:nx, 0:ny, 0:nz))
        do k = 0, nz; do j = 0, ny; do i = 0, nx
            blk%node(:, i, j, k) = [i*H, j*H, k*H]
        enddo; enddo; enddo
        allocate(blk%center(3, 1:nx, 1:ny, 1:nz))
        call blk%compute_geometry                 ! allocates face(f)%cell, fills centre + normal
    end subroutine buildGeo

    !> Every face cell tagged with code.
    subroutine tagAll(geo, code)
        type(obj_block), intent(inout) :: geo
        integer,         intent(in)    :: code
        integer :: f
        do f = 1, 6
            geo%face(f)%cell(:,:)%bcdef = code
        enddo
    end subroutine tagAll

    subroutine tag(geo, f, code)
        type(obj_block), intent(inout) :: geo
        integer,         intent(in)    :: f, code
        geo%face(f)%cell(:,:)%bcdef = code
    end subroutine tag

    !> Uniform gas with the G1 scalar set and a linear T profile in y.
    subroutine setUniform(sol, v)
        type(obj_flowblock), intent(inout) :: sol
        real(real64),        intent(in)    :: v(3)
        integer :: j
        sol%density  = 1.2_real64
        sol%velocity(1,:,:,:) = v(1)
        sol%velocity(2,:,:,:) = v(2)
        sol%velocity(3,:,:,:) = v(3)
        sol%mil = 1.8e-5_real64
        sol%kl  = 0.026_real64
        sol%gam = 1.4_real64
        sol%R   = 287.0_real64
        sol%temperature = 300.0_real64
        do j = 1, sol%Ny
            sol%temperature(:, j, :) = tprof(j)
        enddo
    end subroutine setUniform

    !> Gas sampled at p inside dual cell (i,1,1), through the production 2D path.
    subroutine plane2D(sol, i, p, gs)
        type(obj_flowblock), intent(in)  :: sol
        integer,             intent(in)  :: i
        real(real64),        intent(in)  :: p(3)
        real(real64),        intent(out) :: gs(NSP)
        real(real64) :: gn(NSP,4), vert(3,4)
        call sol%gasProperties(gn, [i, 1, 1])
        call sol%getVertices([i, 1, 1], vert, .true.)
        call sampleGas2D(vert, gn, p, NSP, gs)
    end subroutine plane2D

    !> Gas sampled at p inside dual cell si, through the production 3D path.
    subroutine plane3D(sol, si, p, gs)
        type(obj_flowblock), intent(in)  :: sol
        integer,             intent(in)  :: si(3)
        real(real64),        intent(in)  :: p(3)
        real(real64),        intent(out) :: gs(NSP)
        real(real64) :: gn(NSP,8), vert(3,8), xi(3)
        call sol%gasProperties(gn, si)
        call sol%getVertices(si, vert, .false.)
        xi = 0.5_real64
        call interp2ndOrder(vert, gn, p, NSP, xi, gs)
    end subroutine plane3D

    subroutine freeAll(geo, sol)
        type(obj_block),     intent(inout) :: geo
        type(obj_flowblock), intent(inout) :: sol
        call geo%freeBlock
        call sol%freeBlock
    end subroutine freeAll

    pure function fieldT(p) result(t)
        real(real64), intent(in) :: p(3)
        real(real64) :: t
        t = fa + fb*p(1) + fc*p(2) + fd*p(3) + fq*p(1)**2
    end function fieldT

    pure function fieldU(p) result(u)
        real(real64), intent(in) :: p(3)
        real(real64) :: u
        u = 1.0_real64 + (fieldT(p) - fa)/100.0_real64
    end function fieldU

    !> A block on the node lines xn, yn, zn, x scaled by (1 + alpha*y) and y by (1 + alpha*x)
    !  (alpha = 0: a box; else every cell a non-parallelogram bilinear image of the index cube);
    !  planar = one layer (mesh2D). Every face tagged 400.
    subroutine buildLines(geo, sol, xn, yn, zn, alpha, planar)
        type(obj_block),     intent(out) :: geo
        type(obj_flowblock), intent(out) :: sol
        real(real64),        intent(in)  :: xn(0:), yn(0:), zn(0:), alpha
        logical,             intent(in)  :: planar
        integer :: i, j, k, nx, ny, nz
        nx = size(xn) - 1; ny = size(yn) - 1; nz = size(zn) - 1
        ord2 = .true.; mesh2D = planar
        geo%Nx = nx; geo%Ny = ny; geo%Nz = nz
        allocate(geo%node(3, 0:nx, 0:ny, 0:nz))
        do k = 0, nz; do j = 0, ny; do i = 0, nx
            geo%node(:, i, j, k) = [xn(i)*(1.0_real64 + alpha*yn(j)), yn(j)*(1.0_real64 + alpha*xn(i)), zn(k)]
        enddo; enddo; enddo
        allocate(geo%center(3, 1:nx, 1:ny, 1:nz))
        call geo%compute_geometry
        sol%Nx = nx; sol%Ny = ny; sol%Nz = nz
        if (planar) then
            allocate(sol%node(3, 0:nx+1, 0:ny+1, 1:1))
        else
            allocate(sol%node(3, 0:nx+1, 0:ny+1, 0:nz+1))
        endif
        call sol%allocate(1, nx, ny, nz, .true.)
        call fill_dual_nodes(sol, geo, planar)
        call tagAll(geo, 400)
    end subroutine buildLines

    !> Two blocks along x, nx1 cells of dx1 then nx2 of dx2, joined by a 101 face (block 1's face 2,
    !  block 2's face 1); swap puts the second block first in the arrays. The field at the centres.
    subroutine buildPair(geo, sol, nx1, dx1, nx2, dx2, ny, nz, alpha, planar, swap)
        type(obj_block),     intent(out) :: geo(2)
        type(obj_flowblock), intent(out) :: sol(2)
        integer,             intent(in)  :: nx1, nx2, ny, nz
        real(real64),        intent(in)  :: dx1, dx2, alpha
        logical,             intent(in)  :: planar
        logical, optional,   intent(in)  :: swap
        real(real64) :: x1(0:nx1), x2(0:nx2), yn(0:ny), zn(0:nz)
        integer :: i, a, b
        x1 = [(i*dx1, i = 0, nx1)]
        x2 = [(x1(nx1) + i*dx2, i = 0, nx2)]
        yn = [(i*H, i = 0, ny)]
        zn = [(i*H, i = 0, nz)]
        a = 1; b = 2
        if (present(swap)) then
            if (swap) then; a = 2; b = 1; endif
        endif
        call buildLines(geo(a), sol(a), x1, yn, zn, alpha, planar)
        call buildLines(geo(b), sol(b), x2, yn, zn, alpha, planar)
        call link101(geo, a, b)
        call setField(geo(1), sol(1)); call setField(geo(2), sol(2))
    end subroutine buildPair

    !> Block a's face 2 joined to block b's face 1, cell for cell.
    subroutine link101(geo, a, b)
        type(obj_block), intent(inout) :: geo(:)
        integer,         intent(in)    :: a, b
        integer :: j, k
        do k = 1, geo(a)%Nz; do j = 1, geo(a)%Ny
            geo(a)%face(2)%cell(j,k)%bcdef = 101
            geo(a)%face(2)%cell(j,k)%connection = [b, 1, j, k]
            geo(a)%face(2)%cell(j,k)%connectionFace = 1
            geo(b)%face(1)%cell(j,k)%bcdef = 101
            geo(b)%face(1)%cell(j,k)%connection = [a, geo(a)%Nx, j, k]
            geo(b)%face(1)%cell(j,k)%connectionFace = 2
        enddo; enddo
    end subroutine link101

    !> Two one-layer blocks of 6 cells of 0.3 along x overlapping by 1.5 cells (block 2 from x = 1.35),
    !  joined by 102 faces whose layer-1 donors are the two cells of the other block across the
    !  receiver's ghost cell, weight 1/2 each.
    subroutine buildOverlap(geo, sol)
        type(obj_block),     intent(out) :: geo(2)
        type(obj_flowblock), intent(out) :: sol(2)
        real(real64) :: x1(0:6), x2(0:6), yn(0:3), zn(0:1)
        integer :: i, j
        x1 = [(i*0.3_real64, i = 0, 6)]
        x2 = [(1.35_real64 + i*0.3_real64, i = 0, 6)]
        yn = [(i*H, i = 0, 3)]
        zn = [0.0_real64, H]
        call buildLines(geo(1), sol(1), x1, yn, zn, 0.0_real64, .true.)
        call buildLines(geo(2), sol(2), x2, yn, zn, 0.0_real64, .true.)
        do j = 1, 3
            call chimera(geo(1)%face(2)%cell(j,1), 2, [2, 3], j)
            call chimera(geo(2)%face(1)%cell(j,1), 1, [4, 5], j)
        enddo
        call setField(geo(1), sol(1)); call setField(geo(2), sol(2))
    end subroutine buildOverlap

    subroutine chimera(cell, db, di, j)
        type(obj_bc_cell), intent(inout) :: cell
        integer,           intent(in)    :: db, di(2), j
        cell%bcdef = 102
        cell%chimN = [2, 0]
        allocate(cell%chimDonor(2, 4), cell%chimW(2))
        cell%chimDonor(1,:) = [db, di(1), j, 1]
        cell%chimDonor(2,:) = [db, di(2), j, 1]
        cell%chimW = 0.5_real64
    end subroutine chimera

    !> Uniform gas but T and u: the field at the geo cell centres.
    subroutine setField(geo, sol)
        type(obj_block),     intent(in)    :: geo
        type(obj_flowblock), intent(inout) :: sol
        integer :: i, j, k
        sol%density = 1.2_real64; sol%mil = 1.8e-5_real64; sol%kl = 0.026_real64
        sol%gam = 1.4_real64; sol%R = 287.0_real64
        sol%velocity = 0.0_real64; sol%temperature = 0.0_real64
        do k = 1, sol%Nz; do j = 1, sol%Ny; do i = 1, sol%Nx
            sol%temperature(i,j,k) = fieldT(geo%center(:,i,j,k))
            sol%velocity(1,i,j,k)  = fieldU(geo%center(:,i,j,k))
        enddo; enddo; enddo
    end subroutine setField

    !> Stage 1 on both blocks, then the connected-face override.
    subroutine fillPair(geo, sol)
        type(obj_block),     intent(in)    :: geo(2)
        type(obj_flowblock), intent(inout) :: sol(2)
        call sol(1)%fillGhostGradient(geo(1))
        call sol(2)%fillGhostGradient(geo(2))
        call fillGhostPartners(sol, geo)
    end subroutine fillPair

    !> Worst |ghost - field at the node| over the ghost plane i = gi (edges and corners included).
    subroutine planeErr(sol, gi, eT, eU)
        type(obj_flowblock), intent(in)  :: sol
        integer,             intent(in)  :: gi
        real(real64),        intent(out) :: eT, eU
        integer :: j, k, k0, k1
        k0 = merge(1, 0, mesh2D); k1 = merge(1, sol%Nz + 1, mesh2D)
        eT = 0.0_real64; eU = 0.0_real64
        do k = k0, k1; do j = 0, sol%Ny + 1
            eT = max(eT, abs(sol%temperature(gi,j,k) - fieldT(sol%node(:,gi,j,k))))
            eU = max(eU, abs(sol%velocity(1,gi,j,k) - fieldU(sol%node(:,gi,j,k))))
        enddo; enddo
    end subroutine planeErr

    subroutine freePair(geo, sol)
        type(obj_block),     intent(inout) :: geo(2)
        type(obj_flowblock), intent(inout) :: sol(2)
        call freeAll(geo(1), sol(1)); call freeAll(geo(2), sol(2))
    end subroutine freePair


    ! ------------------------------------------------------------------ report

    subroutine expectLE(label, got, tol)
        character(*), intent(in) :: label
        real(real64), intent(in) :: got, tol
        if (abs(got) <= tol) then
            print '(a,a,es13.6)', '  [ok] ', label//' = ', got
        else
            print '(a,a,es13.6,a,es9.2)', '  [FAIL] ', label//' = ', got, ' expected |.| <= ', tol
            nfail = nfail + 1
        endif
    end subroutine expectLE

    subroutine expectTrue(label, ok)
        character(*), intent(in) :: label
        logical,      intent(in) :: ok
        if (ok) then
            print '(a,a)', '  [ok] ', label
        else
            print '(a,a)', '  [FAIL] ', label
            nfail = nfail + 1
        endif
    end subroutine expectTrue

    subroutine expectEq(label, got, want, tol)
        character(*), intent(in) :: label
        real(real64), intent(in) :: got, want, tol
        if (abs(got - want) <= tol*max(1.0_real64, abs(want))) then
            print '(a,a,es13.6)', '  [ok] ', label//' = ', got
        else
            print '(a,a,es13.6,a,es13.6)', '  [FAIL] ', label//' = ', got, ' expected ', want
            nfail = nfail + 1
        endif
    end subroutine expectEq

end program test_ghost_bc
