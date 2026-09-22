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
    use IGLOO_data_block, only: obj_block, obj_flowblock, fillGhostPartners
    use IGLOO_allocation, only: fill_dual_nodes
    use IGLOO_variables,  only: ord2, mesh2D
    use Lib_Equations,    only: sampleGas2D, interp2ndOrder
    implicit none

    real(real64), parameter :: H = 0.25_real64
    integer,      parameter :: NSP = 9              ! rho u v w T mu gamma R k
    integer :: nfail

    nfail = 0

    call g1_uniform_inflow
    call g2_quadratic
    call g3_kplane
    call g4_corner2D
    call g5_partner
    call g6_guard
    call g7_corner3D

    if (nfail == 0) then
        print '(a)', ''
        print '(a)', '[PASS] test_ghost_bc: the ord2 ghost ring honours 300/301 and 201'
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

    !> bcdef has no default initializer -- every face cell must be tagged or the select case
    !  reads garbage.
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
