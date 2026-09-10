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
    ! no e2e case in the suite reaches (every case is mesh2D). D3/D4 run the 2D branch, which
    ! shipped clipped -- they are the regression half, and they also pin that adding the k
    ! direction did not disturb i/j.
    !
    ! The fixture is a uniform Cartesian box built through the PRODUCTION geometry path
    ! (compute_geometry + precomputeMetric), with the dual assembled exactly as
    ! allocation.f90 does it, so the ghost placement under test is the shipped rule and not a
    ! restatement of it.
    !
    use, intrinsic :: iso_fortran_env, only: real64
    use IGLOO_data_block, only: obj_block           ! PRODUCTION type + precomputeDualMetric
    implicit none

    real(real64), parameter :: H = 0.25_real64      ! uniform spacing, all directions
    integer :: nfail

    nfail = 0

    call runCase('D1 3D  4x3x2', 4, 3, 2, .false.)
    call runCase('D2 3D  2x2x2', 2, 2, 2, .false.)
    call runCase('D3 2D  4x3',   4, 3, 1, .true. )
    call runCase('D4 2D  2x2',   2, 2, 1, .true. )

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
    !  must exist. ⚠ The i/j corner block runs in BOTH modes -- allocation.f90 puts it outside
    !  the `.not.mesh2D` guard. Omitting it in 2D leaves node(0,0) undefined and the boundary
    !  quad area is then garbage; this test caught exactly that in its own fixture.
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

end program test_dual_clip
