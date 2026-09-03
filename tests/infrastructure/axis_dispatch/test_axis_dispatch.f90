program test_axis_dispatch
    !
    ! Orientation independence of the bcdef-200 face dispatch.
    !
    ! ATLAS emits 200 for BOTH the wedge faces and the axis face of an axisymmetric wedge
    ! block. They need opposite treatment -- the wedge folds by a rotation about x, the axis
    ! must reflect, because a rotation about x cannot change a particle's radius -- and the
    ! old code told them apart by face INDEX (`f==5 .or. f==6`). Nothing guarantees that
    ! ordering: the block may be built with the azimuthal direction along i or j instead.
    !
    ! Production `faceAzimuth` decides from the face's own geometry: the azimuthal projection
    ! of its outward normal. This test pins that contract at EVERY face index, for three
    ! different block orientations:
    !
    !   D1  wedge on the k-faces (5/6)  -- the conventional ordering
    !   D2  wedge on the i-faces (1/2)
    !   D3  wedge on the j-faces (3/4)
    !
    ! For each: the two wedge faces must classify as wedge (|azim| > wedgeAzimTol) with the
    ! correct SIGN (+ on the +theta side, which selects the -delthe fold), and the other four
    ! -- axis face, outer radial face, and the two axial faces -- must classify as reflect.
    !
    ! D4 additionally pins the equivalence the refactor must preserve: on the D1 ordering the
    ! geometric rule must agree with the retired `f==6 -> -delthe, f==5 -> +delthe` rule.
    ! No e2e case covers this: axisymFold re-sectors the particle every outer step, so bcDef
    ! sees zero wedge-face 200 calls on both db-2daxi and JPL (measured). This is the coverage.
    !
    ! Oracle is analytic throughout: face normals are constructed from the cylindrical frame,
    ! never read back from the routine under test. Tolerances are machine-precision because
    ! every quantity here is an exact dot product of unit vectors.
    !
    use, intrinsic :: iso_fortran_env, only: real64
    use IGLOO_bcBox, only: faceAzimuth, wedgeAzimTol       ! PRODUCTION under test
    implicit none

    real(real64), parameter :: EPS  = epsilon(1.0_real64)
    real(real64), parameter :: TOL  = 64.0_real64*EPS      ! exact dot products of unit vectors
    real(real64), parameter :: RIN  = 1.0e-8_real64        ! GRIB's axis line
    real(real64), parameter :: ROUT = 5.0e-3_real64        ! first radial cell, JPL-like
    real(real64), parameter :: X0 = 0.10_real64, X1 = 0.11_real64
    real(real64), parameter :: DTH = 0.0174533_real64      ! 1 deg wedge, as GRIB writes

    integer      :: nfail
    real(real64) :: vert(3,8), nrm(3,6), azim
    integer      :: wedgeLo, wedgeHi, f

    nfail = 0

    ! ---- D1: azimuthal along k -> wedge on faces 5/6 -------------------------------------
    call buildWedge(1, 2, 3, vert, nrm, wedgeLo, wedgeHi)   ! axial=i, radial=j, azim=k
    call checkOrientation('D1 wedge on k-faces', vert, nrm, wedgeLo, wedgeHi)

    ! ---- D2: azimuthal along i -> wedge on faces 1/2 -------------------------------------
    call buildWedge(2, 3, 1, vert, nrm, wedgeLo, wedgeHi)   ! axial=j, radial=k, azim=i
    call checkOrientation('D2 wedge on i-faces', vert, nrm, wedgeLo, wedgeHi)

    ! ---- D3: azimuthal along j -> wedge on faces 3/4 -------------------------------------
    call buildWedge(3, 1, 2, vert, nrm, wedgeLo, wedgeHi)   ! axial=k, radial=i, azim=j
    call checkOrientation('D3 wedge on j-faces', vert, nrm, wedgeLo, wedgeHi)

    ! ---- D4: on the conventional ordering, reproduce the retired index rule ---------------
    call buildWedge(1, 2, 3, vert, nrm, wedgeLo, wedgeHi)
    do f = 5, 6
        azim = faceAzimuth(vert, f, nrm(:,f))
        ! retired rule: f==6 -> -delthe, f==5 -> +delthe.  new rule: -sign(delthe, azim)
        call expect('D4 fold sign matches the retired index rule', f,                        &
                    -sign(DTH, azim) - merge(-DTH, DTH, f==6), TOL*DTH)
    enddo

    if (nfail == 0) then
        write(*,'(a)') '[PASS] test_axis_dispatch: bcdef-200 dispatch is orientation-independent.'
        stop 0
    else
        write(*,'(a,i0,a)') '[FAIL] test_axis_dispatch: ', nfail, ' assertion(s) failed.'
        stop 1
    endif

contains

    !> Build a hex spanning [X0,X1] axially, [RIN,ROUT] radially and [-DTH/2,+DTH/2]
    !  azimuthally about the x-axis, with the three logical directions mapped freely onto
    !  (axial, radial, azimuthal). `la/lr/lt` name which logical direction (1=i, 2=j, 3=k)
    !  carries the axial / radial / azimuthal extent. Also returns the exact outward unit
    !  normal of each of the six faces, and which face indices are the wedge pair.
    !
    !  Vertex convention, read off the production `guide` table:
    !    v1=(0,0,0) v2=(0,0,1) v3=(0,1,1) v4=(0,1,0) v5=(1,0,0) v6=(1,0,1) v7=(1,1,1) v8=(1,1,0)
    !  in logical (i,j,k). Faces: 1:i=0 2:i=1 3:j=0 4:j=1 5:k=0 6:k=1, i.e. direction d owns
    !  faces (2d-1, 2d) low/high.
    subroutine buildWedge(la, lr, lt, vert, nrm, wLo, wHi)
        integer,      intent(in)  :: la, lr, lt
        real(real64), intent(out) :: vert(3,8), nrm(3,6)
        integer,      intent(out) :: wLo, wHi
        integer      :: c, lg(3,8)
        real(real64) :: xx, rr, th
        lg = reshape([0,0,0,  0,0,1,  0,1,1,  0,1,0,                &
                      1,0,0,  1,0,1,  1,1,1,  1,1,0], [3,8])
        do c = 1, 8
            xx = X0  + real(lg(la,c),real64)*(X1 - X0)
            rr = RIN + real(lg(lr,c),real64)*(ROUT - RIN)
            th = -0.5_real64*DTH + real(lg(lt,c),real64)*DTH
            vert(:,c) = [xx, rr*cos(th), rr*sin(th)]
        enddo
        nrm = 0.0_real64
        nrm(:,2*la-1) = [-1.0_real64, 0.0_real64, 0.0_real64]   ! -axial
        nrm(:,2*la  ) = [ 1.0_real64, 0.0_real64, 0.0_real64]   ! +axial
        nrm(:,2*lr-1) = -rhat(0.0_real64)                       ! -radial: the AXIS face
        nrm(:,2*lr  ) =  rhat(0.0_real64)                       ! +radial: the outer face
        nrm(:,2*lt-1) = -that(-0.5_real64*DTH)                  ! -theta wedge face
        nrm(:,2*lt  ) =  that( 0.5_real64*DTH)                  ! +theta wedge face
        wLo = 2*lt - 1
        wHi = 2*lt
    end subroutine buildWedge

    pure function rhat(th) result(u)
        real(real64), intent(in) :: th
        real(real64) :: u(3)
        u = [0.0_real64, cos(th), sin(th)]
    end function rhat

    pure function that(th) result(u)
        real(real64), intent(in) :: th
        real(real64) :: u(3)
        u = [0.0_real64, -sin(th), cos(th)]
    end function that

    subroutine checkOrientation(tag, vert, nrm, wLo, wHi)
        character(len=*), intent(in) :: tag
        real(real64),     intent(in) :: vert(3,8), nrm(3,6)
        integer,          intent(in) :: wLo, wHi
        integer      :: f
        real(real64) :: a
        do f = 1, 6
            a = faceAzimuth(vert, f, nrm(:,f))
            if (f == wHi) then
                call expect(tag//' +theta wedge face azim=+1', f, a - 1.0_real64, TOL)
            else if (f == wLo) then
                call expect(tag//' -theta wedge face azim=-1', f, a + 1.0_real64, TOL)
            else if (abs(a) > wedgeAzimTol) then
                write(*,'(a,a,a,i0,a,es12.5)') '[FAIL] ', tag, ' f=', f,                     &
                    ' misclassified as a wedge face, |azim|=', abs(a)
                nfail = nfail + 1
            endif
        enddo
    end subroutine checkOrientation

    subroutine expect(what, f, residual, tol)
        character(len=*), intent(in) :: what
        integer,          intent(in) :: f
        real(real64),     intent(in) :: residual, tol
        if (abs(residual) > tol) then
            write(*,'(a,a,a,i0,a,es12.5,a,es12.5)') '[FAIL] ', what, ' f=', f,               &
                ' residual=', abs(residual), ' > tol=', tol
            nfail = nfail + 1
        endif
    end subroutine expect

end program test_axis_dispatch
