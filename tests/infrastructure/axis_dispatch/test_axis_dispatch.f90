program test_axis_dispatch
    !
    ! The bcdef-200 face dispatch must depend on neither the face INDEX nor the direction of
    ! the symmetry AXIS.
    !
    ! ATLAS emits 200 for BOTH the wedge faces and the axis face of an axisymmetric wedge
    ! block. They need opposite treatment -- the wedge folds by a rotation about the axis, the
    ! axis face must reflect, because a rotation about the axis cannot change a particle's
    ! radius -- and the old code told them apart by face index (`f==5 .or. f==6`) while
    ! rotating about a hardcoded x. Neither is guaranteed: the block may be built with the
    ! azimuthal direction along i or j, and the symmetry axis need not be x.
    !
    ! Production `faceAzimuth` decides from the face's own geometry, projected on the
    ! azimuthal direction of the `axisDir` frame. This test pins that contract at EVERY face
    ! index, for three block orientations, under two different symmetry axes:
    !
    !   D1/D2/D3   wedge on the k-, i-, j-faces, symmetry axis = x   (the shipped frame)
    !   D5/D6/D7   the same three, symmetry axis = (1,1,1)/sqrt(3)   (a tilted axis)
    !
    ! For each: the two wedge faces must classify as wedge (|azim| > wedgeAzimTol) with the
    ! correct SIGN (+ on the +theta side, which selects the -delthe fold), and the other four
    ! -- axis face, outer radial face, and the two axial faces -- must classify as reflect.
    !
    ! D4 pins the equivalence the refactor must preserve: on the conventional k-ordering about
    ! x, the geometric rule must agree with the retired `f==6 -> -delthe, f==5 -> +delthe`.
    ! No e2e case covers the fold branch: axisymFold re-sectors the particle every outer step,
    ! so bcDef sees zero wedge-face 200 calls on both db-2daxi and JPL (measured). This is it.
    !
    ! Oracle is analytic throughout: face normals are built from the cylindrical frame, never
    ! read back from the routine under test. Tolerances are machine-precision because every
    ! quantity is an exact dot product of unit vectors.
    !
    use, intrinsic :: iso_fortran_env, only: real64
    use IGLOO_bcBox,     only: faceAzimuth, wedgeAzimTol     ! PRODUCTION under test
    use IGLOO_variables, only: axisDir                       ! PRODUCTION axis frame
    implicit none

    real(real64), parameter :: EPS  = epsilon(1.0_real64)
    real(real64), parameter :: TOL  = 256.0_real64*EPS     ! exact dot products of unit vectors
    real(real64), parameter :: RIN  = 1.0e-8_real64        ! GRIB's axis line
    real(real64), parameter :: ROUT = 5.0e-3_real64        ! first radial cell, JPL-like
    real(real64), parameter :: X0 = 0.10_real64, X1 = 0.11_real64
    real(real64), parameter :: DTH = 0.0174533_real64      ! 1 deg wedge, as GRIB writes
    real(real64), parameter :: S3  = 1.0_real64/sqrt(3.0_real64)

    integer      :: nfail
    real(real64) :: vert(3,8), nrm(3,6), azim
    real(real64) :: aX(3), eX(3), aT(3), eT(3)
    integer      :: wedgeLo, wedgeHi, f

    nfail = 0

    aX = [1.0_real64, 0.0_real64, 0.0_real64]              ! shipped axis
    eX = [0.0_real64, 1.0_real64, 0.0_real64]              ! shipped azimuth origin
    aT = [S3, S3, S3]                                      ! tilted axis
    eT = unitPerp(aT)                                      ! any unit vector normal to it

    ! ---- axis = x -------------------------------------------------------------------------
    call runFrame('D1 axis=x  wedge on k', 1, 2, 3, aX, eX)
    call runFrame('D2 axis=x  wedge on i', 2, 3, 1, aX, eX)
    call runFrame('D3 axis=x  wedge on j', 3, 1, 2, aX, eX)

    ! ---- D4: reproduce the retired index rule on the conventional ordering -----------------
    axisDir = aX
    call buildWedge(1, 2, 3, aX, eX, vert, nrm, wedgeLo, wedgeHi)
    do f = 5, 6
        azim = faceAzimuth(vert, f, nrm(:,f))
        ! retired rule: f==6 -> -delthe, f==5 -> +delthe.  new rule: -sign(delthe, azim)
        call expect('D4 fold sign matches the retired index rule', f,                        &
                    -sign(DTH, azim) - merge(-DTH, DTH, f==6), TOL*DTH)
    enddo

    ! ---- tilted axis ----------------------------------------------------------------------
    call runFrame('D5 axis=(1,1,1) wedge on k', 1, 2, 3, aT, eT)
    call runFrame('D6 axis=(1,1,1) wedge on i', 2, 3, 1, aT, eT)
    call runFrame('D7 axis=(1,1,1) wedge on j', 3, 1, 2, aT, eT)

    axisDir = aX                                            ! leave the module as we found it

    if (nfail == 0) then
        write(*,'(a)') '[PASS] test_axis_dispatch: bcdef-200 dispatch is independent of both '// &
                       'face index and axis direction.'
        stop 0
    else
        write(*,'(a,i0,a)') '[FAIL] test_axis_dispatch: ', nfail, ' assertion(s) failed.'
        stop 1
    endif

contains

    !> Point the production frame at `a`, build a wedge with the azimuthal extent on logical
    !  direction `lt`, and check every face classifies from geometry alone.
    subroutine runFrame(tag, la, lr, lt, a, e1)
        character(len=*), intent(in) :: tag
        integer,          intent(in) :: la, lr, lt
        real(real64),     intent(in) :: a(3), e1(3)
        real(real64) :: vv(3,8), nn(3,6)
        integer      :: wLo, wHi
        axisDir = a                                         ! production reads this
        call buildWedge(la, lr, lt, a, e1, vv, nn, wLo, wHi)
        call checkOrientation(tag, vv, nn, wLo, wHi)
    end subroutine runFrame

    !> Build a hex spanning [X0,X1] along the axis `a`, [RIN,ROUT] radially and [-DTH/2,+DTH/2]
    !  azimuthally about `a`, with the three logical directions mapped freely onto
    !  (axial, radial, azimuthal). `la/lr/lt` name which logical direction (1=i, 2=j, 3=k)
    !  carries each extent. Also returns the exact outward unit normal of all six faces and
    !  which face indices form the wedge pair.
    !
    !  Vertex convention, read off the production `guide` table:
    !    v1=(0,0,0) v2=(0,0,1) v3=(0,1,1) v4=(0,1,0) v5=(1,0,0) v6=(1,0,1) v7=(1,1,1) v8=(1,1,0)
    !  in logical (i,j,k). Faces: 1:i=0 2:i=1 3:j=0 4:j=1 5:k=0 6:k=1, i.e. logical direction d
    !  owns faces (2d-1, 2d) as its low/high sides.
    subroutine buildWedge(la, lr, lt, a, e1, vert, nrm, wLo, wHi)
        integer,      intent(in)  :: la, lr, lt
        real(real64), intent(in)  :: a(3), e1(3)
        real(real64), intent(out) :: vert(3,8), nrm(3,6)
        integer,      intent(out) :: wLo, wHi
        integer      :: c, lg(3,8)
        real(real64) :: xx, rr, th, e2(3)
        lg = reshape([0,0,0,  0,0,1,  0,1,1,  0,1,0,                &
                      1,0,0,  1,0,1,  1,1,1,  1,1,0], [3,8])
        e2 = xprod(a, e1)
        do c = 1, 8
            xx = X0  + real(lg(la,c),real64)*(X1 - X0)
            rr = RIN + real(lg(lr,c),real64)*(ROUT - RIN)
            th = -0.5_real64*DTH + real(lg(lt,c),real64)*DTH
            vert(:,c) = xx*a + rr*(cos(th)*e1 + sin(th)*e2)
        enddo
        nrm = 0.0_real64
        nrm(:,2*la-1) = -a                                   ! -axial
        nrm(:,2*la  ) =  a                                   ! +axial
        nrm(:,2*lr-1) = -rdir(0.0_real64, e1, e2)            ! -radial: the AXIS face
        nrm(:,2*lr  ) =  rdir(0.0_real64, e1, e2)            ! +radial: the outer face
        nrm(:,2*lt-1) = -tdir(-0.5_real64*DTH, e1, e2)       ! -theta wedge face
        nrm(:,2*lt  ) =  tdir( 0.5_real64*DTH, e1, e2)       ! +theta wedge face
        wLo = 2*lt - 1
        wHi = 2*lt
    end subroutine buildWedge

    pure function rdir(th, e1, e2) result(u)
        real(real64), intent(in) :: th, e1(3), e2(3)
        real(real64) :: u(3)
        u = cos(th)*e1 + sin(th)*e2
    end function rdir

    pure function tdir(th, e1, e2) result(u)
        real(real64), intent(in) :: th, e1(3), e2(3)
        real(real64) :: u(3)
        u = -sin(th)*e1 + cos(th)*e2
    end function tdir

    pure function xprod(p, q) result(u)
        real(real64), intent(in) :: p(3), q(3)
        real(real64) :: u(3)
        u = [p(2)*q(3)-p(3)*q(2), p(3)*q(1)-p(1)*q(3), p(1)*q(2)-p(2)*q(1)]
    end function xprod

    !> Any unit vector orthogonal to `a`: cross with whichever cardinal axis is least aligned.
    pure function unitPerp(a) result(u)
        real(real64), intent(in) :: a(3)
        real(real64) :: u(3), c(3)
        integer :: k
        k = minloc(abs(a), dim=1)
        c = 0.0_real64; c(k) = 1.0_real64
        u = xprod(a, c)
        u = u/norm2(u)
    end function unitPerp

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
