module IGLOO_Lib_Statistics
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    private
    !> Sampler functions
    public :: RosinRammler, LogNormal, Normal, NormalStandard, ChiSquare

    !> PDFs
    public :: RosinRammlerPDF

    !> Moments calculation
    public :: RonsinRammlerMoments

    !> Injection-diameter sampling
    public :: sampleDiameter, lawCode, initRandomSeed
    public :: DiracDistr, NormalDistr, LogNorDistr, RosRamDistr

    !> Distribution-law codes (stored as col-8 of cell%properties; see IO.f90)
    integer, parameter :: DiracDistr  = 0
    integer, parameter :: NormalDistr = 1
    integer, parameter :: LogNorDistr = 2
    integer, parameter :: RosRamDistr = 3

    logical,  save :: hasBuffer=.false.
    real(R8), save :: z2buffer

contains

    !> standard normal distribution N(0, 1)
    function NormalStandard() result(z)
        use IGLOO_variables, only: pi
        implicit none
        real(R8) :: u1, u2, v1, v2, w, y
        real(R8) :: z
        
        if (hasBuffer) then
            z = z2buffer
            hasBuffer = .false.
        else
            w = 0._R8
            do while (w > 1._R8 .or. w == 0._R8)
                call random_number(u1)
                call random_number(u2)
                v1 = 2._R8*u1 - 1._R8
                v2 = 2._R8*u2 - 1._R8
                w  = v1*v1 + v2*v2
            enddo
            y = sqrt((-2._R8*log(w))/w)
            z = v1*y
            z2buffer  = v2*y
            hasBuffer = .true.
        end if

    end function NormalStandard

    !> Normal distribution N(mean, sigma)
    function Normal(mean, sigma) result(x)
        implicit none
        real(R8), intent(in) :: mean, sigma
        real(R8) :: x
        x = mean + sigma*NormalStandard()
    end function Normal

    !> Rosin-Rammler ditribution wit size parameter x0 and spread parameter n
    function RosinRammler(x0, n) result(x)
        implicit none
        real(R8), intent(in) :: x0, n
        real(R8) :: x

        if (x0<=0._R8 .or. n<=0._R8) error stop ( ' [ERROR] Rosin-Rammler distribution: x0 and n must be positive.' )

        x = 0._R8
        do while (x == 0._R8)
            call random_number(x)
        enddo
        x = x0*(-log(x))**(1._R8/n)

    end function RosinRammler

    pure function RosinRammlerPDF(r,x0,n) result(p)
        implicit none
        real(R8), intent(in) :: r, x0, n
        real(R8) :: p
        p = (n/x0)*(r/x0)**(n-1.0_R8)*exp(-(r/x0)**n)
    end function RosinRammlerPDF

    !> Log-Normal distribution with mean value mu and standard deviation sigma
    function LogNormal(mu,sigma) result(x)
        implicit none
        real(R8), intent(in) :: mu, sigma
        real(R8) :: y, z
        real(R8) :: x

        if (sigma <= 0._R8) error stop ( ' [ERROR] Log-Normal distribution: sigma must be non-negative.' )
        
        z = NormalStandard()
        y = mu + sigma*z
        x = exp(y)

    end function LogNormal

    !> Chi-squared distribution with k degrees of freedom
    function ChiSquare(k) result(x)
        implicit none
        integer, intent(in) :: k
        integer  :: i 
        real(R8) :: x, z

        x = 0._R8
        if (k<1) error stop ( ' [ERROR] Chi-Squared distribution: number of degrees of freedom < 1.' )

        do i = 1, k
            z = NormalStandard()
            x = x + z*z
        enddo

    end function ChiSquare

    !> Truncated integral of r^j * f(r) dr by Simpson's rule over [xMin, xMax].
    pure function RonsinRammlerMoments(order, x0, n, xMin, xMax) result(integral)
        implicit none
        integer,  intent(in) :: order
        real(R8), intent(in) :: x0, n, xMin, xMax
        integer,  parameter  :: Nbin=1000
        integer  :: i
        real(R8) :: h, r, f, coeff
        real(R8) :: integral

        ! Ampiezza del segmento
        h = (xMax - xMin) / real(Nbin)
        if (h <= 0.0_R8) then; integral = 0.0_R8; return; endif

        integral = xMin**order * RosinRammlerPDF(xMin, x0, n)
        do i = 1, Nbin-1
            r = xMin + real(i) * h
            f = r**order * RosinRammlerPDF(r, x0, n)            
            if (mod(i, 2)==1) then; coeff = 4._R8; else; coeff = 2._R8; endif
            integral = integral + coeff*f
        enddo
        f = xMax**order * RosinRammlerPDF(xMax, x0, n)
        integral = (integral + f) * (h/3.0_R8)

    end function RonsinRammlerMoments

    !> Sample an injection diameter from the distribution law selected per BC inlet.
    !>   p1   = mean / characteristic diameter (= dp = 2*rp), always > 0
    !>   p2   = dispersion width (sigmap): std for Normal/LogNormal, spread n for Rosin
    !>   code = *Distr selector (stored in cell%properties(:,8))
    !> p2 <= 0 returns p1 verbatim with NO random_number call, so a deterministic
    !> (Dirac) inlet is bit-identical to the legacy d = 2*rp behaviour.
    function sampleDiameter(p1, p2, code) result(d)
        implicit none
        real(R8), intent(in) :: p1, p2
        integer,  intent(in) :: code
        integer, parameter   :: MAXTRY = 1000
        real(R8) :: d, sig2
        integer  :: i

        if (p2 <= 0._R8) then          ! Dirac / opt-out: no RNG consumed
            d = p1
            return
        endif

        select case (code)
        case (DiracDistr)
            d = p1
        case (NormalDistr)
            d = -1._R8                 ! truncate to physical (d > 0) by rejection
            do i = 1, MAXTRY
                d = Normal(p1, p2)
                if (d > 0._R8) exit
            enddo
            if (d <= 0._R8) error stop ( ' [ERROR] sampleDiameter: Normal truncation failed (mean << sigma?).' )
        case (LogNorDistr)
            ! p1,p2 are the statistical mean & std of the diameter; convert to the
            ! log-space (mu, sigma) the LogNormal sampler expects.
            if (p1 <= 0._R8) error stop ( ' [ERROR] sampleDiameter: LogNormal mean must be positive.' )
            sig2 = log(1._R8 + (p2/p1)**2)
            d = LogNormal(log(p1) - 0.5_R8*sig2, sqrt(sig2))
        case (RosRamDistr)
            d = RosinRammler(p1, p2)   ! p1 = x0, p2 = n
        case default
            error stop ( ' [ERROR] sampleDiameter: unknown distribution.' )
        end select

    end function sampleDiameter

    !> Map an ATLAS distribution-law name (col 8 of the BC file) to a DIST_* code.
    !> Empty / 'none' / 'dirac' -> deterministic. A file path or any other token is
    !> rejected (tabulated-file sampling is not implemented in IGLOO yet).
    function lawCode(name) result(code)
        implicit none
        character(len=*), intent(in) :: name
        integer :: code
        character(len=len(name)) :: nm
        nm = to_lower(adjustl(name))
        select case (trim(nm))
        case ('', 'none', 'dirac')
            code = DiracDistr
        case ('normal')
            code = NormalDistr
        case ('lognormal')
            code = LogNorDistr
        case ('rosinrammler')
            code = RosRamDistr
        case default
            write(*,*) ' [IGLOO::read_cdp_bc_file] unsupported distribution law: ', trim(name)
            error stop 1
        end select
    end function lawCode

    !> Seed the intrinsic RNG deterministically from a single integer seed, so a
    !> given seed reproduces the same stochastic injection across runs.
    subroutine initRandomSeed(seed)
        implicit none
        integer, intent(in) :: seed
        integer :: n, i
        integer, allocatable :: s(:)
        call random_seed(size=n)
        allocate(s(n))
        do i = 1, n
            s(i) = seed + 37*(i-1)
        enddo
        call random_seed(put=s)
        deallocate(s)
    end subroutine initRandomSeed

    pure function to_lower(s) result(t)
        implicit none
        character(len=*), intent(in) :: s
        character(len=len(s)) :: t
        integer :: i, ic
        t = s
        do i = 1, len(s)
            ic = iachar(s(i:i))
            if (ic >= iachar('A') .and. ic <= iachar('Z')) t(i:i) = achar(ic + 32)
        enddo
    end function to_lower

end module IGLOO_Lib_Statistics