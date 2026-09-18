!> Random samplers (intrinsic and per-particle streamed), Rosin-Rammler moments and
!  injection-diameter sampling.
module IGLOO_Lib_Statistics
    use, intrinsic :: iso_fortran_env, only : I8 => int64, R8 => real64
    implicit none
    private
    !> Bare samplers draw from the intrinsic random_number (serial pin-time path only); the *S
    !  twins take a per-particle stream state and are the only ones callable inside the OMP region.
    public :: RosinRammler, LogNormal, Normal, NormalStandard, ChiSquare
    public :: RosinRammlerS, NormalStandardS, ChiSquareS
    public :: rngSeedFor, rngNext

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

    !> splitmix64 constants as signed decimal (two's-complement images of the published hex values).
    integer(I8), parameter :: SM_GOLDEN = -7046029254386353131_I8
    integer(I8), parameter :: SM_IDMIX  = -4417276706812531889_I8
    integer(I8), parameter :: SM_MIX_A  = -4658895280553007687_I8
    integer(I8), parameter :: SM_MIX_B  = -7723592293110705685_I8

contains

    !> Per-particle deterministic RNG stream (splitmix64, Steele-Lea-Flood 2014): seeded from
    !  (rng_seed, famID, ID) and advanced only by that particle's own draws.

    !> Deterministic stream seed for one particle, mixing famID and the within-group ID.
    pure function rngSeedFor(famID, ID) result(state)
        use IGLOO_variables, only: rng_seed
        implicit none
        integer, intent(in) :: famID, ID
        integer(I8) :: state

        state = mix64( int(rng_seed,I8)                                &
                     + int(famID,I8) * SM_GOLDEN           &
                     + int(ID,   I8) * SM_IDMIX )

    end function rngSeedFor

    !> Next uniform in (0,1) from a particle's own stream; advances `state`.
    function rngNext(state) result(u)
        implicit none
        integer(I8), intent(inout) :: state
        real(R8) :: u
        integer(I8) :: z
        real(R8), parameter :: twoM53 = 1._R8/9007199254740992._R8   ! 2**-53

        do
            state = state + SM_GOLDEN      ! splitmix64 Weyl increment
            z     = mix64(state)
            !> Top 53 bits -> [0,1), taken with a logical shift.
            u = real(ibits(z, 11, 53), R8) * twoM53
            if (u > 0._R8) exit
        enddo

    end function rngNext

    !> splitmix64 finalizer; the multiplies wrap mod 2**64 by design.
    pure function mix64(x) result(z)
        implicit none
        integer(I8), intent(in) :: x
        integer(I8) :: z

        z = x
        z = ieor(z, ishft(z, -30)) * SM_MIX_A
        z = ieor(z, ishft(z, -27)) * SM_MIX_B
        z = ieor(z, ishft(z, -31))

    end function mix64

    !> Streamed standard normal N(0,1), Marsaglia polar method.
    function NormalStandardS(state) result(z)
        implicit none
        integer(I8), intent(inout) :: state
        real(R8) :: z, u1, u2, s, y

        do
            u1 = 2._R8*rngNext(state) - 1._R8
            u2 = 2._R8*rngNext(state) - 1._R8
            s  = u1*u1 + u2*u2
            if (s > 0._R8 .and. s < 1._R8) exit
        enddo
        y = sqrt(-2._R8*log(s)/s)
        z = u1*y                    !> second deviate discarded; see NormalStandard

    end function NormalStandardS

    !> Streamed Rosin-Rammler sample with size parameter x0 and spread parameter n.
    function RosinRammlerS(x0, n, state) result(x)
        implicit none
        real(R8),    intent(in)    :: x0, n
        integer(I8), intent(inout) :: state
        real(R8) :: x

        if (x0<=0._R8 .or. n<=0._R8) &
            error stop ( ' [ERROR] Rosin-Rammler distribution: x0 and n must be positive.' )
        x = x0*(-log(rngNext(state)))**(1._R8/n)

    end function RosinRammlerS

    !> Streamed chi-squared sample with k degrees of freedom.
    function ChiSquareS(k, state) result(x)
        implicit none
        integer,     intent(in)    :: k
        integer(I8), intent(inout) :: state
        real(R8) :: x, z
        integer  :: i

        if (k < 1) error stop ( ' [ERROR] ChiSquare: degrees of freedom must be >= 1.' )
        x = 0._R8
        do i = 1, k
            z = NormalStandardS(state)
            x = x + z*z
        enddo

    end function ChiSquareS


    !> Standard normal N(0,1), Marsaglia polar method; stateless (the second deviate is discarded).
    function NormalStandard() result(z)
        implicit none
        real(R8) :: u1, u2, v1, v2, w, y
        real(R8) :: z

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

    end function NormalStandard

    !> Normal distribution N(mean, sigma)
    function Normal(mean, sigma) result(x)
        implicit none
        real(R8), intent(in) :: mean, sigma
        real(R8) :: x
        x = mean + sigma*NormalStandard()
    end function Normal

    !> Rosin-Rammler distribution with size parameter x0 and spread parameter n
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

    !> Rosin-Rammler probability density at r.
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

    !> Samples an injection diameter from the per-inlet distribution law `code`: p1 = mean /
    !  characteristic diameter, p2 = width (std, or Rosin-Rammler n); p2 <= 0 returns p1 with no draw.
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
            ! convert the diameter mean/std to log-space (mu, sigma)
            if (p1 <= 0._R8) error stop ( ' [ERROR] sampleDiameter: LogNormal mean must be positive.' )
            sig2 = log(1._R8 + (p2/p1)**2)
            d = LogNormal(log(p1) - 0.5_R8*sig2, sqrt(sig2))
        case (RosRamDistr)
            d = RosinRammler(p1, p2)   ! p1 = x0, p2 = n
        case default
            error stop ( ' [ERROR] sampleDiameter: unknown distribution.' )
        end select

    end function sampleDiameter

    !> Maps a distribution-law name (empty/none/dirac, normal, lognormal, rosinrammler) to its
    !  *Distr code; any other token is rejected.
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

    !> Seeds the intrinsic RNG deterministically from one integer.
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

    !> ASCII lower-case copy of s.
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