program test_psd_table
    !
    ! Injection size distributions, sampled through the production sampleDiameter and tested against
    ! their prescribed CDFs, typed here and never read back from the solver.
    !
    !   T1  knots of table6 (six radii 2.5..30 um, bimodal mass fractions): the class edges (diameters)
    !       and the cumulative mass fraction at each edge against hand values, 1e-15
    !   T2  one entry per token: loading the same token again returns the same index
    !   T3  a token naming no file as written resolves to INPUT/<basename>
    !   T4  a one-row table is a Dirac at 2*r1 and draws nothing from the random stream
    !   T5  Kolmogorov-Smirnov and Anderson-Darling (case 0) of 1e5 draws, seed 42, for table6, Normal
    !       (20 um, 10 um, truncated to d > 0), LogNormal (mean 20 um, std 10 um) and Rosin-Rammler
    !       (x0 = 30 um, n = 2.5), alpha = 1e-3: D_n < 1.94947/sqrt(n), A^2 < 5.969. Seeds 43..61 are
    !       printed, not gated.
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64
    use IGLOO_Lib_Statistics, only: sampleDiameter, loadPSDTable, psdTab, initRandomSeed, &
                                    TableDistrBase, NormalDistr, LogNorDistr, RosRamDistr
    implicit none

    integer,  parameter :: N = 100000, NSEED = 20, SEED0 = 42, NLAW = 4
    real(R8), parameter :: KCRIT = 1.94947_R8, A2CRIT = 5.969_R8
    real(R8), parameter :: R6(6) = [2.5e-6_R8, 4.0e-6_R8, 6.0e-6_R8, 1.5e-5_R8, 2.25e-5_R8, 3.0e-5_R8]
    real(R8), parameter :: M6(6) = [0.05_R8, 0.10_R8, 0.15_R8, 0.35_R8, 0.25_R8, 0.10_R8]
    !> edges = 2 x (r1 - (r2-r1)/2, the mid-points, r6 + (r6-r5)/2); cdf = cumulative mass fraction
    real(R8), parameter :: E6(0:6) = [3.5e-6_R8, 6.5e-6_R8, 10.0e-6_R8, 21.0e-6_R8, 37.5e-6_R8, 52.5e-6_R8, &
                                      67.5e-6_R8]
    real(R8), parameter :: C6(0:6) = [0.0_R8, 0.05_R8, 0.15_R8, 0.30_R8, 0.65_R8, 0.90_R8, 1.0_R8]
    real(R8), parameter :: MU = 20.0e-6_R8, SIG = 10.0e-6_R8, X0 = 30.0e-6_R8, NRR = 2.5_R8
    character(len=*), parameter :: NAMES(NLAW) = [character(len=12) :: 'table6', 'Normal', 'LogNormal', &
                                                  'RosinRammler']
    real(R8), allocatable :: d(:)
    integer,  allocatable :: s1(:), s2(:)
    real(R8) :: Dn, A2, Dmax(NLAW), A2max(NLAW), dirac
    integer  :: nbad, idx6, idx6b, idxIn, idx1, law, iseed, nsz, nAbove(NLAW), u

    nbad = 0
    call execute_command_line('mkdir -p INPUT', wait=.true.)
    open(newunit=u, file='psdtest-table6.dat', status='replace', action='write')
    write(u,'(A)') '# particle radius [m]   mass fraction'
    write(u,'(A)') '2.5e-6    0.05', '4.0e-6    0.10', '6.0e-6    0.15', &
                   '1.5e-5    0.35', '2.25e-5   0.25', '3.0e-5    0.10'
    close(u)
    open(newunit=u, file='INPUT/psdtest-fallback.dat', status='replace', action='write')
    write(u,'(A)') '; a table found under INPUT/', '1.0e-5  1.0', '2.0e-5  1.0'
    close(u)
    open(newunit=u, file='psdtest-one-row.dat', status='replace', action='write')
    write(u,'(A)') '! one row: a Dirac', '1.2e-5   3.0'
    close(u)

    ! ---- T1: knots -----------------------------------------------------------------------------
    idx6 = loadPSDTable('psdtest-table6.dat')
    call check(psdTab(idx6)%K == 6, 'T1 table6 holds 6 classes', nbad)
    call check(all(abs(psdTab(idx6)%edge - E6) <= 1.0e-15_R8*E6), &
               'T1 class edges 3.5, 6.5, 10, 21, 37.5, 52.5, 67.5 um (diameter), 1e-15 relative', nbad)
    call check(all(abs(psdTab(idx6)%cdf - C6) <= 1.0e-15_R8), &
               'T1 cumulative mass fraction 0, .05, .15, .30, .65, .90, 1 at the edges, 1e-15', nbad)
    write(*,'(A,7ES24.16)') '        edges ', psdTab(idx6)%edge
    write(*,'(A,7ES24.16)') '        cdf   ', psdTab(idx6)%cdf

    ! ---- T2: one entry per token ------------------------------------------------------------------
    nsz = size(psdTab)
    idx6b = loadPSDTable('  psdtest-table6.dat')
    call check(idx6b == idx6 .and. size(psdTab) == nsz, 'T2 the same token loads once', nbad)

    ! ---- T3: INPUT/ + basename ---------------------------------------------------------------------
    idxIn = loadPSDTable('nowhere/psdtest-fallback.dat')
    call check(psdTab(idxIn)%file == 'INPUT/psdtest-fallback.dat' .and. psdTab(idxIn)%K == 2, &
               'T3 nowhere/psdtest-fallback.dat resolves to INPUT/psdtest-fallback.dat', nbad)

    ! ---- T4: one row = Dirac, no draw ---------------------------------------------------------------
    idx1 = loadPSDTable('psdtest-one-row.dat')
    call initRandomSeed(SEED0)
    call random_seed(size=nsz)
    allocate(s1(nsz), s2(nsz))
    call random_seed(get=s1)
    dirac = sampleDiameter(0._R8, 0._R8, TableDistrBase + idx1)
    call random_seed(get=s2)
    call check(dirac == 2._R8*1.2e-5_R8, 'T4 one row: d = 2 r1 = 24 um exactly', nbad)
    call check(all(s1 == s2), 'T4 one row: the random stream is not advanced', nbad)

    ! ---- T5: K-S and A^2 ----------------------------------------------------------------------------
    allocate(d(N))
    Dmax = 0._R8; A2max = 0._R8; nAbove = 0
    write(*,'(A,I0,A,ES10.3,A,F6.3)') '  T5 n = ', N, ': D_crit = ', KCRIT/sqrt(real(N,R8)), ', A2_crit = ', A2CRIT
    do law = 1, NLAW
        do iseed = SEED0, SEED0 + NSEED - 1
            call draw(law, iseed, d)
            call stats(law, d, Dn, A2)
            if (iseed == SEED0) then
                call check(Dn < KCRIT/sqrt(real(N,R8)), 'T5 '//trim(NAMES(law))//' K-S D_n < D_crit (seed 42)', nbad)
                call check(A2 < A2CRIT, 'T5 '//trim(NAMES(law))//' Anderson-Darling A^2 < 5.969 (seed 42)', nbad)
                write(*,'(A,ES11.4,A,ES11.4)') '        seed 42: D_n = ', Dn, '   A^2 = ', A2
            endif
            Dmax(law)  = max(Dmax(law), Dn)
            A2max(law) = max(A2max(law), A2)
            if (Dn >= KCRIT/sqrt(real(N,R8)) .or. A2 >= A2CRIT) nAbove(law) = nAbove(law) + 1
        enddo
        write(*,'(A,I0,A,ES11.4,A,ES11.4,A,I0)') '        seeds 42..', SEED0 + NSEED - 1, ': max D_n = ', Dmax(law), &
            '   max A^2 = ', A2max(law), '   above a threshold: ', nAbove(law)
    enddo

    write(*,*)
    if (nbad == 0) then
        write(*,'(A)') '[PASS] test_psd_table: the tabulated and parametric injection laws follow their CDFs'
    else
        write(*,'(A,I0,A)') '[FAIL] test_psd_table: ', nbad, ' check(s) failed'
        error stop 1
    endif

contains

    !> n diameters of law `law` through sampleDiameter, the intrinsic stream seeded with `seed`.
    subroutine draw(law, seed, x)
        integer,  intent(in)  :: law, seed
        real(R8), intent(out) :: x(:)
        integer :: i
        call initRandomSeed(seed)
        do i = 1, size(x)
            select case (law)
            case (1); x(i) = sampleDiameter(0._R8, 0._R8, TableDistrBase + idx6)
            case (2); x(i) = sampleDiameter(MU, SIG, NormalDistr)
            case (3); x(i) = sampleDiameter(MU, SIG, LogNorDistr)
            case (4); x(i) = sampleDiameter(X0, NRR, RosRamDistr)
            end select
        enddo
    end subroutine draw

    !> The prescribed CDF of law `law` at diameter x.
    pure function cdf(law, x) result(F)
        integer,  intent(in) :: law
        real(R8), intent(in) :: x
        real(R8) :: F, p0, s2
        integer  :: k
        select case (law)
        case (1)                                 ! histogram: linear between the hand-typed knots
            F = 0._R8
            if (x >= E6(6)) F = 1._R8
            do k = 1, 6
                if (x >= E6(k-1) .and. x < E6(k)) F = C6(k-1) + (x - E6(k-1))/(E6(k) - E6(k-1))*(C6(k) - C6(k-1))
            enddo
        case (2)                                 ! Normal truncated to d > 0
            p0 = phi(-MU/SIG)
            F = max(0._R8, (phi((x - MU)/SIG) - p0)/(1._R8 - p0))
        case (3)                                 ! LogNormal of mean MU and std SIG
            s2 = log(1._R8 + (SIG/MU)**2)
            F = phi((log(x) - (log(MU) - 0.5_R8*s2))/sqrt(s2))
        case default                             ! Rosin-Rammler
            F = 1._R8 - exp(-(x/X0)**NRR)
        end select
    end function cdf

    pure function phi(z) result(p)
        real(R8), intent(in) :: z
        real(R8) :: p
        p = 0.5_R8*erfc(-z/sqrt(2._R8))
    end function phi

    !> Kolmogorov-Smirnov D_n and Anderson-Darling A^2 (F clamped to [1e-15, 1-1e-15]) of x against cdf.
    subroutine stats(law, x, Dn, A2)
        integer,  intent(in)    :: law
        real(R8), intent(inout) :: x(:)
        real(R8), intent(out)   :: Dn, A2
        real(R8), allocatable   :: F(:)
        real(R8) :: rn, s
        integer  :: i, n
        n = size(x); rn = real(n, R8)
        call sort_r8(x)
        allocate(F(n))
        do i = 1, n
            F(i) = cdf(law, x(i))
        enddo
        Dn = 0._R8
        do i = 1, n
            Dn = max(Dn, real(i,R8)/rn - F(i), F(i) - real(i-1,R8)/rn)
        enddo
        F = min(max(F, 1.0e-15_R8), 1._R8 - 1.0e-15_R8)
        s = 0._R8
        do i = 1, n
            s = s + real(2*i-1,R8)*(log(F(i)) + log(1._R8 - F(n+1-i)))
        enddo
        A2 = -rn - s/rn
    end subroutine stats

    subroutine check(cond, name, nbad)
        logical,          intent(in)    :: cond
        character(len=*), intent(in)    :: name
        integer,          intent(inout) :: nbad
        if (cond) then
            write(*,'(a,a)') '  [PASS] ', name
        else
            write(*,'(a,a)') '  [FAIL] ', name
            nbad = nbad + 1
        end if
    end subroutine check

    !> Heapsort, ascending.
    subroutine sort_r8(a)
        real(R8), intent(inout) :: a(:)
        integer  :: n, s, r
        real(R8) :: t
        n = size(a)
        do s = n/2, 1, -1;  call sift(a, s, n);  enddo
        do r = n, 2, -1
            t = a(1); a(1) = a(r); a(r) = t
            call sift(a, 1, r-1)
        enddo
    end subroutine sort_r8

    subroutine sift(a, lo, hi)
        real(R8), intent(inout) :: a(:)
        integer,  intent(in)    :: lo, hi
        integer  :: root, child
        real(R8) :: t
        root = lo
        do
            child = 2*root
            if (child > hi) exit
            if (child < hi) then
                if (a(child) < a(child+1)) child = child + 1
            endif
            if (a(root) >= a(child)) exit
            t = a(root); a(root) = a(child); a(child) = t
            root = child
        enddo
    end subroutine sift

end program test_psd_table
