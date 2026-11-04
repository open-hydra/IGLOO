program test_fmt
    !
    ! IGLOO_fmt against the compiler runtime: each record of every layout the particle files write (trajectories,
    ! outloc, planes, snapshot, scatter) and of the zone headers' I0 is formatted by fmt_row and by the runtime's
    ! internal write with the same edit descriptors and the same values; the bytes must agree.
    !
    !   P1      fixed: the scatter layout on 2e6 plausible records with negatives that round to zero, -0.0 and zeros,
    !           and 1e6 decimal half-way points of F12.6 with their neighbours
    !   L1-L11  the eleven particle-record layouts, each in two parts:
    !           C  fixed: every field runs through its descriptor's constructed values -- +-0.0; NaN, -NaN, +-Inf;
    !              subnormals and the smallest normals; values that round to zero; exact ties; decimal half-way points;
    !              powers of ten and the values half a last digit below; the fast range's edges, values wider than
    !              the field, exponents beyond two digits (the points with their neighbours up to 8 ulps away, both
    !              signs); integers at every decade and at their width's edges
    !           R  NREC records, every field drawn at random: every decade of the double range, the decades around
    !              the field, the fast range, and the classes of C at random points
    !   L12     the zone headers' integers (I0) up to 2**53, both parts, through fmt_row and through fmt_int
    !   S       an integer column that holds no integer below 2**53 prints asterisks; fmt_int beyond the fast path
    !
    ! NREC (default 300000) is the first argument; P1 and the C parts do not depend on it. Mismatched fields are counted
    ! per descriptor, the fixed parts apart from the random ones; the throughput of both routes is measured on records
    ! of in-range values (one thread).
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64, I8 => int64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
    use IGLOO_fmt, only: fmtSpec, fmt_compile, fmt_row, fmt_int
    implicit none

    integer, parameter :: NL = 12, LREC = 512, NTPUT = 200000, ND = 10, NCLS = 9, NCMAX = 8192
    character(len=56), parameter :: LAY(NL) = [character(len=56) :: &
        '(7F12.6,2E13.6E2,I8)', &
        '(7F12.6,2E13.6E2,I8,ES16.8E2)', &
        '(7F12.6,2E13.6E2,I8,F10.6)', &
        '(7F12.6,2E13.6E2,I8,ES16.8E2,F10.6)', &
        '(7F12.6,2E13.6E2,I8,ES16.8E2,F12.6)', &
        '(7F12.6,2E13.6E2,I8,F12.6)', &
        '(6F12.6,2E13.6E2,I8,ES16.8E2,I5,E13.6E2,F10.6)', &
        '(6F12.6,2E13.6E2,I8,ES16.8E2,I5,E13.6E2)', &
        '(6F12.6,2E13.6E2,I8,I5,E13.6E2,F10.6)', &
        '(6F12.6,2E13.6E2,I8,I5,E13.6E2)', &
        '(I3,7F12.6,2E13.6E2,2ES16.8E2,F10.6,2I10,ES16.8E2,I3)', &
        '(I0,I0,I0,I0)']
    character(len=40), parameter :: WHAT(NL) = [character(len=40) :: &
        'trajectories, scatter', 'trajectories + t, snapshot', 'trajectories + f, scatter + f', &
        'trajectories + t + f, snapshot + f', 'trajectories P2T + t', 'trajectories P2T', 'outloc + t + f', &
        'outloc + t', 'outloc + f', 'outloc', 'planes', 'zone header integers']
    !> the descriptors, the surface-temperature column of the P2T layouts apart
    character(len=9), parameter :: DNAME(ND) = [character(len=9) :: 'F12.6', 'E13.6E2', 'I8', 'ES16.8E2', 'F10.6', &
        'F12.6 T_s', 'I3', 'I5', 'I10', 'I0']
    !> the constructed classes
    character(len=13), parameter :: CNAME(NCLS) = [character(len=13) :: '+-0', 'NaN,Inf', 'subnormal', &
        'round-to-0', 'exact-tie', 'half-way', 'power-of-10', 'edge,beyond', 'integer']

    integer :: nrec, l, k, nfail, nshow
    integer(I8) :: nfield(ND, 2), nbadf(ND, 2), nbadrec
    logical :: shown(ND)
    character(len=32) :: arg
    type(fmtSpec) :: spec

    nrec = 300000
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read(arg, *) nrec
    endif
    nfail = 0
    nfield = 0; nbadf = 0; shown = .false.

    call prototype_population(2000000)
    print '(a)', 'set   fixed: records mismatched   random: records mismatched   runtime-fields runtime-records  layout'
    do l = 1, NL
        call fmt_compile(spec, trim(LAY(l)))
        call layout_population(l, nrec)
    enddo
    call spot_checks()

    print '(a)', ''
    print '(a)', 'mismatched fields by descriptor   fixed parts (fields tested)   random parts (fields tested)'
    do k = 1, ND
        print '(2x,a9,i14,a,i0,a,i14,a,i0,a)', DNAME(k), nbadf(k, 1), '  (', nfield(k, 1), ')', nbadf(k, 2), '  (', &
            nfield(k, 2), ')'
        if (any(nbadf(k, :) > 0)) nfail = nfail + 1
    enddo
    call throughput()

    if (nfail == 0) then
        print '(a)', ''
        print '(a)', '[PASS] test_fmt: every record equals the runtime''s internal write'
    else
        print '(a)', ''
        print '(a,i0,a)', '[FAIL] test_fmt: ', nfail, ' assertion(s) failed'
        error stop 1
    endif

contains

    !> The runtime's record: the production writers' format and item list for layout l.
    subroutine oracle(l, v, r)
        integer, intent(in) :: l
        real(R8), intent(in) :: v(:)
        character(len=*), intent(out) :: r
        integer :: iv(size(v))
        integer(I8) :: jv(size(v))
        integer :: k
        do k = 1, size(v)
            iv(k) = 0; jv(k) = 0
            if (v(k) >= -2147483648._R8 .and. v(k) < 2147483648._R8) iv(k) = int(v(k))
            if (abs(v(k)) < 9007199254740992._R8) jv(k) = int(v(k), I8)
        enddo
        r = ' '
        select case (l)
        case (1);  write(r, '(7F12.6,2E13.6E2,I8)') v(1:9), iv(10)
        case (2);  write(r, '(7F12.6,2E13.6E2,I8,ES16.8E2)') v(1:9), iv(10), v(11)
        case (3);  write(r, '(7F12.6,2E13.6E2,I8,F10.6)') v(1:9), iv(10), v(11)
        case (4);  write(r, '(7F12.6,2E13.6E2,I8,ES16.8E2,F10.6)') v(1:9), iv(10), v(11), v(12)
        case (5);  write(r, '(7F12.6,2E13.6E2,I8,ES16.8E2,F12.6)') v(1:9), iv(10), v(11), v(12)
        case (6);  write(r, '(7F12.6,2E13.6E2,I8,F12.6)') v(1:9), iv(10), v(11)
        case (7);  write(r, '(6F12.6,2E13.6E2,I8,ES16.8E2,I5,E13.6E2,F10.6)') v(1:8), iv(9), v(10), iv(11), v(12), v(13)
        case (8);  write(r, '(6F12.6,2E13.6E2,I8,ES16.8E2,I5,E13.6E2)') v(1:8), iv(9), v(10), iv(11), v(12)
        case (9);  write(r, '(6F12.6,2E13.6E2,I8,I5,E13.6E2,F10.6)') v(1:8), iv(9), iv(10), v(11), v(12)
        case (10); write(r, '(6F12.6,2E13.6E2,I8,I5,E13.6E2)') v(1:8), iv(9), iv(10), v(11)
        case (11); write(r, '(I3,7F12.6,2E13.6E2,2ES16.8E2,F10.6,2I10,ES16.8E2,I3)') iv(1), v(2:13), iv(14), iv(15), &
                                                                                       v(16), iv(17)
        case (12); write(r, '(I0,I0,I0,I0)') jv(1:4)
        end select
    end subroutine oracle

    !> Index in DNAME of field k of layout l.
    integer function dindex(l, k)
        integer, intent(in) :: l, k
        character(len=16) :: f
        f = spec%rt(k)
        dindex = 0
        if ((l == 5 .and. k == 12) .or. (l == 6 .and. k == 11)) then
            dindex = 6
            return
        endif
        select case (trim(f))
        case ('(F12.6)');    dindex = 1
        case ('(E13.6E2)');  dindex = 2
        case ('(I8)');       dindex = 3
        case ('(ES16.8E2)'); dindex = 4
        case ('(F10.6)');    dindex = 5
        case ('(I3)');       dindex = 7
        case ('(I5)');       dindex = 8
        case ('(I10)');      dindex = 9
        case ('(I0)');       dindex = 10
        end select
    end function dindex

    !> Compare one record field by field (layouts of fixed width) or whole (I0) and count the mismatches of part
    !  (1 fixed, 2 random); show the first ones.
    subroutine compare(l, part, v, got, ngot, ref)
        integer, intent(in) :: l, part, ngot
        character(len=*), intent(in) :: got, ref
        real(R8), intent(in) :: v(:)
        integer :: k, o, w, di, nref
        logical :: bad
        nref = len_trim(ref)
        bad = ngot /= nref
        if (l == 12) then
            nfield(10, part) = nfield(10, part) + 4
            if (.not. bad) bad = got(1:ngot) /= ref(1:nref)
            if (bad) nbadf(10, part) = nbadf(10, part) + 1
        else
            o = 0
            do k = 1, spec%nf
                w = spec%w(k); di = dindex(l, k)
                nfield(di, part) = nfield(di, part) + 1
                if (got(o+1:o+w) /= ref(o+1:o+w)) then
                    nbadf(di, part) = nbadf(di, part) + 1
                    bad = .true.
                endif
                o = o + w
            enddo
        endif
        if (.not. bad) return
        nbadrec = nbadrec + 1
        if (nshow < 3) then
            nshow = nshow + 1
            print '(2a)', '  RTL ', ref(1:nref)
            print '(2a)', '  FMT ', got(1:max(ngot, 1))
            print '(a,*(es24.16e3))', '  values ', v
        endif
    end subroutine compare

    !> One record of layout l in part (1 fixed, 2 random), the runtime's bytes against fmt_row's; nsrec counts the
    !  records with a field the runtime wrote.
    subroutine check_record(l, part, v, ns, nsrec)
        integer, intent(in) :: l, part
        real(R8), intent(in) :: v(:)
        integer, intent(inout) :: ns
        integer(I8), intent(inout) :: nsrec
        character(len=LREC) :: got, ref
        integer :: ngot, ns0, f
        call oracle(l, v, ref)
        ngot = 0; ns0 = ns
        call fmt_row(spec, v, got, ngot, ns)
        if (ns > ns0) nsrec = nsrec + 1
        call compare(l, part, v, got, ngot, ref)
        if (l == 12) then
            !> the same integers through fmt_int
            ngot = 0
            do f = 1, 4
                call fmt_int(int(v(f), I8), 0, got, ngot)
            enddo
            call compare(l, part, v, got, ngot, ref)
        endif
    end subroutine check_record

    !> P1: the scatter layout on n records of plausible values (positions, velocities, temperature, diameter, mass, ID)
    !  and n/2 decimal half-way points of F12.6 with their neighbours.
    subroutine prototype_population(n)
        integer, intent(in) :: n
        integer(I8) :: s
        integer :: k, m, ii, idv, ns
        integer(I8) :: nsrec
        real(R8) :: w(10), x
        call fmt_compile(spec, trim(LAY(1)))
        nbadrec = 0; nshow = 0; ns = 0; nsrec = 0
        s = 12345_I8
        do k = 1, n
            call genrec(s, w(1:9), idv)
            if (mod(k, 7) == 0) w(4) = -w(4)*1.e-7_R8
            if (mod(k, 11) == 0) w(5) = real(nint(w(5)*128._R8), R8)/128._R8 + 1._R8/256._R8
            if (mod(k, 13) == 0) w(8) = real(100000 + mod(k, 900000), R8) + 0.5_R8
            if (mod(k, 17) == 0) w(1) = negzero()
            if (mod(k, 19) == 0) w(9) = 0._R8
            w(10) = real(idv, R8)
            call check_record(1, 1, w, ns, nsrec)
        enddo
        do k = 1, n/10
            x = (real(k, R8) + 0.5_R8)*1.e-6_R8
            do m = -2, 2
                w = 0.1_R8; w(1) = x
                do ii = 1, abs(m)
                    if (m > 0) w(1) = nearest(w(1), 1._R8)
                    if (m < 0) w(1) = nearest(w(1), -1._R8)
                enddo
                w(8) = x*1.e-3_R8; w(9) = -x; w(10) = real(k, R8)
                call check_record(1, 1, w, ns, nsrec)
            enddo
        enddo
        print '(a,i10,a,i9,a,i9,a,i9,2x,a)', 'P1    fixed records ', n + 5*(n/10), '  mismatched ', nbadrec, &
            '  runtime-fields ', ns, '  runtime-records ', nsrec, trim(LAY(1))//'  plausible records'
        if (nbadrec > 0) nfail = nfail + 1
    end subroutine prototype_population

    !> L<l>: the constructed records (C, fixed), then n records of random values (R).
    subroutine layout_population(l, n)
        integer, intent(in) :: l, n
        integer(I8) :: s
        integer :: k, f, ns, ncr, nc(spec%nf)
        integer(I8) :: nsrec, nbadc
        real(R8) :: v(spec%nf)
        real(R8), allocatable :: cv(:,:)
        integer :: cls(NCMAX)
        allocate(cv(NCMAX, spec%nf))
        do f = 1, spec%nf
            call constructed(f, l == 12, cv(:, f), cls, nc(f))
            call show_classes(dindex(l, f), cls, nc(f))
        enddo
        ncr = maxval(nc)
        nbadrec = 0; nshow = 0; ns = 0; nsrec = 0
        do k = 1, ncr
            do f = 1, spec%nf
                v(f) = cv(1 + mod(k - 1 + 37*(f - 1), nc(f)), f)
            enddo
            call check_record(l, 1, v, ns, nsrec)
        enddo
        nbadc = nbadrec
        nbadrec = 0; nshow = 0
        s = 88172645463325252_I8 + 7919_I8*l
        do k = 1, n
            do f = 1, spec%nf
                v(f) = draw(s, f, l == 12)
            enddo
            call check_record(l, 2, v, ns, nsrec)
        enddo
        print '(a,i2.2,i15,i11,i17,i11,i17,i16,2x,a)', 'L', l, ncr, nbadc, n, nbadrec, ns, nsrec, &
            trim(LAY(l))//'  '//trim(WHAT(l))
        if (nbadc + nbadrec > 0) nfail = nfail + 1
    end subroutine layout_population

    !> The constructed values for field f of the current layout (i0: the zone headers' integers), fixed whatever NREC,
    !  with their class (CNAME): both signs; ties, half-way points, powers of ten, the values below them and the edges
    !  with their neighbours up to 8 ulps away.
    subroutine constructed(f, i0, c, cls, n)
        integer, intent(in) :: f
        logical, intent(in) :: i0
        real(R8), intent(out) :: c(:)
        integer, intent(out) :: cls(:)
        integer, intent(out) :: n
        character(len=16) :: rt
        real(R8) :: lim, h, x, top, mm, nn(3)
        integer :: w, d, nd, lo, hi, q, i, k
        integer(I8), parameter :: SBIT = shiftl(1_I8, 63)
        integer(I8), parameter :: SUB(8) = [1_I8, 2_I8, 3_I8, 1000_I8, 1048576_I8, 1099511627776_I8, &
                                            2251799813685248_I8, 4503599627370495_I8]
        rt = spec%rt(f); w = spec%w(f); d = spec%d(f)
        n = 0
        if (rt(2:2) == 'I') then
            top = 2147483647._R8
            if (i0) top = 9007199254740991._R8
            !> every decade, so every width's edges: 10**k - 1, 10**k, 10**k + 1
            call put(0._R8, 9, c, cls, n)
            do k = 0, 16
                x = 10._R8**k
                if (x > top) exit
                call put2(x, 9, c, cls, n); call put2(x - 1._R8, 9, c, cls, n); call put2(x + 1._R8, 9, c, cls, n)
            enddo
            call put2(top, 9, c, cls, n); call put2(top - 1._R8, 9, c, cls, n)
            if (i0) then
                call put2(2._R8**31, 9, c, cls, n); call put2(2._R8**32, 9, c, cls, n)
                call put2(2._R8**52, 9, c, cls, n)
            else
                call put(-2147483648._R8, 9, c, cls, n)
            endif
            return
        endif
        !> zeros, NaN and infinities, subnormals and the smallest normals
        call put(0._R8, 1, c, cls, n); call put(negzero(), 1, c, cls, n)
        x = ieee_value(1._R8, ieee_quiet_nan)
        call put(x, 2, c, cls, n); call put(transfer(ior(transfer(x, 0_I8), SBIT), 0._R8), 2, c, cls, n)
        call put(ieee_value(1._R8, ieee_positive_inf), 2, c, cls, n)
        call put(-ieee_value(1._R8, ieee_positive_inf), 2, c, cls, n)
        do i = 1, size(SUB)
            call put(transfer(SUB(i), 0._R8), 3, c, cls, n); call put(transfer(ior(SUB(i), SBIT), 0._R8), 3, c, cls, n)
        enddo
        do k = 0, 2
            call put2(ulps(tiny(1._R8), k), 3, c, cls, n)
        enddo
        if (rt(2:2) == 'F') then
            lim = 10._R8**(w - d - 2) - 0.5_R8
            h = 0.5_R8*10._R8**(-d)
            !> rounding to zero: half the last digit and its neighbours, fractions of it, tiny values
            call nbr(h, 4, c, cls, n)
            do k = 1, 9
                call put2(h*real(k, R8)/10._R8, 4, c, cls, n)
            enddo
            do q = -300, -d - 1, 10
                call put2(10._R8**q, 4, c, cls, n)
            enddo
            !> exact ties: odd multiples of 2**-(d+1), from the smallest to the fast range's bound
            mm = 0._R8
            do while ((2._R8*mm + 1._R8)*2._R8**(-d - 1) < lim)
                call nbr((2._R8*mm + 1._R8)*2._R8**(-d - 1), 5, c, cls, n)
                mm = aint(1.5_R8*mm) + 1._R8
            enddo
            call nbr((2._R8*aint((lim*2._R8**(d + 1) - 1._R8)/2._R8) + 1._R8)*2._R8**(-d - 1), 5, c, cls, n)
            !> decimal half-way points of the last digit, up to the bound
            mm = 1._R8
            do while ((mm + 0.5_R8)*10._R8**(-d) < lim)
                call nbr((mm + 0.5_R8)*10._R8**(-d), 6, c, cls, n)
                mm = aint(1.7_R8*mm) + 1._R8
            enddo
            call nbr((aint(lim*10._R8**d - 0.5_R8) + 0.5_R8)*10._R8**(-d), 6, c, cls, n)
            !> powers of ten over the field and the values half a last digit below them
            do q = -d - 2, w - d
                call nbr(10._R8**q, 7, c, cls, n)
            enddo
            do q = -d + 1, w - d - 1
                call nbr(10._R8**q - h, 7, c, cls, n)
            enddo
            !> the fast range's bound, the largest value that fits, values wider than the field
            call nbr(lim, 8, c, cls, n); call nbr(10._R8**(w - d - 1) - h, 8, c, cls, n)
            do q = 10, 300, 10
                call put2(10._R8**q, 8, c, cls, n)
            enddo
            call put2(huge(1._R8), 8, c, cls, n)
        else
            nd = d; if (rt(2:3) == 'ES') nd = d + 1
            lo = nd - 23; hi = nd
            nn = [10._R8**(nd - 1), aint(3.7_R8*10._R8**(nd - 1)) + 7._R8, 10._R8**nd - 1._R8]
            !> exact ties: halves of nd-digit integers times exact powers of ten
            do i = 1, 3
                do q = 0, 22
                    x = (nn(i) + 0.5_R8)*10._R8**q
                    if (x >= 2._R8**53) exit
                    call nbr(x, 5, c, cls, n)
                enddo
            enddo
            !> decimal half-way points of the last digit at every fast decade and one beyond
            do i = 1, 3
                do q = lo - nd, hi - nd + 1
                    if (q >= 0) then
                        x = (nn(i) + 0.5_R8)*10._R8**q
                    else
                        x = (nn(i) + 0.5_R8)/10._R8**(-q)
                    endif
                    call nbr(x, 6, c, cls, n)
                enddo
            enddo
            !> powers of ten and the values half a last digit below them, the fast decades and two beyond
            do q = lo - 2, hi + 2
                call nbr(10._R8**q, 7, c, cls, n); call nbr(10._R8**q*(1._R8 - 0.5_R8*10._R8**(-nd)), 7, c, cls, n)
            enddo
            !> exponents at two digits and beyond
            do q = 97, 101
                call nbr(10._R8**q, 8, c, cls, n); call nbr(10._R8**(-q), 8, c, cls, n)
            enddo
            call nbr(10._R8**99*(1._R8 - 0.5_R8*10._R8**(-nd)), 8, c, cls, n)
            call put2(huge(1._R8), 8, c, cls, n)
        endif
    end subroutine constructed

    !> Append x of class k to the constructed list.
    subroutine put(x, k, c, cls, n)
        real(R8), intent(in) :: x
        integer, intent(in) :: k
        real(R8), intent(inout) :: c(:)
        integer, intent(inout) :: cls(:), n
        if (n >= size(c)) error stop 'test_fmt: constructed list too long'
        n = n + 1; c(n) = x; cls(n) = k
    end subroutine put

    !> Append x and -x.
    subroutine put2(x, k, c, cls, n)
        real(R8), intent(in) :: x
        integer, intent(in) :: k
        real(R8), intent(inout) :: c(:)
        integer, intent(inout) :: cls(:), n
        call put(x, k, c, cls, n); call put(-x, k, c, cls, n)
    end subroutine put2

    !> Append x and its neighbours up to 8 ulps away, both signs.
    subroutine nbr(x, k, c, cls, n)
        real(R8), intent(in) :: x
        integer, intent(in) :: k
        real(R8), intent(inout) :: c(:)
        integer, intent(inout) :: cls(:), n
        integer :: j
        do j = -8, 8
            call put2(ulps(x, j), k, c, cls, n)
        enddo
    end subroutine nbr

    !> The constructed values of a descriptor by class, the first time the descriptor is met.
    subroutine show_classes(di, cls, n)
        integer, intent(in) :: di, n
        integer, intent(in) :: cls(:)
        integer :: k, cnt(NCLS)
        if (shown(di)) return
        shown(di) = .true.
        cnt = 0
        do k = 1, n
            cnt(cls(k)) = cnt(cls(k)) + 1
        enddo
        print '(a,a9,a,i5,a,*(1x,a,1x,i0))', 'C ', DNAME(di), ' constructed ', n, ':', &
            (trim(CNAME(k)), cnt(k), k = 1, NCLS)
    end subroutine show_classes

    !> A value for field f of the current layout; zone integers (i0) range up to 2**53.
    real(R8) function draw(s, f, i0)
        integer(I8), intent(inout) :: s
        integer, intent(in) :: f
        logical, intent(in) :: i0
        character(len=16) :: rt
        integer :: w, d, nd, lo, hi
        rt = spec%rt(f); w = spec%w(f); d = spec%d(f)
        if (rt(2:2) == 'I') then
            draw = draw_int(s, w, i0)
        elseif (rt(2:3) == 'ES' .or. rt(2:2) == 'E') then
            nd = d; if (rt(2:3) == 'ES') nd = d + 1
            !> the fast decades: 10**lo <= |x| < 10**hi
            lo = nd - 23; hi = nd
            draw = draw_e(s, nd, lo, hi)
        else
            draw = draw_f(s, w, d)
        endif
    end function draw

    !> F(w,d) values; lim is the fast range's bound, h half the last digit.
    real(R8) function draw_f(s, w, d) result(x)
        integer(I8), intent(inout) :: s
        integer, intent(in) :: w, d
        real(R8) :: lim, h, u
        integer :: c, kk
        lim = 10._R8**(w - d - 2) - 0.5_R8
        h = 0.5_R8*10._R8**(-d)
        u = urand(s)
        c = int(urand(s)*100._R8)
        if (c < 6) then                                    !> every decade of the double range
            x = (1._R8 + 9._R8*u)*10._R8**(-307 + int(urand(s)*615._R8))
        elseif (c < 24) then                                !> the decades around the field
            x = (1._R8 + 9._R8*u)*10._R8**(-d - 6 + int(urand(s)*real(w + 6, R8)))
        elseif (c < 44) then                                !> uniform in the fast range
            x = u*lim
        elseif (c < 58) then                                !> decimal half-way points of the last digit
            x = (aint(10._R8**(urand(s)*real(w - 2, R8))) + 0.5_R8)*10._R8**(-d)
            x = ulps(x, int(urand(s)*17._R8) - 8)
        elseif (c < 66) then                                !> exact ties: odd multiples of 2**-(d+1)
            x = (2._R8*aint(u*lim*2._R8**d) + 1._R8)*2._R8**(-d - 1)
            x = ulps(x, int(urand(s)*17._R8) - 8)
        elseif (c < 74) then                                !> rounding to zero, the half quantum's neighbours
            if (urand(s) < 0.5_R8) then
                x = u*h
            else
                x = ulps(h, int(urand(s)*17._R8) - 8)
            endif
        elseif (c < 78) then
            x = 0._R8
        elseif (c < 90) then                                !> edges: the bound, powers of ten, values just below
            kk = int(urand(s)*3._R8)
            if (kk == 0) then
                x = lim
            elseif (kk == 1) then
                x = 10._R8**(-d + int(urand(s)*real(w - 1, R8)))
            else
                x = 10._R8**(-d + 1 + int(urand(s)*real(w - 2, R8))) - h
            endif
            x = ulps(x, int(urand(s)*17._R8) - 8)
        elseif (c < 93) then
            x = special(s)
        elseif (c < 96) then
            x = subnormal(s)
        else                                                !> wider than the field
            x = (1._R8 + 9._R8*u)*10._R8**(w - d - 3 + int(urand(s)*20._R8))
        endif
        if (urand(s) < 0.5_R8) x = -x
        if (x == 0._R8 .and. urand(s) < 0.5_R8) x = negzero()
    end function draw_f

    !> E/ES values with nd significant digits; the fast range is 10**lo <= |x| < 10**hi.
    real(R8) function draw_e(s, nd, lo, hi) result(x)
        integer(I8), intent(inout) :: s
        integer, intent(in) :: nd, lo, hi
        real(R8) :: u
        integer :: c, q, kk
        u = urand(s)
        c = int(urand(s)*100._R8)
        if (c < 6) then                                    !> every decade of the double range
            x = (1._R8 + 9._R8*u)*10._R8**(-307 + int(urand(s)*615._R8))
        elseif (c < 36) then                                !> the fast decades and two beyond
            x = (1._R8 + 9._R8*u)*10._R8**(lo - 2 + int(urand(s)*real(hi - lo + 4, R8)))
        elseif (c < 52) then                                !> decimal half-way points of the last digit
            q = lo - nd + int(urand(s)*real(hi - lo + 1, R8))
            x = aint(10._R8**(nd - 1) + u*9._R8*10._R8**(nd - 1)) + 0.5_R8
            if (q >= 0) then
                x = x*10._R8**q
            else
                x = x/10._R8**(-q)
            endif
            x = ulps(x, int(urand(s)*17._R8) - 8)
        elseif (c < 60) then                                !> small integers and halves (exact ties among them)
            x = aint(u*10._R8**(1 + int(urand(s)*real(nd + 2, R8))))
            if (urand(s) < 0.5_R8) x = x + 0.5_R8
        elseif (c < 64) then
            x = 0._R8
        elseif (c < 86) then                                !> powers of ten, the values just below, the range edges
            kk = int(urand(s)*3._R8)
            q = lo - 2 + int(urand(s)*real(hi - lo + 4, R8))
            if (kk == 0) then
                x = 10._R8**q
            elseif (kk == 1) then
                x = 10._R8**q*(1._R8 - 0.5_R8*10._R8**(-nd))
            else
                x = 10._R8**lo
                if (urand(s) < 0.5_R8) x = 10._R8**hi*(1._R8 - 0.5_R8*10._R8**(-nd))
            endif
            x = ulps(x, int(urand(s)*17._R8) - 8)
        elseif (c < 89) then
            x = special(s)
        elseif (c < 93) then
            x = subnormal(s)
        else                                                !> exponents beyond two digits
            x = (1._R8 + 9._R8*u)*10._R8**(100 + int(urand(s)*208._R8))
            if (urand(s) < 0.5_R8) x = (1._R8 + 9._R8*u)*10._R8**(-307 + int(urand(s)*207._R8))
        endif
        if (urand(s) < 0.5_R8) x = -x
        if (x == 0._R8 .and. urand(s) < 0.5_R8) x = negzero()
    end function draw_e

    !> Integers for I(w): every decade, the width's edges, the int32 extremes; i0: up to 2**53.
    real(R8) function draw_int(s, w, i0) result(x)
        integer(I8), intent(inout) :: s
        integer, intent(in) :: w
        logical, intent(in) :: i0
        real(R8) :: top
        integer :: c, wl
        top = 2147483647._R8
        if (i0) top = 9007199254740991._R8
        wl = w; if (w == 0) wl = 16
        c = int(urand(s)*10._R8)
        if (c < 5) then
            x = aint(min(10._R8**(urand(s)*log10(top)), top))
        elseif (c < 8) then
            select case (int(urand(s)*6._R8))
            case (0); x = 10._R8**(wl - 1) - 1._R8
            case (1); x = 10._R8**wl - 1._R8
            case (2); x = 10._R8**wl
            case (3); x = 10._R8**(wl - 1)
            case (4); x = 0._R8
            case default; x = 1._R8
            end select
            x = min(x, top)
        elseif (c < 9) then
            x = top
        else
            x = aint(urand(s)*1000._R8)
        endif
        if (urand(s) < 0.5_R8) x = -x
        if (x == -top .and. .not. i0 .and. urand(s) < 0.5_R8) x = -2147483648._R8
    end function draw_int

    !> NaN, +Inf or -Inf.
    real(R8) function special(s)
        integer(I8), intent(inout) :: s
        select case (int(urand(s)*3._R8))
        case (0);     special = ieee_value(1._R8, ieee_quiet_nan)
        case (1);     special = ieee_value(1._R8, ieee_positive_inf)
        case default; special = -ieee_value(1._R8, ieee_positive_inf)
        end select
    end function special

    !> A positive subnormal.
    real(R8) function subnormal(s)
        integer(I8), intent(inout) :: s
        subnormal = transfer(max(1_I8, iand(xorshift(s), 4503599627370495_I8)), 0._R8)
    end function subnormal

    !> -0.0, built from its bits.
    real(R8) function negzero()
        negzero = transfer(shiftl(1_I8, 63), 0._R8)
    end function negzero

    !> x moved by k ulps away from zero (k < 0: towards it), its sign kept.
    real(R8) function ulps(x, k)
        real(R8), intent(in) :: x
        integer, intent(in) :: k
        integer(I8) :: b
        b = transfer(abs(x), 0_I8) + int(k, I8)
        b = max(0_I8, min(b, 9218868437227405311_I8))
        ulps = transfer(b, 0._R8)
        if (x < 0._R8) ulps = -ulps
    end function ulps

    !> Marsaglia's xorshift64 (shifts only).
    integer(I8) function xorshift(s)
        integer(I8), intent(inout) :: s
        s = ieor(s, shiftl(s, 13)); s = ieor(s, shiftr(s, 7)); s = ieor(s, shiftl(s, 17))
        xorshift = s
    end function xorshift

    !> Uniform in [0, 1).
    real(R8) function urand(s)
        integer(I8), intent(inout) :: s
        urand = real(shiftr(xorshift(s), 11), R8)*2._R8**(-53)
    end function urand

    !> xorshift64* deviate in [0, 1) (P1's stream).
    real(R8) function urand_star(s)
        integer(I8), intent(inout) :: s
        s = ieor(s, ishft(s, -12)); s = ieor(s, ishft(s, 25)); s = ieor(s, ishft(s, -27))
        urand_star = real(ishft(s*2685821657736338717_I8, -11), R8)*2._R8**(-53)
    end function urand_star

    !> P1's plausible record: positions, velocities, temperature, diameter, mass, ID.
    subroutine genrec(s, v, id)
        integer(I8), intent(inout) :: s
        real(R8), intent(out) :: v(9)
        integer, intent(out) :: id
        v(1) = 0.15_R8*urand_star(s); v(2) = 0.05_R8*urand_star(s); v(3) = 0.05_R8*urand_star(s)
        v(4) = 300._R8*urand_star(s) - 50._R8; v(5) = 20._R8*urand_star(s) - 10._R8
        v(6) = 20._R8*urand_star(s) - 10._R8
        v(7) = 280._R8 + 40._R8*urand_star(s)
        v(8) = 10._R8**(-6._R8 + 2.6_R8*urand_star(s)); v(9) = 10._R8**(-16._R8 + 8.5_R8*urand_star(s))
        id = 1 + int(99999998._R8*urand_star(s))
    end subroutine genrec

    !> S: asterisks for an integer column without an integer below 2**53; fmt_int's runtime write; the layouts' widths.
    subroutine spot_checks()
        type(fmtSpec) :: sp
        character(len=64) :: got, ref
        real(R8) :: bad(4)
        integer :: k, n, ns
        integer, parameter :: WREC(NL) = [118, 134, 128, 144, 146, 130, 150, 140, 134, 124, 194, 80]
        bad = [ieee_value(1._R8, ieee_quiet_nan), ieee_value(1._R8, ieee_positive_inf), 9007199254740992._R8, &
               -1.e300_R8]
        call fmt_compile(sp, '(I8)')
        do k = 1, 4
            n = 0
            call fmt_row(sp, bad(k:k), got, n)
            if (n /= 8 .or. got(1:8) /= '********') then
                print '(a,es12.4,2a)', '  S1 FAIL: I8 of ', bad(k), ' gives ', got(1:max(n, 1))
                nfail = nfail + 1
            endif
        enddo
        n = 0; ns = 0
        call fmt_int(-huge(1_I8) - 1_I8, 0, got, n, ns)
        write(ref, '(I0)') -huge(1_I8) - 1_I8
        if (got(1:n) /= trim(ref) .or. ns /= 1) then
            print '(4a)', '  S2 FAIL: I0 of -2**63 gives ', got(1:max(n, 1)), ' runtime ', trim(ref)
            nfail = nfail + 1
        endif
        n = 0; ns = 0
        call fmt_int(123456789_I8, 8, got, n, ns)
        write(ref, '(I8)') 123456789_I8
        if (got(1:n) /= ref(1:8) .or. ns /= 1) then
            print '(4a)', '  S3 FAIL: I8 of 123456789 gives ', got(1:max(n, 1)), ' runtime ', ref(1:8)
            nfail = nfail + 1
        endif
        do k = 1, NL
            call fmt_compile(sp, trim(LAY(k)))
            if (k < NL .and. sp%maxlen /= WREC(k)) then
                print '(a,i0,a,i0,a,i0)', '  S4 FAIL: layout ', k, ' width ', sp%maxlen, ' expected ', WREC(k)
                nfail = nfail + 1
            endif
        enddo
        print '(a)', 'S     spot checks done'
    end subroutine spot_checks

    !> Records/s of both routes on records of in-range values, one thread.
    subroutine throughput()
        real(R8), allocatable :: v(:,:)
        character(len=LREC) :: got
        integer(I8) :: s, t0, t1, rate
        integer :: l, k, f, ngot, ns
        real(R8) :: tf, tr, x, chk
        character(len=16) :: rt
        print '(a)', ''
        print '(a,i0,a)', 'throughput on ', NTPUT, ' records of in-range values [records/s]: fmt_row  runtime  ratio  '// &
            'runtime-fields'
        do l = 1, NL
            call fmt_compile(spec, trim(LAY(l)))
            allocate(v(spec%nf, NTPUT))
            s = 1234567_I8 + l
            do k = 1, NTPUT
                do f = 1, spec%nf
                    rt = spec%rt(f)
                    x = urand(s)
                    if (rt(2:2) == 'I') then
                        v(f, k) = aint(1._R8 + x*min(10._R8**max(spec%w(f) - 1, 1), 99999._R8))
                    elseif (rt(2:3) == 'ES') then
                        v(f, k) = 10._R8**(-9._R8 + 12._R8*x)
                    elseif (rt(2:2) == 'E') then
                        v(f, k) = 10._R8**(-16._R8 + 13._R8*x)
                    elseif (spec%w(f) == 10) then
                        v(f, k) = x
                    else
                        v(f, k) = 1100._R8*x - 100._R8
                    endif
                enddo
            enddo
            ns = 0; chk = 0._R8
            call system_clock(t0, rate)
            do k = 1, NTPUT
                ngot = 0
                call fmt_row(spec, v(:, k), got, ngot, ns)
                chk = chk + iachar(got(ngot:ngot))
            enddo
            call system_clock(t1)
            tf = real(t1 - t0, R8)/real(rate, R8)
            call system_clock(t0)
            do k = 1, NTPUT
                call oracle(l, v(:, k), got)
                chk = chk + iachar(got(1:1))
            enddo
            call system_clock(t1)
            tr = real(t1 - t0, R8)/real(rate, R8)
            print '(a,i2.2,2es11.3,f7.2,i10,2x,a)', 'L', l, NTPUT/tf, NTPUT/tr, tr/tf, ns, trim(LAY(l))
            if (chk < 0._R8) print *, chk
            deallocate(v)
        enddo
    end subroutine throughput

end program test_fmt
