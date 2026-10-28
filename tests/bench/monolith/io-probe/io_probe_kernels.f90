!> Kernels of the I/O probe (compiled with -fp-model precise under Intel): the data regions and tokens of an ASCII
!  Tecplot file as ORION's reader sees them (blank, tab, comma and line end separate the values), ORION's per-token
!  list-directed READ, one list-directed READ per run of plain tokens, the exact fast conversion (Clinger) with the
!  runtime for every other token, glibc's strtod, and E23.15E3 by the runtime and by a lock-free formatter.
module io_probe_kern
  use, intrinsic :: iso_fortran_env, only: R8 => real64, I8 => int64
  use, intrinsic :: iso_c_binding, only: c_char, c_double, c_int, c_int64_t
  implicit none
  private
  public :: R8, I8, LF
  public :: is_sep, numeric_start, zone_line, find_regions, count_tokens, fill_tokens, next_line_end
  public :: plain_token, fast_token, fastg_token, x87_token, rtl_token, strtod_token, batch_read
  public :: fmt_e23_fast, fmt_e23_rtl
  public :: probe_fma, probe_x87, probe_strtod, probe_times, probe_cyc_open, probe_cyc_read, probe_cpu, probe_evict

  interface
    real(c_double) function probe_fma(a, b, c) bind(C, name='probe_fma')
      import :: c_double
      real(c_double), value :: a, b, c
    end function probe_fma
    real(c_double) function probe_x87(m, e) bind(C, name='probe_x87')
      import :: c_double, c_int64_t, c_int
      integer(c_int64_t), value :: m
      integer(c_int), value :: e
    end function probe_x87
    real(c_double) function probe_strtod(s, pos, endpos) bind(C, name='probe_strtod')
      import :: c_char, c_double, c_int64_t
      character(kind=c_char), intent(in) :: s(*)
      integer(c_int64_t), value :: pos
      integer(c_int64_t), intent(out) :: endpos
    end function probe_strtod
    subroutine probe_times(user, sys) bind(C, name='probe_times')
      import :: c_double
      real(c_double), intent(out) :: user, sys
    end subroutine probe_times
    integer(c_int) function probe_cyc_open() bind(C, name='probe_cyc_open')
      import :: c_int
    end function probe_cyc_open
    integer(c_int64_t) function probe_cyc_read() bind(C, name='probe_cyc_read')
      import :: c_int64_t
    end function probe_cyc_read
    integer(c_int) function probe_cpu() bind(C, name='probe_cpu')
      import :: c_int
    end function probe_cpu
    integer(c_int) function probe_evict(path) bind(C, name='probe_evict')
      import :: c_char, c_int
      character(kind=c_char), intent(in) :: path(*)
    end function probe_evict
  end interface

  character, parameter :: LF = achar(10)
  real(R8), parameter :: P10(0:22) = [1.e0_R8, 1.e1_R8, 1.e2_R8, 1.e3_R8, 1.e4_R8, 1.e5_R8, 1.e6_R8, 1.e7_R8, &
    1.e8_R8, 1.e9_R8, 1.e10_R8, 1.e11_R8, 1.e12_R8, 1.e13_R8, 1.e14_R8, 1.e15_R8, 1.e16_R8, 1.e17_R8, 1.e18_R8, &
    1.e19_R8, 1.e20_R8, 1.e21_R8, 1.e22_R8]
  integer(I8), parameter :: I10(0:18) = [1_I8, 10_I8, 100_I8, 1000_I8, 10000_I8, 100000_I8, 1000000_I8, &
    10000000_I8, 100000000_I8, 1000000000_I8, 10000000000_I8, 100000000000_I8, 1000000000000_I8, &
    10000000000000_I8, 100000000000000_I8, 1000000000000000_I8, 10000000000000000_I8, 100000000000000000_I8, &
    1000000000000000000_I8]
  integer(I8), parameter :: TWO53 = 9007199254740992_I8

contains

  !> A separator of ORION's tokenizer (blank, tab, comma) or a line end, by character code (a comparison of
  !  one-character strings costs a runtime call per byte under ifx).
  pure logical function is_sep(c)
    integer, intent(in) :: c
    is_sep = (c == 32 .or. c == 9 .or. c == 44 .or. c == 10)
  end function is_sep

  !> Position of the next LF at or after p (n + 1 when none).
  pure integer(I8) function next_line_end(buf, p, n) result(q)
    character(len=*), intent(in) :: buf
    integer(I8), intent(in) :: p, n
    q = p
    do while (q <= n)
      if (iachar(buf(q:q)) == 10) return
      q = q + 1
    enddo
  end function next_line_end

  !> The line's first token starts like a number (its first character is a digit, a sign or a point).
  pure logical function numeric_start(line)
    character(len=*), intent(in) :: line
    integer :: i
    integer :: c
    numeric_start = .false.
    do i = 1, len(line)
      c = iachar(line(i:i))
      if (c == 32 .or. c == 9 .or. c == 44) cycle
      numeric_start = ((c >= 48 .and. c <= 57) .or. c == 43 .or. c == 45 .or. c == 46)
      return
    enddo
  end function numeric_start

  !> ORION's zone record: the keyword ZONE, in any case, as the first word of the line.
  pure logical function zone_line(line)
    character(len=*), intent(in) :: line
    integer :: i, c, k
    integer, parameter :: Z(4) = [90, 79, 78, 69]
    zone_line = .false.
    do i = 1, len(line)
      if (iachar(line(i:i)) /= 32) exit
    enddo
    if (i + 3 > len(line)) return
    do c = 0, 3
      k = iachar(line(i+c:i+c))
      if (k >= 97 .and. k <= 122) k = k - 32
      if (k /= Z(c+1)) return
    enddo
    zone_line = .true.
    if (i + 4 <= len(line)) then
      k = iachar(line(i+4:i+4))
      zone_line = .not. ((k >= 65 .and. k <= 90) .or. (k >= 97 .and. k <= 122) .or. (k >= 48 .and. k <= 57) .or. k == 95)
    endif
  end function zone_line

  !> The data regions of buf(1:n): the runs of lines that start like a number (blank lines inside a run kept);
  !  a zone record or any other line ends a run (the header lines of a zone start otherwise). Byte ranges [rs, re].
  subroutine find_regions(buf, n, rs, re, nreg, nzone)
    character(len=*), intent(in) :: buf
    integer(I8), intent(in) :: n
    integer(I8), allocatable, intent(out) :: rs(:), re(:)
    integer, intent(out) :: nreg, nzone
    integer(I8) :: p, q
    logical :: inrun
    integer(I8), allocatable :: ts(:), te(:)
    allocate(rs(16), re(16))
    nreg = 0; nzone = 0; inrun = .false.
    p = 1
    do while (p <= n)
      q = next_line_end(buf, p, n)
      if (zone_line(buf(p:q-1))) then
        inrun = .false.; nzone = nzone + 1
      elseif (numeric_start(buf(p:q-1))) then
        if (.not. inrun) then
          if (nreg == size(rs)) then
            allocate(ts(2*nreg), te(2*nreg)); ts(1:nreg) = rs; te(1:nreg) = re
            call move_alloc(ts, rs); call move_alloc(te, re)
          endif
          nreg = nreg + 1; rs(nreg) = p; inrun = .true.
        endif
        re(nreg) = q - 1
      elseif (len_trim(buf(p:q-1)) > 0) then
        if (inrun) inrun = .false.
      endif
      p = q + 1
    enddo
  end subroutine find_regions

  !> Tokens in buf(a:b).
  pure integer(I8) function count_tokens(buf, a, b) result(nt)
    character(len=*), intent(in) :: buf
    integer(I8), intent(in) :: a, b
    integer(I8) :: p
    logical :: sep, prev
    nt = 0; prev = .true.
    do p = a, b
      sep = is_sep(iachar(buf(p:p)))
      if (prev .and. .not. sep) nt = nt + 1
      prev = sep
    enddo
  end function count_tokens

  !> The tokens of buf(a:b) into ts(k0+1:), tl(k0+1:) (start, length); k0 advances; longest token length.
  subroutine fill_tokens(buf, a, b, ts, tl, k0, longest)
    character(len=*), intent(in) :: buf
    integer(I8), intent(in) :: a, b
    integer(I8), intent(inout) :: ts(:)
    integer, intent(inout) :: tl(:)
    integer(I8), intent(inout) :: k0
    integer, intent(inout) :: longest
    integer(I8) :: p
    logical :: intok
    intok = .false.
    do p = a, b
      if (is_sep(iachar(buf(p:p)))) then
        if (intok) then
          tl(k0) = int(p - ts(k0)); longest = max(longest, tl(k0)); intok = .false.
        endif
      elseif (.not. intok) then
        intok = .true.; k0 = k0 + 1; ts(k0) = p
      endif
    enddo
    if (intok) then
      tl(k0) = int(b + 1 - ts(k0)); longest = max(longest, tl(k0))
    endif
  end subroutine fill_tokens

  !> A plain number token: [+-]?(d+[.d*]|.d+)([eEdD][+-]?d+)? -- the tokens one list-directed READ of many values
  !  converts as it converts each alone (no repeat count, slash, null value or other form).
  pure logical function plain_token(s)
    character(len=*), intent(in) :: s
    integer :: i, n, nd, c
    logical :: dot
    plain_token = .false.
    n = len(s); i = 1
    if (n == 0) return
    c = iachar(s(1:1))
    if (c == 43 .or. c == 45) i = 2
    nd = 0; dot = .false.
    do while (i <= n)
      c = iachar(s(i:i))
      if (c >= 48 .and. c <= 57) then
        nd = nd + 1
      elseif (c == 46 .and. .not. dot) then
        dot = .true.
      else
        exit
      endif
      i = i + 1
    enddo
    if (nd == 0) return
    if (i <= n) then
      if (.not. is_exp(iachar(s(i:i)))) return
      i = i + 1
      if (i <= n) then
        c = iachar(s(i:i))
        if (c == 43 .or. c == 45) i = i + 1
      endif
      if (i > n) return
      do while (i <= n)
        c = iachar(s(i:i))
        if (c < 48 .or. c > 57) return
        i = i + 1
      enddo
    endif
    plain_token = .true.
  end function plain_token

  !> An exponent letter: e, E, d or D.
  pure logical function is_exp(c)
    integer, intent(in) :: c
    is_exp = (c == 101 .or. c == 69 .or. c == 100 .or. c == 68)
  end function is_exp

  !> A plain token as m * 10**ex with the trailing zeros of its digits in ex: ok = .false. for any other token,
  !  for more than 18 significant digits and for a minus sign on a zero mantissa.
  pure subroutine parse_token(s, m, ex, neg, ok)
    character(len=*), intent(in) :: s
    integer(I8), intent(out) :: m
    integer, intent(out) :: ex
    logical, intent(out) :: neg, ok
    integer :: i, n, nd, nz, frac, d, c
    logical :: dot, dig, eneg
    ok = .false.; m = 0_I8; ex = 0; neg = .false.
    n = len(s); i = 1
    if (n == 0) return
    c = iachar(s(1:1))
    if (c == 45) then
      neg = .true.; i = 2
    elseif (c == 43) then
      i = 2
    endif
    nd = 0; nz = 0; frac = 0; dot = .false.; dig = .false.
    do while (i <= n)
      c = iachar(s(i:i))
      d = c - 48
      if (d >= 0 .and. d <= 9) then
        dig = .true.
        if (dot) frac = frac + 1
        if (d == 0) then
          if (nd > 0) nz = nz + 1
        else
          if (nd + nz + 1 > 18) return
          m = m*I10(nz+1) + d
          nd = nd + nz + 1; nz = 0
        endif
      elseif (c == 46 .and. .not. dot) then
        dot = .true.
      else
        exit
      endif
      i = i + 1
    enddo
    if (.not. dig) return
    if (i <= n) then
      if (.not. is_exp(iachar(s(i:i)))) return
      i = i + 1
      eneg = .false.
      if (i <= n) then
        c = iachar(s(i:i))
        if (c == 45) then
          eneg = .true.; i = i + 1
        elseif (c == 43) then
          i = i + 1
        endif
      endif
      if (i > n) return
      do while (i <= n)
        d = iachar(s(i:i)) - 48
        if (d < 0 .or. d > 9) return
        if (ex < 100000) ex = ex*10 + d
        i = i + 1
      enddo
      if (eneg) ex = -ex
    endif
    if (m == 0_I8 .and. neg) return
    ex = ex + nz - frac
    ok = .true.
  end subroutine parse_token

  !> The exact fast path (Clinger): m <= 2**53 and |ex| <= 22, one IEEE multiply or divide of exact operands
  !  (correctly rounded); ok = .false. for every other token.
  pure subroutine fast_token(s, x, ok)
    character(len=*), intent(in) :: s
    real(R8), intent(out) :: x
    logical, intent(out) :: ok
    integer(I8) :: m
    integer :: ex
    logical :: neg
    x = 0._R8
    call parse_token(s, m, ex, neg, ok)
    if (.not. ok) return
    if (m == 0_I8) return
    ok = .false.
    if (m > TWO53 .or. abs(ex) > 22) return
    x = real(m, R8)
    if (ex >= 0) then
      x = x*P10(ex)
    else
      x = x/P10(-ex)
    endif
    if (neg) x = -x
    ok = .true.
  end subroutine fast_token

  !> The guarded fast path: Clinger's conversion, and the runtime's (ok = .false.) for the tokens whose exact value
  !  lies within 2**-10 ulp of a midpoint between two doubles (where a runtime converting through extended
  !  precision can round twice) or is a power of two, besides those off Clinger's range.
  subroutine fastg_token(s, x, ok)
    character(len=*), intent(in) :: s
    real(R8), intent(out) :: x
    logical, intent(out) :: ok
    real(R8), parameter :: G = 2._R8**(-10)
    integer(I8) :: m
    integer :: ex
    logical :: neg
    real(R8) :: xm, p, e, u, d
    x = 0._R8
    call parse_token(s, m, ex, neg, ok)
    if (.not. ok) return
    if (m == 0_I8) return
    ok = .false.
    if (m > TWO53 .or. abs(ex) > 22) return
    xm = real(m, R8)
    if (ex >= 0) then
      p = xm*P10(ex)
      if (fraction(p) == 0.5_R8) return
      e = probe_fma(xm, P10(ex), -p)
      u = spacing(p)
      if (abs(abs(e) - 0.5_R8*u) <= G*u) return
    else
      d = P10(-ex)
      p = xm/d
      if (fraction(p) == 0.5_R8) return
      e = probe_fma(-p, d, xm)
      u = spacing(p)
      if (abs(abs(e) - 0.5_R8*u*d) <= G*u*d) return
    endif
    x = p
    if (neg) x = -x
    ok = .true.
  end subroutine fastg_token

  !> The conversion in x87 extended precision, rounded to double: m (18 digits at most, exact in 64 bits) times or
  !  divided by 10**|ex| (exact for |ex| <= 27); ok = .false. for every other token.
  subroutine x87_token(s, x, ok)
    character(len=*), intent(in) :: s
    real(R8), intent(out) :: x
    logical, intent(out) :: ok
    integer(I8) :: m
    integer :: ex
    logical :: neg
    x = 0._R8
    call parse_token(s, m, ex, neg, ok)
    if (.not. ok) return
    if (m == 0_I8) return
    ok = .false.
    if (abs(ex) > 27) return
    x = probe_x87(m, ex)
    if (neg) x = -x
    ok = .true.
  end subroutine x87_token

  !> ORION's conversion of one token: a list-directed READ from a blank-padded character(100) (ios -999: longer).
  subroutine rtl_token(s, x, ios)
    character(len=*), intent(in) :: s
    real(R8), intent(out) :: x
    integer, intent(out) :: ios
    character(len=100) :: tok
    x = 0._R8
    if (len(s) > len(tok)) then
      ios = -999
      return
    endif
    tok = ' '
    tok(1:len(s)) = s
    read(tok, *, iostat=ios) x
  end subroutine rtl_token

  !> glibc's strtod on the token at buf(p:p+l-1); ok = .false. when it does not take the whole token.
  subroutine strtod_token(buf, p, l, x, ok)
    character(len=*), intent(in) :: buf
    integer(I8), intent(in) :: p
    integer, intent(in) :: l
    real(R8), intent(out) :: x
    logical, intent(out) :: ok
    integer(c_int64_t) :: e
    x = probe_strtod(buf, int(p - 1, c_int64_t), e)
    ok = (e == p - 1 + l)
  end subroutine strtod_token

  !> Tokens k1..k2 (all plain, starts ts, lengths tl) by one list-directed READ of their text, line ends, tabs and
  !  commas made blanks in the copy cp; ios of the READ.
  subroutine batch_read(buf, ts, tl, k1, k2, x, cp, ios)
    character(len=*), intent(in) :: buf
    integer(I8), intent(in) :: ts(:)
    integer, intent(in) :: tl(:)
    integer(I8), intent(in) :: k1, k2
    real(R8), intent(inout) :: x(:)
    character(len=:), allocatable, intent(inout) :: cp
    integer, intent(out) :: ios
    integer(I8) :: a, b, L, p
    a = ts(k1); b = ts(k2) + tl(k2) - 1; L = b - a + 1
    if (.not. allocated(cp)) allocate(character(len=max(L, 1048576_I8)) :: cp)
    if (len(cp, kind=I8) < L) then
      deallocate(cp); allocate(character(len=L) :: cp)
    endif
    do p = 1, L
      if (is_sep(iachar(buf(a+p-1:a+p-1)))) then
        cp(p:p) = ' '
      else
        cp(p:p) = buf(a+p-1:a+p-1)
      endif
    enddo
    read(cp(1:L), *, iostat=ios) x(k1:k2)
  end subroutine batch_read

  !> The runtime's E23.15E3 of v(:) into the records rec(:) (24 characters, the last a line feed): one WRITE.
  subroutine fmt_e23_rtl(v, rec)
    real(R8), intent(in) :: v(:)
    character(len=24), intent(inout) :: rec(:)
    integer :: i
    write(rec, '(E23.15E3)') v
    do i = 1, size(rec)
      rec(i)(24:24) = LF
    enddo
  end subroutine fmt_e23_rtl

  !> The floor nf and the nearest integer n of a*10**j by an exact decision (the product's or the quotient's error
  !  by fma); ok = .false. on an exact tie, out of the scaled range or in the two-step scale's error band.
  subroutine round_scaled(a, j, n, nf, ok)
    real(R8), intent(in) :: a
    integer, intent(in) :: j
    integer(I8), intent(out) :: n, nf
    logical, intent(out) :: ok
    real(R8) :: p, e, fr, d, dv, q, r, p1, e1, p2, e2, t
    logical :: below
    ok = .false.; n = 0_I8; nf = 0_I8
    if (j >= 0 .and. j <= 22) then
      p = a*P10(j)
      if (.not. (p < 4.e15_R8)) return
      e = probe_fma(a, P10(j), -p)
      n = int(p, I8)
      fr = p - real(n, R8)
      d = (fr - 0.5_R8) + e
      below = (fr == 0._R8 .and. e < 0._R8)
    elseif (j < 0 .and. j >= -22) then
      dv = P10(-j)
      q = a/dv
      if (.not. (q < 4.e15_R8)) return
      r = probe_fma(-q, dv, a)
      n = int(q, I8)
      fr = q - real(n, R8)
      d = probe_fma(fr - 0.5_R8, dv, r)
      below = (fr == 0._R8 .and. r < 0._R8)
    elseif (j > 22 .and. j <= 44) then
      p1 = a*P10(22)
      e1 = probe_fma(a, P10(22), -p1)
      p2 = p1*P10(j-22)
      if (.not. (p2 < 4.e15_R8)) return
      e2 = probe_fma(p1, P10(j-22), -p2)
      t = e2 + e1*P10(j-22)
      n = int(p2, I8)
      fr = p2 - real(n, R8)
      d = (fr - 0.5_R8) + t
      if (abs(d) <= 1.e-12_R8) return
      if (fr == 0._R8 .and. abs(t) <= 1.e-12_R8) return
      below = (fr == 0._R8 .and. t < 0._R8)
    else
      return
    endif
    nf = n
    if (below) nf = n - 1_I8
    if (d > 0._R8) then
      n = n + 1_I8
    elseif (.not. (d < 0._R8)) then
      return
    endif
    ok = .true.
  end subroutine round_scaled

  !> E23.15E3 of x into s(1:23) without the runtime: ok = .false. for NaN, Inf, a subnormal, |x| outside
  !  [1e-30, 1e37) and the values round_scaled cannot decide, which the caller formats with the runtime. The decimal
  !  exponent k is the one with 10**14 <= floor(|x|*10**(15-k)) < 10**15, decided on the exact floor.
  subroutine fmt_e23_fast(x, s, ok)
    real(R8), intent(in) :: x
    character(len=*), intent(inout) :: s
    logical, intent(out) :: ok
    real(R8) :: a
    integer(I8) :: n, nf, bits
    integer :: k, it, ae, q
    logical :: done
    ok = .false.
    ! zero and subnormals by their bits: with denormals-are-zero (ifx's default) a subnormal compares equal to 0
    bits = iand(transfer(x, 0_I8), huge(0_I8))
    if (bits == 0_I8) then
      s(1:23) = ' 0.000000000000000E+000'
      ok = .true.
      return
    endif
    if (ishft(bits, -52) == 0_I8) return
    a = abs(x)
    if (.not. (a >= 1.e-30_R8 .and. a < 1.e37_R8)) return
    k = floor(log10(a)) + 1
    done = .false.
    do it = 1, 4
      call round_scaled(a, 15 - k, n, nf, ok)
      if (.not. ok) return
      if (nf >= 1000000000000000_I8) then
        k = k + 1
      elseif (nf < 100000000000000_I8) then
        k = k - 1
      else
        if (n == 1000000000000000_I8) then
          n = 100000000000000_I8; k = k + 1
        endif
        done = .true.
        exit
      endif
    enddo
    ok = done
    if (.not. ok) return
    if (x < 0._R8) then
      s(1:1) = '-'
    else
      s(1:1) = ' '
    endif
    s(2:3) = '0.'
    do q = 18, 4, -1
      s(q:q) = achar(48 + int(mod(n, 10_I8))); n = n/10_I8
    enddo
    s(19:19) = 'E'
    if (k < 0) then
      s(20:20) = '-'
    else
      s(20:20) = '+'
    endif
    ae = abs(k)
    s(21:21) = achar(48 + ae/100); s(22:22) = achar(48 + mod(ae/10, 10)); s(23:23) = achar(48 + mod(ae, 10))
  end subroutine fmt_e23_fast

end module io_probe_kern
