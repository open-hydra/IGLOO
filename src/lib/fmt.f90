!> The particle files' edit descriptors written without the runtime: F(w,d), E(w,d)Ee, ES(w,d)Ee, I(w) and I0, byte for
!  byte the compiler runtime's internal write. A field the fast path cannot decide exactly (a near-tie, a value outside
!  the fast range, a subnormal or non-finite value, an integer wider than its field) is written by that runtime write.
module IGLOO_fmt
  use, intrinsic :: iso_fortran_env, only: R8 => real64, I8 => int64
  implicit none
  private
  public :: fmtSpec, fmt_compile, fmt_row, fmt_int

  integer, parameter :: KF = 1, KE = 2, KES = 3, KI = 4
  !> 10**j, exact for 0 <= j <= 22
  real(R8), parameter :: P10(0:22) = [1._R8, 1.e1_R8, 1.e2_R8, 1.e3_R8, 1.e4_R8, 1.e5_R8, 1.e6_R8, 1.e7_R8, 1.e8_R8, &
    1.e9_R8, 1.e10_R8, 1.e11_R8, 1.e12_R8, 1.e13_R8, 1.e14_R8, 1.e15_R8, 1.e16_R8, 1.e17_R8, 1.e18_R8, 1.e19_R8, &
    1.e20_R8, 1.e21_R8, 1.e22_R8]
  integer(I8), parameter :: I10(0:18) = [1_I8, 10_I8, 100_I8, 1000_I8, 10000_I8, 100000_I8, 1000000_I8, 10000000_I8, &
    100000000_I8, 1000000000_I8, 10000000000_I8, 100000000000_I8, 1000000000000_I8, 10000000000000_I8, &
    100000000000000_I8, 1000000000000000_I8, 10000000000000000_I8, 100000000000000000_I8, 1000000000000000000_I8]
  !> Half-integer window of the scaled value, in its ulps: inside it the runtime decides
  real(R8), parameter :: TIEULP = 4._R8
  !> Integer columns arrive as reals, exact below 2**53
  real(R8), parameter :: IMAX = 9007199254740992._R8

  !> A record layout compiled from a format such as '(7F12.6,2E13.6E2,I8)'.
  type :: fmtSpec
    integer :: nf = 0                                   !> fields
    integer :: maxlen = 0                               !> longest record [characters]; an I0 field counts 20
    integer, allocatable :: kind(:), w(:), d(:), e(:)
    character(len=16), allocatable :: rt(:)             !> the field's own format, for the runtime write
  end type fmtSpec

contains

  !> Compile '([r]Fw.d, [r]Ew.dEe, [r]ESw.dEe, [r]Iw, ...)', I0 included; any other item stops.
  subroutine fmt_compile(spec, form)
    type(fmtSpec), intent(out) :: spec
    character(len=*), intent(in) :: form
    character(len=len(form)) :: s
    character :: c
    integer :: ls, i, pass, nf, r, k, kd, w, d, e
    logical :: ok
    ls = 0; s = ' '
    do i = 1, len(form)
      c = form(i:i)
      if (c == ' ') cycle
      if (c >= 'a' .and. c <= 'z') c = achar(iachar(c) - 32)
      ls = ls + 1; s(ls:ls) = c
    enddo
    if (ls < 3) call refuse(form)
    if (s(1:1) /= '(' .or. s(ls:ls) /= ')') call refuse(form)
    do pass = 1, 2
      nf = 0; i = 2
      do
        call scan_item(s(1:ls), i, r, kd, w, d, e, ok)
        if (.not. ok) call refuse(form)
        do k = 1, r
          nf = nf + 1
          if (pass == 2) then
            spec%kind(nf) = kd; spec%w(nf) = w; spec%d(nf) = d; spec%e(nf) = e
          endif
        enddo
        if (i == ls) exit
        if (s(i:i) /= ',') call refuse(form)
        i = i + 1
      enddo
      if (pass == 1) allocate(spec%kind(nf), spec%w(nf), spec%d(nf), spec%e(nf), spec%rt(nf))
    enddo
    spec%nf = nf
    spec%maxlen = 0
    do k = 1, nf
      select case (spec%kind(k))
      case (KF);  write(spec%rt(k), '(A,I0,A,I0,A)') '(F', spec%w(k), '.', spec%d(k), ')'
      case (KE);  write(spec%rt(k), '(A,I0,A,I0,A,I0,A)') '(E', spec%w(k), '.', spec%d(k), 'E', spec%e(k), ')'
      case (KES); write(spec%rt(k), '(A,I0,A,I0,A,I0,A)') '(ES', spec%w(k), '.', spec%d(k), 'E', spec%e(k), ')'
      case (KI);  write(spec%rt(k), '(A,I0,A)') '(I', spec%w(k), ')'
      end select
      if (spec%kind(k) == KI .and. spec%w(k) == 0) then
        spec%maxlen = spec%maxlen + 20
      else
        spec%maxlen = spec%maxlen + spec%w(k)
      endif
    enddo
  end subroutine fmt_compile


  !> Stop on a format the module does not write.
  subroutine refuse(form)
    character(len=*), intent(in) :: form
    write(*,'(A)') 'IGLOO_fmt: unsupported format '//trim(form)
    error stop 'IGLOO_fmt: unsupported format'
  end subroutine refuse


  !> One item at s(i:): repeat count, descriptor and its limits; i ends past it. ok = .false. on anything else.
  subroutine scan_item(s, i, r, kd, w, d, e, ok)
    character(len=*), intent(in) :: s
    integer, intent(inout) :: i
    integer, intent(out) :: r, kd, w, d, e
    logical, intent(out) :: ok
    integer :: nd
    ok = .false.; kd = 0; w = -1; d = 0; e = 0
    r = number(s, i)
    if (r == -1) r = 1
    if (r < 1) return
    if (i + 1 <= len(s)) then
      if (s(i:i+1) == 'ES') then
        kd = KES; i = i + 2
      endif
    endif
    if (kd == 0) then
      select case (s(i:i))
      case ('F'); kd = KF
      case ('E'); kd = KE
      case ('I'); kd = KI
      case default; return
      end select
      i = i + 1
    endif
    w = number(s, i)
    if (w < 0) return
    if (kd /= KI) then
      if (s(i:i) /= '.') return
      i = i + 1
      d = number(s, i)
      if (d < 0) return
    endif
    if (kd == KE .or. kd == KES) then
      if (s(i:i) /= 'E') return
      i = i + 1
      e = number(s, i)
      if (e < 1 .or. e > 9) return
    endif
    select case (kd)
    case (KF)
      ok = d >= 1 .and. w >= d + 3 .and. w <= 17
    case (KE, KES)
      nd = d; if (kd == KES) nd = d + 1
      ok = nd >= 1 .and. nd <= 15 .and. w >= d + e + 5
    case (KI)
      ok = w <= 20
    end select
  end subroutine scan_item


  !> The unsigned decimal at s(i:), i past it; -1 when none, -2 when too long.
  integer function number(s, i)
    character(len=*), intent(in) :: s
    integer, intent(inout) :: i
    number = -1
    do while (i <= len(s))
      if (s(i:i) < '0' .or. s(i:i) > '9') exit
      if (number < 0) number = 0
      if (number > 99999) then
        number = -2; return
      endif
      number = 10*number + (iachar(s(i:i)) - 48)
      i = i + 1
    enddo
  end function number


  !> One record of spec's layout from v(1:spec%nf) at buf(n+1:), n past it; nslow counts the fields the runtime wrote.
  !  The caller keeps len(buf) - n >= spec%maxlen; an integer column holds an exact integer.
  pure subroutine fmt_row(spec, v, buf, n, nslow)
    type(fmtSpec), intent(in) :: spec
    real(R8), intent(in) :: v(:)
    character(len=*), intent(inout) :: buf
    integer, intent(inout) :: n
    integer, intent(inout), optional :: nslow
    integer :: k, w
    logical :: slow
    do k = 1, spec%nf
      w = spec%w(k)
      select case (spec%kind(k))
      case (KF)
        call put_f(v(k), w, spec%d(k), buf, n, slow)
      case (KE)
        call put_e(v(k), w, spec%d(k), spec%e(k), .false., buf, n, slow)
      case (KES)
        call put_e(v(k), w, spec%d(k), spec%e(k), .true., buf, n, slow)
      case default
        if (fclass(v(k)) <= 1 .and. abs(v(k)) < IMAX) then
          call fmt_int(int(v(k), I8), w, buf, n, nslow)
        else
          !> not an integer below 2**53: asterisks, as the runtime fills an integer that does not fit
          buf(n+1:n+max(w, 1)) = repeat('*', max(w, 1))
          n = n + max(w, 1)
        endif
        cycle
      end select
      if (slow) then
        write(buf(n+1:n+w), spec%rt(k)) v(k)
        n = n + w
        if (present(nslow)) nslow = nslow + 1
      endif
    enddo
  end subroutine fmt_row


  !> I(w) of iv at buf(n+1:), n past it; w = 0 is I0. nslow counts a runtime write.
  pure subroutine fmt_int(iv, w, buf, n, nslow)
    integer(I8), intent(in) :: iv
    integer, intent(in) :: w
    character(len=*), intent(inout) :: buf
    integer, intent(inout) :: n
    integer, intent(inout), optional :: nslow
    character(len=8) :: f
    logical :: slow
    call put_i(iv, w, buf, n, slow)
    if (.not. slow) return
    write(f, '(A,I0,A)') '(I', w, ')'
    if (w == 0) then
      write(buf(n+1:n+20), f) iv
      n = n + len_trim(buf(n+1:n+20))
    else
      write(buf(n+1:n+w), f) iv
      n = n + w
    endif
    if (present(nslow)) nslow = nslow + 1
  end subroutine fmt_int


  !> The runtime's sign: gfortran signs -0.0, the Intel runtime does not; both sign a negative that rounds to zero.
  pure logical function negative(x)
    real(R8), intent(in) :: x
#ifdef __GFORTRAN__
    negative = btest(transfer(x, 0_I8), 63)
#else
    negative = btest(transfer(x, 0_I8), 63) .and. iand(transfer(x, 0_I8), huge(0_I8)) /= 0_I8
#endif
  end function negative


  !> From the bits, whatever flush-to-zero and fast-math assume: 0 zero, 1 normal, 2 subnormal, 3 NaN or Inf.
  pure integer function fclass(x)
    real(R8), intent(in) :: x
    integer(I8) :: b, ex
    b  = iand(transfer(x, 0_I8), huge(0_I8))
    ex = shiftr(b, 52)
    if (ex == 2047_I8) then
      fclass = 3
    elseif (ex /= 0_I8) then
      fclass = 1
    elseif (b == 0_I8) then
      fclass = 0
    else
      fclass = 2
    endif
  end function fclass


  !> Nearest integer m to p = a*10**j (a >= 0, 0 <= j <= 22, p < 2**52); tie when p is within TIEULP ulps of a half.
  pure subroutine round_scaled(a, j, m, tie)
    real(R8), intent(in) :: a
    integer, intent(in) :: j
    integer(I8), intent(out) :: m
    logical, intent(out) :: tie
    real(R8) :: p, fr
    p   = a*P10(j)
    m   = int(p, I8)
    fr  = p - real(m, R8)
    tie = abs(fr - 0.5_R8) <= TIEULP*spacing(p)
    if (fr > 0.5_R8) m = m + 1
  end subroutine round_scaled


  !> F(w,d) at buf(n+1:n+w); fast while |x| < 10**(w-d-2) - 1/2 (the sign keeps a column).
  pure subroutine put_f(x, w, d, buf, n, slow)
    real(R8), intent(in) :: x
    integer, intent(in) :: w, d
    character(len=*), intent(inout) :: buf
    integer, intent(inout) :: n
    logical, intent(out) :: slow
    real(R8) :: ax
    integer(I8) :: m
    integer :: q, k
    slow = .true.
    if (fclass(x) > 1) return
    ax = abs(x)
    if (.not. ax < P10(w-d-2) - 0.5_R8) return
    call round_scaled(ax, d, m, slow)
    if (slow) return
    q = n + w
    do k = 1, d
      buf(q:q) = achar(48 + int(mod(m, 10_I8))); m = m/10_I8; q = q - 1
    enddo
    buf(q:q) = '.'; q = q - 1
    do
      buf(q:q) = achar(48 + int(mod(m, 10_I8))); m = m/10_I8; q = q - 1
      if (m == 0) exit
    enddo
    if (negative(x)) then
      buf(q:q) = '-'; q = q - 1
    endif
    if (q > n) buf(n+1:q) = ' '
    n = n + w
  end subroutine put_f


  !> E(w,d)Ee (es false: 0.d1..dd) or ES(w,d)Ee (es true: d0.d1..dd) at buf(n+1:n+w); fast while the significand's
  !  scale 10**j has 0 <= j <= 22.
  pure subroutine put_e(x, w, d, e, es, buf, n, slow)
    real(R8), intent(in) :: x
    integer, intent(in) :: w, d, e
    logical, intent(in) :: es
    character(len=*), intent(inout) :: buf
    integer, intent(inout) :: n
    logical, intent(out) :: slow
    real(R8) :: ax
    integer(I8) :: m
    integer :: nd, j, it, ex, aex, q, k, c
    slow = .true.
    c = fclass(x)
    if (c > 1) return
    nd = d; if (es) nd = d + 1
    ax = abs(x)
    if (c == 0) then
      m = 0_I8; ex = 0
    else
      j = nd - 1 - floor(log10(ax))
      do it = 1, 3
        slow = .true.
        if (j < 0 .or. j > 22) return
        call round_scaled(ax, j, m, slow)
        if (slow) return
        if (m >= I10(nd)) then
          j = j - 1
        elseif (m < I10(nd-1)) then
          j = j + 1
        else
          exit
        endif
      enddo
      slow = .true.
      if (m >= I10(nd) .or. m < I10(nd-1)) return
      ex = nd - j
      if (es) ex = ex - 1
    endif
    aex = abs(ex)
    if (aex >= I10(e)) return
    slow = .false.
    q = n + w
    do k = 1, e
      buf(q:q) = achar(48 + mod(aex, 10)); aex = aex/10; q = q - 1
    enddo
    if (ex < 0) then
      buf(q:q) = '-'
    else
      buf(q:q) = '+'
    endif
    buf(q-1:q-1) = 'E'; q = q - 2
    do k = 1, d
      buf(q:q) = achar(48 + int(mod(m, 10_I8))); m = m/10_I8; q = q - 1
    enddo
    buf(q:q) = '.'; q = q - 1
    buf(q:q) = achar(48 + int(m)); q = q - 1
    if (negative(x)) then
      buf(q:q) = '-'; q = q - 1
    endif
    if (q > n) buf(n+1:q) = ' '
    n = n + w
  end subroutine put_e


  !> I(w) at buf(n+1:n+w), right-justified, or I0 (w = 0) at its own width; slow when wider than w.
  pure subroutine put_i(iv, w, buf, n, slow)
    integer(I8), intent(in) :: iv
    integer, intent(in) :: w
    character(len=*), intent(inout) :: buf
    integer, intent(inout) :: n
    logical, intent(out) :: slow
    integer(I8) :: a
    integer :: nd, nc, q, k
    slow = .true.
    if (iv == -huge(iv) - 1_I8) return
    a = abs(iv)
    nd = 1
    do while (nd < 19)
      if (a < I10(nd)) exit
      nd = nd + 1
    enddo
    nc = nd
    if (iv < 0_I8) nc = nc + 1
    if (w > 0 .and. nc > w) return
    slow = .false.
    q = n + max(w, nc)
    do k = 1, nd
      buf(q:q) = achar(48 + int(mod(a, 10_I8))); a = a/10_I8; q = q - 1
    enddo
    if (iv < 0_I8) then
      buf(q:q) = '-'; q = q - 1
    endif
    if (q > n) buf(n+1:q) = ' '
    n = n + max(w, nc)
  end subroutine put_i

end module IGLOO_fmt
