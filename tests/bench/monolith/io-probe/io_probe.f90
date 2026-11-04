!> The I/O probe of the timing harness (tests/bench/monolith/README.md): where the gas read and the grid write
!  spend their time, kernel by kernel. A measurement tool, never a test (CMake option IGLOO_IO_PROBE). One IPROBE
!  line of key=value fields per measurement on standard output.
!    io-probe read-split <file.tec>   one thread: stream reads (cold, warm), ORION's reader, one formatted line read
!                                     per line (with and without ORION's zone test), the data regions and tokens,
!                                     the conversion kernels rtl-token, rtl-batch, fast, fastg, x87 and strtod, the
!                                     store into ORION's layout, and E23.15E3 by the runtime (fmt-batch) and by
!                                     fmt-fast
!    io-probe randsets <dir>          the same kernels on make_values.py's sets (values.bin, v15E.txt, v16E.txt,
!                                     v17E.txt) and on vE23.txt, written here from values.bin by the runtime
!    io-probe census <list>           per .tec file of the list: tokens, the fast paths' shares, the mismatches
!    io-probe scale <file.tec> <ntok> <T>c|<T>s ...   the read pipeline (pieces counted, then converted) of fast,
!                                     fastg and rtl-batch, the store, fmt-batch and fmt-fast at T threads, proc_bind
!                                     close (c) or spread (s); outputs compared with the first spec's (one thread)
!    io-probe swrite <file> <bytes> <block>            a buffer written by stream access in blocks, to close
!    io-probe bcast <I> <J> <K> <nvar> <reps>          IGLOO's gas broadcast of a one-block field (under mpirun)
module io_probe_run
#ifdef USE_MPI
  use mpi
#endif
  use omp_lib
  use, intrinsic :: iso_c_binding, only: c_null_char
  use io_probe_kern
  use Lib_ORION_data, only: orion_data
  use Lib_Tecplot, only: tec_read_structured_multiblock
  use IGLOO_Mod_MPI, only: mpi_init_env, mpi_finalize_env, mpi_bcast_orion, mpi_size_, mpi_is_root
  implicit none
  private
  public :: emit, ival8, read_split, randsets, census, scale, swrite, bcast

  type :: meter
    real(R8) :: w = 0._R8, u = 0._R8, s = 0._R8
    integer(I8) :: c = -1_I8
  end type meter

  integer, parameter :: SLICE = 4096
  character(len=8), parameter :: FKN(3) = ['fast    ', 'fastg   ', 'x87     ']
  integer(I8), parameter :: PIECE_BYTES = 1048576_I8

contains

  ! ---------------------------------------------------------------- reporting
  subroutine emit(line)
    character(len=*), intent(in) :: line
    write(*,'(A)') trim(line)
    flush(6)
  end subroutine emit

  function fr(key, x) result(s)
    character(len=*), intent(in) :: key
    real(R8), intent(in) :: x
    character(len=:), allocatable :: s
    character(len=32) :: t
    write(t,'(ES14.6E3)') x
    s = ' '//key//'='//trim(adjustl(t))
  end function fr

  function fi(key, i) result(s)
    character(len=*), intent(in) :: key
    integer(I8), intent(in) :: i
    character(len=:), allocatable :: s
    character(len=32) :: t
    write(t,'(I0)') i
    s = ' '//key//'='//trim(t)
  end function fi

  function fs(key, v) result(s)
    character(len=*), intent(in) :: key, v
    character(len=:), allocatable :: s
    s = ' '//key//'='//trim(v)
  end function fs

  integer(I8) function ival8(t)
    character(len=*), intent(in) :: t
    read(t, *) ival8
  end function ival8

  function now() result(m)
    type(meter) :: m
    m%w = omp_get_wtime()
    call probe_times(m%u, m%s)
    m%c = probe_cyc_read()
  end function now

  !> wall, user and system seconds, and user cycles of the calling thread, from m0 to m1; per value when n > 0.
  function cost(m0, m1, n) result(s)
    type(meter), intent(in) :: m0, m1
    integer(I8), intent(in) :: n
    character(len=:), allocatable :: s
    s = fr('wall_s', m1%w - m0%w)//fr('user_s', m1%u - m0%u)//fr('sys_s', m1%s - m0%s)
    if (m0%c >= 0 .and. m1%c >= 0) then
      s = s//fi('cycles', m1%c - m0%c)
      if (n > 0) s = s//fr('cyc_per_value', real(m1%c - m0%c, R8)/real(n, R8))
    endif
    if (n > 0) s = s//fr('ns_per_value', 1.e9_R8*(m1%w - m0%w)/real(n, R8))
  end function cost

  ! ---------------------------------------------------------------- comparisons
  !> Elements whose bits differ (-0 against +0 and NaN payloads count).
  integer(I8) function ndiff(a, b, n)
    integer(I8), intent(in) :: n
    real(R8), intent(in) :: a(*), b(*)
    integer(I8) :: i
    ndiff = 0
    do i = 1, n
      if (transfer(a(i), 0_I8) /= transfer(b(i), 0_I8)) ndiff = ndiff + 1
    enddo
  end function ndiff

  integer(I8) function rdiff(a, b)
    character(len=24), intent(in) :: a(:), b(:)
    integer(I8) :: i
    rdiff = 0
    do i = 1, size(a, kind=I8)
      if (a(i) /= b(i)) rdiff = rdiff + 1
    enddo
  end function rdiff

  ! ---------------------------------------------------------------- files
  !> The whole file into buf(1:n) by stream reads of 256 MiB; buf(n+1:n+1) is a NUL (strtod's end).
  subroutine slurp(path, buf, n)
    character(len=*), intent(in) :: path
    character(len=:), allocatable, intent(inout) :: buf
    integer(I8), intent(out) :: n
    integer :: u, ios
    integer(I8) :: p, m
    integer(I8), parameter :: CH = 268435456_I8
    open(newunit=u, file=path, access='stream', form='unformatted', action='read', status='old', iostat=ios)
    if (ios /= 0) then
      call emit('IPROBE error=cannot-open file='//path)
      error stop 1
    endif
    inquire(unit=u, size=n)
    if (allocated(buf)) deallocate(buf)
    allocate(character(len=n+1) :: buf)
    p = 1
    do while (p <= n)
      m = min(CH, n - p + 1)
      read(u, pos=p, iostat=ios) buf(p:p+m-1)
      if (ios /= 0) then
        call emit('IPROBE error=read file='//path)
        error stop 1
      endif
      p = p + m
    enddo
    close(u)
    buf(n+1:n+1) = achar(0)
  end subroutine slurp

  !> Tokens of the data regions into ts, tl; the longest token and the longest line; the seconds of the region
  !  scan, the token count, the token fill and the line-length scan in tsub.
  subroutine tokenize(buf, n, ts, tl, ntok, longest, lline, nreg, nzone, tsub)
    character(len=*), intent(in) :: buf
    integer(I8), intent(in) :: n
    integer(I8), allocatable, intent(out) :: ts(:)
    integer, allocatable, intent(out) :: tl(:)
    integer(I8), intent(out) :: ntok
    integer, intent(out) :: longest, nreg, nzone
    integer(I8), intent(out) :: lline
    real(R8), intent(out), optional :: tsub(4)
    integer(I8), allocatable :: rs(:), re(:)
    integer(I8) :: k, p, q
    integer :: r
    real(R8) :: t(5)
    t(1) = omp_get_wtime()
    call find_regions(buf, n, rs, re, nreg, nzone)
    t(2) = omp_get_wtime()
    ntok = 0
    do r = 1, nreg
      ntok = ntok + count_tokens(buf, rs(r), re(r))
    enddo
    t(3) = omp_get_wtime()
    allocate(ts(max(ntok, 1_I8)), tl(max(ntok, 1_I8)))
    k = 0; longest = 0
    do r = 1, nreg
      call fill_tokens(buf, rs(r), re(r), ts, tl, k, longest)
    enddo
    t(4) = omp_get_wtime()
    lline = 0; p = 1
    do while (p <= n)
      q = next_line_end(buf, p, n)
      lline = max(lline, q - p)
      p = q + 1
    enddo
    t(5) = omp_get_wtime()
    if (present(tsub)) tsub = t(2:5) - t(1:4)
  end subroutine tokenize

  ! ---------------------------------------------------------------- serial kernels on tokens
  subroutine k_rtl_token(buf, ts, tl, x, nfail)
    character(len=*), intent(in) :: buf
    integer(I8), intent(in) :: ts(:)
    integer, intent(in) :: tl(:)
    real(R8), intent(out) :: x(:)
    integer(I8), intent(out) :: nfail
    integer(I8) :: k
    integer :: ios
    nfail = 0
    do k = 1, size(x, kind=I8)
      call rtl_token(buf(ts(k):ts(k)+tl(k)-1), x(k), ios)
      if (ios /= 0) nfail = nfail + 1
    enddo
  end subroutine k_rtl_token

  !> One list-directed READ per run of at most 65536 plain tokens; every other token alone.
  subroutine k_rtl_batch(buf, ts, tl, x, nfail, nonplain)
    character(len=*), intent(in) :: buf
    integer(I8), intent(in) :: ts(:)
    integer, intent(in) :: tl(:)
    real(R8), intent(inout) :: x(:)
    integer(I8), intent(out) :: nfail, nonplain
    integer(I8) :: k, k1, nt
    integer :: ios
    character(len=:), allocatable :: cp
    nfail = 0; nonplain = 0
    nt = size(x, kind=I8)
    k1 = 1
    do k = 1, nt
      if (.not. plain_token(buf(ts(k):ts(k)+tl(k)-1))) then
        if (k > k1) call run(k1, k - 1)
        nonplain = nonplain + 1
        call rtl_token(buf(ts(k):ts(k)+tl(k)-1), x(k), ios)
        if (ios /= 0) nfail = nfail + 1
        k1 = k + 1
      elseif (k - k1 + 1 == 65536) then
        call run(k1, k)
        k1 = k + 1
      endif
    enddo
    if (nt >= k1) call run(k1, nt)
  contains
    subroutine run(a, b)
      integer(I8), intent(in) :: a, b
      call batch_read(buf, ts, tl, a, b, x, cp, ios)
      if (ios /= 0) nfail = nfail + (b - a + 1)
    end subroutine run
  end subroutine k_rtl_batch

  !> kind 1: Clinger's fast path; 2: the guarded fast path; 3: x87; the runtime's READ for the other tokens.
  subroutine k_fast(buf, ts, tl, x, nfb, nfail, kind)
    character(len=*), intent(in) :: buf
    integer(I8), intent(in) :: ts(:)
    integer, intent(in) :: tl(:)
    real(R8), intent(out) :: x(:)
    integer(I8), intent(out) :: nfb, nfail
    integer, intent(in) :: kind
    integer(I8) :: k
    integer :: ios
    logical :: ok
    nfb = 0; nfail = 0
    do k = 1, size(x, kind=I8)
      select case (kind)
      case (1)
        call fast_token(buf(ts(k):ts(k)+tl(k)-1), x(k), ok)
      case (2)
        call fastg_token(buf(ts(k):ts(k)+tl(k)-1), x(k), ok)
      case default
        call x87_token(buf(ts(k):ts(k)+tl(k)-1), x(k), ok)
      end select
      if (.not. ok) then
        nfb = nfb + 1
        call rtl_token(buf(ts(k):ts(k)+tl(k)-1), x(k), ios)
        if (ios /= 0) nfail = nfail + 1
      endif
    enddo
  end subroutine k_fast

  subroutine k_strtod(buf, ts, tl, x, npart)
    character(len=*), intent(in) :: buf
    integer(I8), intent(in) :: ts(:)
    integer, intent(in) :: tl(:)
    real(R8), intent(out) :: x(:)
    integer(I8), intent(out) :: npart
    integer(I8) :: k
    logical :: ok
    npart = 0
    do k = 1, size(x, kind=I8)
      call strtod_token(buf, ts(k), tl(k), x(k), ok)
      if (.not. ok) npart = npart + 1
    enddo
  end subroutine k_strtod

  subroutine k_fmt_batch(x, rec)
    real(R8), intent(in) :: x(:)
    character(len=24), intent(inout) :: rec(:)
    integer(I8) :: k, nv
    nv = size(x, kind=I8)
    do k = 1, nv, SLICE
      call fmt_e23_rtl(x(k:min(nv, k+SLICE-1)), rec(k:min(nv, k+SLICE-1)))
    enddo
  end subroutine k_fmt_batch

  subroutine k_fmt_fast(x, rec, nfb)
    real(R8), intent(in) :: x(:)
    character(len=24), intent(inout) :: rec(:)
    integer(I8), intent(out) :: nfb
    integer(I8) :: k
    logical :: ok
    nfb = 0
    do k = 1, size(x, kind=I8)
      call fmt_e23_fast(x(k), rec(k), ok)
      if (.not. ok) then
        write(rec(k)(1:23), '(E23.15E3)') x(k)
        nfb = nfb + 1
      endif
      rec(k)(24:24) = LF
    enddo
  end subroutine k_fmt_fast

  !> The field-major value stream into ORION's layout, ORION's loop order (BLOCK packing).
  subroutine store_orion(x, ni, nj, nk, ndir, nvar, start, mesh, vars)
    real(R8), intent(in) :: x(:)
    integer, intent(in) :: ni, nj, nk, ndir, nvar, start
    real(R8), intent(inout) :: mesh(1:, 0:, 0:, 0:), vars(1:, start:, start:, start:)
    integer(I8) :: t
    integer :: d, s, i, j, k
    t = 0
    do d = 1, ndir
      do k = 0, nk-1; do j = 0, nj-1; do i = 0, ni-1
        t = t + 1; mesh(d,i,j,k) = x(t)
      enddo; enddo; enddo
    enddo
    do s = 1, nvar
      do k = start, max(start, nk-1); do j = start, nj-1; do i = start, ni-1
        t = t + 1; vars(s,i,j,k) = x(t)
      enddo; enddo; enddo
    enddo
  end subroutine store_orion

  ! ---------------------------------------------------------------- read-split
  subroutine read_split(path)
    character(len=*), intent(in) :: path
    character(len=:), allocatable :: buf
    type(orion_data) :: orion
    type(meter) :: m0, m1
    integer(I8) :: n, ntok, nval, nfail, nfb, nonplain, npart, lline, nl, nz, mm, mv
    integer(I8), allocatable :: ts(:)
    integer, allocatable :: tl(:)
    integer :: ev, err, longest, nreg, nzone, ni, nj, nk, ndir, nvar, start, kk
    real(R8), allocatable :: xr(:), xb(:), xf(:), xs(:), mesh(:,:,:,:), vars(:,:,:,:)
    character(len=24), allocatable :: ra(:), rb(:)
    character(len=:), allocatable :: tag
    real(R8) :: tsub(4)

    tag = 'IPROBE step=read'//fs('file', path)
    ev = probe_evict(path//c_null_char)
    m0 = now(); call slurp(path, buf, n); m1 = now()
    call emit(tag//' mode=slurp-cold'//fi('bytes', n)//fi('evict_rc', int(ev, I8))//cost(m0, m1, 0_I8)// &
              fr('MBps', real(n, R8)/1.e6_R8/(m1%w - m0%w)))
    m0 = now(); call slurp(path, buf, n); m1 = now()
    call emit(tag//' mode=slurp-warm'//fi('bytes', n)//cost(m0, m1, 0_I8)//fr('MBps', real(n, R8)/1.e6_R8/(m1%w - m0%w)))

    orion%tec%node = .false.; orion%tec%bc = .false.; orion%tec%format = 'ascii'
    m0 = now(); err = tec_read_structured_multiblock(orion=orion, filename=path); m1 = now()
    if (err /= 0 .or. .not. allocated(orion%block)) then
      call emit(tag//' mode=orion-read'//fi('err', int(err, I8))//' stop=1')
      return
    endif
    ndir = size(orion%block(1)%mesh, 1)
    ni = size(orion%block(1)%mesh, 2); nj = size(orion%block(1)%mesh, 3); nk = size(orion%block(1)%mesh, 4)
    nvar = 0; start = 0
    if (allocated(orion%block(1)%vars)) then
      nvar = size(orion%block(1)%vars, 1); start = lbound(orion%block(1)%vars, 2)
    endif
    mm = size(orion%block(1)%mesh, kind=I8); mv = 0
    if (nvar > 0) mv = size(orion%block(1)%vars, kind=I8)
    nval = mm + mv
    call emit(tag//' mode=orion-read'//fi('err', int(err, I8))//fi('blocks', size(orion%block, kind=I8))// &
              fi('values', nval)//cost(m0, m1, nval))

    m0 = now(); call lines_pass(path, .true., nl, nz); m1 = now()
    call emit(tag//' mode=lines'//fi('lines', nl)//fi('zones', nz)//cost(m0, m1, nl))
    m0 = now(); call lines_pass(path, .false., nl, nz); m1 = now()
    call emit(tag//' mode=lines-raw'//fi('lines', nl)//cost(m0, m1, nl))

    m0 = now(); call tokenize(buf, n, ts, tl, ntok, longest, lline, nreg, nzone, tsub); m1 = now()
    call emit(tag//' mode=tokenize'//fi('tokens', ntok)//fi('expected', nval)//fi('regions', int(nreg, I8))// &
              fi('zones', int(nzone, I8))//fi('longest_token', int(longest, I8))//fi('longest_line', lline)// &
              fr('regions_s', tsub(1))//fr('count_s', tsub(2))//fr('fill_s', tsub(3))//fr('linescan_s', tsub(4))// &
              cost(m0, m1, ntok))
    if (ntok /= nval .or. size(orion%block) /= 1) then
      call emit(tag//' stop=token-count-or-blocks (one block, the expected count, needed below)')
      return
    endif

    allocate(xr(ntok), xb(ntok), xf(ntok), xs(ntok))
    m0 = now(); call k_rtl_token(buf, ts, tl, xr, nfail); m1 = now()
    call emit(tag//' mode=rtl-token'//fi('tokens', ntok)//fi('fail', nfail)//cost(m0, m1, ntok))

    allocate(mesh(1:ndir, 0:ni-1, 0:nj-1, 0:nk-1))
    allocate(vars(1:max(nvar, 1), start:max(start, ni-1), start:max(start, nj-1), start:max(start, nk-1)))
    m0 = now(); call store_orion(xr, ni, nj, nk, ndir, nvar, start, mesh, vars); m1 = now()
    nfail = ndiff(orion%block(1)%mesh, mesh, mm)
    if (nvar > 0) nfail = nfail + ndiff(orion%block(1)%vars, vars, mv)
    call emit(tag//' mode=store'//fi('values', ntok)//fi('mism_vs_orion', nfail)//cost(m0, m1, ntok))
    deallocate(orion%block)

    m0 = now(); call k_rtl_batch(buf, ts, tl, xb, nfail, nonplain); m1 = now()
    call emit(tag//' mode=rtl-batch'//fi('tokens', ntok)//fi('fail', nfail)//fi('nonplain', nonplain)// &
              fi('mism_vs_rtl_token', ndiff(xb, xr, ntok))//cost(m0, m1, ntok))
    do kk = 1, 3
      m0 = now(); call k_fast(buf, ts, tl, xf, nfb, nfail, kk); m1 = now()
      call emit(tag//' mode='//trim(FKN(kk))//fi('tokens', ntok)//fi('fallback', nfb)// &
                fr('fallback_share', real(nfb, R8)/real(ntok, R8))//fi('fail', nfail)// &
                fi('mism_vs_rtl_token', ndiff(xf, xr, ntok))//cost(m0, m1, ntok))
    enddo
    m0 = now(); call k_strtod(buf, ts, tl, xs, npart); m1 = now()
    call emit(tag//' mode=strtod'//fi('tokens', ntok)//fi('partial', npart)// &
              fi('mism_vs_rtl_token', ndiff(xs, xr, ntok))//cost(m0, m1, ntok))
    deallocate(xb, xf, xs, ts, tl)

    allocate(ra(ntok), rb(ntok))
    m0 = now(); call k_fmt_batch(xr, ra); m1 = now()
    call emit(tag//' mode=fmt-batch'//fi('values', ntok)//cost(m0, m1, ntok))
    m0 = now(); call k_fmt_fast(xr, rb, nfb); m1 = now()
    call emit(tag//' mode=fmt-fast'//fi('values', ntok)//fi('fallback', nfb)// &
              fi('mism_vs_fmt_batch', rdiff(rb, ra))//cost(m0, m1, ntok))
  end subroutine read_split

  !> One formatted line read per line into character(1000), ORION's zone test on each when ztest (ORION's
  !  passes 1a and 1b); the zone test is a copy of ORION's is_zone_header and upper_case.
  subroutine lines_pass(path, ztest, nl, nz)
    character(len=*), intent(in) :: path
    logical, intent(in) :: ztest
    integer(I8), intent(out) :: nl, nz
    character(1000) :: line
    integer :: u, ios
    nl = 0; nz = 0
    open(newunit=u, file=path, status='old', action='read')
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      nl = nl + 1
      if (ztest) then
        if (orion_is_zone_header(line)) nz = nz + 1
      endif
    enddo
    close(u)
  end subroutine lines_pass

  logical function orion_is_zone_header(text)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: u
    u = adjustl(orion_upper_case(text))
    orion_is_zone_header = (index(u,'ZONE')==1)
    if (orion_is_zone_header .and. len(u)>4) &
      orion_is_zone_header = .not.(scan(u(5:5),'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_')>0)
  end function orion_is_zone_header

  function orion_upper_case(text) result(out)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: out
    integer :: q, code
    out = text
    do q = 1, len(text)
      code = iachar(out(q:q))
      if (code>=iachar('a') .and. code<=iachar('z')) out(q:q) = achar(code-iachar('a')+iachar('A'))
    enddo
  end function orion_upper_case

  ! ---------------------------------------------------------------- randsets
  subroutine randsets(dir)
    character(len=*), intent(in) :: dir
    character(len=:), allocatable :: buf, tag
    real(R8), allocatable :: v(:), xr(:), xb(:), xf(:), xs(:)
    character(len=24), allocatable :: ra(:), rb(:)
    integer(I8), allocatable :: ts(:)
    integer, allocatable :: tl(:)
    type(meter) :: m0, m1
    integer(I8) :: nv, n, ntok, nfail, nfb, nonplain, npart, lline
    integer :: u, f, longest, nreg, nzone, kk
    character(len=8), parameter :: SETS(4) = ['v15E    ', 'v16E    ', 'v17E    ', 'vE23    ']

    open(newunit=u, file=dir//'/values.bin', access='stream', form='unformatted', action='read', status='old')
    inquire(unit=u, size=nv)
    nv = nv/8
    allocate(v(nv))
    read(u) v
    close(u)
    allocate(ra(nv), rb(nv))
    tag = 'IPROBE step=random set=values'
    m0 = now(); call k_fmt_batch(v, ra); m1 = now()
    call emit(tag//' mode=fmt-batch'//fi('values', nv)//cost(m0, m1, nv))
    m0 = now(); call k_fmt_fast(v, rb, nfb); m1 = now()
    call emit(tag//' mode=fmt-fast'//fi('values', nv)//fi('fallback', nfb)// &
              fi('mism_vs_fmt_batch', rdiff(rb, ra))//cost(m0, m1, nv))
    open(newunit=u, file=dir//'/vE23.txt', access='stream', form='unformatted', action='write', status='replace')
    write(u) ra
    close(u)
    deallocate(ra, rb)

    do f = 1, size(SETS)
      tag = 'IPROBE step=random set='//trim(SETS(f))
      call slurp(dir//'/'//trim(SETS(f))//'.txt', buf, n)
      call tokenize(buf, n, ts, tl, ntok, longest, lline, nreg, nzone)
      call emit(tag//' mode=tokenize'//fi('bytes', n)//fi('tokens', ntok)//fi('expected', nv)// &
                fi('longest_token', int(longest, I8)))
      if (ntok /= nv) cycle
      allocate(xr(ntok), xb(ntok), xf(ntok), xs(ntok))
      m0 = now(); call k_rtl_token(buf, ts, tl, xr, nfail); m1 = now()
      call emit(tag//' mode=rtl-token'//fi('tokens', ntok)//fi('fail', nfail)// &
                fi('mism_vs_values', ndiff(xr, v, ntok))//cost(m0, m1, ntok))
      m0 = now(); call k_rtl_batch(buf, ts, tl, xb, nfail, nonplain); m1 = now()
      call emit(tag//' mode=rtl-batch'//fi('tokens', ntok)//fi('fail', nfail)//fi('nonplain', nonplain)// &
                fi('mism_vs_rtl_token', ndiff(xb, xr, ntok))//fi('mism_vs_values', ndiff(xb, v, ntok))// &
                cost(m0, m1, ntok))
      do kk = 1, 3
        m0 = now(); call k_fast(buf, ts, tl, xf, nfb, nfail, kk); m1 = now()
        call emit(tag//' mode='//trim(FKN(kk))//fi('tokens', ntok)//fi('fallback', nfb)// &
                  fr('fallback_share', real(nfb, R8)/real(ntok, R8))//fi('fail', nfail)// &
                  fi('mism_vs_rtl_token', ndiff(xf, xr, ntok))//fi('mism_vs_values', ndiff(xf, v, ntok))// &
                  cost(m0, m1, ntok))
      enddo
      m0 = now(); call k_strtod(buf, ts, tl, xs, npart); m1 = now()
      call emit(tag//' mode=strtod'//fi('tokens', ntok)//fi('partial', npart)// &
                fi('mism_vs_rtl_token', ndiff(xs, xr, ntok))//fi('mism_vs_values', ndiff(xs, v, ntok))// &
                cost(m0, m1, ntok))
      deallocate(xr, xb, xf, xs)
    enddo
  end subroutine randsets

  ! ---------------------------------------------------------------- census
  subroutine census(list)
    character(len=*), intent(in) :: list
    character(len=4096) :: path
    character(len=:), allocatable :: buf
    real(R8), allocatable :: xr(:), xf(:), xs(:)
    integer(I8), allocatable :: ts(:)
    integer, allocatable :: tl(:)
    integer(I8) :: n, ntok, k, nfast, nonplain, nfail, mfast, mstr, npart, ncr, lline, n87, m87, ng, mg, tt(13)
    integer :: u, ios, longest, nreg, nzone, nfile
    logical :: ok
    tt = 0; nfile = 0
    open(newunit=u, file=list, status='old', action='read')
    do
      read(u, '(A)', iostat=ios) path
      if (ios /= 0) exit
      if (len_trim(path) == 0) cycle
      call slurp(trim(path), buf, n)
      call tokenize(buf, n, ts, tl, ntok, longest, lline, nreg, nzone)
      allocate(xr(max(ntok, 1_I8)), xf(max(ntok, 1_I8)), xs(max(ntok, 1_I8)))
      nfast = 0; nonplain = 0; nfail = 0; mfast = 0; mstr = 0; npart = 0; n87 = 0; m87 = 0; ng = 0; mg = 0
      do k = 1, ntok
        call rtl_token(buf(ts(k):ts(k)+tl(k)-1), xr(k), ios)
        if (ios /= 0) nfail = nfail + 1
        if (.not. plain_token(buf(ts(k):ts(k)+tl(k)-1))) nonplain = nonplain + 1
        call fast_token(buf(ts(k):ts(k)+tl(k)-1), xf(k), ok)
        if (ok) then
          nfast = nfast + 1
          if (transfer(xf(k), 0_I8) /= transfer(xr(k), 0_I8)) mfast = mfast + 1
        endif
        call fastg_token(buf(ts(k):ts(k)+tl(k)-1), xf(k), ok)
        if (ok) then
          ng = ng + 1
          if (transfer(xf(k), 0_I8) /= transfer(xr(k), 0_I8)) mg = mg + 1
        endif
        call x87_token(buf(ts(k):ts(k)+tl(k)-1), xf(k), ok)
        if (ok) then
          n87 = n87 + 1
          if (transfer(xf(k), 0_I8) /= transfer(xr(k), 0_I8)) m87 = m87 + 1
        endif
        call strtod_token(buf, ts(k), tl(k), xs(k), ok)
        if (.not. ok) npart = npart + 1
        if (ok .and. transfer(xs(k), 0_I8) /= transfer(xr(k), 0_I8)) mstr = mstr + 1
      enddo
      ncr = count_char(buf, n, achar(13))
      call emit('IPROBE step=census'//fs('file', path)//fi('bytes', n)//fi('zones', int(nzone, I8))// &
                fi('regions', int(nreg, I8))//fi('tokens', ntok)//fi('fast', nfast)//fi('fallback', ntok - nfast)// &
                fr('fallback_share', real(ntok - nfast, R8)/real(max(ntok, 1_I8), R8))//fi('nonplain', nonplain)// &
                fi('rtl_fail', nfail)//fi('mism_fast', mfast)//fi('fastg', ng)//fi('mism_fastg', mg)// &
                fi('x87', n87)//fi('mism_x87', m87)// &
                fi('strtod_partial', npart)//fi('mism_strtod', mstr)// &
                fi('longest_token', int(longest, I8))//fi('longest_line', lline)//fi('cr', ncr))
      nfile = nfile + 1
      tt = tt + [n, ntok, nfast, ntok - nfast, nonplain, nfail, mfast, npart, mstr, n87, m87, ng, mg]
      deallocate(xr, xf, xs)
    enddo
    close(u)
    call emit('IPROBE step=census file=ALL'//fi('files', int(nfile, I8))//fi('bytes', tt(1))//fi('tokens', tt(2))// &
              fi('fast', tt(3))//fi('fallback', tt(4))//fr('fallback_share', real(tt(4), R8)/real(max(tt(2), 1_I8), R8))// &
              fi('nonplain', tt(5))//fi('rtl_fail', tt(6))//fi('mism_fast', tt(7))//fi('fastg', tt(12))// &
              fi('mism_fastg', tt(13))//fi('x87', tt(10))// &
              fi('mism_x87', tt(11))//fi('strtod_partial', tt(8))//fi('mism_strtod', tt(9)))
  end subroutine census

  integer(I8) function count_char(buf, n, c)
    character(len=*), intent(in) :: buf
    integer(I8), intent(in) :: n
    character, intent(in) :: c
    integer(I8) :: p
    integer :: ic
    count_char = 0
    ic = iachar(c)
    do p = 1, n
      if (iachar(buf(p:p)) == ic) count_char = count_char + 1
    enddo
  end function count_char

  ! ---------------------------------------------------------------- scale
  subroutine scale(path, ntokmax)
    character(len=*), intent(in) :: path
    integer(I8), intent(in) :: ntokmax
    character(len=:), allocatable :: buf, tag
    character(len=16) :: spec
    integer(I8), allocatable :: ts(:), rs(:), re(:), pa(:), pb(:), cnt(:), first(:), nfbp(:), fbi(:), fbs(:)
    integer, allocatable :: tl(:), fbl(:), cpus(:)
    real(R8), allocatable :: x(:), xref(:), xbr(:), xgr(:), mesh(:,:,:,:), vars(:,:,:,:), mref(:,:,:,:), vref(:,:,:,:)
    character(len=24), allocatable :: rec(:), rref(:)
    integer(I8) :: n, ntok, lline, npc, ip, p, q, k, nfb, mism, nfail
    integer :: longest, nreg, nzone, r, isp, T, ni, nj, nk, ndir, nvar, start, ios, kf
    logical :: spread, havelayout
    type(meter) :: m0, m1, m2, m3
    real(R8) :: t1(7)

    call slurp(path, buf, n)
    call find_regions(buf, n, rs, re, nreg, nzone)
    ! pieces of about PIECE_BYTES at line ends, over the data regions
    allocate(pa(16), pb(16)); npc = 0
    do r = 1, nreg
      p = rs(r)
      do while (p <= re(r))
        q = min(re(r), p + PIECE_BYTES - 1)
        if (q < re(r)) q = min(re(r), next_line_end(buf, q, re(r)))
        call push(p, q)
        p = q + 1
      enddo
    enddo
    allocate(cnt(npc), first(npc), nfbp(npc))
    do ip = 1, npc
      cnt(ip) = count_tokens(buf, pa(ip), pb(ip))
    enddo
    ntok = sum(cnt)
    if (ntokmax > 0 .and. ntok > ntokmax) then
      k = 0
      do ip = 1, npc
        if (k + cnt(ip) > ntokmax) exit
        k = k + cnt(ip)
      enddo
      npc = ip - 1; ntok = k
    endif
    first(1) = 0
    do ip = 2, npc
      first(ip) = first(ip-1) + cnt(ip-1)
    enddo
    ! the box's layout (one zone, BLOCK, cell-centred) when the token count is complete
    call layout_of(buf, n, ni, nj, nk, ndir, nvar, start, havelayout)
    if (havelayout) havelayout = (ntok == int(ndir, I8)*ni*nj*nk + &
      int(nvar, I8)*(ni - start)*(nj - start)*(max(start, nk - 1) - start + 1))
    tag = 'IPROBE step=threads'//fs('file', path)
    call emit(tag//' mode=setup'//fi('tokens', ntok)//fi('pieces', npc)//fi('piece_bytes', PIECE_BYTES)// &
              fs('layout', merge('box ', 'none', havelayout)))
    allocate(x(ntok), xref(ntok), xbr(ntok), xgr(ntok), fbi(ntok), fbs(ntok), fbl(ntok), cpus(0:omp_get_max_threads()-1))
    allocate(rec(ntok), rref(ntok))
    if (havelayout) then
      allocate(mesh(1:ndir, 0:ni-1, 0:nj-1, 0:nk-1), mref(1:ndir, 0:ni-1, 0:nj-1, 0:nk-1))
      allocate(vars(1:nvar, start:max(start, ni-1), start:max(start, nj-1), start:max(start, nk-1)))
      allocate(vref(1:nvar, start:max(start, ni-1), start:max(start, nj-1), start:max(start, nk-1)))
    endif

    do isp = 4, command_argument_count()
      call get_command_argument(isp, spec)
      spread = (index(spec, 's') > 0)
      read(spec(1:verify(spec, '0123456789') - 1), *) T
      ! fast and fastg: pass A (count), pass B (tokenize and convert into the slots), the fallbacks after the region
      do kf = 1, 2
        cpus = -1
        m0 = now()
        call par_count(T, spread)
        m1 = now()
        call par_fast(T, spread, kf)
        m2 = now()
        nfb = 0; nfail = 0
        do ip = 1, npc
          do k = first(ip) + 1, first(ip) + nfbp(ip)
            call rtl_token(buf(fbs(k):fbs(k)+fbl(k)-1), x(fbi(k)), ios)
            if (ios /= 0) nfail = nfail + 1
          enddo
          nfb = nfb + nfbp(ip)
        enddo
        m3 = now()
        if (kf == 1 .and. isp == 4) xref = x
        if (kf == 2 .and. isp == 4) xgr = x
        if (kf == 1) then
          mism = ndiff(x, xref, ntok)
        else
          mism = ndiff(x, xgr, ntok)
        endif
        call emit(tag//' kernel='//trim(FKN(kf))//fi('threads', int(T, I8))// &
                  fs('bind', merge('spread', 'close ', spread))// &
                  fr('passA_s', m1%w - m0%w)//fr('passB_s', m2%w - m1%w)//fr('fallback_s', m3%w - m2%w)// &
                  fi('fallback', nfb)//fi('fail', nfail)//rate(m3%w - m0%w, isp, 5 + kf)// &
                  fi('mism_vs_1', mism)//cpuset())
      enddo
      ! rtl-batch: pass A, then one list-directed READ per piece
      cpus = -1
      m0 = now()
      call par_count(T, spread)
      call par_batch(T, spread, nfail)
      m1 = now()
      if (isp == 4) xbr = x
      call emit(tag//' kernel=rtlbatch'//fi('threads', int(T, I8))//fs('bind', merge('spread', 'close ', spread))// &
                fi('fail', nfail)//rate(m1%w - m0%w, isp, 2)//fi('mism_vs_1', ndiff(x, xbr, ntok))// &
                fi('mism_vs_fast', ndiff(x, xref, ntok))//fi('mism_vs_fastg', ndiff(x, xgr, ntok))//cpuset())
      ! store into ORION's layout
      if (havelayout) then
        cpus = -1
        m0 = now(); call par_store(T, spread); m1 = now()
        if (isp == 4) then
          mref = mesh; vref = vars
        endif
        mism = ndiff(mesh, mref, size(mesh, kind=I8)) + ndiff(vars, vref, size(vars, kind=I8))
        call emit(tag//' kernel=store'//fi('threads', int(T, I8))//fs('bind', merge('spread', 'close ', spread))// &
                  rate(m1%w - m0%w, isp, 3)//fi('mism_vs_1', mism)//cpuset())
      endif
      ! the runtime's E23.15E3, one WRITE per SLICE values
      cpus = -1
      m0 = now(); call par_fmt(T, spread, .false., nfb); m1 = now()
      if (isp == 4) rref = rec
      call emit(tag//' kernel=fmtbatch'//fi('threads', int(T, I8))//fs('bind', merge('spread', 'close ', spread))// &
                rate(m1%w - m0%w, isp, 4)//fi('mism_vs_1', rdiff(rec, rref))//cpuset())
      ! the lock-free E23.15E3
      cpus = -1
      m0 = now(); call par_fmt(T, spread, .true., nfb); m1 = now()
      call emit(tag//' kernel=fmtfast'//fi('threads', int(T, I8))//fs('bind', merge('spread', 'close ', spread))// &
                fi('fallback', nfb)//rate(m1%w - m0%w, isp, 5)//fi('mism_vs_fmtbatch_1', rdiff(rec, rref))//cpuset())
    enddo

  contains

    subroutine push(a, b)
      integer(I8), intent(in) :: a, b
      integer(I8), allocatable :: t(:)
      if (npc == size(pa)) then
        allocate(t(2*npc)); t(1:npc) = pa; call move_alloc(t, pa)
        allocate(t(2*npc)); t(1:npc) = pb; call move_alloc(t, pb)
      endif
      npc = npc + 1; pa(npc) = a; pb(npc) = b
    end subroutine push

    !> Mvalues/s, and the speed-up against the first spec's time of kernel kk.
    function rate(w, i, kk) result(s)
      real(R8), intent(in) :: w
      integer, intent(in) :: i, kk
      character(len=:), allocatable :: s
      if (i == 4) t1(kk) = w
      s = fr('wall_s', w)//fr('Mvalues_per_s', real(ntok, R8)/1.e6_R8/w)//fr('speedup', t1(kk)/w)
    end function rate

    !> The distinct CPUs the region's threads started on, and their count.
    function cpuset() result(s)
      character(len=:), allocatable :: s
      integer :: c, nd, lo
      logical :: seen(0:4095)
      seen = .false.
      do c = 0, size(cpus) - 1
        if (cpus(c) >= 0 .and. cpus(c) <= 4095) seen(cpus(c)) = .true.
      enddo
      nd = count(seen)
      s = fi('cpus', int(nd, I8))//' cpuset='
      c = 0
      do while (c <= 4095)
        if (seen(c)) then
          lo = c
          do while (c < 4095)
            if (.not. seen(c+1)) exit
            c = c + 1
          enddo
          if (lo == c) then
            s = s//trim(itoa(lo))//','
          else
            s = s//trim(itoa(lo))//'-'//trim(itoa(c))//','
          endif
        endif
        c = c + 1
      enddo
    end function cpuset

    subroutine par_count(T, spread)
      integer, intent(in) :: T
      logical, intent(in) :: spread
      integer(I8) :: i
      if (spread) then
        !$omp parallel num_threads(T) proc_bind(spread) default(shared) private(i)
        cpus(omp_get_thread_num()) = probe_cpu()
        !$omp do schedule(static)
        do i = 1, npc
          cnt(i) = count_tokens(buf, pa(i), pb(i))
        enddo
        !$omp end do
        !$omp end parallel
      else
        !$omp parallel num_threads(T) proc_bind(close) default(shared) private(i)
        cpus(omp_get_thread_num()) = probe_cpu()
        !$omp do schedule(static)
        do i = 1, npc
          cnt(i) = count_tokens(buf, pa(i), pb(i))
        enddo
        !$omp end do
        !$omp end parallel
      endif
      first(1) = 0
      do i = 2, npc
        first(i) = first(i-1) + cnt(i-1)
      enddo
    end subroutine par_count

    subroutine par_fast(T, spread, kind)
      integer, intent(in) :: T, kind
      logical, intent(in) :: spread
      integer(I8) :: i
      if (spread) then
        !$omp parallel num_threads(T) proc_bind(spread) default(shared) private(i)
        cpus(omp_get_thread_num()) = probe_cpu()
        !$omp do schedule(static)
        do i = 1, npc
          call fast_piece(buf, pa(i), pb(i), first(i), x, fbi, fbs, fbl, nfbp(i), kind)
        enddo
        !$omp end do
        !$omp end parallel
      else
        !$omp parallel num_threads(T) proc_bind(close) default(shared) private(i)
        cpus(omp_get_thread_num()) = probe_cpu()
        !$omp do schedule(static)
        do i = 1, npc
          call fast_piece(buf, pa(i), pb(i), first(i), x, fbi, fbs, fbl, nfbp(i), kind)
        enddo
        !$omp end do
        !$omp end parallel
      endif
    end subroutine par_fast

    subroutine par_batch(T, spread, nfail)
      integer, intent(in) :: T
      logical, intent(in) :: spread
      integer(I8), intent(out) :: nfail
      integer(I8) :: i, nf
      nfail = 0
      if (spread) then
        !$omp parallel num_threads(T) proc_bind(spread) default(shared) private(i, nf) reduction(+:nfail)
        cpus(omp_get_thread_num()) = probe_cpu()
        !$omp do schedule(static)
        do i = 1, npc
          call batch_piece(buf, pa(i), pb(i), first(i), cnt(i), x, nf)
          nfail = nfail + nf
        enddo
        !$omp end do
        !$omp end parallel
      else
        !$omp parallel num_threads(T) proc_bind(close) default(shared) private(i, nf) reduction(+:nfail)
        cpus(omp_get_thread_num()) = probe_cpu()
        !$omp do schedule(static)
        do i = 1, npc
          call batch_piece(buf, pa(i), pb(i), first(i), cnt(i), x, nf)
          nfail = nfail + nf
        enddo
        !$omp end do
        !$omp end parallel
      endif
    end subroutine par_batch

    subroutine par_store(T, spread)
      integer, intent(in) :: T
      logical, intent(in) :: spread
      integer :: kk
      if (spread) then
        !$omp parallel num_threads(T) proc_bind(spread) default(shared) private(kk)
        cpus(omp_get_thread_num()) = probe_cpu()
        !$omp do schedule(static)
        do kk = 0, nk - 1
          call store_plane(kk)
        enddo
        !$omp end do
        !$omp end parallel
      else
        !$omp parallel num_threads(T) proc_bind(close) default(shared) private(kk)
        cpus(omp_get_thread_num()) = probe_cpu()
        !$omp do schedule(static)
        do kk = 0, nk - 1
          call store_plane(kk)
        enddo
        !$omp end do
        !$omp end parallel
      endif
    end subroutine par_store

    !> Plane kk of the mesh and, when it is a cell plane, of the variables, from xref (ORION's loop order).
    subroutine store_plane(kk)
      integer, intent(in) :: kk
      integer :: d, s, i, j
      integer(I8) :: nn, nc, b0, nic, njc
      nn = int(ni, I8)*nj*nk
      do d = 1, ndir
        b0 = (d - 1)*nn + int(kk, I8)*ni*nj
        do j = 0, nj - 1; do i = 0, ni - 1
          mesh(d, i, j, kk) = xref(b0 + int(j, I8)*ni + i + 1)
        enddo; enddo
      enddo
      if (kk < start .or. kk > max(start, nk - 1)) return
      nic = ni - start; njc = nj - start
      nc = nic*njc*(max(start, nk - 1) - start + 1)
      do s = 1, nvar
        b0 = ndir*nn + (s - 1)*nc + int(kk - start, I8)*nic*njc
        do j = start, nj - 1; do i = start, ni - 1
          vars(s, i, j, kk) = xref(b0 + int(j - start, I8)*nic + (i - start) + 1)
        enddo; enddo
      enddo
    end subroutine store_plane

    subroutine par_fmt(T, spread, fast, nfb)
      integer, intent(in) :: T
      logical, intent(in) :: spread, fast
      integer(I8), intent(out) :: nfb
      integer(I8) :: i, lo, hi, nsl, nf
      nsl = (ntok + SLICE - 1)/SLICE
      nfb = 0
      if (spread) then
        !$omp parallel num_threads(T) proc_bind(spread) default(shared) private(i, lo, hi, nf) reduction(+:nfb)
        cpus(omp_get_thread_num()) = probe_cpu()
        !$omp do schedule(static)
        do i = 1, nsl
          lo = (i - 1)*SLICE + 1; hi = min(ntok, i*SLICE)
          if (fast) then
            call k_fmt_fast(xref(lo:hi), rec(lo:hi), nf); nfb = nfb + nf
          else
            call fmt_e23_rtl(xref(lo:hi), rec(lo:hi))
          endif
        enddo
        !$omp end do
        !$omp end parallel
      else
        !$omp parallel num_threads(T) proc_bind(close) default(shared) private(i, lo, hi, nf) reduction(+:nfb)
        cpus(omp_get_thread_num()) = probe_cpu()
        !$omp do schedule(static)
        do i = 1, nsl
          lo = (i - 1)*SLICE + 1; hi = min(ntok, i*SLICE)
          if (fast) then
            call k_fmt_fast(xref(lo:hi), rec(lo:hi), nf); nfb = nfb + nf
          else
            call fmt_e23_rtl(xref(lo:hi), rec(lo:hi))
          endif
        enddo
        !$omp end do
        !$omp end parallel
      endif
    end subroutine par_fmt

  end subroutine scale

  !> Pass B of the fast reader on one piece: its tokens converted into x(first+1:); the tokens off the fast path
  !  recorded at fbi/fbs/fbl(first+1:first+nfb) for the serial runtime conversion after the region.
  subroutine fast_piece(buf, a, b, first, x, fbi, fbs, fbl, nfb, kind)
    character(len=*), intent(in) :: buf
    integer(I8), intent(in) :: a, b, first
    real(R8), intent(inout) :: x(:)
    integer(I8), intent(inout) :: fbi(:), fbs(:)
    integer, intent(inout) :: fbl(:)
    integer(I8), intent(out) :: nfb
    integer, intent(in) :: kind
    integer(I8) :: p, s0, k
    logical :: ok
    nfb = 0; k = first; p = a
    do while (p <= b)
      if (is_sep(iachar(buf(p:p)))) then
        p = p + 1
        cycle
      endif
      s0 = p
      do while (p <= b)
        if (is_sep(iachar(buf(p:p)))) exit
        p = p + 1
      enddo
      k = k + 1
      if (kind == 1) then
        call fast_token(buf(s0:p-1), x(k), ok)
      else
        call fastg_token(buf(s0:p-1), x(k), ok)
      endif
      if (.not. ok) then
        nfb = nfb + 1
        fbi(first + nfb) = k; fbs(first + nfb) = s0; fbl(first + nfb) = int(p - s0)
      endif
    enddo
  end subroutine fast_piece

  !> Pass B of the batched reader on one piece: its cnt tokens by one list-directed READ of a copy whose line ends,
  !  tabs and commas are blanks, when every token is plain; else runs of plain tokens and the others alone.
  subroutine batch_piece(buf, a, b, first, cnt, x, nfail)
    character(len=*), intent(in) :: buf
    integer(I8), intent(in) :: a, b, first, cnt
    real(R8), intent(inout) :: x(:)
    integer(I8), intent(out) :: nfail
    character(len=:), allocatable :: cp
    integer(I8), allocatable :: ts(:)
    integer, allocatable :: tl(:)
    integer(I8) :: p, k, k1, L
    integer :: ios, longest
    logical :: allplain
    nfail = 0
    if (cnt == 0) return
    allocate(ts(cnt), tl(cnt))
    k = 0; longest = 0
    call fill_tokens(buf, a, b, ts, tl, k, longest)
    allplain = .true.
    do k = 1, cnt
      if (.not. plain_token(buf(ts(k):ts(k)+tl(k)-1))) then
        allplain = .false.
        exit
      endif
    enddo
    if (allplain) then
      L = b - a + 1
      allocate(character(len=L) :: cp)
      do p = 1, L
        if (is_sep(iachar(buf(a+p-1:a+p-1)))) then
          cp(p:p) = ' '
        else
          cp(p:p) = buf(a+p-1:a+p-1)
        endif
      enddo
      read(cp, *, iostat=ios) x(first+1:first+cnt)
      if (ios /= 0) nfail = cnt
      return
    endif
    k1 = 1
    do k = 1, cnt
      if (.not. plain_token(buf(ts(k):ts(k)+tl(k)-1))) then
        if (k > k1) call run(k1, k - 1)
        call rtl_token(buf(ts(k):ts(k)+tl(k)-1), x(first + k), ios)
        if (ios /= 0) nfail = nfail + 1
        k1 = k + 1
      endif
    enddo
    if (cnt >= k1) call run(k1, cnt)
  contains
    subroutine run(r1, r2)
      integer(I8), intent(in) :: r1, r2
      integer(I8) :: aa, bb, q
      aa = ts(r1); bb = ts(r2) + tl(r2) - 1
      if (allocated(cp)) deallocate(cp)
      allocate(character(len=bb-aa+1) :: cp)
      do q = 1, bb - aa + 1
        if (is_sep(iachar(buf(aa+q-1:aa+q-1)))) then
          cp(q:q) = ' '
        else
          cp(q:q) = buf(aa+q-1:aa+q-1)
        endif
      enddo
      read(cp, *, iostat=ios) x(first+r1:first+r2)
      if (ios /= 0) nfail = nfail + (r2 - r1 + 1)
    end subroutine run
  end subroutine batch_piece

  !> The first zone's I, J, K, the coordinate count (the plane rule aside) and the centring of the variables, from
  !  the VARIABLES list and the zone header as ORION reads them (one block, BLOCK packing only).
  subroutine layout_of(buf, n, ni, nj, nk, ndir, nvar, start, ok)
    character(len=*), intent(in) :: buf
    integer(I8), intent(in) :: n
    integer, intent(out) :: ni, nj, nk, ndir, nvar, start
    logical, intent(out) :: ok
    ok = .false.; ni = 0; nj = 0; nk = 0; ndir = 0; nvar = 0; start = 0
    call header_dims(buf, n, ni, nj, nk, nvar, start, ok)
    if (.not. ok) return
    ndir = 3
    if (nk == 1) ndir = 2
    if (nj == 1 .and. nk == 1) ndir = 1
    nvar = nvar - ndir
    ok = (nvar >= 0)
  end subroutine layout_of

  !> I, J, K of the first zone record and its continuation lines, the number of quoted names after the keyword
  !  VARIABLES, and start = 1 when the zone header lists CELLCENTERED; ok = .false. otherwise (unquoted names).
  subroutine header_dims(buf, n, ni, nj, nk, nnames, start, ok)
    character(len=*), intent(in) :: buf
    integer(I8), intent(in) :: n
    integer, intent(out) :: ni, nj, nk, nnames, start
    logical, intent(out) :: ok
    integer(I8) :: p, q, zp
    character(len=:), allocatable :: hdr, vl
    integer :: i
    ok = .false.; ni = 0; nj = 0; nk = 1; nnames = 0; start = 0
    p = 1; zp = 0; vl = ''
    do while (p <= n)
      q = next_line_end(buf, p, n)
      if (zone_line(buf(p:q-1))) then
        zp = p
        exit
      endif
      vl = vl//' '//buf(p:q-1)
      p = q + 1
    enddo
    if (zp == 0) return
    hdr = buf(p:q-1)
    p = q + 1
    do while (p <= n)
      q = next_line_end(buf, p, n)
      if (numeric_start(buf(p:q-1))) exit
      hdr = hdr//' '//buf(p:q-1)
      p = q + 1
    enddo
    call upcase(hdr)
    ni = int_key(hdr, 'I'); nj = int_key(hdr, 'J'); nk = int_key(hdr, 'K')
    if (nk <= 0) nk = 1
    if (index(hdr, 'CELLCENTERED') > 0) start = 1
    call upcase(vl)
    i = index(vl, 'VARIABLES')
    if (i == 0) return
    nnames = count_quotes(vl(i:))/2
    ok = (ni > 0 .and. nj > 0 .and. nk > 0 .and. nnames > 0)
  end subroutine header_dims

  subroutine upcase(s)
    character(len=*), intent(inout) :: s
    integer :: i
    do i = 1, len(s)
      if (s(i:i) >= 'a' .and. s(i:i) <= 'z') s(i:i) = achar(iachar(s(i:i)) - 32)
    enddo
  end subroutine upcase

  integer function count_quotes(s)
    character(len=*), intent(in) :: s
    integer :: i
    count_quotes = 0
    do i = 1, len(s)
      if (s(i:i) == '"') count_quotes = count_quotes + 1
    enddo
  end function count_quotes

  !> The integer after "key=" (key a whole word outside quotes, blanks allowed around '='); 0 when absent.
  integer function int_key(s, key)
    character(len=*), intent(in) :: s, key
    integer :: p, q, ios
    logical :: inq
    int_key = 0; inq = .false.
    do p = 1, len(s) - len(key) + 1
      if (s(p:p) == '"') inq = .not. inq
      if (inq .or. s(p:p+len(key)-1) /= key) cycle
      if (p > 1) then
        if (scan(s(p-1:p-1), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_') > 0) cycle
      endif
      q = p + len(key)
      do while (q <= len(s))
        if (s(q:q) /= ' ') exit
        q = q + 1
      enddo
      if (q > len(s)) return
      if (s(q:q) /= '=') cycle
      read(s(q+1:min(len(s), q+20)), *, iostat=ios) int_key
      if (ios /= 0) int_key = 0
      return
    enddo
  end function int_key

  function itoa(i) result(s)
    integer, intent(in) :: i
    character(len=16) :: s
    write(s, '(I0)') i
  end function itoa

  ! ---------------------------------------------------------------- swrite
  subroutine swrite(path, nbytes, blk)
    character(len=*), intent(in) :: path
    integer(I8), intent(in) :: nbytes, blk
    character(len=:), allocatable :: big
    character(len=24), parameter :: PAT = ' 0.123456789012345E+001'//achar(10)
    integer(I8) :: p, q
    integer :: u
    type(meter) :: m0, m1, m2
    allocate(character(len=nbytes) :: big)
    do p = 1, nbytes, 24
      q = min(nbytes, p + 23)
      big(p:q) = PAT(1:q-p+1)
    enddo
    m0 = now()
    open(newunit=u, file=path, access='stream', form='unformatted', status='replace', action='write')
    p = 1
    do while (p <= nbytes)
      q = min(nbytes, p + blk - 1)
      write(u) big(p:q)
      p = q + 1
    enddo
    m1 = now()
    close(u)
    m2 = now()
    call emit('IPROBE step=stream'//fs('file', path)//fi('bytes', nbytes)//fi('block', blk)// &
              fr('write_s', m1%w - m0%w)//fr('close_s', m2%w - m1%w)//cost(m0, m2, 0_I8)// &
              fr('GBps', real(nbytes, R8)/1.e9_R8/(m2%w - m0%w)))
  end subroutine swrite

  ! ---------------------------------------------------------------- bcast
  subroutine bcast(ni, nj, nk, nvar, reps)
    integer, intent(in) :: ni, nj, nk, nvar, reps
#ifdef USE_MPI
    type(orion_data) :: orion
    real(R8) :: t0, dt
    real(R8), allocatable :: dts(:)
    integer(I8) :: h, nbytes
    integer(I8), allocatable :: hs(:)
    integer :: r, ierr, i, j, k, s
    call mpi_init_env()
    allocate(dts(0:mpi_size_-1), hs(0:mpi_size_-1))
    if (mpi_is_root) then
      allocate(orion%varnames(3 + nvar))
      orion%varnames = 'v'
      allocate(orion%block(1))
      orion%block(1)%name = 'Block1'
      orion%block(1)%Ni = ni - 1; orion%block(1)%Nj = nj - 1; orion%block(1)%Nk = nk - 1
      allocate(orion%block(1)%mesh(1:3, 0:ni-1, 0:nj-1, 0:nk-1))
      allocate(orion%block(1)%vars(1:nvar, 1:ni-1, 1:nj-1, 1:nk-1))
      do k = 0, nk-1; do j = 0, nj-1; do i = 0, ni-1
        orion%block(1)%mesh(:, i, j, k) = [real(i, R8), real(j, R8), real(k, R8)]*1.e-3_R8
      enddo; enddo; enddo
      do k = 1, nk-1; do j = 1, nj-1; do i = 1, ni-1; do s = 1, nvar
        orion%block(1)%vars(s, i, j, k) = real(s, R8) + 1.e-6_R8*real(i + 7*j + 13*k, R8)
      enddo; enddo; enddo; enddo
    endif
    do r = 1, reps
      if (.not. mpi_is_root) then
        if (allocated(orion%block)) deallocate(orion%block)
        if (allocated(orion%varnames)) deallocate(orion%varnames)
      endif
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      t0 = MPI_Wtime()
      call mpi_bcast_orion(orion)
      dt = MPI_Wtime() - t0
      h = field_hash(orion)
      call MPI_Gather(dt, 1, MPI_DOUBLE_PRECISION, dts, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_Gather(h, 1, MPI_INTEGER8, hs, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, ierr)
      if (mpi_is_root) then
        nbytes = 8*(size(orion%block(1)%mesh, kind=I8) + size(orion%block(1)%vars, kind=I8))
        call emit('IPROBE step=bcast'//fi('ranks', int(mpi_size_, I8))//fi('rep', int(r, I8))//fi('bytes', nbytes)// &
                  fr('root_s', dts(0))//fr('max_s', maxval(dts))//fr('min_receiver_s', minval(dts(1:)))// &
                  fi('last_receiver', int(maxloc(dts(1:), 1), I8))//fi('hash_ok', merge(1_I8, 0_I8, all(hs == hs(0)))))
      endif
    enddo
    call mpi_finalize_env()
#else
    call emit('IPROBE step=bcast error=not-built-with-MPI'//fi('ni', int(ni + nj + nk + nvar + reps, I8)))
#endif
  end subroutine bcast

  !> An order-dependent 64-bit digest of every value of the field (bits, in storage order; rotate and xor).
  integer(I8) function field_hash(o)
    type(orion_data), intent(in) :: o
    integer(I8) :: i
    field_hash = 1469598103934665603_I8
    if (.not. allocated(o%block)) return
    call mix(o%block(1)%mesh, size(o%block(1)%mesh, kind=I8))
    if (allocated(o%block(1)%vars)) call mix(o%block(1)%vars, size(o%block(1)%vars, kind=I8))
  contains
    subroutine mix(a, n)
      integer(I8), intent(in) :: n
      real(R8), intent(in) :: a(*)
      do i = 1, n
        field_hash = ieor(ishftc(field_hash, 13), transfer(a(i), 0_I8))
      enddo
    end subroutine mix
  end function field_hash

end module io_probe_run

program io_probe
  use io_probe_run
  use io_probe_kern, only: R8, probe_cyc_open, probe_fma
  implicit none


  character(len=4096) :: mode, a1, a2, a3, a4, a5
  integer :: ier

  call get_command_argument(1, mode)
  call get_command_argument(2, a1)
  call get_command_argument(3, a2)
  call get_command_argument(4, a3)
  call get_command_argument(5, a4)
  call get_command_argument(6, a5)
  ier = probe_cyc_open()
  if (probe_fma(0.1_R8, 1.e15_R8, -1.e14_R8) /= 6103515625._R8*2._R8**(-40)) &
    call emit('IPROBE step=selftest fma=WRONG (the lock-free formatter is not exact)')
  select case (trim(mode))
  case ('read-split')
    call read_split(trim(a1))
  case ('randsets')
    call randsets(trim(a1))
  case ('census')
    call census(trim(a1))
  case ('scale')
    call scale(trim(a1), ival8(a2))
  case ('swrite')
    call swrite(trim(a1), ival8(a2), ival8(a3))
  case ('bcast')
    call bcast(int(ival8(a1)), int(ival8(a2)), int(ival8(a3)), int(ival8(a4)), int(ival8(a5)))
  case default
    write(*,'(A)') 'usage: io-probe read-split|randsets|census|scale|swrite|bcast ... (see the source header)'
    error stop 2
  end select
end program io_probe
