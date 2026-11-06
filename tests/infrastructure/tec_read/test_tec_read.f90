!> The gas file's one-pass reader (IGLOO_Mod_TecASCII) against ORION's reader in the same executable: every tracked
!  .tec input of the suite and a synthetic set written here, every value compared bit for bit, at every chunk size and
!  thread count.
!    test_tec_read --list <file>                the suite's files (the list CMake writes at configure time), then the
!                                               synthetic set
!    test_tec_read [--threads n,m,...] <file>   the given files only, read once by ORION and once by the reader at each
!                                               thread count (default: the environment's), timed
program test_tec_read
  use, intrinsic :: iso_fortran_env, only: I8 => int64, R8 => real64
  use, intrinsic :: iso_c_binding, only: c_int
  use Lib_ORION_data, only: orion_data
  use Lib_Tecplot, only: tec_read_structured_multiblock, tec_write_structured_multiblock
  use IGLOO_Mod_TecASCII, only: tec_read_ascii_fast
  !$ use omp_lib, only: omp_get_max_threads
  implicit none

  interface
    !> The runtime's C entry, the count passed by value.
    subroutine set_num_threads_c(n) bind(c, name='omp_set_num_threads')
      import :: c_int
      integer(c_int), value :: n
    end subroutine set_num_threads_c
  end interface

  integer, parameter :: SUITE_FILES = 72      !> the tracked .tec inputs of tests/: a new one moves this number
  integer(I8), parameter :: CHUNKS(4) = [64_I8, 4096_I8, 1048576_I8, 0_I8]   !> 0: the reader's default
  integer, parameter :: THREADS(4) = [1, 2, 3, 7]
  integer, parameter :: REPEATS = 5           !> the synthetic set's reads at 7 threads, per chunk size
  character, parameter :: LF = achar(10)
  integer :: nfail = 0, nread = 0, nrep7 = 1, nthr = 1
  integer, allocatable :: tlist(:)
  character(len=4096) :: arg
  character(len=:), allocatable :: txt
  integer :: i

  !$ nthr = size(THREADS)
  call get_command_argument(1, arg)
  if (trim(arg) == '--list') then
    call get_command_argument(2, arg)
    call suite(trim(arg))
    nrep7 = REPEATS
    call synthetic()
  else
    i = 1
    tlist = [0]
    if (trim(arg) == '--threads') then
      call get_command_argument(2, arg)
      call parse_threads(trim(arg))
      i = 3
    endif
    do while (i <= command_argument_count())
      call get_command_argument(i, arg)
      call timed(trim(arg))
      i = i + 1
    enddo
  endif
  write(*,'(A,I0,A,I0,A)') ' test_tec_read: ', nread, ' reads compared, ', nfail, ' failure(s)'
  if (nfail > 0) then
    write(*,'(A)') ' [FAIL] the reader differs from ORION'
    error stop 1
  endif
  write(*,'(A)') ' [PASS] the reader returns ORION''s values on every file'

contains

  subroutine fail(msg)
    character(len=*), intent(in) :: msg
    nfail = nfail + 1
    write(*,'(A)') '   FAIL '//trim(msg)
  end subroutine fail

  !> ORION's reader as IGLOO calls it (read_TECsolfile).
  integer function orion_read(f, o) result(rc)
    character(len=*), intent(in) :: f
    type(orion_data), intent(inout) :: o
    o%tec%node = .false.; o%tec%bc = .false.; o%tec%format = 'ascii'
    rc = tec_read_structured_multiblock(orion=o, filename=f)
  end function orion_read

  integer function fast_read(f, o, chunk, why) result(rc)
    character(len=*), intent(in) :: f
    type(orion_data), intent(inout) :: o
    integer(I8), intent(in) :: chunk
    character(len=*), intent(out) :: why
    o%tec%node = .false.; o%tec%bc = .false.; o%tec%format = 'ascii'
    if (chunk > 0) then
      rc = tec_read_ascii_fast(o, f, chunk_bytes=chunk, why=why)
    else
      rc = tec_read_ascii_fast(o, f, why=why)
    endif
  end function fast_read

  !> Values whose bits differ, over n elements (-0 against +0 and NaN payloads count).
  integer(I8) function ndiff(a, b, n)
    integer(I8), intent(in) :: n
    real(R8), intent(in) :: a(*), b(*)
    integer(I8) :: k
    ndiff = 0
    do k = 1, n
      if (transfer(a(k), 0_I8) /= transfer(b(k), 0_I8)) ndiff = ndiff + 1
    enddo
  end function ndiff

  !> b against a: the names, the solution time, the strand, the centring, the blocks, their bounds and values.
  subroutine compare(tag, a, b)
    character(len=*), intent(in) :: tag
    type(orion_data), intent(in) :: a, b
    integer :: k
    integer(I8) :: nd
    character(len=:), allocatable :: m
    nread = nread + 1
    if (allocated(a%varnames) .neqv. allocated(b%varnames)) then
      call fail(tag//': variable names allocated in one only'); return
    endif
    if (allocated(a%varnames)) then
      if (size(a%varnames) /= size(b%varnames)) then
        call fail(tag//': a different number of variable names'); return
      endif
      if (any(a%varnames /= b%varnames)) then
        call fail(tag//': different variable names'); return
      endif
    endif
    if (transfer(a%solutiontime, 0_I8) /= transfer(b%solutiontime, 0_I8)) call fail(tag//': solution time')
    if (a%strandid /= b%strandid) call fail(tag//': STRANDID')
    if (a%tec%node .neqv. b%tec%node) call fail(tag//': the centring')
    if (allocated(a%block) .neqv. allocated(b%block)) then
      call fail(tag//': blocks allocated in one only'); return
    endif
    if (.not. allocated(a%block)) return
    if (size(a%block) /= size(b%block)) then
      call fail(tag//': a different number of blocks'); return
    endif
    do k = 1, size(a%block)
      m = tag//' block '//itoa(k)
      if (a%block(k)%Ni /= b%block(k)%Ni .or. a%block(k)%Nj /= b%block(k)%Nj .or. a%block(k)%Nk /= b%block(k)%Nk) &
        call fail(m//': Ni/Nj/Nk')
      if (allocated(a%block(k)%mesh) .neqv. allocated(b%block(k)%mesh)) then
        call fail(m//': mesh allocated in one only'); cycle
      endif
      if (allocated(a%block(k)%mesh)) then
        if (any(lbound(a%block(k)%mesh) /= lbound(b%block(k)%mesh)) .or. &
            any(ubound(a%block(k)%mesh) /= ubound(b%block(k)%mesh))) then
          call fail(m//': mesh bounds'); cycle
        endif
        nd = ndiff(a%block(k)%mesh, b%block(k)%mesh, size(a%block(k)%mesh, kind=I8))
        if (nd > 0) call fail(m//': '//itoa(int(nd))//' mesh values differ')
      endif
      if (allocated(a%block(k)%vars) .neqv. allocated(b%block(k)%vars)) then
        call fail(m//': variables allocated in one only'); cycle
      endif
      if (allocated(a%block(k)%vars)) then
        if (any(lbound(a%block(k)%vars) /= lbound(b%block(k)%vars)) .or. &
            any(ubound(a%block(k)%vars) /= ubound(b%block(k)%vars))) then
          call fail(m//': variable bounds'); cycle
        endif
        nd = ndiff(a%block(k)%vars, b%block(k)%vars, size(a%block(k)%vars, kind=I8))
        if (nd > 0) call fail(m//': '//itoa(int(nd))//' variable values differ')
      endif
    enddo
  end subroutine compare

  !> expect 'agree': ORION reads it -> the reader reads the same values; ORION refuses it -> so does the reader.
  !  'decline': the reader leaves it to ORION (1). 'written': the reader returns ref (ORION is not the oracle).
  subroutine check_file(f, expect, ref)
    character(len=*), intent(in) :: f, expect
    type(orion_data), intent(in), optional :: ref
    type(orion_data) :: a
    integer :: rco, c, t, r, nr
    rco = 0
    if (expect == 'agree') rco = orion_read(f, a)
    do t = 1, nthr
      call set_threads(THREADS(t))
      nr = 1
      if (THREADS(t) == 7) nr = nrep7
      do r = 1, nr
        do c = 1, size(CHUNKS)
          call check_chunk(f, expect, CHUNKS(c), THREADS(t), rco, a, ref)
        enddo
      enddo
    enddo
    call set_threads(1)
  end subroutine check_file

  subroutine set_threads(n)
    integer, intent(in) :: n
    !$ call set_num_threads_c(int(n, c_int))
    !$ if (omp_get_max_threads() /= n) call fail('the thread count is not '//itoa(n))
  end subroutine set_threads

  subroutine check_chunk(f, expect, chunk, nt, rco, a, ref)
    character(len=*), intent(in) :: f, expect
    integer(I8), intent(in) :: chunk
    integer, intent(in) :: nt, rco
    type(orion_data), intent(in) :: a
    type(orion_data), intent(in), optional :: ref
    type(orion_data) :: b
    integer :: rcf
    character(len=96) :: why
    character(len=:), allocatable :: tag
    tag = short(f)//' chunk '//itoa(int(chunk))//' threads '//itoa(nt)
    rcf = fast_read(f, b, chunk, why)
    select case (expect)
    case ('agree')
      if (rco == 0) then
        if (rcf /= 0) then
          call fail(tag//': ORION reads it, the reader returns '//itoa(rcf)//' ('//trim(why)//')')
        else
          call compare(tag, a, b)
        endif
      else
        nread = nread + 1
        if (rcf == 0) call fail(tag//': ORION refuses it, the reader reads it')
      endif
    case ('decline')
      nread = nread + 1
      if (rcf /= 1) call fail(tag//': not declined (returns '//itoa(rcf)//')')
    case ('written')
      if (rcf /= 0) then
        call fail(tag//': the reader returns '//itoa(rcf)//' ('//trim(why)//')')
      else
        call compare(tag, ref, b)
      endif
    end select
  end subroutine check_chunk

  !> The path from its tests/ directory on (the suite's files), else as given.
  function short(f) result(s)
    character(len=*), intent(in) :: f
    character(len=:), allocatable :: s
    integer :: k
    k = index(f, '/tests/', back=.true.)
    if (k > 0) then
      s = f(k+1:)
    else
      s = f
    endif
  end function short

  function itoa(k) result(s)
    integer, intent(in) :: k
    character(len=:), allocatable :: s
    character(len=16) :: t
    write(t,'(I0)') k
    s = trim(t)
  end function itoa

  ! ---------------------------------------------------------------- the suite's files
  subroutine suite(list)
    character(len=*), intent(in) :: list
    character(len=4096) :: f
    integer :: u, ios, n
    open(newunit=u, file=list, status='old', action='read')
    n = 0
    do
      read(u,'(A)', iostat=ios) f
      if (ios /= 0) exit
      if (len_trim(f) == 0) cycle
      n = n + 1
      call scan_bytes(trim(f))
      call check_file(trim(f), 'agree')
    enddo
    close(u)
    write(*,'(A,I0,A,I0)') ' suite files: ', n, ', expected ', SUITE_FILES
    if (n /= SUITE_FILES) call fail('the number of tracked .tec inputs changed: update SUITE_FILES (and INFO.md)')
  end subroutine suite

  !> A suite file must be one ORION reads right: no line of 1000 characters or more, no CR (a POINT zone, which the
  !  reader leaves to ORION, fails check_file).
  subroutine scan_bytes(f)
    character(len=*), intent(in) :: f
    character(len=:), allocatable :: b
    integer(I8) :: n, p, a
    integer :: u
    open(newunit=u, file=f, access='stream', form='unformatted', action='read', status='old')
    inquire(unit=u, size=n)
    allocate(character(len=n) :: b)
    read(u) b
    close(u)
    a = 1
    do p = 1, n
      if (iachar(b(p:p)) == 13) then
        call fail(f//': a CR byte'); exit
      endif
      if (iachar(b(p:p)) == 10) then
        if (p - a >= 1000) call fail(f//': a line of 1000 characters or more')
        a = p + 1
      endif
    enddo
  end subroutine scan_bytes

  ! ---------------------------------------------------------------- the synthetic set
  subroutine add(line)
    character(len=*), intent(in) :: line
    txt = txt//line//LF
  end subroutine add

  subroutine save_file(f, final_lf, crlf)
    character(len=*), intent(in) :: f
    logical, intent(in), optional :: final_lf, crlf
    character(len=:), allocatable :: out
    integer :: u, k
    out = txt
    if (present(crlf)) then
      if (crlf) then
        out = ''
        do k = 1, len(txt)
          if (txt(k:k) == LF) then
            out = out//achar(13)//LF
          else
            out = out//txt(k:k)
          endif
        enddo
      endif
    endif
    if (present(final_lf)) then
      if (.not. final_lf .and. len(out) > 0) out = out(1:len(out)-1)
    endif
    open(newunit=u, file=f, access='stream', form='unformatted', action='write', status='replace')
    write(u) out
    close(u)
    txt = ''
  end subroutine save_file

  !> The k-th value of a deterministic sequence of both signs over 1e-19 .. 1e18: an integer times a power of two,
  !  exact (no rounding a compiler could contract differently at two call sites).
  real(R8) function val(k)
    integer, intent(in) :: k
    val = scale(real(1000003 + mod(37*k, 999983), R8), mod(13*k, 121) - 83)
    if (mod(k, 3) == 0) val = -val
  end function val

  function es(x) result(s)
    real(R8), intent(in) :: x
    character(len=:), allocatable :: s
    character(len=32) :: t
    write(t,'(ES25.17E3)') x
    s = trim(adjustl(t))
  end function es

  !> n values one per line, the format chosen by k (ES, F, list-directed, G).
  subroutine add_values(k0, n)
    integer, intent(in) :: k0, n
    integer :: k
    character(len=40) :: t
    do k = k0, k0 + n - 1
      select case (mod(k, 4))
      case (0)
        write(t,'(ES25.17E3)') val(k)
      case (1)
        write(t,'(F24.12)') 10._R8*fraction(val(k))
      case (2)
        write(t,*) val(k)
      case default
        write(t,'(ES15.7)') val(k)
      end select
      call add('  '//trim(t))
    enddo
  end subroutine add_values

  subroutine synthetic()
    type(orion_data) :: ref
    integer :: k, n
    character(len=2048) :: long
    character(len=:), allocatable :: s
    character(len=40), allocatable :: tok(:)

    txt = ''
    ! three zones of different sizes
    call add(' VARIABLES = "X" "Y" "Z" "A" "B"')
    call add(' ZONE T="z1", I=3, J=2, K=2, DATAPACKING=BLOCK, VARLOCATION=([4-5]=CELLCENTERED)')
    call add_values(1, 3*12 + 2*2)
    call add(' ZONE T="z2", I=4, J=3, K=2, DATAPACKING=BLOCK, VARLOCATION=([4-5]=CELLCENTERED)')
    call add_values(100, 3*24 + 2*6)
    call add(' ZONE T="z3", I=2, J=2, K=3, DATAPACKING=BLOCK, VARLOCATION=([4-5]=CELLCENTERED)')
    call add_values(300, 3*12 + 2*2)
    call save_file('s01-three-zones.tec')
    call check_file('s01-three-zones.tec', 'agree')

    ! MOSE's layout: a title, names with commas, the zone header over three lines, list-directed values
    call add(' TITLE="Soluzione tempo =    3.31800141305266       sec."')
    call add(' VARIABLES = "X", "Y", "Z","Roi( 1 )","U"')
    call add(' Zone T="block: 1",I=3,J=3,K=  2')
    call add(' ,DATAPACKING=BLOCK,')
    call add('  VARLOCATION=([4-5]=CELLCENTERED)')
    do k = 1, 3*18 + 2*4
      write(long,*) val(k)
      call add(trim(long)//'     ')
    enddo
    call save_file('s02-mose-header.tec')
    call check_file('s02-mose-header.tec', 'agree')

    ! separators: tabs, commas, several values per line, blank lines
    call add(' VARIABLES = "x","y","z","p"')
    call add(' ZONE I=2,J=2,K=2, VARLOCATION=([4]=CELLCENTERED)')
    call add('1.0,2.0'//achar(9)//'3.0   4.0')
    call add('')
    call add(' 5.0 ,6.0 , 7.0'//achar(9)//achar(9)//'8.0,')
    call add(',9.0 10.0 11.0 12.0')
    call add('')
    call add('')
    call add('13.0 14.0 15.0 16.0 17.0 18.0 19.0 20.0 21.0 22.0 23.0 24.0 25.0')
    call save_file('s03-separators.tec')
    call check_file('s03-separators.tec', 'agree')

    ! token forms: exponent letters, signs, points, zeros, subnormals, NaN, mantissas around 2**53, powers of ten
    allocate(tok(0))
    tok = [character(len=40) :: '1.5D+02', '2.0d-3', '1.0D0', '12', '-7', '+3.5', '+.5e+1', '.5', '5.', '-0.0', &
           '0.0', '-0', '0', '1e-310', '2.2250738585072014E-308', '4.9E-324', 'NaN', '12345678901234567', &
           '9.0071992547409931E+15', '9007199254740992', '9007199254740993', '123456789012345678', &
           '1.7976931348623157E+308', '-1.7976931348623157E+308', '3*1.5', '7e5', '7E+05', '0.000123', '1200', &
           '1.000000000000000E+01', '1.487282881060780E-008', '2.816068501015889E-008', '-.25e-0', '1.E2']
    do k = 0, 22
      write(long,'(A,I0)') '1e', k
      tok = [character(len=40) :: tok, trim(long)]
      write(long,'(A,I0)') '3.25e-', k
      tok = [character(len=40) :: tok, trim(long)]
    enddo
    n = 4*2*2*int((size(tok) + 15)/16)
    call add(' VARIABLES = "x" "y" "z" "p"')
    write(long,'(A,I0,A)') ' ZONE I=2, J=2, K=', n/16, ', DATAPACKING=BLOCK'
    call add(trim(long))
    do k = 1, n
      if (k <= size(tok)) then
        call add(' '//trim(tok(k)))
      else
        call add(' '//es(val(k)))
      endif
    enddo
    call save_file('s04-token-forms.tec')
    call check_file('s04-token-forms.tec', 'agree')

    ! an overflowing token
    call add(' VARIABLES = "x" "y" "z"')
    call add(' ZONE I=2, J=2, K=2')
    call add_values(1, 11)
    call add(' 1e400')
    call add_values(20, 12)
    call save_file('s05-overflow.tec')
    call check_file('s05-overflow.tec', 'agree')

    ! a mesh-only file
    call add(' VARIABLES = "x" "y" "z"')
    call add(' ZONE I=3, J=2, K=2')
    call add_values(1, 36)
    call save_file('s06-mesh-only.tec')
    call check_file('s06-mesh-only.tec', 'agree')

    ! a 2-D file (K = 1: two coordinates)
    call add(' VARIABLES = "x" "y" "p" "q"')
    call add(' ZONE I=4, J=3, VARLOCATION=([3-4]=CELLCENTERED)')
    call add_values(1, 2*12 + 2*6)
    call save_file('s07-2d.tec')
    call check_file('s07-2d.tec', 'agree')

    ! a plane of a 3-D field: x, y, z with z nodal
    call add(' VARIABLES = "x" "y" "z" "p"')
    call add(' ZONE I=4, J=3, K=1, VARLOCATION=([4]=CELLCENTERED)')
    call add_values(1, 3*12 + 6)
    call save_file('s08-plane-xyz.tec')
    call check_file('s08-plane-xyz.tec', 'agree')

    ! POINT packing: left to ORION
    call add(' VARIABLES = "x" "y" "z" "p"')
    call add(' ZONE I=2, J=2, K=2, DATAPACKING=POINT')
    call add_values(1, 32)
    call save_file('s09-point.tec')
    call check_file('s09-point.tec', 'decline')

    ! extra tokens after a zone's values, an extra line and an auxiliary line between zones
    call add(' VARIABLES = "x" "y" "z" "p"')
    call add(' ZONE I=2, J=2, K=2, VARLOCATION=([4]=CELLCENTERED)')
    call add_values(1, 24)
    call add(' '//es(val(25))//' 98.5 99.5')
    call add(' 1.0 2.0 3.0')
    call add(' DATASETAUXDATA Common.Time="1"')
    call add(' ZONE I=3, J=2, K=2, VARLOCATION=([4]=CELLCENTERED)')
    call add_values(50, 36 + 2)
    call save_file('s10-extra-aux.tec')
    call check_file('s10-extra-aux.tec', 'agree')

    ! a data line of about 2000 characters: ORION cuts it at 1000, the reader reads it whole (the written values)
    call add(' VARIABLES = "x" "y" "z" "p"')
    call add(' ZONE I=4, J=4, K=5')
    s = ''
    do k = 1, 80
      s = s//' '//es(val(k))
    enddo
    call add(s)
    do k = 81, 320
      call add('  '//es(val(k)))
    enddo
    call save_file('s11-long-line.tec')
    allocate(ref%varnames(4))
    ref%varnames = [character(len=32) :: 'x', 'y', 'z', 'p']
    ref%solutiontime = -10._R8
    ref%strandid = -1
    ref%tec%node = .true.
    allocate(ref%block(1))
    ref%block(1)%Ni = 3; ref%block(1)%Nj = 3; ref%block(1)%Nk = 4
    allocate(ref%block(1)%mesh(1:3, 0:3, 0:3, 0:4), ref%block(1)%vars(1:1, 0:3, 0:3, 0:4))
    call fill_written(ref)
    call check_file('s11-long-line.tec', 'written', ref)

    ! a VARIABLES line of 1 500 characters, and a VARIABLES list over three lines
    s = ' VARIABLES = "x" "y" "z"'
    do k = 1, 60
      write(long,'(A,I2.2,A)') ' "a_long_variable_name_', k, '"'
      s = s//trim(long)
    enddo
    call add(s)
    call add(' ZONE I=2, J=2, K=2, VARLOCATION=([4-63]=CELLCENTERED)')
    call add_values(1, 24 + 60)
    call save_file('s12-long-variables.tec')
    call check_file('s12-long-variables.tec', 'agree')
    call add(' VARIABLES = "x"')
    call add(' "y"')
    call add('  "z", "p"')
    call add(' ZONE I=2, J=2, K=2, VARLOCATION=([4]=CELLCENTERED)')
    call add_values(1, 25)
    call save_file('s13-variables-lines.tec')
    call check_file('s13-variables-lines.tec', 'agree')

    ! a file written by ORION's own ASCII writer, two blocks
    call orion_file('s14-orion-written.tec')
    call check_file('s14-orion-written.tec', 'agree')

    ! a last line without a line feed; lower case keywords
    call add(' variables = "x" "y" "z"')
    call add(' zone i=3, j=2, k=2')
    call add_values(1, 36)
    call save_file('s15-no-final-lf.tec', final_lf=.false.)
    call check_file('s15-no-final-lf.tec', 'agree')

    ! null values and CR line ends: left to ORION
    call add(' VARIABLES = "x" "y" "z"')
    call add(' ZONE I=2, J=2, K=2')
    call add_values(1, 10)
    call add(' /')
    call add_values(20, 13)
    call save_file('s16-null-slash.tec')
    call check_file('s16-null-slash.tec', 'decline')
    call add(' VARIABLES = "x" "y" "z"')
    call add(' ZONE I=2, J=2, K=2')
    call add_values(1, 10)
    call add(' 3*')
    call add_values(20, 13)
    call save_file('s17-null-repeat.tec')
    call check_file('s17-null-repeat.tec', 'decline')
    call add(' VARIABLES = "x" "y" "z"')
    call add(' ZONE I=2, J=2, K=2')
    call add_values(1, 24)
    call save_file('s18-crlf.tec', crlf=.true.)
    call check_file('s18-crlf.tec', 'decline')

    ! a zone that ends early (ORION refuses), a token of 120 characters (ORION refuses)
    call add(' VARIABLES = "x" "y" "z"')
    call add(' ZONE I=2, J=2, K=2')
    call add_values(1, 19)
    call add(' ZONE I=2, J=2, K=2')
    call add_values(30, 24)
    call save_file('s19-short-zone.tec')
    call check_file('s19-short-zone.tec', 'agree')
    call add(' VARIABLES = "x" "y" "z"')
    call add(' ZONE I=2, J=2, K=2')
    call add_values(1, 10)
    call add(' 1'//repeat('0', 119))
    call add_values(20, 13)
    call save_file('s20-long-token.tec')
    call check_file('s20-long-token.tec', 'agree')
  end subroutine synthetic

  !> The values of s11 as ORION converts each written token alone (the runtime's list-directed READ, which is not
  !  correctly rounded on every 18-digit token): the coordinates, then p, over the nodes, i fastest.
  subroutine fill_written(o)
    type(orion_data), intent(inout) :: o
    integer :: d, i, j, k, t
    t = 0
    do d = 1, 3
      do k = 0, 4; do j = 0, 3; do i = 0, 3
        t = t + 1; o%block(1)%mesh(d, i, j, k) = as_read(t)
      enddo; enddo; enddo
    enddo
    do k = 0, 4; do j = 0, 3; do i = 0, 3
      t = t + 1; o%block(1)%vars(1, i, j, k) = as_read(t)
    enddo; enddo; enddo
  end subroutine fill_written

  real(R8) function as_read(t)
    integer, intent(in) :: t
    character(len=100) :: tok
    tok = es(val(t))
    read(tok, *) as_read
  end function as_read

  !> Two cell-centred blocks written by ORION's writer as write_outfield calls it (E23.15E3, STRANDID = 0).
  subroutine orion_file(f)
    character(len=*), intent(in) :: f
    type(orion_data) :: o
    integer :: b, i, j, k, s, t, err
    o%tec%node = .false.; o%tec%bc = .false.; o%tec%format = 'ascii'
    o%solutiontime = -10._R8
    allocate(o%block(2))
    t = 0
    do b = 1, 2
      o%block(b)%name = 'Block'//itoa(b)
      o%block(b)%Ni = 2 + b; o%block(b)%Nj = 2; o%block(b)%Nk = b
      allocate(o%block(b)%mesh(1:3, 0:2+b, 0:2, 0:b), o%block(b)%vars(1:2, 1:2+b, 1:2, 1:b))
      do k = 0, b; do j = 0, 2; do i = 0, 2 + b
        t = t + 1; o%block(b)%mesh(:, i, j, k) = [val(t), val(t+1000), val(t+2000)]
      enddo; enddo; enddo
      do k = 1, b; do j = 1, 2; do i = 1, 2 + b; do s = 1, 2
        t = t + 1; o%block(b)%vars(s, i, j, k) = val(3*t + s)
      enddo; enddo; enddo; enddo
    enddo
    err = tec_write_structured_multiblock(orion=o, varnames='"a" "b"', filename=f)
    if (err /= 0) call fail(f//': ORION''s writer')
  end subroutine orion_file

  ! ---------------------------------------------------------------- one large file, timed
  !> n,m,... into tlist.
  subroutine parse_threads(t)
    character(len=*), intent(in) :: t
    integer :: p, q, n
    tlist = [integer ::]
    p = 1
    do while (p <= len(t))
      q = index(t(p:), ',')
      if (q == 0) then
        q = len(t) + 1
      else
        q = p + q - 1
      endif
      read(t(p:q-1), *) n
      tlist = [tlist, n]
      p = q + 1
    enddo
  end subroutine parse_threads

  subroutine timed(f)
    character(len=*), intent(in) :: f
    type(orion_data) :: a
    integer :: rco, k
    integer(I8) :: c0, c1, rate
    call system_clock(c0, rate)
    rco = orion_read(f, a)
    call system_clock(c1)
    write(*,'(A,A,A,I0,A,F10.3,A)') ' ', f, ': ORION rc ', rco, ' in ', real(c1 - c0, R8)/rate, ' s'
    do k = 1, size(tlist)
      call timed_fast(f, a, rco, tlist(k))
    enddo
  end subroutine timed

  !> One read by the reader at nt threads (0: the environment's), timed, against ORION's a.
  subroutine timed_fast(f, a, rco, nt)
    character(len=*), intent(in) :: f
    type(orion_data), intent(in) :: a
    integer, intent(in) :: rco, nt
    type(orion_data) :: b
    integer :: rcf, n
    integer(I8) :: c1, c2, rate
    character(len=96) :: why
    n = 1
    if (nt > 0) call set_threads(nt)
    !$ n = omp_get_max_threads()
    call system_clock(c1, rate)
    rcf = fast_read(f, b, 0_I8, why)
    call system_clock(c2)
    write(*,'(A,A,A,I0,A,I0,A,F10.3,A)') ' ', f, ': the reader at ', n, ' thread(s) rc ', rcf, ' in ', &
      real(c2 - c1, R8)/rate, ' s'
    if (rco /= 0) then
      nread = nread + 1
      if (rcf == 0) call fail(f//': ORION refuses it, the reader reads it')
    elseif (rcf /= 0) then
      call fail(f//': the reader returns '//itoa(rcf)//' ('//trim(why)//')')
    else
      call compare(f//' at '//itoa(n)//' thread(s)', a, b)
    endif
  end subroutine timed_fast

end program test_tec_read
