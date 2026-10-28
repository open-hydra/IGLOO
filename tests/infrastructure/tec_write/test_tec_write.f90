!> The grid files' stream writer (IGLOO_Mod_TecASCII's tec_write_ascii_fast) against ORION's ASCII writer in the same
!  executable: the orion_data write_outfield makes, written by both and compared byte for byte at every chunk size, a
!  longer file at the same name before every write of the stream writer.
program test_tec_write
  use, intrinsic :: iso_fortran_env, only: I8 => int64, R8 => real64
  use Lib_ORION_data, only: orion_data
  use Lib_Tecplot, only: tec_write_structured_multiblock
  use IGLOO_Mod_TecASCII, only: tec_write_ascii_fast
  implicit none

  integer, parameter :: LLEN = 200                             !> write_outfield's varnames length (IGLOO_variables llen)
  integer(I8), parameter :: CHUNKS(3) = [1_I8, 4096_I8, 0_I8]   !> values per stream WRITE; 0: the writer's default
  integer, parameter :: EXTRA = 4096                           !> bytes the junk file holds beyond ORION's file
  character, parameter :: LF = achar(10)
  integer :: nfail = 0, ncmp = 0
  type(orion_data) :: o
  character(len=LLEN) :: vn

  ! source.tec of one material, one block, every field above 4096 values
  call make(o, [40], [30], [25], 5, -10._R8)
  vn = source_names([character(len=16) :: 'A'])
  call check('f1-source', o, vn)
  call declines(o, vn)
  ! euler<fam>.tec over two blocks, a solution time of a time-accurate run
  call make(o, [9, 4], [7, 6], [5, 3], 6, 2.5_R8)
  vn = euler_names()
  call check('f2-euler-2blk', o, vn)
  ! source.tec of six materials, one cell along z: a VARIABLES line of 168 characters, ten variables
  call make(o, [30], [20], [1], 10, -10._R8)
  vn = source_names([character(len=16) :: 'Al2O3_particles', 'H2O_droplets', 'n-decane_drops', 'n-heptane_drops', &
                     'kerosene_jp10', 'aluminium_bulk'])
  call check('f3-source-nz1', o, vn)

  write(*,'(A,I0,A,I0,A)') ' test_tec_write: ', ncmp, ' writes compared, ', nfail, ' failure(s)'
  if (nfail > 0) then
    write(*,'(A)') ' [FAIL] the stream writer differs from ORION''s writer'
    error stop 1
  endif
  write(*,'(A)') ' [PASS] the stream writer writes ORION''s bytes'

contains

  subroutine fail(msg)
    character(len=*), intent(in) :: msg
    nfail = nfail + 1
    write(*,'(A)') '   FAIL '//trim(msg)
  end subroutine fail

  function itoa(k) result(s)
    integer(I8), intent(in) :: k
    character(len=:), allocatable :: s
    character(len=24) :: t
    write(t,'(I0)') k
    s = trim(t)
  end function itoa

  !> The varnames write_outfield builds for source.tec (its expressions, its length).
  function source_names(mat) result(v)
    character(len=*), intent(in) :: mat(:)
    character(len=LLEN) :: v
    integer :: i
    v = ''
    do i = 1, size(mat)
      v = trim(v)//'"'//trim('wdot(')//trim(mat(i))//')"'
    enddo
    v = trim(v)//trim('"Fx" "Fy" "Fz" "E"')
  end function source_names

  !> The varnames write_outfield builds for euler<fam>.tec.
  function euler_names() result(v)
    character(len=LLEN) :: v
    v = ''
    v = trim(v)//trim('"rho<sub>p" "u<sub>p" "v<sub>p" "w<sub>p" "T<sub>p" "n<sub>p"')
  end function euler_names

  !> An orion_data as write_outfield makes it: blocks of nx x ny x nz cells, the mesh over the nodes from index 0, nvar
  !  cell-centred variables, the solution time; every value from val.
  subroutine make(o, nx, ny, nz, nvar, time)
    type(orion_data), intent(out) :: o
    integer, intent(in) :: nx(:), ny(:), nz(:), nvar
    real(R8), intent(in) :: time
    integer :: b, i, j, k, d
    integer(I8) :: t
    o%tec%node = .false.
    o%tec%bc = .false.
    o%tec%format = 'ascii'
    o%solutiontime = time
    allocate(o%block(size(nx)))
    t = 0
    do b = 1, size(nx)
      associate(blk => o%block(b))
        blk%name = 'Block'//itoa(int(b, I8))
        blk%Ni = nx(b); blk%Nj = ny(b); blk%Nk = nz(b)
        allocate(blk%mesh(1:3, 0:nx(b), 0:ny(b), 0:nz(b)))
        allocate(blk%vars(1:nvar, 1:nx(b), 1:ny(b), 1:nz(b)))
        do k = 0, nz(b); do j = 0, ny(b); do i = 0, nx(b); do d = 1, 3
          t = t + 1; blk%mesh(d, i, j, k) = val(t)
        enddo; enddo; enddo; enddo
        do k = 1, nz(b); do j = 1, ny(b); do i = 1, nx(b); do d = 1, nvar
          t = t + 1; blk%vars(d, i, j, k) = val(t)
        enddo; enddo; enddo; enddo
      end associate
    enddo
  end subroutine make

  !> Fixture value t: the bits of two 31-bit draws of t (every exponent, both signs, NaN and Inf among them); every
  !  seventh value one of the specials.
  real(R8) function val(t)
    integer(I8), intent(in) :: t
    integer(I8) :: a, c
    real(R8) :: sp(21)
    if (mod(t, 7_I8) == 0) then
      call specials(sp)
      val = sp(int(mod(t/7, 21_I8)) + 1)
      return
    endif
    a = mod(1103515245_I8*mod(t, 2147483647_I8) + 12345_I8, 2147483648_I8)
    c = mod(69069_I8*mod(7919_I8*t + a, 2147483647_I8) + 1_I8, 2147483648_I8)
    val = transfer(ior(ishft(a, 33), ior(ishft(c, 2), mod(t, 4_I8))), 1._R8)
  end function val

  !> Signed zeros, the subnormal ends, tiny, huge, infinities, quiet NaNs of both signs, a carry to the next decade.
  subroutine specials(x)
    real(R8), intent(out) :: x(21)
    integer(I8), parameter :: INF = int(z'7FF0000000000000', I8), QNAN = int(z'7FF8000000000000', I8)
    integer(I8), parameter :: SUBMAX = int(z'000FFFFFFFFFFFFF', I8)
    x(1) = transfer(0_I8, 1._R8)
    x(2) = transfer(ibset(0_I8, 63), 1._R8)
    x(3) = transfer(1_I8, 1._R8)
    x(4) = transfer(ibset(1_I8, 63), 1._R8)
    x(5) = transfer(SUBMAX, 1._R8)
    x(6) = transfer(ibset(SUBMAX, 63), 1._R8)
    x(7) = tiny(1._R8)
    x(8) = huge(1._R8)
    x(9) = -huge(1._R8)
    x(10) = transfer(INF, 1._R8)
    x(11) = transfer(ibset(INF, 63), 1._R8)
    x(12) = transfer(QNAN, 1._R8)
    x(13) = transfer(ibset(QNAN, 63), 1._R8)
    x(14) = 1._R8
    x(15) = -1._R8
    x(16) = 0.1_R8
    x(17) = nearest(1._R8, -1._R8)
    x(18) = 9.9999999999999999e299_R8
    x(19) = 1.e-300_R8
    x(20) = -1.e300_R8
    x(21) = -1.e-310_R8
  end subroutine specials

  !> The bytes of file f.
  function slurp(f) result(b)
    character(len=*), intent(in) :: f
    character(len=:), allocatable :: b
    integer(I8) :: n
    integer :: u
    open(newunit=u, file=f, access='stream', form='unformatted', action='read', status='old')
    inquire(unit=u, size=n)
    allocate(character(len=n) :: b)
    if (n > 0) read(u) b
    close(u)
  end function slurp

  !> n bytes of junk lines at f.
  subroutine junk(f, n)
    character(len=*), intent(in) :: f
    integer(I8), intent(in) :: n
    character(len=:), allocatable :: b
    integer(I8) :: k
    integer :: u
    allocate(character(len=n) :: b)
    do k = 1, n
      b(k:k) = 'J'
      if (mod(k, 64_I8) == 0) b(k:k) = LF
    enddo
    open(newunit=u, file=f, access='stream', form='unformatted', action='write', status='replace')
    write(u) b
    close(u)
  end subroutine junk

  !> got against ref: their lengths, then the first byte that differs and its line.
  subroutine same_bytes(tag, ref, got)
    character(len=*), intent(in) :: tag, ref, got
    integer(I8) :: k, n
    n = min(len(ref, kind=I8), len(got, kind=I8))
    do k = 1, n
      if (ref(k:k) /= got(k:k)) then
        call fail(tag//': byte '//itoa(k)//' (line '//itoa(count_lf(ref(1:k)) + 1)//') differs')
        return
      endif
    enddo
    if (len(ref) /= len(got)) call fail(tag//': '//itoa(len(got, kind=I8))//' bytes against ORION''s '// &
                                        itoa(len(ref, kind=I8)))
  end subroutine same_bytes

  integer(I8) function count_lf(s)
    character(len=*), intent(in) :: s
    integer(I8) :: k
    count_lf = 0
    do k = 1, len(s, kind=I8)
      if (s(k:k) == LF) count_lf = count_lf + 1
    enddo
  end function count_lf

  !> ORION's writer on o into orion-<tag>.tec, the stream writer on o at every chunk size into fast-<tag>.tec over a
  !  longer junk file, every fast file against ORION's byte for byte.
  subroutine check(tag, o, vn)
    character(len=*), intent(in) :: tag
    type(orion_data), intent(in) :: o
    character(len=*), intent(in) :: vn
    character(len=:), allocatable :: ref, got, t
    integer :: c, rc
    rc = tec_write_structured_multiblock(orion=o, varnames=vn, filename='orion-'//tag//'.tec')
    ref = slurp('orion-'//tag//'.tec')
    do c = 1, size(CHUNKS)
      t = tag//' chunk '//itoa(CHUNKS(c))
      call junk('fast-'//tag//'.tec', len(ref, kind=I8) + EXTRA)
      if (CHUNKS(c) > 0) then
        rc = tec_write_ascii_fast(o, vn, 'fast-'//tag//'.tec', chunk_values=CHUNKS(c))
      else
        rc = tec_write_ascii_fast(o, vn, 'fast-'//tag//'.tec')
      endif
      ncmp = ncmp + 1
      if (rc /= 0) then
        call fail(t//': the stream writer returns '//itoa(int(rc, I8)))
        cycle
      endif
      got = slurp('fast-'//tag//'.tec')
      call same_bytes(t, ref, got)
    enddo
  end subroutine check

  !> What write_outfield never makes is declined: nothing opened, the file at the name left as it was.
  subroutine declines(o, vn)
    type(orion_data), intent(in) :: o
    character(len=*), intent(in) :: vn
    type(orion_data) :: d
    d = o
    d%tec%node = .true.
    call decline_case('nodal variables', d, vn)
    d = o
    d%tec%bc = .true.
    call decline_case('ghost cells', d, vn)
    d = o
    d%tec%format = 'binary'
    call decline_case('binary', d, vn)
    d = o
    deallocate(d%block(1)%vars)
    call decline_case('mesh only', d, vn)
  end subroutine declines

  subroutine decline_case(tag, d, vn)
    character(len=*), intent(in) :: tag, vn
    type(orion_data), intent(in) :: d
    character(len=:), allocatable :: got
    integer :: rc
    call junk('decline.tec', 1000_I8)
    rc = tec_write_ascii_fast(d, vn, 'decline.tec')
    ncmp = ncmp + 1
    if (rc /= 1) call fail(tag//': not declined (returns '//itoa(int(rc, I8))//')')
    got = slurp('decline.tec')
    if (len(got) /= 1000 .or. got(1:1) /= 'J') call fail(tag//': the file at the name was touched')
  end subroutine decline_case

end program test_tec_write
