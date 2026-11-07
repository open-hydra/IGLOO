!> The gas file's one-pass reader and the grid files' stream writer. The reader reads every ASCII Tecplot file it
!  accepts with ORION's semantics, its values converted by the run-time library's list-directed READ as ORION converts
!  them, one READ per run of plain number tokens, on every OpenMP thread; a file it declines is left to ORION's reader.
!  The writer writes the bytes of ORION's ASCII writer in four steps (open, zone, field, close), its values formatted by
!  the run-time library's internal WRITE many at a time, on every OpenMP thread, and written by stream access.
!  IGLOO_TEC_FAST=0 in the environment brings back ORION's reader and writer. The header routines and expressions are
!  copies of ORION's (Lib_Tecplot, GPL-3.0); the reader's messages on the paths that end in a refusal are left to
!  ORION's reader, which reads every declined file.
module IGLOO_Mod_TecASCII
  use, intrinsic :: iso_fortran_env, only: I8 => int64, R8 => real64, stderr => error_unit
  use Lib_ORION_data, only: orion_data
  use IR_Precision, only: str, FR_P, DR_P
  !$ use omp_lib, only: omp_get_max_threads, omp_get_thread_num
  implicit none
  private
  public :: tec_read_ascii_fast, tec_write_ascii_fast, tec_fast_enabled
  public :: tec_grid_open, tec_grid_zone, tec_grid_field, tec_grid_close

  integer(I8), parameter :: CHUNK_DEFAULT = 268435456_I8   !> stream buffer [bytes]
  integer, parameter :: RUN_MAX = 65536                    !> plain tokens per list-directed READ
  integer, parameter :: TOK_MAX = 100                      !> ORION's token buffer
  integer, parameter :: HDR = 1000, VHDR = 32768           !> ORION's line buffers: zone and data lines, VARIABLES
  integer(I8), parameter :: NULL_BITS = int(z'7FF8DEAD0BEEF001', I8)   !> a NaN no text reads as: a null value
  integer(I8), parameter :: WCHUNK_DEFAULT = 4194304_I8    !> the writer's values per stream WRITE
  integer, parameter :: WSLICE = 4096                      !> the writer's values per internal WRITE
  integer, parameter :: WREC = DR_P + 1                    !> the writer's bytes per value: ORION's FR_P field, a line feed
  character, parameter :: LF = achar(10)

  !> The file seen through a buffer: b(lo:hi) unread, next the file position of the byte after hi.
  type :: stream_t
    integer :: u = -1
    integer(I8) :: size = 0, next = 1, lo = 1, hi = 0
    character(len=:), allocatable :: b
    logical :: eof = .false., ioerr = .false.
  end type stream_t

  !> A grid file being written: its unit (open or not), the I/O status, the values per stream WRITE and the chunk's text.
  type, public :: tec_grid_t
    integer :: u = 0, ios = 0, nvar = 0
    logical :: opened = .false.
    integer(I8) :: chunk = WCHUNK_DEFAULT
    real(R8) :: time = 0
    logical :: time_accurate = .true.
    character(len=500) :: varloc = ''
    character(len=WREC), allocatable :: rec(:)
  end type tec_grid_t

  !> The zone being read: its dimensions, ORION's count of its values and flat views of its arrays.
  type :: zone_t
    integer :: ni = 0, nj = 0, nk = 0, ndir = 0, nvar = 0, start = 0
    integer(I8) :: nn = 0, nc = 0, need = 0
    real(R8), pointer, contiguous :: m(:) => null(), v(:) => null()
  end type zone_t

contains

  !> .false. when IGLOO_TEC_FAST is 0: ORION's reader and writer as they are (read once).
  logical function tec_fast_enabled()
    logical, save :: known = .false., on = .true.
    character(len=16) :: val
    integer :: n, st
    if (.not. known) then
      call get_environment_variable('IGLOO_TEC_FAST', val, n, st)
      on = .not. (st == 0 .and. trim(adjustl(val)) == '0')
      known = .true.
    endif
    tec_fast_enabled = on
  end function tec_fast_enabled

  !> Reads filename into orion as ORION's tec_read_ascii does (the BLOCK zones of an ASCII Tecplot file; data lines of
  !  any length read whole): 0 read, 1 declined (the caller reads the file with ORION), 2 an I/O error. orion is
  !  untouched unless the read succeeds. chunk_bytes sizes the stream buffer; why tells a decline's reason.
  integer function tec_read_ascii_fast(orion, filename, chunk_bytes, why) result(stat)
    type(orion_data), intent(inout) :: orion
    character(len=*), intent(in) :: filename
    integer(I8), intent(in), optional :: chunk_bytes
    character(len=*), intent(out), optional :: why
    type(orion_data), target :: o
    type(stream_t) :: s
    character(len=96) :: reason
    integer(I8) :: cap
    integer :: ios

    stat = 2
    reason = 'the file cannot be opened or read'
    open(newunit=s%u, file=filename, access='stream', form='unformatted', action='read', status='old', iostat=ios)
    if (ios == 0) then
      inquire(unit=s%u, size=s%size)
      cap = CHUNK_DEFAULT
      if (present(chunk_bytes)) cap = max(chunk_bytes, 64_I8)
      allocate(character(len=cap) :: s%b)
      call read_file(s, o, stat, reason)
      if (s%ioerr) then
        stat = 2
        reason = 'the file cannot be read'
      endif
      close(s%u)
    endif
    if (present(why)) why = ''
    if (stat /= 0) then
      if (present(why)) why = reason
      return
    endif
    call move_alloc(o%varnames, orion%varnames)
    call move_alloc(o%block, orion%block)
    orion%solutiontime = o%solutiontime
    orion%strandid = o%strandid
    orion%tec%node = o%tec%node
  end function tec_read_ascii_fast

  !> Writes orion to filename as ORION's tec_write_structured_multiblock writes it in ASCII, byte for byte, for an
  !  orion_data writable accepts, through tec_grid_open, tec_grid_zone, tec_grid_field and tec_grid_close. 0 written,
  !  1 declined (nothing opened: ORION's writer takes it), 2 an I/O error (the unit closed).
  integer function tec_write_ascii_fast(orion, varnames, filename, chunk_values) result(stat)
    type(orion_data), intent(in) :: orion
    character(len=*), intent(in) :: varnames, filename
    integer(I8), intent(in), optional :: chunk_values
    type(tec_grid_t) :: g
    integer(I8) :: nmax
    integer :: b, d, k, ndir, nvar, ni1, ni2, nj1, nj2, nk1, nk2

    stat = 1
    if (.not. writable(orion)) return
    ndir = 3
    nvar = size(orion%block(1)%vars, 1)
    nmax = 1
    do b = 1, size(orion%block)
      nmax = max(nmax, size(orion%block(b)%mesh, kind=I8)/ndir)
    enddo
    stat = tec_grid_open(g, filename, varnames, nvar, orion%solutiontime, nmax, chunk_values)
    if (stat /= 0) return
    do b = 1, size(orion%block)
      associate(blk => orion%block(b))
        ni1 = 0; ni2 = size(blk%mesh, 2) - 1
        nj1 = 0; nj2 = size(blk%mesh, 3) - 1
        nk1 = 0; nk2 = size(blk%mesh, 4) - 1
        call tec_grid_zone(g, blk%name, ni2-ni1+1, nj2-nj1+1, nk2-nk1+1)
        do d = 1, ndir
          call tec_grid_field(g, blk%mesh(d, ni1:ni2, nj1:nj2, nk1:nk2))
        enddo
        if (ni2 == 0) ni2 = 1
        if (nj2 == 0) nj2 = 1
        if (nk2 == 0) nk2 = 1
        do k = 1, nvar
          call tec_grid_field(g, blk%vars(k, ni1+1:ni2, nj1+1:nj2, nk1+1:nk2))
        enddo
      end associate
    enddo
    stat = tec_grid_close(g)
  end function tec_write_ascii_fast

  !> Opens filename for a grid file of three coordinates and nvar cell-centred variables and writes ORION's VARIABLES
  !  line; nmax is the largest field's value count, solutiontime ORION's (negative: steady, STRANDID = 0). 0 opened,
  !  2 an I/O error (nothing left open).
  integer function tec_grid_open(g, filename, varnames, nvar, solutiontime, nmax, chunk_values) result(stat)
    type(tec_grid_t), intent(out) :: g
    character(len=*), intent(in) :: filename, varnames
    integer, intent(in) :: nvar
    real(R8), intent(in) :: solutiontime
    integer(I8), intent(in) :: nmax
    integer(I8), intent(in), optional :: chunk_values
    character(len=32768) :: tecvarname
    integer :: ndir, nvartot, ios

    if (present(chunk_values)) g%chunk = max(1_I8, chunk_values)
    ndir = 3
    g%nvar = nvar
    nvartot = ndir + nvar
    g%time_accurate = .true.
    g%time = solutiontime
    if (g%time < 0) then
      g%time = abs(g%time)
      g%time_accurate = .false.
    endif
    ! ORION's header text (tec_init, tec_blk_data: the ASCII branch)
    tecvarname = ' VARIABLES = x'
    tecvarname = trim(tecvarname)//' y'
    tecvarname = trim(tecvarname)//' z'
    tecvarname = trim(tecvarname)//' '//unquote(trim(varnames))
    if (nvar == 1) then
      g%varloc = ', VARLOCATION=([1-'//trim(str(.true.,ndir))//']=NODAL,['//trim(str(.true.,ndir+1))// &
                 ']=CELLCENTERED)'
    else
      g%varloc = ', VARLOCATION=([1-'//trim(str(.true.,ndir))//']=NODAL,['//trim(str(.true.,ndir+1))//'-'// &
                 trim(str(.true.,nvartot))//']=CELLCENTERED)'
    endif
    allocate(g%rec(min(g%chunk, nmax)))
    stat = 2
    g%ios = 2
    open(newunit=g%u, file=trim(filename), access='stream', form='unformatted', status='replace', action='write', &
         iostat=ios)
    if (ios /= 0) return
    write(g%u, iostat=g%ios) trim(tecvarname)//LF
    if (g%ios /= 0) then
      close(g%u, iostat=ios)
      return
    endif
    g%opened = .true.
    stat = 0
  end function tec_grid_open

  !> ORION's zone line of a block of ni x nj x nk nodes; nothing after an I/O error.
  subroutine tec_grid_zone(g, name, ni, nj, nk)
    type(tec_grid_t), intent(inout) :: g
    character(len=*), intent(in) :: name
    integer, intent(in) :: ni, nj, nk
    character(len=500) :: teczoneheader
    if (g%ios /= 0) return
    teczoneheader = ' ZONE  T = '//trim(name)//                     &
                    ', I='//trim(str(no_sign=.true.,n=ni))//        &
                    ', J='//trim(str(no_sign=.true.,n=nj))//        &
                    ', K='//trim(str(no_sign=.true.,n=nk))//        &
                    ', DATAPACKING=BLOCK'//adjustl(trim(g%varloc))
    teczoneheader = trim(teczoneheader)//', SOLUTIONTIME='//trim(str(no_sign=.true.,n=g%time))
    if (.not. g%time_accurate) teczoneheader = trim(teczoneheader)//', STRANDID = 0'
    write(g%u, iostat=g%ios) trim(teczoneheader)//LF
  end subroutine tec_grid_zone

  !> Field f (i fastest, at most the nmax values tec_grid_open was given) after its zone line, one FR_P value per
  !  line; nothing after an I/O error.
  subroutine tec_grid_field(g, f)
    type(tec_grid_t), intent(inout) :: g
    real(R8), intent(in) :: f(:,:,:)
    if (g%ios /= 0) return
    call write_field(g%u, f, g%chunk, g%rec, g%ios)
  end subroutine tec_grid_field

  !> Closes the grid file: 0 written, 2 an I/O error (the unit closed either way).
  integer function tec_grid_close(g) result(stat)
    type(tec_grid_t), intent(inout) :: g
    integer :: ios
    stat = 2
    if (.not. g%opened) return
    close(g%u, iostat=ios)
    g%opened = .false.
    if (allocated(g%rec)) deallocate(g%rec)
    if (g%ios == 0 .and. ios == 0) stat = 0
  end function tec_grid_close

  !> The orion_data tec_write_ascii_fast writes: ASCII, no ghost cells, cell-centred variables, three coordinates over
  !  nodes from index 0, the same variables in every block over cells from index 1 (one cell along a direction without
  !  cells, as ORION's writer loops).
  logical function writable(orion) result(ok)
    type(orion_data), intent(in) :: orion
    integer :: b, nvar
    ok = .false.
    if (trim(orion%tec%format) /= 'ascii' .or. orion%tec%node .or. orion%tec%bc) return
    if (.not. allocated(orion%block)) return
    if (size(orion%block) < 1) return
    if (.not. allocated(orion%block(1)%vars)) return
    nvar = size(orion%block(1)%vars, 1)
    if (nvar < 1) return
    do b = 1, size(orion%block)
      associate(blk => orion%block(b))
        if (.not. (allocated(blk%mesh) .and. allocated(blk%vars))) return
        if (size(blk%mesh, 1) /= 3 .or. any(lbound(blk%mesh) /= [1, 0, 0, 0])) return
        if (any(lbound(blk%vars) /= 1)) return
        if (any(shape(blk%vars) /= [nvar, max(1, size(blk%mesh, 2) - 1), max(1, size(blk%mesh, 3) - 1), &
                                    max(1, size(blk%mesh, 4) - 1)])) return
      end associate
    enddo
    ok = .true.
  end function writable

  !> Field f (i fastest, as ORION loops) to unit u, chunk values per stream WRITE, rec the chunk's text: each chunk cut
  !  into one piece per thread, each piece formatted on its thread; the pieces' errors joined after the region.
  subroutine write_field(u, f, chunk, rec, ios)
    integer, intent(in) :: u
    real(R8), intent(in) :: f(:,:,:)
    integer(I8), intent(in) :: chunk
    character(len=WREC), intent(inout) :: rec(:)
    integer, intent(out) :: ios
    integer(I8) :: n, c0, m, lo, hi, q
    integer :: nt, np, p, st, nerr
    n = size(f, kind=I8)
    nt = 1
    !$ nt = omp_get_max_threads()
    ios = 0
    c0 = 0
    do while (c0 < n)
      m = min(chunk, n - c0)
      np = int(min(int(nt, I8), m))
      nerr = 0
      !$omp parallel do schedule(static) default(shared) private(p, lo, hi, q, st) reduction(+:nerr) if(np > 1)
      do p = 1, np
        lo = (m*(p - 1))/np
        hi = (m*p)/np
        do q = lo, hi - 1, WSLICE
          call format_slice(f, c0 + q, min(int(WSLICE, I8), hi - q), rec(q+1:), st)
          if (st /= 0) nerr = nerr + 1
        enddo
      enddo
      !$omp end parallel do
      if (nerr > 0) then
        ios = 1
        return
      endif
      write(u, iostat=ios) rec(1:m)
      if (ios /= 0) return
      c0 = c0 + m
    enddo
  end subroutine write_field

  !> Values t0 .. t0+n-1 (0-based, i fastest) of f into rec(1:n): one internal WRITE of FR_P, then a line feed each.
  subroutine format_slice(f, t0, n, rec, ios)
    real(R8), intent(in) :: f(:,:,:)
    integer(I8), intent(in) :: t0, n
    character(len=WREC), intent(inout) :: rec(:)
    integer, intent(out) :: ios
    real(R8) :: v(WSLICE)
    integer(I8) :: ni, nj, nij, r
    integer :: i, j, k, q
    ni = size(f, 1, kind=I8); nj = size(f, 2, kind=I8); nij = ni*nj
    k = int(t0/nij); r = t0 - k*nij
    j = int(r/ni); i = int(r - j*ni)
    do q = 1, int(n)
      v(q) = f(i+1, j+1, k+1)
      i = i + 1
      if (i == ni) then
        i = 0; j = j + 1
        if (j == nj) then
          j = 0; k = k + 1
        endif
      endif
    enddo
    write(rec(1:n), FR_P, iostat=ios) v(1:n)
    do q = 1, int(n)
      rec(q)(WREC:WREC) = LF
    enddo
  end subroutine format_slice

  !> `"a""b"` -> ` a  b `: ORION's writer's quotes made blanks (Lib_Tecplot unquote).
  pure function unquote(sin) result(sout)
    character(len=*), intent(in) :: sin
    character(len=len(sin)) :: sout
    integer :: i
    do i = 1, len(sin)
      if (sin(i:i) == '"') then
        sout(i:i) = ' '
      else
        sout(i:i) = sin(i:i)
      endif
    enddo
  end function unquote

  !> The VARIABLES text, then every zone in file order: its header, ORION's checks, its arrays, its values.
  subroutine read_file(s, o, stat, reason)
    type(stream_t), intent(inout) :: s
    type(orion_data), target, intent(inout) :: o
    integer, intent(out) :: stat
    character(len=*), intent(inout) :: reason
    character(len=VHDR) :: vline, variables_header
    character(len=HDR) :: line, header
    integer(I8) :: a, e
    integer :: nb, ni, nj, nk, ndir, nvar, start, strandid, sid
    logical :: found, zone_point, zone_node, node1, z_cell, plane_xyz, have
    real(R8) :: solutiontime

    stat = 1
    ! VARIABLES: the first line naming it, joined with the non-blank lines that follow up to the first zone line
    found = .false.
    variables_header = ''
    do while (next_line(s, a, e))
      if (has_cr(s%b, a, e)) then
        reason = 'a CR byte'; return
      endif
      vline = s%b(a:e)
      if (index(upper_case(vline), 'VARIABLES') > 0) then
        variables_header = trim(vline)
        found = .true.
        exit
      endif
      if (is_zone_header(vline)) then
        reason = 'a zone record before the VARIABLES line'; return
      endif
    enddo
    if (.not. found) then
      reason = 'no VARIABLES line'; return
    endif
    have = .false.
    do while (next_line(s, a, e))
      if (has_cr(s%b, a, e)) then
        reason = 'a CR byte'; return
      endif
      vline = s%b(a:e)
      if (is_zone_header(vline)) then
        have = .true.
        exit
      endif
      if (len_trim(vline) > 0) variables_header = trim(variables_header)//' '//trim(vline)
    enddo
    call read_variables(variables_header, o%varnames)
    if (.not. allocated(o%varnames)) then
      reason = 'no variables in the VARIABLES text'; return
    endif
    if (.not. have) then
      reason = 'no zone record'; return
    endif

    solutiontime = -10._R8
    strandid = -1
    nb = 0
    allocate(o%block(4))
    zone: do
      ! the zone line in hand: its continuation lines up to the first line that starts with a number
      line = s%b(a:e)
      header = trim(line)
      found = .false.
      do while (next_line(s, a, e))
        if (has_cr(s%b, a, e)) then
          reason = 'a CR byte'; return
        endif
        line = s%b(a:e)
        if (line_is_numeric_start(line)) then
          found = .true.
          exit
        endif
        header = trim(header)//' '//trim(line)
      enddo
      call parse_zone_header(header, ni, nj, nk, zone_point, zone_node, solutiontime, sid)
      strandid = sid
      z_cell = var_cellcentered(upper_case(header), 3)
      if (ni <= 0 .or. nj <= 0 .or. nk <= 0) then
        reason = 'a zone without positive dimensions'; return
      endif
      if (.not. found) then
        reason = 'a zone without data'; return
      endif
      nb = nb + 1
      if (nb == 1) then
        node1 = zone_node
        ndir = 3
        if (nk == 1) ndir = 2
        if (nj == 1 .and. nk == 1) ndir = 1
        plane_xyz = .false.
        if (ndir == 2 .and. size(o%varnames) >= 3) then
          plane_xyz = is_coordinate_name(o%varnames(1), 'X') .and. is_coordinate_name(o%varnames(2), 'Y') .and. &
                      is_coordinate_name(o%varnames(3), 'Z') .and. .not. z_cell
          if (plane_xyz) ndir = 3
        endif
        if (size(o%varnames) < ndir) then
          reason = 'fewer variables than coordinates'; return
        endif
        nvar = size(o%varnames) - ndir
        start = 0
        if (.not. zone_node) start = 1
      endif
      if (zone_node .neqv. node1) then
        reason = 'zones of different centring'; return
      endif
      if (plane_xyz .and. z_cell) then
        reason = 'z cell-centred in a zone after a plane of x, y, z'; return
      endif
      if ((nk == 1 .and. ndir == 3 .and. .not. plane_xyz) .or. (nk > 1 .and. (ndir < 3 .or. plane_xyz)) .or. &
          (nj == 1 .and. ndir > 1)) then
        reason = 'a zone that does not fit the coordinates'; return
      endif
      if (zone_point) then
        reason = 'POINT packing'; return
      endif
      if (start == 1 .and. nvar > 0 .and. (ni < 2 .or. nj < 2)) then
        reason = 'a cell-centred zone one node wide'; return
      endif
      call grow_blocks(o, nb)
      call open_zone(o, nb, ni, nj, nk, ndir, nvar, start)
      ! the first data line in hand: the zone's values from its start
      s%lo = a
      call read_zone_values(s, o, nb, stat, reason)
      if (stat /= 0) return
      stat = 1
      ! the lines after the zone's values, up to the next zone line, are ignored
      have = .false.
      do while (next_line(s, a, e))
        line = s%b(a:e)
        if (is_zone_header(line)) then
          have = .true.
          exit
        endif
      enddo
      if (.not. have) exit zone
    enddo zone
    call shrink_blocks(o, nb)
    o%solutiontime = solutiontime
    o%strandid = strandid
    o%tec%node = node1
    stat = 0
  end subroutine read_file

  !> The next line b(a:e) without its line feed (the stream moves past it); .false. at the end of the file.
  logical function next_line(s, a, e) result(got)
    type(stream_t), intent(inout) :: s
    integer(I8), intent(out) :: a, e
    integer(I8) :: q
    got = .false.
    a = 0; e = -1
    do
      q = s%lo
      do while (q <= s%hi)
        if (iachar(s%b(q:q)) == 10) exit
        q = q + 1
      enddo
      if (q <= s%hi) then
        a = s%lo; e = q - 1; s%lo = q + 1
        got = .true.
        return
      endif
      if (s%eof) then
        if (s%hi >= s%lo) then
          a = s%lo; e = s%hi; s%lo = s%hi + 1
          got = .true.
        endif
        return
      endif
      call refill(s)
    enddo
  end function next_line

  !> Moves b(lo:hi) to the buffer's start and fills the rest from the file; a full buffer doubles first.
  subroutine refill(s)
    type(stream_t), intent(inout) :: s
    character(len=:), allocatable :: nb
    integer(I8) :: keep, m, i
    integer :: ios
    if (s%eof) return
    keep = s%hi - s%lo + 1
    if (keep > 0 .and. s%lo > 1) then
      do i = 1, keep
        s%b(i:i) = s%b(s%lo+i-1:s%lo+i-1)
      enddo
    endif
    s%lo = 1
    s%hi = max(keep, 0_I8)
    if (s%hi == len(s%b, kind=I8)) then
      allocate(character(len=2*len(s%b, kind=I8)) :: nb)
      nb(1:s%hi) = s%b(1:s%hi)
      call move_alloc(nb, s%b)
    endif
    m = min(len(s%b, kind=I8) - s%hi, s%size - s%next + 1)
    if (m > 0) then
      read(s%u, pos=s%next, iostat=ios) s%b(s%hi+1:s%hi+m)
      if (ios /= 0) then
        s%ioerr = .true.; s%eof = .true.
        return
      endif
      s%next = s%next + m
      s%hi = s%hi + m
    endif
    if (s%next > s%size) s%eof = .true.
  end subroutine refill

  logical function has_cr(b, a, e)
    character(len=*), intent(in) :: b
    integer(I8), intent(in) :: a, e
    integer(I8) :: p
    has_cr = .false.
    do p = a, e
      if (iachar(b(p:p)) == 13) then
        has_cr = .true.
        return
      endif
    enddo
  end function has_cr

  !> Room for block nb (the components of the blocks read so far move, no copy).
  subroutine grow_blocks(o, nb)
    type(orion_data), intent(inout) :: o
    integer, intent(in) :: nb
    type(orion_data) :: t
    integer :: i
    if (nb <= size(o%block)) return
    allocate(t%block(2*size(o%block)))
    do i = 1, nb - 1
      call move_block(o, i, t, i)
    enddo
    call move_alloc(t%block, o%block)
  end subroutine grow_blocks

  subroutine shrink_blocks(o, nb)
    type(orion_data), intent(inout) :: o
    integer, intent(in) :: nb
    type(orion_data) :: t
    integer :: i
    if (nb == size(o%block)) return
    allocate(t%block(nb))
    do i = 1, nb
      call move_block(o, i, t, i)
    enddo
    call move_alloc(t%block, o%block)
  end subroutine shrink_blocks

  subroutine move_block(o, i, t, j)
    type(orion_data), intent(inout) :: o, t
    integer, intent(in) :: i, j
    t%block(j)%Ni = o%block(i)%Ni; t%block(j)%Nj = o%block(i)%Nj; t%block(j)%Nk = o%block(i)%Nk
    if (allocated(o%block(i)%mesh)) call move_alloc(o%block(i)%mesh, t%block(j)%mesh)
    if (allocated(o%block(i)%vars)) call move_alloc(o%block(i)%vars, t%block(j)%vars)
  end subroutine move_block

  !> Block b's arrays with ORION's bounds, and the count of values ORION reads into them.
  subroutine open_zone(o, b, ni, nj, nk, ndir, nvar, start)
    type(orion_data), target, intent(inout) :: o
    integer, intent(in) :: b, ni, nj, nk, ndir, nvar, start
    o%block(b)%Ni = ni - 1; o%block(b)%Nj = nj - 1; o%block(b)%Nk = nk - 1
    allocate(o%block(b)%mesh(1:ndir, 0:ni-1, 0:nj-1, 0:nk-1))
    if (nvar > 0) allocate(o%block(b)%vars(1:nvar, start:max(start, ni-1), start:max(start, nj-1), &
                                             start:max(start, nk-1)))
  end subroutine open_zone

  !> Block b's values: the first need tokens from the stream's position, across every line, as ORION's next_value
  !  takes them; the stream stops after the line of the last one.
  subroutine read_zone_values(s, o, b, stat, reason)
    type(stream_t), intent(inout) :: s
    type(orion_data), target, intent(inout) :: o
    integer, intent(in) :: b
    integer, intent(out) :: stat
    character(len=*), intent(inout) :: reason
    type(zone_t) :: z
    integer(I8) :: got, re, pend, q
    real(R8), allocatable :: tmp(:,:)
    integer :: st, nt

    stat = 1
    z%ndir = size(o%block(b)%mesh, 1)
    z%ni = size(o%block(b)%mesh, 2); z%nj = size(o%block(b)%mesh, 3); z%nk = size(o%block(b)%mesh, 4)
    z%nn = int(z%ni, I8)*z%nj*z%nk
    z%m(1:size(o%block(b)%mesh, kind=I8)) => o%block(b)%mesh
    z%nvar = 0; z%nc = 0
    if (allocated(o%block(b)%vars)) then
      z%nvar = size(o%block(b)%vars, 1)
      z%start = lbound(o%block(b)%vars, 2)
      z%nc = int(z%ni - z%start, I8)*(z%nj - z%start)*(max(z%start, z%nk - 1) - z%start + 1)
      z%v(1:size(o%block(b)%vars, kind=I8)) => o%block(b)%vars
    endif
    z%need = z%ndir*z%nn + z%nvar*z%nc
    nt = 1
    !$ nt = omp_get_max_threads()
    allocate(tmp(RUN_MAX, nt))
    got = 0
    do while (got < z%need)
      ! the complete lines in the buffer (or its every byte at the end of the file)
      re = s%hi
      if (.not. s%eof) then
        do while (re >= s%lo)
          if (iachar(s%b(re:re)) == 10) exit
          re = re - 1
        enddo
      endif
      if (re < s%lo) then
        if (s%eof) then
          reason = 'fewer values than a zone holds'; return
        endif
        call refill(s)
        if (s%ioerr) return
        cycle
      endif
      call convert_region(s%b, s%lo, re, z, got, pend, tmp, st)
      if (st /= 0) then
        reason = decline_text(st); return
      endif
      if (got == z%need) then
        q = pend + 1
        do while (q <= s%hi)
          if (iachar(s%b(q:q)) == 10) exit
          q = q + 1
        enddo
        s%lo = q + 1
      else
        s%lo = re + 1
        if (.not. s%eof) call refill(s)
        if (s%ioerr) return
      endif
    enddo
    stat = 0
  end subroutine read_zone_values

  !> The tokens of b(a:e) into the zone's slots from got on, until the zone is full, on every thread: the span cut into
  !  one piece per thread at line ends, pass A counts each piece's tokens, their prefix sum gives each piece its first
  !  slot, pass B converts each piece's tokens (convert_span), the piece that fills the zone up to the zone's count.
  !  pend: the last byte of the last token taken; st /= 0: the decline of the first piece that declines; tmp(:,t): the
  !  scratch of thread t.
  subroutine convert_region(b, a, e, z, got, pend, tmp, st)
    character(len=*), intent(inout) :: b
    integer(I8), intent(in) :: a, e
    type(zone_t), intent(in) :: z
    integer(I8), intent(inout) :: got
    integer(I8), intent(out) :: pend
    real(R8), intent(inout) :: tmp(:,:)
    integer, intent(out) :: st
    integer(I8), allocatable :: pa(:), pb(:), cnt(:), t0(:), lim(:), nk(:), pe(:)
    integer, allocatable :: sp(:)
    integer(I8) :: q, nom, n1, e1
    integer :: np, k, m, tid

    np = size(tmp, 2)
    if (np == 1) then
      ! one thread: the tokens in order until the zone is full or the span ends
      call convert_span(b, a, e, z, got, z%need - got, n1, e1, tmp(:,1), st)
      pend = e1
      got = got + n1
      return
    endif
    allocate(pa(np), pb(np))
    ! the cuts: piece k ends at the line feed at or after a + k*len/np - 1
    m = 0
    q = a
    do k = 1, np
      if (q > e) exit
      m = m + 1
      pa(m) = q
      nom = max(q, a + ((e - a + 1)*k)/np - 1)
      if (k == np .or. nom >= e) then
        pb(m) = e
      else
        do while (nom < e)
          if (iachar(b(nom:nom)) == 10) exit
          nom = nom + 1
        enddo
        pb(m) = nom
      endif
      q = pb(m) + 1
    enddo
    np = m
    allocate(cnt(np), t0(np), lim(np), nk(np), pe(np), sp(np))
    !$omp parallel do schedule(static) default(shared) private(k) if(np > 1)
    do k = 1, np
      cnt(k) = count_tokens(b, pa(k), pb(k))
    enddo
    !$omp end parallel do
    t0(1) = got
    do k = 2, np
      t0(k) = t0(k-1) + cnt(k-1)
    enddo
    do k = 1, np
      lim(k) = max(0_I8, min(cnt(k), z%need - t0(k)))
    enddo
    sp = 0
    nk = 0
    pe = 0
    !$omp parallel do schedule(static) default(shared) private(k, tid) if(np > 1)
    do k = 1, np
      tid = 1
      !$ tid = omp_get_thread_num() + 1
      if (lim(k) > 0) call convert_span(b, pa(k), pb(k), z, t0(k), lim(k), nk(k), pe(k), tmp(:,tid), sp(k))
    enddo
    !$omp end parallel do
    st = 0
    pend = a - 1
    do k = 1, np
      if (sp(k) == 0 .and. nk(k) /= lim(k)) sp(k) = 5
      if (sp(k) /= 0) then
        st = sp(k)
        return
      endif
      if (lim(k) > 0) then
        got = got + lim(k)
        pend = pe(k)
      endif
    enddo
  end subroutine convert_region

  !> The number of tokens of b(a:e): the runs of bytes other than the separators and CR.
  pure function count_tokens(b, a, e) result(n)
    character(len=*), intent(in) :: b
    integer(I8), intent(in) :: a, e
    integer(I8) :: n, p
    integer :: c
    logical :: gap, cut

    n = 0
    gap = .true.
    do p = a, e
      c = iachar(b(p:p))
      cut = is_sep(c) .or. c == 13
      if (gap .and. .not. cut) n = n + 1
      gap = cut
    enddo
  end function count_tokens

  !> The tokens of b(a:e), at most nmax, into the zone's slots t0, t0+1, ...: runs of plain tokens by one list-directed
  !  READ each (their separators made blanks in place), any other token by ORION's own READ. ntok: the tokens taken;
  !  pend: the last byte of the last one; st /= 0: a decline (decline_text).
  subroutine convert_span(b, a, e, z, t0, nmax, ntok, pend, tmp, st)
    character(len=*), intent(inout) :: b
    integer(I8), intent(in) :: a, e, t0, nmax
    type(zone_t), intent(in) :: z
    integer(I8), intent(out) :: ntok, pend
    real(R8), intent(inout) :: tmp(:)
    integer, intent(out) :: st
    integer(I8) :: p, ts, ra, rb, got
    integer :: c, n, ios
    real(R8) :: x

    st = 0
    pend = a - 1
    got = t0
    n = 0; ra = 0; rb = 0
    p = a
    do while (p <= e .and. got + n < t0 + nmax)
      c = iachar(b(p:p))
      if (c == 13) then
        st = 1; exit
      endif
      if (is_sep(c)) then
        p = p + 1
        cycle
      endif
      ts = p
      do while (p <= e)
        c = iachar(b(p:p))
        if (is_sep(c) .or. c == 13) exit
        p = p + 1
      enddo
      if (p - ts > TOK_MAX) then
        st = 2; exit
      endif
      if (plain_token(b(ts:p-1))) then
        if (n == 0) ra = ts
        rb = p - 1
        n = n + 1
        if (n == size(tmp)) then
          call run(st)
          if (st /= 0) exit
        endif
      else
        call run(st)
        if (st /= 0) exit
        call token_read(b(ts:p-1), x, st)
        if (st /= 0) exit
        tmp(1) = x
        call put_values(z, got, tmp, 1)
        got = got + 1
        pend = p - 1
      endif
    enddo
    if (st == 0) call run(st)
    ntok = got - t0
  contains
    subroutine run(st)
      integer, intent(out) :: st
      integer(I8) :: q, q0, q1
      integer :: k
      st = 0
      if (n == 0) return
      do q = ra, rb
        if (is_sep(iachar(b(q:q)))) b(q:q) = ' '
      enddo
      read(b(ra:rb), *, iostat=ios) tmp(1:n)
      if (ios /= 0) then
        ! the run token by token, as ORION reads it
        k = 0; q = ra
        do while (q <= rb)
          if (b(q:q) == ' ') then
            q = q + 1
            cycle
          endif
          q0 = q
          do while (q <= rb)
            if (b(q:q) == ' ') exit
            q = q + 1
          enddo
          q1 = q - 1
          k = k + 1
          call token_read(b(q0:q1), tmp(k), st)
          if (st /= 0) return
        enddo
      endif
      call put_values(z, got, tmp, n)
      got = got + n
      pend = rb
      n = 0
    end subroutine run
  end subroutine convert_span

  !> One token as ORION reads it: a list-directed READ from a blank-padded character(100). st 3: the runtime refuses
  !  it, 4: it reads as a null value (ORION's value would be left unassigned).
  subroutine token_read(t, x, st)
    character(len=*), intent(in) :: t
    real(R8), intent(out) :: x
    integer, intent(out) :: st
    character(len=TOK_MAX) :: tok
    integer :: ios
    st = 0
    tok = ' '
    tok(1:len(t)) = t
    x = transfer(NULL_BITS, 1._R8)
    read(tok, *, iostat=ios) x
    if (ios /= 0) then
      st = 3
    elseif (transfer(x, 0_I8) == NULL_BITS) then
      st = 4
    endif
  end subroutine token_read

  function decline_text(st) result(t)
    integer, intent(in) :: st
    character(len=96) :: t
    select case (st)
    case (1)
      t = 'a CR byte'
    case (2)
      t = 'a token longer than 100 characters'
    case (3)
      t = 'a token the runtime cannot read'
    case (4)
      t = 'a null value (/, r*)'
    case (5)
      t = 'fewer tokens than counted'
    case default
      t = 'an internal state'
    end select
  end function decline_text

  !> Values v(1:n) into the slots of the zone's tokens t0 .. t0+n-1 (0-based): the coordinates field after field over
  !  the nodes, then the variables field after field over the cells, i fastest, as ORION's loops store them.
  subroutine put_values(z, t0, v, n)
    type(zone_t), intent(in) :: z
    integer(I8), intent(in) :: t0
    real(R8), intent(in) :: v(:)
    integer, intent(in) :: n
    integer(I8) :: t, f, d, nm
    integer :: k
    nm = z%ndir*z%nn
    t = t0
    k = 1
    do while (k <= n)
      if (t < nm) then
        d = t/z%nn; f = t - d*z%nn
        do while (k <= n .and. f < z%nn)
          z%m(z%ndir*f + d + 1) = v(k)
          k = k + 1; f = f + 1; t = t + 1
        enddo
      else
        d = (t - nm)/z%nc; f = t - nm - d*z%nc
        do while (k <= n .and. f < z%nc)
          z%v(z%nvar*f + d + 1) = v(k)
          k = k + 1; f = f + 1; t = t + 1
        enddo
      endif
    enddo
  end subroutine put_values

  pure logical function is_sep(c)
    integer, intent(in) :: c
    is_sep = (c == 32 .or. c == 9 .or. c == 44 .or. c == 10)
  end function is_sep

  !> [+-]?(d+[.d*]|.d+)([eEdD][+-]?d+)?: a token one list-directed READ of many values converts as it converts it
  !  alone (no repeat count, slash, null value or other form).
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
      c = iachar(s(i:i))
      if (.not. (c == 101 .or. c == 69 .or. c == 100 .or. c == 68)) return
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

  ! ------------------------------------------------------------------------------------------------------------------
  ! Copies of ORION's header routines (Lib_Tecplot, tec_read_ascii and read_variables), with their messages on the
  ! refusal paths left to ORION's reader.
  ! ------------------------------------------------------------------------------------------------------------------

  ! A zone record starts with the keyword ZONE as a whole word: a name that merely begins with
  ! it (zone_velocity, in a continuation line of the VARIABLES list) is not one.
  logical function is_zone_header(text)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: u
    u = adjustl(upper_case(text))
    is_zone_header = (index(u,'ZONE')==1)
    if (is_zone_header .and. len(u)>4) &
      is_zone_header = .not.(scan(u(5:5),'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_')>0)
  end function is_zone_header

  logical function line_is_numeric_start(text)
    character(len=*), intent(in) :: text
    character(1000) :: t
    character(100) :: tok
    integer :: p, q, ios_
    real(R8) :: value
    t = adjustl(text)
    if (len_trim(t)==0) then
      line_is_numeric_start = .false.
      return
    endif
    p = 1
    do while (p<=len_trim(t))
      if (t(p:p)/=' ' .and. t(p:p)/=char(9) .and. t(p:p)/=',') exit
      p = p+1
    enddo
    if (p>len_trim(t)) then
      line_is_numeric_start = .false.
      return
    endif
    q = p
    do while (q<=len_trim(t))
      if (t(q:q)==' ' .or. t(q:q)==char(9) .or. t(q:q)==',') exit
      q = q+1
    enddo
    tok = ' '
    tok(1:min(len(tok),q-p)) = t(p:q-1)
    read(tok,*,iostat=ios_) value
    line_is_numeric_start = (ios_==0)
  end function line_is_numeric_start

  subroutine parse_zone_header(text,I_,J_,K_,point_,node_,time_,strand_)
    character(len=*), intent(in) :: text
    integer, intent(out) :: I_,J_,K_
    logical, intent(out) :: point_,node_
    real(R8), intent(inout) :: time_
    integer, intent(out) :: strand_
    character(1000) :: work
    character(100) :: token
    integer :: ios_
    real(R8) :: t_
    logical :: found

    I_ = 0; J_ = 0; K_ = 1
    point_ = .false.
    node_ = .true.

    work = upper_case(text)

    call get_integer_keyword(work,'I',I_)
    call get_integer_keyword(work,'J',J_)
    call get_integer_keyword(work,'K',K_,found)
    if (.not.found) K_ = 1

    call get_word_keyword(work,'DATAPACKING',token,found)
    if (.not.found) call get_word_keyword(work,'F',token,found)
    point_ = found .and. trim(token)=='POINT'

    node_ = (index(varlocation_value(work),'CELLCENTERED')==0)

    call get_word_keyword(work,'SOLUTIONTIME',token,found)
    if (found) then
      read(token,*,iostat=ios_) t_
      if (ios_==0) time_ = t_
    endif

    call get_integer_keyword(work,'STRANDID',strand_,found)
    if (.not.found) strand_ = -1
  end subroutine parse_zone_header

  integer function keyword_value_pos(text,key) result(pos)
    character(len=*), intent(in) :: text,key
    integer :: p, q, lt, lk
    logical :: in_quotes

    pos = 0
    lt = len_trim(text)
    lk = len_trim(key)
    in_quotes = .false.

    do p = 1, lt-lk+1
      if (text(p:p)=='"') then
        in_quotes = .not.in_quotes
        cycle
      endif
      if (in_quotes .or. text(p:p+lk-1)/=key(1:lk)) cycle
      if (p>1) then
        if (is_word_char(text(p-1:p-1))) cycle
      endif
      q = skip_blanks(text,p+lk,lt)
      if (q>lt) return
      if (text(q:q)/='=') cycle
      q = skip_blanks(text,q+1,lt)
      if (q<=lt) pos = q
      return
    enddo
  end function keyword_value_pos

  subroutine get_word_keyword(text,key,out,found)
    character(len=*), intent(in) :: text,key
    character(len=*), intent(out) :: out
    logical, intent(out) :: found
    integer :: q, r, n, lt

    out = ' '
    q = keyword_value_pos(text,key)
    found = (q>0)
    if (.not.found) return

    lt = len_trim(text)
    n = 0
    do r = q, lt
      if (index(', )'//char(9),text(r:r))>0 .or. n>=len(out)) exit
      n = n + 1
      out(n:n) = text(r:r)
    enddo
    found = (n>0)
  end subroutine get_word_keyword

  subroutine get_integer_keyword(text,key,out,found)
    character(len=*), intent(in) :: text,key
    integer, intent(out) :: out
    logical, intent(out), optional :: found
    character(100) :: token
    integer :: ios_
    logical :: ok

    out = 0
    call get_word_keyword(text,key,token,ok)
    if (ok) then
      read(token,*,iostat=ios_) out
      ok = (ios_==0)
      if (.not.ok) out = 0
    endif
    if (present(found)) found = ok
  end subroutine get_integer_keyword

  function varlocation_value(text) result(out)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: out
    integer :: q, r

    out = ' '
    q = keyword_value_pos(text,'VARLOCATION')
    if (q<=0) return
    if (text(q:q)/='(') return
    r = index(text(q:),')')
    if (r<=0) then
      r = len_trim(text)
    else
      r = q + r - 1
    endif
    out = text(q:r)
  end function varlocation_value

  pure integer function skip_blanks(text,start,last) result(q)
    character(len=*), intent(in) :: text
    integer, intent(in) :: start, last
    q = start
    do while (q<=last)
      if (text(q:q)/=' ' .and. text(q:q)/=char(9)) exit
      q = q + 1
    enddo
  end function skip_blanks

  pure logical function is_word_char(c)
    character(1), intent(in) :: c
    is_word_char = (c>='A' .and. c<='Z') .or. (c>='a' .and. c<='z') .or. &
                   (c>='0' .and. c<='9') .or. c=='_'
  end function is_word_char

  function upper_case(text) result(out)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: out
    integer :: q, code
    out = text
    do q = 1, len(text)
      code = iachar(out(q:q))
      if (code>=iachar('a') .and. code<=iachar('z')) &
        out(q:q) = achar(code-iachar('a')+iachar('A'))
    enddo
  end function upper_case

  logical function is_coordinate_name(text,c)
    character(len=*), intent(in) :: text, c
    is_coordinate_name = (trim(adjustl(upper_case(text)))==c)
  end function is_coordinate_name

  logical function var_cellcentered(text,iv)
    character(len=*), intent(in) :: text
    integer,          intent(in) :: iv
    integer :: p, q, r, c
    character(len=len(text)) :: vl
    var_cellcentered = .false.
    vl = varlocation_value(text)
    p = 1
    do
      q = index(vl(p:),'[')
      if (q<=0) return
      q = p + q - 1
      r = index(vl(q:),']')
      if (r<=0) return
      r = q + r - 1
      p = r + 1
      if (p>len(vl)) return
      c = verify(vl(p:),' =')
      if (c<=0) return
      c = p + c - 1
      if (c+11<=len(vl)) then
        if (vl(c:c+11)=='CELLCENTERED') then
          if (in_var_list(vl(q+1:r-1),iv)) then
            var_cellcentered = .true.
            return
          endif
        endif
      endif
    enddo
  end function var_cellcentered

  logical function in_var_list(list,iv)
    character(len=*), intent(in) :: list
    integer,          intent(in) :: iv
    integer :: s, e, d, lo, hi, ios1, ios2
    in_var_list = .false.
    s = 1
    do while (s<=len(list))
      e = index(list(s:),',')
      if (e==0) then
        e = len(list) + 1
      else
        e = s + e - 1
      endif
      if (e>s) then
        d = index(list(s:e-1),'-')
        if (d>0) then
          read(list(s:s+d-2),*,iostat=ios1) lo
          read(list(s+d:e-1),*,iostat=ios2) hi
        else
          read(list(s:e-1),*,iostat=ios1) lo
          ios2 = ios1
          hi = lo
        endif
        if (ios1==0 .and. ios2==0) then
          if (iv>=lo .and. iv<=hi) then
            in_var_list = .true.
            return
          endif
        endif
      endif
      s = e + 1
    enddo
  end function in_var_list

  subroutine read_variables(line,variables)
    character(len=*), intent(inout) :: line
    character(len=32), allocatable, intent(out) :: variables(:)

    integer, parameter :: MAX_VARIABLES = 512
    character(len=32) :: names(MAX_VARIABLES)
    character(len=len(line)) :: name
    character(len=len(line)) :: work
    character(len=len(line)) :: upper
    integer :: nvar
    integer :: i, p, q, L, n, c
    character :: quote
    logical :: quoted
    logical :: quote_closed

    nvar = 0
    names = ' '

    work  = adjustl(line)
    upper = work
    do c = 1, len_trim(upper)
      if (iachar(upper(c:c)) >= iachar('a') .and. &
          iachar(upper(c:c)) <= iachar('z')) then
        upper(c:c) = achar(iachar(upper(c:c)) - iachar('a') + iachar('A'))
      endif
    enddo

    if (index(adjustl(upper),'VARIABLES') /= 1) return

    p = index(work,'=')
    if (p <= 0) return

    L = len_trim(work)
    i = p + 1

    do while (i <= L)

      do while (i <= L)
        if (work(i:i) /= ' ' .and. work(i:i) /= char(9) .and. work(i:i) /= ',') exit
        i = i + 1
      enddo
      if (i > L) exit

      if (work(i:i) == '"' .or. work(i:i) == "'") then
        quoted = .true.
        quote = work(i:i)
        i = i + 1
        name = ' '
        n = 0
        quote_closed = .false.

        do while (i <= L)
          if (work(i:i) == quote) then
            i = i + 1
            quote_closed = .true.
            exit
          endif

          if (n < len(name)) then
            n = n + 1
            name(n:n) = work(i:i)
          endif
          i = i + 1
        enddo

        if (.not.quote_closed) return

        if (n > 0) then
          call append_variable_name(names,nvar,MAX_VARIABLES,trim(name))
          if (nvar >= MAX_VARIABLES .and. i <= L) return
        endif

      else
        quoted = .false.
        q = i
        do while (q <= L)
          if (work(q:q) == ' ' .or. work(q:q) == char(9) .or. work(q:q) == ',') exit
          q = q + 1
        enddo

        if (q > i) then
          name = ' '
          name(1:min(len(name),q-i)) = work(i:q-1)
          call append_variable_name(names,nvar,MAX_VARIABLES,trim(name))
          if (nvar >= MAX_VARIABLES .and. q < L) return
        endif
        i = q
      endif

      if (quoted .and. i <= L) cycle
    enddo

    if (nvar <= 0) return

    allocate(character(32)::variables(1:nvar))
    variables = names(1:nvar)

  end subroutine read_variables

  subroutine append_variable_name(names,nvar,max_variables,value)
    character(len=32), intent(inout) :: names(:)
    integer, intent(inout) :: nvar
    integer, intent(in) :: max_variables
    character(len=*), intent(in) :: value

    if (len_trim(value) == 0) return

    if (nvar >= max_variables) then
      write(stderr,'(A,I0)') &
        'TECPLOT ASCII READ ERROR: more than MAX_VARIABLES=',max_variables
      return
    endif

    nvar = nvar + 1
    if (len_trim(value) > len(names(nvar))) then
      write(stderr,'(A)') &
        'TECPLOT ASCII READ WARNING: variable name exceeds 32 characters; truncating.'
      write(stderr,'(A)') 'Variable: '//trim(value)
    endif
    names(nvar) = value(1:min(len(value),len(names(nvar))))
  end subroutine append_variable_name

end module IGLOO_Mod_TecASCII
