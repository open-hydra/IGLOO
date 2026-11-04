!> A parcel's particle records: the rows of each stream held as numbers while the parcel is integrated, formatted by
!  IGLOO_fmt at its end and handed as one byte block per stream to that stream's writer, which writes the blocks of a
!  parcel loop in slot (parcel) order: a block that arrives before its turn waits in memory, within a budget over every
!  stream, or else in the stream's spill file. Under MPI each file (a rank's shard) is indexed, one row per header line
!  and parcel block, and the index is gathered on the root rank for the merge.
module IGLOO_records
  use, intrinsic :: iso_fortran_env, only: R8 => real64, I8 => int64
  use IGLOO_fmt, only: fmtSpec, fmt_compile, fmt_row, fmt_int
  use IGLOO_Mod_Timing, only: tm_window
  use IGLOO_Mod_MPI, only: mpi_size_, mpi_gather_i8, mpi_abort_all, region_fail
  implicit none
  private
  public :: recList, rec_setup, rec_emit, rec_emit_traj, rec_push, rec_break, ps_open, ps_text, ps_zone, ps_close, &
            ps_begin, ps_end

  !> The streams (one file each per material)
  integer, parameter, public :: STR_TRAJ = 1, STR_EXIT = 2, STR_PLANE = 3, STR_SNAP = 4, STR_SCAT = 5
  !> The record layouts: trajectory, snapshot and scatter rows (t: time, f: frozen fraction, s: surface temperature),
  !  exit records, plane crossings
  integer, parameter, public :: L_ROW = 1, L_ROW_T = 2, L_ROW_F = 3, L_ROW_TF = 4, L_ROW_S = 5, L_ROW_TS = 6, &
                                L_EXIT = 7, L_EXIT_T = 8, L_EXIT_F = 9, L_EXIT_TF = 10, L_PLANE = 11
  integer, parameter :: NLAY = 11
  character(len=56), parameter :: LAYFMT(NLAY) = [character(len=56) :: &
    '(7F12.6,2E13.6E2,I8)', '(7F12.6,2E13.6E2,I8,ES16.8E2)', '(7F12.6,2E13.6E2,I8,F10.6)', &
    '(7F12.6,2E13.6E2,I8,ES16.8E2,F10.6)', '(7F12.6,2E13.6E2,I8,F12.6)', '(7F12.6,2E13.6E2,I8,ES16.8E2,F12.6)', &
    '(6F12.6,2E13.6E2,I8,I5,E13.6E2)', '(6F12.6,2E13.6E2,I8,ES16.8E2,I5,E13.6E2)', &
    '(6F12.6,2E13.6E2,I8,I5,E13.6E2,F10.6)', '(6F12.6,2E13.6E2,I8,ES16.8E2,I5,E13.6E2,F10.6)', &
    '(I3,7F12.6,2E13.6E2,2ES16.8E2,F10.6,2I10,ES16.8E2,I3)']
  character, parameter :: LF = achar(10)

  !> The compiled layouts: written once by rec_setup, outside the parcel loops
  type(fmtSpec) :: spec(NLAY)
  logical :: ready = .false.

  !> One stream's rows of a parcel, v(:, k) row k in the layout lay (integer columns hold exact integers); a trajectory's
  !  images after the first: brk(:, j) = (its first row, its number)
  type :: recList
    real(R8), allocatable :: v(:,:)
    integer :: n = 0
    integer :: lay = 0
    integer, allocatable :: brk(:,:)
    integer :: nbrk = 0
  end type recList

  !> A window slot: no block yet, nothing (left) to write, a block held in memory, a block in the spill file
  integer, parameter :: ST_WAIT = 0, ST_NONE = 1, ST_MEM = 2, ST_SPILL = 3

  !> A block held in memory
  type :: heldBlock
    character(len=:), allocatable :: b
  end type heldBlock

  !> A particle file's writer: its unit and the window of the current parcel loop (slots first..last of the group's
  !  loop, next = the first slot not written yet), each slot's state, parcel ID, held block, or spill offset and length;
  !  under MPI the file's index, ix(:, j) = (group, key, offset, length) of item j in file order, and the next offset
  type :: partStream
    logical :: on = .false.
    integer :: unit = 0
    integer :: first = 1, last = 0, next = 1, group = 0
    integer, allocatable :: st(:), nb(:), id(:)
    integer(I8), allocatable :: off(:)
    type(heldBlock), allocatable :: mem(:)
    logical :: spillOpen = .false.
    integer :: spill = 0
    integer(I8) :: spillEnd = 1
    integer(I8), allocatable :: ix(:,:)
    integer :: nix = 0
    integer(I8) :: pos = 1
    !> trajectories: the current group's zone title, `ZONE T="Mat <mat> Group<g, I3>`
    character(len=:), allocatable :: zone
  end type partStream

  !> A particle file's index gathered on the root rank: rows(4*j-3:4*j) = (group, key, offset, length) of item j, rank
  !  by rank in file order; cnt(r) = the values from rank r. Key 0, or group 0, marks a line every shard holds.
  type, public :: partIndex
    integer(I8), allocatable :: rows(:)
    integer, allocatable :: cnt(:)
  end type partIndex

  type(partStream) :: ps(5)
  !> The files are indexed (decomposed runs)
  logical :: indexOn = .false.
  !> The bytes the windows may hold in memory over every stream (IGLOO_IOBUF_BYTES, a test hook; default 64 MiB)
  integer(I8) :: budget = 67108864_I8
  !> Bytes held now; the current window's high-water of held bytes, spilled bytes and largest block
  integer(I8) :: held = 0, heldMax = 0, spilled = 0, largest = 0

contains

  !> Compile the layouts and read the window budget (serial, before the first parcel loop).
  subroutine rec_setup()
    integer :: k, st
    character(len=40) :: val
    if (ready) return
    do k = 1, NLAY
      call fmt_compile(spec(k), trim(LAYFMT(k)))
    enddo
    call get_environment_variable('IGLOO_IOBUF_BYTES', val, status=st)
    if (st == 0) then
      read(val, *, iostat=st) budget
      if (st /= 0 .or. budget < 0) call mpi_abort_all('IGLOO: IGLOO_IOBUF_BYTES must be an integer >= 0 [bytes]')
    elseif (st == -1) then
      call mpi_abort_all('IGLOO: IGLOO_IOBUF_BYTES must be an integer >= 0 [bytes]')
    endif
    ready = .true.
  end subroutine rec_setup


  !> Append one row of layout lay (a type(recList) dummy: a class one costs each entry of the caller a set-up).
  subroutine rec_push(list, row, lay)
    type(recList), intent(inout) :: list
    real(R8), intent(in) :: row(:)
    integer, intent(in) :: lay
    real(R8), allocatable :: tmp(:,:)
    if (.not. allocated(list%v)) allocate(list%v(size(row), 16))
    if (list%n == size(list%v, 2)) then
      allocate(tmp(size(list%v, 1), 2*size(list%v, 2)))
      tmp(:, 1:list%n) = list%v(:, 1:list%n)
      call move_alloc(tmp, list%v)
    endif
    list%n = list%n + 1
    list%v(:, list%n) = row
    list%lay = lay
  end subroutine rec_push


  !> The rows pushed from now on start image k of the trajectory (a 201 wrap); an image left without rows is dropped.
  subroutine rec_break(list, k)
    type(recList), intent(inout) :: list
    integer, intent(in) :: k
    integer, allocatable :: tmp(:,:)
    if (list%nbrk > 0) then
      if (list%brk(1, list%nbrk) == list%n + 1) then
        list%brk(2, list%nbrk) = k
        return
      endif
    endif
    if (.not. allocated(list%brk)) allocate(list%brk(2, 4))
    if (list%nbrk == size(list%brk, 2)) then
      allocate(tmp(2, 2*size(list%brk, 2)))
      tmp(:, 1:list%nbrk) = list%brk(:, 1:list%nbrk)
      call move_alloc(tmp, list%brk)
    endif
    list%nbrk = list%nbrk + 1
    list%brk(:, list%nbrk) = [list%n + 1, k]
  end subroutine rec_break


  !> The trajectory block of slot (parcel id, shed or burnt out by parent, 0 if injected): one ordered zone per image
  !  with rows, `<zone title> ID <id> parent <parent>[ image <k>]", I=<rows>, ZONETYPE=ORDERED, DATAPACKING=POINT`, then
  !  its rows; the list is emptied.
  subroutine rec_emit_traj(list, slot, id, parent)
    type(recList), intent(inout) :: list
    integer, intent(in) :: slot, id, parent
    character(len=:), allocatable :: buf
    integer :: k, n, b, r0, r1, img
    if (list%n == 0) then
      call ps_hand(STR_TRAJ, slot, id, '')
      list%nbrk = 0
      return
    endif
    allocate(character(len=list%n*(spec(list%lay)%maxlen + 1) + &
                           (list%nbrk + 1)*(len(ps(STR_TRAJ)%zone) + 128)) :: buf)
    n = 0
    r0 = 1
    img = 0
    do b = 0, list%nbrk
      r1 = list%n
      if (b < list%nbrk) r1 = list%brk(1, b+1) - 1
      if (r1 >= r0) then
        call put(ps(STR_TRAJ)%zone//' ID ')
        call fmt_int(int(id, I8), 0, buf, n)
        call put(' parent ')
        call fmt_int(int(parent, I8), 0, buf, n)
        if (img > 0) then
          call put(' image ')
          call fmt_int(int(img, I8), 0, buf, n)
        endif
        call put('", I=')
        call fmt_int(int(r1 - r0 + 1, I8), 0, buf, n)
        call put(', ZONETYPE=ORDERED, DATAPACKING=POINT'//LF)
        do k = r0, r1
          call fmt_row(spec(list%lay), list%v(:, k), buf, n)
          n = n + 1
          buf(n:n) = LF
        enddo
      endif
      if (b < list%nbrk) then
        r0 = list%brk(1, b+1)
        img = list%brk(2, b+1)
      endif
    enddo
    call ps_hand(STR_TRAJ, slot, id, buf(1:n))
    list%n = 0
    list%nbrk = 0
  contains
    subroutine put(t)
      character(len=*), intent(in) :: t
      buf(n+1:n+len(t)) = t
      n = n + len(t)
    end subroutine put
  end subroutine rec_emit_traj


  !> Format the rows of list into one block and hand it to stream's writer as the block of slot (parcel id), an empty
  !  block included; the list is emptied.
  subroutine rec_emit(list, stream, slot, id)
    type(recList), intent(inout) :: list
    integer, intent(in) :: stream, slot, id
    character(len=:), allocatable :: buf
    integer :: k, n
    if (list%n == 0) then
      call ps_hand(stream, slot, id, '')
      return
    endif
    allocate(character(len=list%n*(spec(list%lay)%maxlen + 1)) :: buf)
    n = 0
    do k = 1, list%n
      call fmt_row(spec(list%lay), list%v(:, k), buf, n)
      n = n + 1
      buf(n:n) = LF
    enddo
    call ps_hand(stream, slot, id, buf(1:n))
    list%n = 0
  end subroutine rec_emit


  !> Hand the block of slot (parcel id) to stream's writer, under that stream's critical section.
  subroutine ps_hand(stream, slot, id, blk)
    integer, intent(in) :: stream, slot, id
    character(len=*), intent(in) :: blk
    select case (stream)
    case (STR_TRAJ)
      !$OMP CRITICAL (igloo_rec_traj)
      call ps_put(ps(STR_TRAJ), slot, id, blk)
      !$OMP END CRITICAL (igloo_rec_traj)
    case (STR_EXIT)
      !$OMP CRITICAL (igloo_rec_exit)
      call ps_put(ps(STR_EXIT), slot, id, blk)
      !$OMP END CRITICAL (igloo_rec_exit)
    case (STR_PLANE)
      !$OMP CRITICAL (igloo_rec_plane)
      call ps_put(ps(STR_PLANE), slot, id, blk)
      !$OMP END CRITICAL (igloo_rec_plane)
    case (STR_SNAP)
      !$OMP CRITICAL (igloo_rec_snap)
      call ps_put(ps(STR_SNAP), slot, id, blk)
      !$OMP END CRITICAL (igloo_rec_snap)
    case (STR_SCAT)
      !$OMP CRITICAL (igloo_rec_scat)
      call ps_put(ps(STR_SCAT), slot, id, blk)
      !$OMP END CRITICAL (igloo_rec_scat)
    end select
  end subroutine ps_hand


  !> The block of slot (parcel id; under the stream's critical section): written now when it is the next slot, else
  !  held, in memory while the bytes held over every stream stay within the budget, else in the spill file; then the
  !  ready slots that follow are written.
  subroutine ps_put(s, slot, id, blk)
    type(partStream), intent(inout) :: s
    integer, intent(in) :: slot, id
    character(len=*), intent(in) :: blk
    integer :: i
    integer(I8) :: n, h
    !> inside the parcel loop: a defect fails the loop, which aborts after the join
    if (.not. s%on .or. slot < s%first .or. slot > s%last) then
      call region_fail('IGLOO: a particle block outside an open file''s window')
      return
    endif
    i = slot - s%first + 1
    if (s%st(i) /= ST_WAIT) then
      call region_fail('IGLOO: a parcel handed two blocks to one particle file')
      return
    endif
    n = len(blk, kind=I8)
    !$OMP ATOMIC UPDATE
    largest = max(largest, n)
    s%id(i) = id
    if (slot == s%next) then
      if (n > 0) call ps_out(s, blk, s%group, id)
      s%st(i) = ST_NONE
      s%next = s%next + 1
      call ps_drain(s)
    elseif (n == 0) then
      s%st(i) = ST_NONE
    else
      !$OMP ATOMIC CAPTURE
      held = held + n
      h = held
      !$OMP END ATOMIC
      if (h <= budget) then
        s%mem(i)%b = blk
        s%st(i) = ST_MEM
        !$OMP ATOMIC UPDATE
        heldMax = max(heldMax, h)
      else
        !$OMP ATOMIC UPDATE
        held = held - n
        if (.not. s%spillOpen) then
          open(newunit=s%spill, status='scratch', access='stream', form='unformatted')
          s%spillOpen = .true.
        endif
        write(s%spill, pos=s%spillEnd) blk
        s%off(i) = s%spillEnd
        s%nb(i) = int(n)
        s%spillEnd = s%spillEnd + n
        s%st(i) = ST_SPILL
        !$OMP ATOMIC UPDATE
        spilled = spilled + n
      endif
    endif
  end subroutine ps_put


  !> Write the ready slots from next on, freeing each, up to the first slot without its block.
  subroutine ps_drain(s)
    type(partStream), intent(inout) :: s
    character(len=:), allocatable :: tmp
    integer :: i
    integer(I8) :: n
    do while (s%next <= s%last)
      i = s%next - s%first + 1
      select case (s%st(i))
      case (ST_WAIT)
        exit
      case (ST_MEM)
        n = len(s%mem(i)%b, kind=I8)
        call ps_out(s, s%mem(i)%b, s%group, s%id(i))
        deallocate(s%mem(i)%b)
        !$OMP ATOMIC UPDATE
        held = held - n
      case (ST_SPILL)
        allocate(character(len=s%nb(i)) :: tmp)
        read(s%spill, pos=s%off(i)) tmp
        call ps_out(s, tmp, s%group, s%id(i))
        deallocate(tmp)
      end select
      s%st(i) = ST_NONE
      s%next = s%next + 1
    enddo
  end subroutine ps_drain


  !> Stream k's file is open on unit, empty (master, outside the parcel loops).
  subroutine ps_open(k, unit)
    integer, intent(in) :: k, unit
    ps(k)%on = .true.
    ps(k)%unit = unit
    ps(k)%first = 1; ps(k)%last = 0; ps(k)%next = 1; ps(k)%group = 0
    indexOn = mpi_size_ > 1
    ps(k)%nix = 0
    ps(k)%pos = 1
    if (indexOn .and. .not. allocated(ps(k)%ix)) allocate(ps(k)%ix(4, 1024))
  end subroutine ps_open


  !> The zone title of stream k's next group (master, outside the parcel loops): `ZONE T="Mat <mat> Group<g, I3>`.
  subroutine ps_zone(k, title)
    integer, intent(in) :: k
    character(len=*), intent(in) :: title
    ps(k)%zone = title
  end subroutine ps_zone


  !> A header or zone line of stream k's file (master, outside the parcel loops): group 0 or key 0, a line every rank's
  !  shard holds.
  subroutine ps_text(k, text, group, key)
    integer, intent(in) :: k, group, key
    character(len=*), intent(in) :: text
    call ps_out(ps(k), text, group, key)
  end subroutine ps_text


  !> Write bytes to the stream's file; when decomposed, index them as (group, key, offset, length).
  subroutine ps_out(s, bytes, group, key)
    type(partStream), intent(inout) :: s
    character(len=*), intent(in) :: bytes
    integer, intent(in) :: group, key
    integer(I8), allocatable :: tmp(:,:)
    write(s%unit) bytes
    if (.not. indexOn) return
    if (s%nix == size(s%ix, 2)) then
      allocate(tmp(4, 2*size(s%ix, 2)))
      tmp(:, 1:s%nix) = s%ix(:, 1:s%nix)
      call move_alloc(tmp, s%ix)
    endif
    s%nix = s%nix + 1
    s%ix(:, s%nix) = [int(group, I8), int(key, I8), s%pos, len(bytes, kind=I8)]
    s%pos = s%pos + len(bytes, kind=I8)
  end subroutine ps_out


  !> Close stream k's file and its spill file (a scratch file, deleted). Decomposed: the index must tile the file, and it
  !  is gathered on the root rank into idx (collective: every rank closes the same streams in the same order).
  subroutine ps_close(k, idx)
    integer, intent(in) :: k
    type(partIndex), intent(inout) :: idx
    integer(I8) :: p
    if (.not. ps(k)%on) return
    if (indexOn) then
      inquire(unit=ps(k)%unit, pos=p)
      if (p /= ps(k)%pos) call mpi_abort_all('IGLOO: a particle file holds bytes its block index does not')
      call mpi_gather_i8(ps(k)%ix, 4*ps(k)%nix, idx%rows, idx%cnt)
    endif
    if (ps(k)%spillOpen) close(ps(k)%spill)
    ps(k)%spillOpen = .false.
    close(ps(k)%unit)
    ps(k)%on = .false.
  end subroutine ps_close


  !> Open the window of a parcel loop of group over slots first..last on every open stream; skip marks the slots whose
  !  parcel is not integrated here (master, before the loop).
  subroutine ps_begin(first, last, skip, group)
    integer, intent(in) :: first, last, group
    logical, intent(in) :: skip(first:)
    integer :: k, n
    n = max(last - first + 1, 0)
    do k = 1, size(ps)
      if (.not. ps(k)%on) cycle
      call ps_grow(ps(k), n)
      ps(k)%first = first; ps(k)%last = last; ps(k)%next = first; ps(k)%spillEnd = 1; ps(k)%group = group
      if (n > 0) ps(k)%st(1:n) = merge(ST_NONE, ST_WAIT, skip(first:last))
      call ps_drain(ps(k))
    enddo
    held = 0; heldMax = 0; spilled = 0; largest = 0
  end subroutine ps_begin


  !> Room for n slots.
  subroutine ps_grow(s, n)
    type(partStream), intent(inout) :: s
    integer, intent(in) :: n
    if (allocated(s%st)) then
      if (size(s%st) >= n) return
      deallocate(s%st, s%nb, s%id, s%off, s%mem)
    endif
    allocate(s%st(max(n, 1)), s%nb(max(n, 1)), s%id(max(n, 1)), s%off(max(n, 1)), s%mem(max(n, 1)))
  end subroutine ps_grow


  !> Close the window of the loop just joined (master, after the failure check): write what is left, check that every
  !  slot was written and nothing is held, and hand the window's figures to the timing report.
  subroutine ps_end()
    integer :: k
    do k = 1, size(ps)
      if (.not. ps(k)%on) cycle
      call ps_drain(ps(k))
      if (ps(k)%next <= ps(k)%last) call mpi_abort_all('IGLOO: a parcel of the loop left no block for a particle file')
    enddo
    if (held /= 0) call mpi_abort_all('IGLOO: the particle-file window still holds bytes after its loop')
    call tm_window(budget, heldMax, spilled, largest)
  end subroutine ps_end

end module IGLOO_records
