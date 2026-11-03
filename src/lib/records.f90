!> A parcel's particle records: the rows of each stream held as numbers while the parcel is integrated, formatted by
!  IGLOO_fmt at its end and handed as one byte block per stream to that stream's writer, which writes the blocks of a
!  parcel loop in slot (parcel) order: a block that arrives before its turn waits in memory, within a budget over every
!  stream, or else in the stream's spill file.
module IGLOO_records
  use, intrinsic :: iso_fortran_env, only: R8 => real64, I8 => int64
  use IGLOO_fmt, only: fmtSpec, fmt_compile, fmt_row
  use IGLOO_Mod_Timing, only: tm_window
  implicit none
  private
  public :: recList, rec_setup, rec_emit, rec_push, ps_open, ps_close, ps_begin, ps_end

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

  !> One stream's rows of a parcel, v(:, k) row k in the layout lay (integer columns hold exact integers).
  type :: recList
    real(R8), allocatable :: v(:,:)
    integer :: n = 0
    integer :: lay = 0
  end type recList

  !> A window slot: no block yet, nothing (left) to write, a block held in memory, a block in the spill file
  integer, parameter :: ST_WAIT = 0, ST_NONE = 1, ST_MEM = 2, ST_SPILL = 3

  !> A block held in memory
  type :: heldBlock
    character(len=:), allocatable :: b
  end type heldBlock

  !> A particle file's writer: its unit and the window of the current parcel loop (slots first..last, next = the first
  !  slot not written yet), each slot's state, held block, or spill offset and length
  type :: partStream
    logical :: on = .false.
    integer :: unit = 0
    integer :: first = 1, last = 0, next = 1
    integer, allocatable :: st(:), nb(:)
    integer(I8), allocatable :: off(:)
    type(heldBlock), allocatable :: mem(:)
    logical :: spillOpen = .false.
    integer :: spill = 0
    integer(I8) :: spillEnd = 1
  end type partStream

  type(partStream) :: ps(5)
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
      if (st /= 0 .or. budget < 0) error stop 'IGLOO: IGLOO_IOBUF_BYTES must be an integer >= 0 [bytes]'
    elseif (st == -1) then
      error stop 'IGLOO: IGLOO_IOBUF_BYTES must be an integer >= 0 [bytes]'
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


  !> Format the rows of list into one block and hand it to stream's writer as the block of slot, an empty block
  !  included; the list is emptied.
  subroutine rec_emit(list, stream, slot)
    type(recList), intent(inout) :: list
    integer, intent(in) :: stream, slot
    character(len=:), allocatable :: buf
    integer :: k, n
    if (list%n == 0) then
      call ps_hand(stream, slot, '')
      return
    endif
    allocate(character(len=list%n*(spec(list%lay)%maxlen + 1)) :: buf)
    n = 0
    do k = 1, list%n
      call fmt_row(spec(list%lay), list%v(:, k), buf, n)
      n = n + 1
      buf(n:n) = LF
    enddo
    call ps_hand(stream, slot, buf(1:n))
    list%n = 0
  end subroutine rec_emit


  !> Hand the block of slot to stream's writer, under that stream's critical section.
  subroutine ps_hand(stream, slot, blk)
    integer, intent(in) :: stream, slot
    character(len=*), intent(in) :: blk
    select case (stream)
    case (STR_TRAJ)
      !$OMP CRITICAL (igloo_rec_traj)
      call ps_put(ps(STR_TRAJ), slot, blk)
      !$OMP END CRITICAL (igloo_rec_traj)
    case (STR_EXIT)
      !$OMP CRITICAL (igloo_rec_exit)
      call ps_put(ps(STR_EXIT), slot, blk)
      !$OMP END CRITICAL (igloo_rec_exit)
    case (STR_PLANE)
      !$OMP CRITICAL (igloo_rec_plane)
      call ps_put(ps(STR_PLANE), slot, blk)
      !$OMP END CRITICAL (igloo_rec_plane)
    case (STR_SNAP)
      !$OMP CRITICAL (igloo_rec_snap)
      call ps_put(ps(STR_SNAP), slot, blk)
      !$OMP END CRITICAL (igloo_rec_snap)
    case (STR_SCAT)
      !$OMP CRITICAL (igloo_rec_scat)
      call ps_put(ps(STR_SCAT), slot, blk)
      !$OMP END CRITICAL (igloo_rec_scat)
    end select
  end subroutine ps_hand


  !> The block of slot (under the stream's critical section): written now when it is the next slot, else held, in memory
  !  while the bytes held over every stream stay within the budget, else in the spill file; then the ready slots that
  !  follow are written.
  subroutine ps_put(s, slot, blk)
    type(partStream), intent(inout) :: s
    integer, intent(in) :: slot
    character(len=*), intent(in) :: blk
    integer :: i
    integer(I8) :: n, h
    if (.not. s%on) error stop 'IGLOO: a particle block for a file that is not open'
    if (slot < s%first .or. slot > s%last) error stop 'IGLOO: a particle block outside its parcel loop''s window'
    i = slot - s%first + 1
    if (s%st(i) /= ST_WAIT) error stop 'IGLOO: a parcel handed two blocks to one particle file'
    n = len(blk, kind=I8)
    !$OMP ATOMIC UPDATE
    largest = max(largest, n)
    if (slot == s%next) then
      if (n > 0) write(s%unit) blk
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
        write(s%unit) s%mem(i)%b
        deallocate(s%mem(i)%b)
        !$OMP ATOMIC UPDATE
        held = held - n
      case (ST_SPILL)
        allocate(character(len=s%nb(i)) :: tmp)
        read(s%spill, pos=s%off(i)) tmp
        write(s%unit) tmp
        deallocate(tmp)
      end select
      s%st(i) = ST_NONE
      s%next = s%next + 1
    enddo
  end subroutine ps_drain


  !> Stream k's file is open on unit (master, outside the parcel loops).
  subroutine ps_open(k, unit)
    integer, intent(in) :: k, unit
    ps(k)%on = .true.
    ps(k)%unit = unit
    ps(k)%first = 1; ps(k)%last = 0; ps(k)%next = 1
  end subroutine ps_open


  !> Close stream k's file and its spill file (a scratch file, deleted).
  subroutine ps_close(k)
    integer, intent(in) :: k
    if (.not. ps(k)%on) return
    if (ps(k)%spillOpen) close(ps(k)%spill)
    ps(k)%spillOpen = .false.
    close(ps(k)%unit)
    ps(k)%on = .false.
  end subroutine ps_close


  !> Open the window of a parcel loop over slots first..last on every open stream; skip marks the slots whose parcel is
  !  not integrated here (master, before the loop).
  subroutine ps_begin(first, last, skip)
    integer, intent(in) :: first, last
    logical, intent(in) :: skip(first:)
    integer :: k, n
    n = max(last - first + 1, 0)
    do k = 1, size(ps)
      if (.not. ps(k)%on) cycle
      call ps_grow(ps(k), n)
      ps(k)%first = first; ps(k)%last = last; ps(k)%next = first; ps(k)%spillEnd = 1
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
      deallocate(s%st, s%nb, s%off, s%mem)
    endif
    allocate(s%st(max(n, 1)), s%nb(max(n, 1)), s%off(max(n, 1)), s%mem(max(n, 1)))
  end subroutine ps_grow


  !> Close the window of the loop just joined (master, after the failure check): write what is left, check that every
  !  slot was written and nothing is held, and hand the window's figures to the timing report.
  subroutine ps_end()
    integer :: k
    do k = 1, size(ps)
      if (.not. ps(k)%on) cycle
      call ps_drain(ps(k))
      if (ps(k)%next <= ps(k)%last) error stop 'IGLOO: a parcel of the loop left no block for a particle file'
    enddo
    if (held /= 0) error stop 'IGLOO: the particle-file window still holds bytes after its loop'
    call tm_window(budget, heldMax, spilled, largest)
  end subroutine ps_end

end module IGLOO_records
