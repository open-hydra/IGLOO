!> A parcel's particle records: the rows of each stream held as numbers while the parcel is integrated, formatted by
!  IGLOO_fmt at its end and written as one byte block per stream, under that stream's critical section.
module IGLOO_records
  use, intrinsic :: iso_fortran_env, only: R8 => real64
  use IGLOO_fmt, only: fmtSpec, fmt_compile, fmt_row
  implicit none
  private
  public :: recList, rec_setup, rec_emit, rec_push

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

contains

  !> Compile the layouts (serial, before the first parcel loop).
  subroutine rec_setup()
    integer :: k
    if (ready) return
    do k = 1, NLAY
      call fmt_compile(spec(k), trim(LAYFMT(k)))
    enddo
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


  !> Format the rows of list and write them to unit as one block of stream; the list is emptied.
  subroutine rec_emit(list, unit, stream)
    type(recList), intent(inout) :: list
    integer, intent(in) :: unit, stream
    character(len=:), allocatable :: buf
    integer :: k, n
    if (list%n == 0) return
    allocate(character(len=list%n*(spec(list%lay)%maxlen + 1)) :: buf)
    n = 0
    do k = 1, list%n
      call fmt_row(spec(list%lay), list%v(:, k), buf, n)
      n = n + 1
      buf(n:n) = LF
    enddo
    select case (stream)
    case (STR_TRAJ)
      !$OMP CRITICAL (igloo_rec_traj)
      write(unit) buf(1:n)
      !$OMP END CRITICAL (igloo_rec_traj)
    case (STR_EXIT)
      !$OMP CRITICAL (igloo_rec_exit)
      write(unit) buf(1:n)
      !$OMP END CRITICAL (igloo_rec_exit)
    case (STR_PLANE)
      !$OMP CRITICAL (igloo_rec_plane)
      write(unit) buf(1:n)
      !$OMP END CRITICAL (igloo_rec_plane)
    case (STR_SNAP)
      !$OMP CRITICAL (igloo_rec_snap)
      write(unit) buf(1:n)
      !$OMP END CRITICAL (igloo_rec_snap)
    case (STR_SCAT)
      !$OMP CRITICAL (igloo_rec_scat)
      write(unit) buf(1:n)
      !$OMP END CRITICAL (igloo_rec_scat)
    end select
    list%n = 0
  end subroutine rec_emit

end module IGLOO_records
