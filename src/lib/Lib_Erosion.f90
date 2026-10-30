!> Wall erosion by particle impact ([IGLOO-Erosion]): the model of each material, the boundary codes where a
!  parcel's exit hits the wall, and the wall-impact sums behind OUTPUT/wall.tec.
module IGLOO_Lib_Erosion
  use, intrinsic :: iso_fortran_env, only : R8 => real64
  implicit none
  private
  public :: EROSION_NONE
  public :: erosionOn, eroModel, eroModelWord, rhoWall
  public :: isWallCode, setupWall, zeroWall, erosionRatio, checkErosionInput

  integer, parameter :: EROSION_NONE = 0   !> the material erodes nothing; its incident flow still counts

  !> [IGLOO-Erosion] present: the wall sums are kept and OUTPUT/wall.tec is written.
  logical :: erosionOn = .false.
  !> The model token of each material of the phase file, as read, and its code.
  character(len=64), allocatable :: eroModelWord(:)
  integer,           allocatable :: eroModel(:)
  !> Wall density [kg/m^3], the same for every wall: rdot = sum over the materials of g*er/rhoWall.
  real(R8) :: rhoWall = 0._R8

contains

  !> A face cell where a parcel leaving the domain hits a wall: a dispersed-phase wall (301), or an interface to a
  !  block of another phase, header only in a dispersed-phase file (103 with no partner cell, 104).
  pure logical function isWallCode(cell)
    use IGLOO_data_block, only: obj_bc_cell
    type(obj_bc_cell), intent(in) :: cell
    select case (cell%bcdef)
    case (301, 104); isWallCode = .true.
    case (103);      isWallCode = cell%connection(1) == 0
    case default;    isWallCode = .false.
    end select
  end function isWallCode


  !> One wall-sum entry per block; with [IGLOO-Erosion], the faces holding a wall cell flagged on and their sums
  !  allocated and zeroed. nface: the faces flagged.
  subroutine setupWall(geo, wall, nmat, nface)
    use IGLOO_data_block, only: obj_block, obj_wallblock
    class(obj_block),                 intent(in)  :: geo(:)
    type(obj_wallblock), allocatable, intent(out) :: wall(:)
    integer,                          intent(in)  :: nmat
    integer,                          intent(out) :: nface
    integer :: b, f, m, n

    allocate(wall(size(geo)))
    nface = 0
    if (.not. erosionOn) return
    do b = 1, size(geo)
      do f = 1, 6
        if (.not. allocated(geo(b)%face(f)%cell)) cycle
        do n = 1, geo(b)%face(f)%Nn
          do m = 1, geo(b)%face(f)%Nm
            if (isWallCode(geo(b)%face(f)%cell(m,n))) wall(b)%on(f) = .true.
          enddo
        enddo
        if (.not. wall(b)%on(f)) cycle
        allocate(wall(b)%face(f)%acc(2, nmat, geo(b)%face(f)%Nm, geo(b)%face(f)%Nn))
        wall(b)%face(f)%acc = 0._R8
        nface = nface + 1
      enddo
    enddo

  end subroutine setupWall


  !> Zero the wall sums (start of every solve).
  subroutine zeroWall(wall)
    use IGLOO_data_block, only: obj_wallblock
    type(obj_wallblock), intent(inout) :: wall(:)
    integer :: b, f

    do b = 1, size(wall)
      do f = 1, 6
        if (wall(b)%on(f)) wall(b)%face(f)%acc = 0._R8
      enddo
    enddo

  end subroutine zeroWall


  !> Wall mass eroded per unit mass of impacting particles [kg/kg] by one parcel of material mID that hits the wall
  !  at speed v [m/s] and angle alpha [rad] from the surface, with diameter d [m] and density rhop [kg/m^3].
  pure function erosionRatio(mID, v, alpha, d, rhop) result(er)
    integer,  intent(in) :: mID
    real(R8), intent(in) :: v, alpha, d, rhop
    real(R8) :: er

    select case (eroModel(mID))
    case default
      er = 0._R8
    end select

  end function erosionRatio


  !> The [IGLOO-Erosion] input against the phase file: one model token per material, each a known model.
  subroutine checkErosionInput(nmat)
    integer, intent(in) :: nmat
    integer :: m

    if (.not. erosionOn) return
    if (size(eroModelWord) /= nmat) then
      write(*,'(A,I0,A,I0,A)') ' [ERROR] [IGLOO-Erosion] model: ', size(eroModelWord), ' token(s) for ', nmat, &
        ' material(s) of the phase file'
      error stop 'IGLOO: [IGLOO-Erosion] model needs one token per material'
    endif
    if (allocated(eroModel)) deallocate(eroModel)
    allocate(eroModel(nmat))
    do m = 1, nmat
      select case (lowercase(trim(adjustl(eroModelWord(m)))))
      case ('none')
        eroModel(m) = EROSION_NONE
      case default
        write(*,'(A,I0,A)') ' [ERROR] [IGLOO-Erosion] model of material ', m, ': unknown token "'// &
          trim(adjustl(eroModelWord(m)))//'"'
        error stop 'IGLOO: unknown [IGLOO-Erosion] model'
      end select
    enddo

  end subroutine checkErosionInput


  !> ASCII lower case.
  pure function lowercase(s) result(t)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: t
    integer :: i

    t = s
    do i = 1, len(s)
      if (s(i:i) >= 'A' .and. s(i:i) <= 'Z') t(i:i) = achar(iachar(s(i:i)) + 32)
    enddo

  end function lowercase

end module IGLOO_Lib_Erosion
