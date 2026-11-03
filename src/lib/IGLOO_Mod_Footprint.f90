!> Allocated real values per rank of the gas, geometry and accumulator arrays, for the
!  [IGLOO-timing] sizes line.
module IGLOO_Mod_Footprint
  use, intrinsic :: iso_fortran_env, only: I8 => int64
  use IGLOO_data_block, only: obj_block, obj_flowblock, obj_sourceblock, obj_eulerblock
  use IGLOO_Mod_Timing, only: timingOn, tm_sizes
  implicit none
  private
  public :: record_footprint

contains

  !> Counts every allocated real array of the blocks and hands the totals to the timing report.
  subroutine record_footprint(geoblock, gasblock, source, euler)
    type(obj_block),       intent(in) :: geoblock(:)
    type(obj_flowblock),   intent(in) :: gasblock(:)
    type(obj_sourceblock), intent(in) :: source(:)
    type(obj_eulerblock),  intent(in) :: euler(:,:)
    integer(I8) :: nGas, nGeo, nAcc
    integer     :: b, f

    if (.not. timingOn) return
    nGas = 0; nGeo = 0; nAcc = 0
    do b = 1, size(geoblock)
      nGeo = nGeo + geoValues(geoblock(b))
    enddo
    do b = 1, size(gasblock)
      nGeo = nGeo + geoValues(gasblock(b))
      associate(g => gasblock(b))
      if (allocated(g%density))     nGas = nGas + size(g%density, kind=I8)
      if (allocated(g%velocity))    nGas = nGas + size(g%velocity, kind=I8)
      if (allocated(g%temperature)) nGas = nGas + size(g%temperature, kind=I8)
      if (allocated(g%mit))         nGas = nGas + size(g%mit, kind=I8)
      if (allocated(g%mil))         nGas = nGas + size(g%mil, kind=I8)
      if (allocated(g%kl))          nGas = nGas + size(g%kl, kind=I8)
      if (allocated(g%gam))         nGas = nGas + size(g%gam, kind=I8)
      if (allocated(g%R))           nGas = nGas + size(g%R, kind=I8)
      end associate
    enddo
    do b = 1, size(source)
      nGeo = nGeo + geoValues(source(b))
      associate(s => source(b))
      if (allocated(s%sourceMass)) nAcc = nAcc + size(s%sourceMass, kind=I8)
      if (allocated(s%sourceMom))  nAcc = nAcc + size(s%sourceMom, kind=I8)
      if (allocated(s%sourceEn))   nAcc = nAcc + size(s%sourceEn, kind=I8)
      end associate
    enddo
    do f = 1, size(euler, 2)
      do b = 1, size(euler, 1)
        nGeo = nGeo + geoValues(euler(b, f))
        associate(e => euler(b, f))
        if (allocated(e%density))     nAcc = nAcc + size(e%density, kind=I8)
        if (allocated(e%np))          nAcc = nAcc + size(e%np, kind=I8)
        if (allocated(e%velocity))    nAcc = nAcc + size(e%velocity, kind=I8)
        if (allocated(e%temperature)) nAcc = nAcc + size(e%temperature, kind=I8)
        end associate
      enddo
    enddo
    call tm_sizes(nGas, nGeo, nAcc)
  end subroutine record_footprint


  !> Geometry values held by one block of any kind.
  pure function geoValues(blk) result(n)
    class(obj_block), intent(in) :: blk
    integer(I8) :: n
    n = 0
    if (allocated(blk%node))    n = n + size(blk%node, kind=I8)
    if (allocated(blk%center))  n = n + size(blk%center, kind=I8)
    if (allocated(blk%dl))      n = n + size(blk%dl, kind=I8)
    if (allocated(blk%dlmin))   n = n + size(blk%dlmin, kind=I8)
    if (allocated(blk%subVol))  n = n + size(blk%subVol, kind=I8)
    if (allocated(blk%dualW))   n = n + size(blk%dualW, kind=I8)
    if (allocated(blk%cellVol)) n = n + size(blk%cellVol, kind=I8)
  end function geoValues

end module IGLOO_Mod_Footprint
