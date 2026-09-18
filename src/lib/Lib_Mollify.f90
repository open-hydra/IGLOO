!> Conservative binomial smoothing of cell-centered feedback fields.
module IGLOO_Lib_Mollify
  use, intrinsic :: iso_fortran_env, only : R8 => real64
  implicit none
  private

  public :: binomial_smooth
  public :: DEFAULT_MOLLIFY_PASSES

  !> Per-sweep Laplacian coefficient: 1/4 gives the [1/4, 1/2, 1/4] binomial filter.
  real(R8), parameter :: theta = 0.25_R8

  !> Passes used when mollification is on and `mollify-passes` is absent (~2-cell width).
  integer, parameter :: DEFAULT_MOLLIFY_PASSES = 8

  !> Cell count below which the OpenMP region runs serially.
  integer, parameter :: MOLLIFY_OMP_MIN = 4096

contains

  !> Conservative dimension-split binomial smoother for a cell-centered density q (extensive
  !  quantities are passed as q/V): npass sweeps over x, y and, for ndim = 3, z, in gather form
  !  on ping-pong buffers with zero-flux block boundaries; `work` is an optional scratch buffer.
  subroutine binomial_smooth(q, V, npass, ndim, work)
    implicit none
    real(R8), intent(inout),         target   :: q(:,:,:)
    real(R8), intent(in)                      :: V(:,:,:)
    integer,  intent(in)                      :: npass, ndim
    real(R8), intent(inout), optional, target :: work(:,:,:)
    real(R8), allocatable, target :: buf(:,:,:)
    real(R8), pointer             :: src(:,:,:), dst(:,:,:), swp(:,:,:)
    real(R8) :: phiL, phiR, s0, s1, denom
    integer  :: Nx, Ny, Nz, i, j, k, p
    logical  :: useOMP

    if (npass <= 0) return
    Nx = size(q,1); Ny = size(q,2); Nz = size(q,3)

    if (any(V <= tiny(1._R8))) &
      error stop 'IGLOO_Lib_Mollify::binomial_smooth: non-positive cell volume'

    !> Ping-pong buffer: the caller's scratch when given, else a local one.
    if (present(work)) then
      if (size(work,1) /= Nx .or. size(work,2) /= Ny .or. size(work,3) /= Nz) &
        error stop 'IGLOO_Lib_Mollify::binomial_smooth: work scratch shape mismatch'
      dst => work
    else
      allocate(buf(Nx,Ny,Nz))
      dst => buf
    endif
    src => q

    useOMP = (Nx*Ny*Nz > MOLLIFY_OMP_MIN)
    s0 = sum(q*V)

    do p = 1, npass
      !> x-sweep
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,phiL,phiR) IF(useOMP)
      do k = 1, Nz
        do j = 1, Ny
          do i = 1, Nx
            phiR = 0._R8; phiL = 0._R8
            if (i < Nx) phiR = theta*min(V(i,j,k),  V(i+1,j,k))*(src(i+1,j,k) - src(i,j,k))
            if (i > 1 ) phiL = theta*min(V(i-1,j,k),V(i,j,k)  )*(src(i,j,k)   - src(i-1,j,k))
            dst(i,j,k) = src(i,j,k) + phiR/V(i,j,k) - phiL/V(i,j,k)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      swp => src; src => dst; dst => swp
      !> y-sweep
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,phiL,phiR) IF(useOMP)
      do k = 1, Nz
        do j = 1, Ny
          do i = 1, Nx
            phiR = 0._R8; phiL = 0._R8
            if (j < Ny) phiR = theta*min(V(i,j,k),  V(i,j+1,k))*(src(i,j+1,k) - src(i,j,k))
            if (j > 1 ) phiL = theta*min(V(i,j-1,k),V(i,j,k)  )*(src(i,j,k)   - src(i,j-1,k))
            dst(i,j,k) = src(i,j,k) + phiR/V(i,j,k) - phiL/V(i,j,k)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      swp => src; src => dst; dst => swp
      !> z-sweep (3D only)
      if (ndim >= 3) then
        !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,phiL,phiR) IF(useOMP)
        do k = 1, Nz
          do j = 1, Ny
            do i = 1, Nx
              phiR = 0._R8; phiL = 0._R8
              if (k < Nz) phiR = theta*min(V(i,j,k),  V(i,j,k+1))*(src(i,j,k+1) - src(i,j,k))
              if (k > 1 ) phiL = theta*min(V(i,j,k-1),V(i,j,k)  )*(src(i,j,k)   - src(i,j,k-1))
              dst(i,j,k) = src(i,j,k) + phiR/V(i,j,k) - phiL/V(i,j,k)
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO
        swp => src; src => dst; dst => swp
      endif
    enddo

    !> The result lives in src after the last swap; copy it back into q if needed.
    if (.not. associated(src, q)) q = src
    if (allocated(buf)) deallocate(buf)
    src => null(); dst => null()

    !> Conservation guard, relative to the L1 mass Sum(|q|*V) (not the signed net).
    s1    = sum(q*V)
    denom = max(sum(abs(q)*V), tiny(1._R8))
    if (abs(s1 - s0)/denom > 1.e-10_R8) then
      write(*,'(A,3ES23.15)') &
        ' IGLOO_Lib_Mollify: Sum(q*V) NOT conserved  s0,s1,reldrift =', &
        s0, s1, abs(s1 - s0)/denom
      error stop 'IGLOO_Lib_Mollify::binomial_smooth conservation drift'
    endif

  end subroutine binomial_smooth

end module IGLOO_Lib_Mollify
