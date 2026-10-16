!> Conservative binomial smoothing of cell-centered feedback fields.
module IGLOO_Lib_Mollify
  use, intrinsic :: iso_fortran_env, only : R8 => real64
  implicit none
  private

  public :: binomial_smooth
  public :: mollifyBlocks
  public :: DEFAULT_MOLLIFY_PASSES

  !> A cell of block b on its 101 or 201 face f joined to cell (pi,pj,pk) of block pb on face pf; dir
  !  is the index direction normal to both faces (0 when they are normal to different directions).
  type, public :: cellLink
    integer :: b = 0, f = 0, i = 0, j = 0, k = 0
    integer :: pb = 0, pf = 0, pi = 0, pj = 0, pk = 0
    integer :: dir = 0
  end type cellLink

  !> One block's cell array, for mollifyBlocks.
  type, public :: blockField
    real(R8), allocatable :: q(:,:,:)
  end type blockField

  !> A block face's cells (m,n) mapped to their aligned link (0 = none).
  type :: faceLinkMap
    integer, allocatable :: l(:,:)
  end type faceLinkMap

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


  !> binomial_smooth over several blocks at once: before every sweep direction of every pass a
  !  one-cell halo takes the partner's current value and volume on each aligned 101/201 face cell,
  !  so the flux crosses the face as it crosses an interior face; the conservation guard sums every
  !  block. Without an aligned link each block goes through binomial_smooth on its own.
  subroutine mollifyBlocks(q, V, links, npass, ndim)
    implicit none
    type(blockField), intent(inout) :: q(:)
    type(blockField), intent(in)    :: V(:)
    type(cellLink),   intent(in)    :: links(:)
    integer,          intent(in)    :: npass, ndim
    type(blockField),  allocatable :: w(:)
    type(faceLinkMap), allocatable :: fmap(:,:)
    real(R8), allocatable :: hq(:), hv(:)
    real(R8) :: s0, s1, denom
    integer  :: nbk, b, f, l, p, d
    logical  :: inQ

    if (npass <= 0) return
    nbk = size(q)
    if (.not. any(links(:)%dir > 0)) then
      do b = 1, nbk
        call binomial_smooth(q(b)%q, V(b)%q, npass, ndim)
      enddo
      return
    endif

    do b = 1, nbk
      if (any(V(b)%q <= tiny(1._R8))) error stop 'IGLOO_Lib_Mollify::mollifyBlocks: non-positive cell volume'
    enddo

    !> face maps of the aligned links and the halo volumes
    allocate(fmap(6, nbk))
    do b = 1, nbk
      associate(Nx => size(q(b)%q,1), Ny => size(q(b)%q,2), Nz => size(q(b)%q,3))
      allocate(fmap(1,b)%l(Ny,Nz), fmap(2,b)%l(Ny,Nz), fmap(3,b)%l(Nx,Nz), fmap(4,b)%l(Nx,Nz), &
               fmap(5,b)%l(Nx,Ny), fmap(6,b)%l(Nx,Ny))
      end associate
      do f = 1, 6
        fmap(f,b)%l = 0
      enddo
    enddo
    allocate(hq(size(links)), hv(size(links)))
    hq = 0._R8; hv = 0._R8
    do l = 1, size(links)
      if (links(l)%dir == 0) cycle
      associate(lk => links(l))
      select case (lk%f)
      case (1, 2); fmap(lk%f, lk%b)%l(lk%j, lk%k) = l
      case (3, 4); fmap(lk%f, lk%b)%l(lk%i, lk%k) = l
      case (5, 6); fmap(lk%f, lk%b)%l(lk%i, lk%j) = l
      end select
      hv(l) = V(lk%pb)%q(lk%pi, lk%pj, lk%pk)
      end associate
    enddo

    allocate(w(nbk))
    s0 = 0._R8
    do b = 1, nbk
      allocate(w(b)%q, mold=q(b)%q)
      s0 = s0 + sum(q(b)%q*V(b)%q)
    enddo

    !> ping-pong between q and w; inQ: the current state is in q
    inQ = .true.
    do p = 1, npass
      do d = 1, merge(3, 2, ndim >= 3)
        do l = 1, size(links)
          if (links(l)%dir /= d) cycle
          associate(lk => links(l))
          if (inQ) then
            hq(l) = q(lk%pb)%q(lk%pi, lk%pj, lk%pk)
          else
            hq(l) = w(lk%pb)%q(lk%pi, lk%pj, lk%pk)
          endif
          end associate
        enddo
        do b = 1, nbk
          if (inQ) then
            call sweepLinked(q(b)%q, w(b)%q, V(b)%q, d, fmap(:,b), hq, hv)
          else
            call sweepLinked(w(b)%q, q(b)%q, V(b)%q, d, fmap(:,b), hq, hv)
          endif
        enddo
        inQ = .not. inQ
      enddo
    enddo
    if (.not. inQ) then
      do b = 1, nbk
        q(b)%q = w(b)%q
      enddo
    endif

    !> Conservation guard over every block, relative to the L1 mass Sum(|q|*V).
    s1 = 0._R8; denom = 0._R8
    do b = 1, nbk
      s1    = s1    + sum(q(b)%q*V(b)%q)
      denom = denom + sum(abs(q(b)%q)*V(b)%q)
    enddo
    denom = max(denom, tiny(1._R8))
    if (abs(s1 - s0)/denom > 1.e-10_R8) then
      write(*,'(A,3ES23.15)') &
        ' IGLOO_Lib_Mollify: Sum(q*V) over the blocks NOT conserved  s0,s1,reldrift =', &
        s0, s1, abs(s1 - s0)/denom
      error stop 'IGLOO_Lib_Mollify::mollifyBlocks conservation drift'
    endif

  end subroutine mollifyBlocks


  !> One sweep of binomial_smooth in direction d on one block, with the flux at a linked face cell
  !  taken against the halo value hq and volume hv of its link.
  subroutine sweepLinked(src, dst, V, d, fmap, hq, hv)
    implicit none
    real(R8),          intent(in)  :: src(:,:,:), V(:,:,:)
    real(R8),          intent(out) :: dst(:,:,:)
    integer,           intent(in)  :: d
    type(faceLinkMap), intent(in)  :: fmap(6)
    real(R8),          intent(in)  :: hq(:), hv(:)
    real(R8) :: phiL, phiR
    integer  :: Nx, Ny, Nz, i, j, k, l
    logical  :: useOMP

    Nx = size(src,1); Ny = size(src,2); Nz = size(src,3)
    useOMP = (Nx*Ny*Nz > MOLLIFY_OMP_MIN)

    select case (d)
    case (1)
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,l,phiL,phiR) IF(useOMP)
      do k = 1, Nz
        do j = 1, Ny
          do i = 1, Nx
            phiR = 0._R8; phiL = 0._R8
            if (i < Nx) then
              phiR = theta*min(V(i,j,k),  V(i+1,j,k))*(src(i+1,j,k) - src(i,j,k))
            else
              l = fmap(2)%l(j,k)
              if (l > 0) phiR = theta*min(V(i,j,k), hv(l))*(hq(l) - src(i,j,k))
            endif
            if (i > 1) then
              phiL = theta*min(V(i-1,j,k),V(i,j,k)  )*(src(i,j,k)   - src(i-1,j,k))
            else
              l = fmap(1)%l(j,k)
              if (l > 0) phiL = theta*min(hv(l), V(i,j,k))*(src(i,j,k) - hq(l))
            endif
            dst(i,j,k) = src(i,j,k) + phiR/V(i,j,k) - phiL/V(i,j,k)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
    case (2)
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,l,phiL,phiR) IF(useOMP)
      do k = 1, Nz
        do j = 1, Ny
          do i = 1, Nx
            phiR = 0._R8; phiL = 0._R8
            if (j < Ny) then
              phiR = theta*min(V(i,j,k),  V(i,j+1,k))*(src(i,j+1,k) - src(i,j,k))
            else
              l = fmap(4)%l(i,k)
              if (l > 0) phiR = theta*min(V(i,j,k), hv(l))*(hq(l) - src(i,j,k))
            endif
            if (j > 1) then
              phiL = theta*min(V(i,j-1,k),V(i,j,k)  )*(src(i,j,k)   - src(i,j-1,k))
            else
              l = fmap(3)%l(i,k)
              if (l > 0) phiL = theta*min(hv(l), V(i,j,k))*(src(i,j,k) - hq(l))
            endif
            dst(i,j,k) = src(i,j,k) + phiR/V(i,j,k) - phiL/V(i,j,k)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
    case (3)
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,l,phiL,phiR) IF(useOMP)
      do k = 1, Nz
        do j = 1, Ny
          do i = 1, Nx
            phiR = 0._R8; phiL = 0._R8
            if (k < Nz) then
              phiR = theta*min(V(i,j,k),  V(i,j,k+1))*(src(i,j,k+1) - src(i,j,k))
            else
              l = fmap(6)%l(i,j)
              if (l > 0) phiR = theta*min(V(i,j,k), hv(l))*(hq(l) - src(i,j,k))
            endif
            if (k > 1) then
              phiL = theta*min(V(i,j,k-1),V(i,j,k)  )*(src(i,j,k)   - src(i,j,k-1))
            else
              l = fmap(5)%l(i,j)
              if (l > 0) phiL = theta*min(hv(l), V(i,j,k))*(src(i,j,k) - hq(l))
            endif
            dst(i,j,k) = src(i,j,k) + phiR/V(i,j,k) - phiL/V(i,j,k)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
    end select

  end subroutine sweepLinked

end module IGLOO_Lib_Mollify
