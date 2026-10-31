!> Wall erosion by particle impact ([IGLOO-Erosion]): the model of each material and its constants, the boundary
!  codes where a parcel's exit hits the wall, and the wall-impact sums behind OUTPUT/wall.tec. Laws: Oka, Okamura
!  and Yoshida (2005), Wear 259, 95-101, with Oka and Yoshida (2005), Wear 259, 102-109; Neilson and Gilchrist
!  (1968), Wear 11, 111-122 (docs/theory/erosion.md).
module IGLOO_Lib_Erosion
  use, intrinsic :: iso_fortran_env, only : R8 => real64
  implicit none
  private
  public :: EROSION_NONE, EROSION_OKA, EROSION_NG
  public :: erosionOn, eroModel, eroModelWord, rhoWall
  public :: eroKey, NERO_KEYS, eroKeyName, eroIn, NERO_WALL, eroWallName, eroWall, eroWallGiven
  public :: okaK, okaK1, okaK2, okaK3, okaN1, okaN2, okaVref, okaDref, ngPhi, ngEps, ngK, ngN
  public :: isWallCode, setupWall, zeroWall, erosionRatio, okaVolume, ngRatio, checkErosionInput, describeErosion

  integer, parameter :: EROSION_NONE = 0   !> the material erodes nothing; its incident flow still counts
  integer, parameter :: EROSION_OKA  = 1   !> Oka et al. (2005)
  integer, parameter :: EROSION_NG   = 2   !> Neilson and Gilchrist (1968)

  !> [IGLOO-Erosion] present: the wall sums are kept and OUTPUT/wall.tec is written.
  logical :: erosionOn = .false.
  !> The model token of each material of the phase file, as read, and its code.
  character(len=64), allocatable :: eroModelWord(:)
  integer,           allocatable :: eroModel(:)
  !> Wall density [kg/m^3], the same for every wall: rdot = sum over the materials of g*er/rhoWall.
  real(R8) :: rhoWall = 0._R8

  !> The constants of the materials as read by read_erosion: one value per material of the phase file (checked by
  !  checkErosionInput), each read for the materials whose model has it, every key of one model only: Oka's K k1 k2
  !  k3 n1 n2 v-ref D-ref, Neilson-Gilchrist's phi epsilon n k (k its threshold velocity, the paper's K). k2, n1 and
  !  n2 are given as values or as the s and q of X = s*Hv^q (Oka et al. 2005: Part 1 Table 3, Part 2 Table 1).
  type :: eroKey
    logical               :: given = .false.
    real(R8), allocatable :: v(:)
  end type eroKey
  integer, parameter :: NERO_KEYS = 18
  integer, parameter :: KY_K = 1, KY_K1 = 2, KY_K2 = 3, KY_K2S = 4, KY_K2Q = 5, KY_K3 = 6, KY_N1 = 7, KY_N1S = 8, &
                        KY_N1Q = 9, KY_N2 = 10, KY_N2S = 11, KY_N2Q = 12, KY_VREF = 13, KY_DREF = 14, KY_PHI = 15, &
                        KY_EPS = 16, KY_N = 17, KY_KNG = 18
  character(len=8), parameter :: eroKeyName(NERO_KEYS) = [character(len=8) :: 'K', 'k1', 'k2', 'k2-s', 'k2-q', &
    'k3', 'n1', 'n1-s', 'n1-q', 'n2', 'n2-s', 'n2-q', 'v-ref', 'D-ref', 'phi', 'epsilon', 'n', 'k']
  type(eroKey) :: eroIn(NERO_KEYS)
  !> The wall constants of Oka's law, one value for every wall: Hv [GPa], and the a and b of its load relaxation
  !  (Part 2 eq. 1); a = b = 1 gives Part 1's E90.
  integer, parameter :: NERO_WALL = 3
  character(len=2), parameter :: eroWallName(NERO_WALL) = [character(len=2) :: 'Hv', 'a', 'b']
  real(R8) :: eroWall(NERO_WALL) = 0._R8
  logical  :: eroWallGiven(NERO_WALL) = .false.

  !> The constants of each material, set by checkErosionInput (0 where its model has none): Oka's K [mm^3/kg],
  !  k1, k2, k3, n1, n2, v' [m/s], D' [m]; Neilson-Gilchrist's phi and epsilon [J/kg], k [m/s] (ngK), n.
  real(R8), allocatable :: okaK(:), okaK1(:), okaK2(:), okaK3(:), okaN1(:), okaN2(:), okaVref(:), okaDref(:)
  real(R8), allocatable :: ngPhi(:), ngEps(:), ngK(:), ngN(:)

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
  !  at speed v [m/s] and angle alpha [rad] from the surface, with diameter d [m] and density rhop [kg/m^3]. Oka's
  !  volume per particle mass [mm^3/kg] times 1e-9 m^3/mm^3 and the wall density; Neilson-Gilchrist's mass ratio.
  pure function erosionRatio(mID, v, alpha, d, rhop) result(er)
    integer,  intent(in) :: mID
    real(R8), intent(in) :: v, alpha, d, rhop
    real(R8) :: er

    select case (eroModel(mID))
    case (EROSION_OKA)
      er = 1.0e-9_R8*rhoWall*okaVolume(v, alpha, d, eroWall(1), eroWall(2), eroWall(3), okaK(mID), okaK1(mID), &
                                       okaK2(mID), okaK3(mID), okaN1(mID), okaN2(mID), okaVref(mID), okaDref(mID))
    case (EROSION_NG)
      er = ngRatio(v, alpha, ngPhi(mID), ngEps(mID), ngK(mID), ngN(mID))
    case default
      er = 0._R8
    end select

  end function erosionRatio


  !> Oka et al. (2005): the wall volume eroded per unit mass of particles [mm^3/kg], E(alpha) = g(alpha) E90
  !  (Part 1 eqs. 1-2, Part 2 eqs. 3-4), g(alpha) = sin(alpha)^n1 (1 + Hv (1 - sin(alpha)))^n2,
  !  E90 = K (a Hv)^(k1 b) (v/vRef)^k2 (d/dRef)^k3: Hv the wall's Vickers hardness [GPa], a and b its load relaxation,
  !  v the impact speed and vRef its reference [m/s], d the particle diameter and dRef its reference (one unit).
  !  Zero at a grazing impact (sin(alpha) = 0) and at rest.
  pure real(R8) function okaVolume(v, alpha, d, Hv, a, b, K, k1, k2, k3, n1, n2, vRef, dRef) result(e)
    real(R8), intent(in) :: v, alpha, d, Hv, a, b, K, k1, k2, k3, n1, n2, vRef, dRef
    real(R8) :: s

    s = sin(alpha)
    if (.not. (s > 0._R8 .and. v > 0._R8)) then
      e = 0._R8
      return
    endif
    e = s**n1 * (1._R8 + Hv*(1._R8 - s))**n2 * K * (a*Hv)**(k1*b) * (v/vRef)**k2 * (d/dRef)**k3

  end function okaVolume


  !> Neilson and Gilchrist (1968), eqs. 3-5: the wall mass eroded per unit mass of particles [kg/kg],
  !  V^2 cos^2(alpha) sin(n alpha)/(2 phi) + (V sin(alpha) - K)^2/(2 eps) for alpha < alpha0 = pi/(2n),
  !  V^2 cos^2(alpha)/(2 phi) + (V sin(alpha) - K)^2/(2 eps) above it; phi and eps the energies of cutting and
  !  deformation per unit eroded mass [J/kg], K the normal velocity "below which no erosion takes place" [m/s]: the
  !  deformation term is zero where V sin(alpha) <= K.
  pure real(R8) function ngRatio(v, alpha, phi, eps, K, n) result(e)
    real(R8), intent(in) :: v, alpha, phi, eps, K, n
    real(R8), parameter  :: halfPi = 2._R8*atan(1._R8)
    real(R8) :: vn

    if (alpha < halfPi/n) then
      e = v*v*cos(alpha)**2*sin(n*alpha)/(2._R8*phi)
    else
      e = v*v*cos(alpha)**2/(2._R8*phi)
    endif
    vn = v*sin(alpha) - K
    if (vn > 0._R8) e = e + vn*vn/(2._R8*eps)

  end function ngRatio


  !> The [IGLOO-Erosion] input against the phase file: one model token per material, each a known model (Oka,
  !  Neilson-Gilchrist, none; any case); every constant given with one value per material; the constants of the
  !  models in use all given, k2, n1 and n2 each once (a value, or its s and q), each in its range for the materials
  !  of its model. Sets the constants of every material.
  subroutine checkErosionInput(nmat)
    integer, intent(in) :: nmat
    integer :: m, i
    logical :: anyOka, anyNG

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
      case ('oka')
        eroModel(m) = EROSION_OKA
      case ('neilson-gilchrist')
        eroModel(m) = EROSION_NG
      case default
        write(*,'(A,I0,A)') ' [ERROR] [IGLOO-Erosion] model of material ', m, ': unknown token "'// &
          trim(adjustl(eroModelWord(m)))//'" (Oka, Neilson-Gilchrist or none)'
        error stop 'IGLOO: unknown [IGLOO-Erosion] model'
      end select
    enddo
    anyOka = any(eroModel == EROSION_OKA)
    anyNG  = any(eroModel == EROSION_NG)

    do i = 1, NERO_KEYS
      if (.not. eroIn(i)%given) cycle
      if (size(eroIn(i)%v) /= nmat) then
        write(*,'(A,I0,A,I0,A)') ' [ERROR] [IGLOO-Erosion] '//trim(eroKeyName(i))//': ', size(eroIn(i)%v), &
          ' value(s) for ', nmat, ' material(s) of the phase file'
        error stop 'IGLOO: [IGLOO-Erosion] a constant needs one value per material'
      endif
    enddo

    if (anyOka) then
      do i = 1, NERO_WALL
        if (.not. eroWallGiven(i)) call missing(eroWallName(i), 'Oka')
      enddo
      if (.not. (eroWall(1) > 0._R8)) call outOfRange('Hv', 0, eroWall(1), '> 0')
      if (.not. (eroWall(2) > 0._R8)) call outOfRange('a', 0, eroWall(2), '> 0')
      if (.not. (eroWall(3) >= 0._R8)) call outOfRange('b', 0, eroWall(3), '>= 0')
      do i = 1, NERO_KEYS
        if (any(i == [KY_K, KY_K1, KY_K3, KY_VREF, KY_DREF]) .and. .not. eroIn(i)%given) &
          call missing(eroKeyName(i), 'Oka')
      enddo
      call oneForm(KY_K2, KY_K2S, KY_K2Q)
      call oneForm(KY_N1, KY_N1S, KY_N1Q)
      call oneForm(KY_N2, KY_N2S, KY_N2Q)
    endif
    if (anyNG) then
      do i = 1, NERO_KEYS
        if (any(i == [KY_PHI, KY_EPS, KY_KNG, KY_N]) .and. .not. eroIn(i)%given) &
          call missing(eroKeyName(i), 'Neilson-Gilchrist')
      enddo
    endif

    call setConstants(okaK); call setConstants(okaK1); call setConstants(okaK2); call setConstants(okaK3)
    call setConstants(okaN1); call setConstants(okaN2); call setConstants(okaVref); call setConstants(okaDref)
    call setConstants(ngPhi); call setConstants(ngEps); call setConstants(ngK); call setConstants(ngN)
    do m = 1, nmat
      select case (eroModel(m))
      case (EROSION_OKA)
        okaK(m)    = inRange(KY_K, m, '> 0', eroIn(KY_K)%v(m) > 0._R8)
        okaK1(m)   = inRange(KY_K1, m, 'finite', .true.)
        okaK2(m)   = byForm(KY_K2, KY_K2S, KY_K2Q, m, .true.)
        okaK3(m)   = inRange(KY_K3, m, 'finite', .true.)
        okaN1(m)   = byForm(KY_N1, KY_N1S, KY_N1Q, m, .true.)
        okaN2(m)   = byForm(KY_N2, KY_N2S, KY_N2Q, m, .false.)
        okaVref(m) = inRange(KY_VREF, m, '> 0', eroIn(KY_VREF)%v(m) > 0._R8)
        okaDref(m) = inRange(KY_DREF, m, '> 0', eroIn(KY_DREF)%v(m) > 0._R8)
      case (EROSION_NG)
        ngPhi(m) = inRange(KY_PHI, m, '> 0', eroIn(KY_PHI)%v(m) > 0._R8)
        ngEps(m) = inRange(KY_EPS, m, '> 0', eroIn(KY_EPS)%v(m) > 0._R8)
        ngK(m)   = inRange(KY_KNG, m, '>= 0', eroIn(KY_KNG)%v(m) >= 0._R8)
        ngN(m)   = inRange(KY_N, m, '> 0', eroIn(KY_N)%v(m) > 0._R8)
      end select
    enddo

  contains

    subroutine setConstants(c)
      real(R8), allocatable, intent(inout) :: c(:)
      if (allocated(c)) deallocate(c)
      allocate(c(nmat))
      c = 0._R8
    end subroutine setConstants

    subroutine missing(key, model)
      character(len=*), intent(in) :: key, model
      write(*,'(A)') ' [ERROR] [IGLOO-Erosion] '//trim(key)//' is missing: the '//model//' model needs it'
      error stop 'IGLOO: [IGLOO-Erosion] a constant of the models in use is missing'
    end subroutine missing

    !> X given as a value or as s and q (X = s*Hv^q), once: the value alone, or s and q together.
    subroutine oneForm(kv, ks, kq)
      integer, intent(in) :: kv, ks, kq
      logical :: direct, corr, partial
      direct  = eroIn(kv)%given
      corr    = eroIn(ks)%given .and. eroIn(kq)%given
      partial = eroIn(ks)%given .neqv. eroIn(kq)%given
      if ((direct .neqv. corr) .and. .not. partial) return
      write(*,'(A)') ' [ERROR] [IGLOO-Erosion] Oka: give '//trim(eroKeyName(kv))//', or '//trim(eroKeyName(ks))// &
        ' with '//trim(eroKeyName(kq))//' ('//trim(eroKeyName(kv))//' = s*Hv^q), once'
      error stop 'IGLOO: [IGLOO-Erosion] give k2, n1 and n2 each as a value or as its s and q, once'
    end subroutine oneForm

    !> The value of key ky for material m when ok (its range) and finite.
    real(R8) function inRange(ky, mm, range, ok) result(x)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      integer,          intent(in) :: ky, mm
      character(len=*), intent(in) :: range
      logical,          intent(in) :: ok
      x = eroIn(ky)%v(mm)
      if (.not. (ok .and. ieee_is_finite(x))) call outOfRange(eroKeyName(ky), mm, x, range)
    end function inRange

    !> k2, n1 or n2 of material m: its value, or s*Hv^q; positive when pos.
    real(R8) function byForm(kv, ks, kq, mm, pos) result(x)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      integer, intent(in) :: kv, ks, kq, mm
      logical, intent(in) :: pos
      character(len=8) :: what
      if (eroIn(kv)%given) then
        x = eroIn(kv)%v(mm); what = eroKeyName(kv)
      else
        x = eroIn(ks)%v(mm)*eroWall(1)**eroIn(kq)%v(mm); what = eroKeyName(ks)
      endif
      if (.not. ieee_is_finite(x) .or. (pos .and. .not. (x > 0._R8))) then
        if (pos) then
          call outOfRange(what, mm, x, '> 0 ('//trim(eroKeyName(kv))//')')
        else
          call outOfRange(what, mm, x, 'finite ('//trim(eroKeyName(kv))//')')
        endif
      endif
    end function byForm

    subroutine outOfRange(key, mm, x, range)
      character(len=*), intent(in) :: key, range
      integer,          intent(in) :: mm
      real(R8),         intent(in) :: x
      if (mm > 0) then
        write(*,'(A,I0,A,ES12.5,A)') ' [ERROR] [IGLOO-Erosion] '//trim(key)//' of material ', mm, ': ', x, &
          ', expected '//range
      else
        write(*,'(A,ES12.5,A)') ' [ERROR] [IGLOO-Erosion] '//trim(key)//': ', x, ', expected '//range
      endif
      error stop 'IGLOO: [IGLOO-Erosion] a constant out of its range'
    end subroutine outOfRange

  end subroutine checkErosionInput


  !> The constants of material m as the setup log names them ('' for none).
  function describeErosion(m) result(txt)
    integer, intent(in) :: m
    character(len=:), allocatable :: txt
    character(len=320) :: buf

    select case (eroModel(m))
    case (EROSION_OKA)
      write(buf,'(8(A,ES12.5))') 'K ', okaK(m), ', k1 ', okaK1(m), ', k2 ', okaK2(m), ', k3 ', okaK3(m), &
        ', n1 ', okaN1(m), ', n2 ', okaN2(m), ', v-ref ', okaVref(m), ', D-ref ', okaDref(m)
      txt = trim(buf)
    case (EROSION_NG)
      write(buf,'(4(A,ES12.5))') 'phi ', ngPhi(m), ', epsilon ', ngEps(m), ', k ', ngK(m), ', n ', ngN(m)
      txt = trim(buf)
    case default
      txt = ''
    end select

  end function describeErosion


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
