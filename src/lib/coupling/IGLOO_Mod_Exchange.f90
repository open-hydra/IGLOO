!> What IGLOO exchanges with the solvers around it: the only module a coupler needs besides obj_IGLOO.
!
!  Sending (read after solve, per material and geometry cell, [per second, the cell's]):
!    exchange_rate    the mass the parcels give the gas, the carried material they deposit and give up, the mass
!                     their process takes off (with the state it leaves at, when tracked), what they absorbed and
!                     what a removal took off them
!    gas_rate         the momentum and energy all parcels give the gas in a cell
!    track_processed  ask for the state (velocity, temperature) the processed mass of a material leaves at
!
!  Receiving (set before solve; kept until clear_import_rates):
!    add_absorption_rate  mass offered to the parcels of a material in a cell (from the gas or any other phase); the
!                         parcels absorb it in proportion to the mass they hold there, and exchange_rate reports what
!                         they took, which the coupler removes from the donor
!    absorption_packet    the offer of a mass at a velocity and temperature, in IGLOO's enthalpy datum
!    add_removal_rate     mass taken off the parcels of a material in a cell (for the gas or any other phase): whole
!                         droplets leave, in proportion to the mass the parcels hold there, and exchange_rate reports
!                         what left, which the coupler gives the receiver
!    add_birth_rate       new parcels of a born-only material (a line with no injection groups), born in a cell
!    clear_import_rates   drop every offer, removal and birth
!
!  material_id gives a material's index from its name, gen_diameter the diameter (d-gen of its combustion row) the
!  coupler gives the drops of the mass a burning material generates in another phase. A refusal comes back in plain words in `why` (empty when
!  accepted); nothing is changed then.
module IGLOO_Mod_Exchange
  use, intrinsic :: iso_fortran_env, only : R8 => real64
  use IGLOO_module, only: obj_IGLOO
  implicit none
  private
  public :: exchange_t
  public :: material_id, gen_diameter, exchange_rate, gas_rate, track_processed
  public :: absorption_packet, add_absorption_rate, add_removal_rate, add_birth_rate, clear_import_rates

  !> One material's exchange in one cell [per second].
  type :: exchange_t
    real(R8) :: mass        = 0._R8   !> mass the parcels give the gas, net [kg/s]
    real(R8) :: carried_in  = 0._R8   !> carried material the burning deposits on the parcels [kg/s]
    real(R8) :: carried_out = 0._R8   !> carried material the parcels give the gas when they burn out [kg/s]
    real(R8) :: processed   = 0._R8   !> mass the process (burning, evaporation) takes off the parcels [kg/s]
    logical  :: tracked     = .false. !> u and T below hold the state the processed mass leaves at
    real(R8) :: u(3)        = 0._R8   !> its mean velocity [m/s]
    real(R8) :: T           = 0._R8   !> its mean temperature [K]
    real(R8) :: absorbed(5) = 0._R8   !> what the parcels absorbed: mass [kg/s], momentum (3) [N], energy [W]
    real(R8) :: removed(5)  = 0._R8   !> what a removal took off them: mass [kg/s], momentum (3) [N], energy [W]
  end type exchange_t

contains

  !> The index of the material named `name`, 0 if there is none.
  integer function material_id(igloo, name) result(m)
    type(obj_IGLOO),  intent(in) :: igloo
    character(len=*), intent(in) :: name
    integer :: k

    m = 0
    do k = 1, size(igloo%material)
      if (trim(igloo%material(k)%matName) == trim(name)) m = k
    enddo
  end function material_id


  !> d-gen of material m's combustion row [m]: the diameter the coupler gives the drops of the mass m generates in
  !  another phase; 0 when the row gives none. IGLOO itself does not use it.
  real(R8) function gen_diameter(igloo, m) result(d)
    type(obj_IGLOO), intent(in) :: igloo
    integer,         intent(in) :: m

    d = igloo%material(m)%dGen
  end function gen_diameter


  !> Material m's exchange in geometry cell (i, j, k) of block b, after solve.
  function exchange_rate(igloo, m, b, i, j, k) result(x)
    use IGLOO_variables, only: carryColumns, procColumns, trackProc, rowCarry, rowProc
    type(obj_IGLOO), intent(in) :: igloo
    integer,         intent(in) :: m, b, i, j, k
    type(exchange_t) :: x
    real(R8) :: w

    if (.not. allocated(igloo%source)) return
    associate (src => igloo%source(b))
    if (.not. allocated(src%sourceMass)) return
    x%mass = src%sourceMass(m, i, j, k)
    if (carryColumns) then
      x%carried_in  = src%sourceMass(rowCarry(m, 1), i, j, k)
      x%carried_out = src%sourceMass(rowCarry(m, 2), i, j, k)
    endif
    if (procColumns .and. size(src%sourceMass, 1) >= rowProc(m, 5)) then
      w = src%sourceMass(rowProc(m, 1), i, j, k)
      x%processed = w
      x%tracked   = trackProc(m) .and. w > 0._R8
      if (x%tracked) then
        x%u = src%sourceMass(rowProc(m, 2):rowProc(m, 4), i, j, k)/w
        x%T = src%sourceMass(rowProc(m, 5), i, j, k)/w
      endif
    endif
    if (allocated(src%absorbed)) x%absorbed = src%absorbed(:, m, i, j, k)
    if (allocated(src%removed))  x%removed  = src%removed(:, m, i, j, k)
    end associate
  end function exchange_rate


  !> The momentum [N] (1:3) and energy [W] (4) all parcels give the gas in geometry cell (i, j, k) of block b, after
  !  solve. What they absorbed is not in it: that came from the donor.
  function gas_rate(igloo, b, i, j, k) result(q)
    type(obj_IGLOO), intent(in) :: igloo
    integer,         intent(in) :: b, i, j, k
    real(R8) :: q(4)

    q = 0._R8
    if (.not. allocated(igloo%source)) return
    associate (src => igloo%source(b))
    if (.not. allocated(src%sourceMom)) return
    q(1:3) = src%sourceMom(:, i, j, k)
    q(4)   = src%sourceEn(i, j, k)
    end associate
  end function gas_rate


  !> Ask for the state the processed mass of material m leaves at (exchange_rate's u and T), from the next solve on.
  subroutine track_processed(igloo, m, why)
    use IGLOO_variables, only: procColumns, trackProc
    type(obj_IGLOO),               intent(inout) :: igloo
    integer,                       intent(in)    :: m
    character(len=:), allocatable, intent(out)   :: why

    why = refuse_material(igloo, m)
    if (len(why) > 0) return
    associate (mat => igloo%material(m))
    if (mat%evapSelect == 0 .and. mat%combSelect == 0) then
      why = trim(mat%matName)//' neither burns nor evaporates: its parcels process no mass'
      return
    endif
    end associate
    trackProc(m) = .true.
    procColumns  = .true.
  end subroutine track_processed


  !> The offer of `mass` [kg/s] at velocity u [m/s] and temperature T [K] to material m: [mass, momentum (3),
  !  energy], the energy at m's own enthalpy (its datum) plus the kinetic part.
  function absorption_packet(igloo, m, mass, u, T) result(p)
    use IGLOO_Lib_Properties, only: lookupTab
    type(obj_IGLOO), intent(in) :: igloo
    integer,         intent(in) :: m
    real(R8),        intent(in) :: mass, u(3), T
    real(R8) :: p(5), h

    associate (mat => igloo%material(m))
    if (mat%cpVariable) then
      h = lookupTab(mat%hTab, T)
    else
      h = mat%cp*T + mat%hOff
    endif
    end associate
    p(1)   = mass
    p(2:4) = mass*u
    p(5)   = mass*(h + 0.5_R8*dot_product(u, u))
  end function absorption_packet


  !> Offer `packet` [mass kg/s, momentum (3) N, energy W] (absorption_packet) to the parcels of material m in
  !  geometry cell (i, j, k) of block b. They absorb it in proportion to the mass they held there in the last sweep:
  !  the first sweep after an offer absorbs nothing. Offers add up until clear_import_rates.
  subroutine add_absorption_rate(igloo, m, b, i, j, k, packet, why)
    use IGLOO_variables, only: nm, absColumns, absorbs
    type(obj_IGLOO),               intent(inout) :: igloo
    integer,                       intent(in)    :: m, b, i, j, k
    real(R8),                      intent(in)    :: packet(5)
    character(len=:), allocatable, intent(out)   :: why
    integer :: c

    why = refuse_cell(igloo, m, b, i, j, k)
    if (len(why) > 0) return
    if (.not. packet(1) > 0._R8) then
      why = 'the mass offered must be positive'
      return
    endif
    select case (parcel_model(igloo, m))
    case (3)
      why = trim(igloo%material(m)%matName)//' breaks up without evaporating: its parcels keep a fixed mass'
    case (6, 7)
      why = trim(igloo%material(m)%matName)//' freezes: absorbing would change its frozen fraction'
    case (8)
      why = trim(igloo%material(m)%matName)//' carries a second material: its parcels cannot absorb'
    end select
    if (len(why) > 0) return
    do c = 1, size(igloo%source)
      associate (src => igloo%source(c))
      if (.not. allocated(src%offered)) then
        allocate(src%offered(5, nm, src%Nx, src%Ny, src%Nz))
        src%offered = 0._R8
      endif
      end associate
    enddo
    igloo%source(b)%offered(:, m, i, j, k) = igloo%source(b)%offered(:, m, i, j, k) + packet
    absorbs(m) = .true.
    absColumns = .true.
  end subroutine add_absorption_rate


  !> Take `mass` [kg/s] off the parcels of material m in geometry cell (i, j, k) of block b: whole droplets leave at
  !  their own state, in proportion to the mass the material held there in the last sweep (the first sweep after a
  !  request takes nothing), so the droplets left are unchanged. A parcel whose share is all of it leaves the run.
  !  exchange_rate's removed is what left. Requests add up until clear_import_rates.
  subroutine add_removal_rate(igloo, m, b, i, j, k, mass, why)
    use IGLOO_variables, only: nm, absColumns, removes
    type(obj_IGLOO),               intent(inout) :: igloo
    integer,                       intent(in)    :: m, b, i, j, k
    real(R8),                      intent(in)    :: mass
    character(len=:), allocatable, intent(out)   :: why
    integer :: c

    why = refuse_cell(igloo, m, b, i, j, k)
    if (len(why) > 0) return
    if (.not. mass > 0._R8) then
      why = 'the mass removed must be positive'
      return
    endif
    do c = 1, size(igloo%source)
      associate (src => igloo%source(c))
      if (.not. allocated(src%asked)) then
        allocate(src%asked(nm, src%Nx, src%Ny, src%Nz))
        src%asked = 0._R8
      endif
      end associate
    enddo
    igloo%source(b)%asked(m, i, j, k) = igloo%source(b)%asked(m, i, j, k) + mass
    removes(m) = .true.
    absColumns = .true.
  end subroutine add_removal_rate


  !> New parcels of born-only material m (a line with no injection groups): `mass` [kg/s] born at the centre of
  !  geometry cell (i, j, k) of block b, with velocity u [m/s], temperature T [K] and diameter d [m], one stream.
  !  Births add up until clear_import_rates.
  subroutine add_birth_rate(igloo, m, b, i, j, k, mass, u, T, d, why)
    use IGLOO_data_phases, only: obj_shed
    type(obj_IGLOO),               intent(inout) :: igloo
    integer,                       intent(in)    :: m, b, i, j, k
    real(R8),                      intent(in)    :: mass, u(3), T, d
    character(len=:), allocatable, intent(out)   :: why
    type(obj_shed) :: rec

    why = refuse_cell(igloo, m, b, i, j, k)
    if (len(why) > 0) return
    associate (mat => igloo%material(m))
    if (.not. mat%bornOnly) then
      why = trim(mat%matName)//' is injected: births go to a line with no injection groups'
    elseif (mat%brkupHasChild) then
      why = trim(mat%matName)//' breaks up into children: its births would mix with their generations'
    elseif (.not. (mass > 0._R8 .and. d > 0._R8 .and. T > 0._R8)) then
      why = 'a birth needs a positive mass, diameter and temperature'
    endif
    if (len(why) > 0) return
    rec%ipos  = [b, i, j, k]
    rec%igas  = [i, j, k]          !> its gas cell at gas order 1; at order 2 drainBirths finds the dual cell
    rec%pos   = igloo%geoblock(b)%center(:, i, j, k)
    rec%vel   = u
    rec%temp  = T
    rec%diam  = d
    rec%npdot = 1._R8
    rec%mass  = mass
    rec%time  = 0._R8
    call mat%births%push(rec)
    end associate
  end subroutine add_birth_rate


  !> Drop every offer to absorb, every removal and every birth.
  subroutine clear_import_rates(igloo)
    use IGLOO_variables, only: absorbs, removes
    type(obj_IGLOO), intent(inout) :: igloo
    integer :: c

    if (allocated(igloo%source)) then
      do c = 1, size(igloo%source)
        if (allocated(igloo%source(c)%offered)) igloo%source(c)%offered = 0._R8
        if (allocated(igloo%source(c)%asked))   igloo%source(c)%asked   = 0._R8
      enddo
    endif
    if (allocated(absorbs)) absorbs = .false.
    if (allocated(removes)) removes = .false.
    do c = 1, size(igloo%material)
      igloo%material(c)%births%n = 0
    enddo
  end subroutine clear_import_rates


  !> '' for a material index of the run, else why not.
  function refuse_material(igloo, m) result(why)
    type(obj_IGLOO), intent(in) :: igloo
    integer,         intent(in) :: m
    character(len=:), allocatable :: why

    why = ''
    if (m < 1 .or. m > size(igloo%material)) why = 'no material of this run has that index'
  end function refuse_material


  !> '' for a material index and a geometry cell of the run, else why not.
  function refuse_cell(igloo, m, b, i, j, k) result(why)
    type(obj_IGLOO), intent(in) :: igloo
    integer,         intent(in) :: m, b, i, j, k
    character(len=:), allocatable :: why

    why = refuse_material(igloo, m)
    if (len(why) > 0) return
    if (b < 1 .or. b > size(igloo%geoblock)) then
      why = 'no block of this run has that index'
    elseif (i < 1 .or. i > igloo%geoblock(b)%Nx .or. j < 1 .or. j > igloo%geoblock(b)%Ny .or. &
            k < 1 .or. k > igloo%geoblock(b)%Nz) then
      why = 'the cell lies outside its block'
    endif
  end function refuse_cell


  !> The ODE model material m's parcels integrate (Lib_RHS determineModel).
  integer function parcel_model(igloo, m) result(model)
    use Lib_RHS, only: determineModel
    type(obj_IGLOO), intent(in) :: igloo
    integer,         intent(in) :: m

    associate (mat => igloo%material(m))
    model = determineModel(mat%evapSelect > 0, mat%brkupEqOde, mat%combSelect, mat%solidSelect, &
                           carry=mat%carryShare > 0._R8 .or. mat%carryLaw > 0)
    end associate
  end function parcel_model

end module IGLOO_Mod_Exchange
