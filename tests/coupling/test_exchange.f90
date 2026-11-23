!> The exchange channels a coupler sees (src/lib/coupling/IGLOO_Mod_Exchange.f90), driven the way a coupler drives
!  them: setup_static, then reset_state + solve per sweep, the calls of the module in between. Run in a case
!  directory; the argument picks the channel:
!
!    track   burn-box's deck (Beckstead burning, gas at 10 m/s, 600 K). track_processed on A:
!            X1 the processed mass equals the net mass the parcels give the gas, cell by cell (nothing is carried);
!            X2 it leaves at the parcels' velocity, the gas's (the parcels coast at u_g); X3 at a temperature
!            between the gas's and the burning drops' ceiling; X4 a material index outside the run is refused; X5
!            gas_rate: the drops coast, so the only momentum they give the gas is the one of the mass they give it.
!    absorb  drag-stokes's deck (inert drops, ODE model 1). The same offer (mass, at the gas velocity and
!            temperature) in cells 21-40 of the box for three sweeps, then none:
!            A1 sweep 1 absorbs nothing (nothing was held before it); A2 sweep 2 absorbs the offer within 4 % in
!            the cells the drops cross (the share is last sweep's held mass, which held none of the offer yet) and
!            sweep 3 within 0.5 %; A3 what was absorbed carries the offer's velocity and energy per kg, cell by
!            cell; A4 the drops leave with the injected flow plus what they absorbed; A5 after clear_import_rates
!            the run is sweep 1's again, bit for bit; A6 the refusals (a cell outside the block, a mass that is
!            not positive).
!    births  zero-group-line's deck (A injected, B a 0-group line). Two births of B:
!            B1 an injected material is refused, and a birth without a diameter; B2 the two births are integrated
!            as two parcels of B, IDs 1 and 2; B3 each reaches the outlet with its mass flow (inert drops); B4 a
!            second sweep integrates them again; B5 after clear_import_rates B has no parcel.
!
!  Exit 0 when every gate holds, 1 otherwise.
program test_exchange
  use, intrinsic :: iso_fortran_env, only: R8 => real64
  use IGLOO_module,       only: obj_IGLOO
  use IGLOO_Mod_MPI,      only: mpi_init_env, mpi_finalize_env
  use IGLOO_Mod_Exchange, only: exchange_t, material_id, exchange_rate, gas_rate, track_processed, &
                                absorption_packet, add_absorption_rate, add_birth_rate, clear_import_rates
  implicit none
  type(obj_IGLOO) :: igloo
  character(len=16) :: mode
  integer :: nFail = 0

  call mpi_init_env()
  call get_command_argument(1, mode)
  select case (trim(mode))
  case ('track');  call run_track()
  case ('absorb'); call run_absorb()
  case ('births'); call run_births()
  case default
    write(*,'(a)') ' usage: test_exchange track|absorb|births'
    nFail = 1
  end select
  call mpi_finalize_env()
  if (nFail > 0) then
    write(*,'(/a,i0,a)') ' [FAIL] ', nFail, ' gate(s)'
    error stop 1
  endif
  write(*,'(/a)') ' [PASS] every gate holds'

contains

  subroutine gate(name, ok, detail)
    character(len=*), intent(in) :: name, detail
    logical,          intent(in) :: ok
    if (ok) then
      write(*,'(a)') ' [PASS] '//name//': '//detail
    else
      write(*,'(a)') ' [FAIL] '//name//': '//detail
      nFail = nFail + 1
    endif
  end subroutine gate


  subroutine run_track()
    character(len=:), allocatable :: why
    character(len=160) :: line
    type(exchange_t) :: x
    real(R8) :: errW, errU, sumW, sumNet, Tlo, Thi, q(4), fMax, eSum
    integer :: m, i, j, k, nCell

    call igloo%setup_static()
    m = material_id(igloo, 'A')
    call track_processed(igloo, m, why)
    call gate('X0', len(why) == 0, 'track_processed accepts A ('//why//')')
    call igloo%reset_state()
    call igloo%solve()
    errW = 0._R8; errU = 0._R8; sumW = 0._R8; sumNet = 0._R8; Tlo = huge(1._R8); Thi = -huge(1._R8); nCell = 0
    do k = 1, igloo%geoblock(1)%Nz; do j = 1, igloo%geoblock(1)%Ny; do i = 1, igloo%geoblock(1)%Nx
      x = exchange_rate(igloo, m, 1, i, j, k)
      sumW = sumW + x%processed; sumNet = sumNet + x%mass
      errW = max(errW, abs(x%processed - x%mass))
      if (x%tracked) then
        nCell = nCell + 1
        errU = max(errU, norm2(x%u - [10._R8, 0._R8, 0._R8]))
        Tlo = min(Tlo, x%T); Thi = max(Thi, x%T)
      endif
    enddo; enddo; enddo
    write(line,'(a,es10.3,a,es10.3,a)') 'worst |W - net| ', errW, ' kg/s of a total ', sumNet, ' kg/s'
    call gate('X1', sumNet > 0._R8 .and. errW <= 1.e-12_R8*sumNet, trim(line))
    write(line,'(a,i0,a,es10.3,a)') 'in ', nCell, ' cells, worst |u - u_g| ', errU, ' m/s'
    call gate('X2', nCell > 1000 .and. errU <= 1.e-9_R8, trim(line))
    write(line,'(a,f8.2,a,f8.2,a)') 'T in [', Tlo, ', ', Thi, '] K'
    call gate('X3', Tlo >= 600._R8 - 1.e-6_R8 .and. Thi < 3000._R8, trim(line))
    call track_processed(igloo, 99, why)
    call gate('X4', len(why) > 0, 'material 99 refused: '//why)
    fMax = 0._R8; eSum = 0._R8
    do k = 1, igloo%geoblock(1)%Nz; do j = 1, igloo%geoblock(1)%Ny; do i = 1, igloo%geoblock(1)%Nx
      q = gas_rate(igloo, 1, i, j, k)
      x = exchange_rate(igloo, m, 1, i, j, k)
      fMax = max(fMax, norm2(q(1:3) - x%mass*[10._R8, 0._R8, 0._R8])); eSum = max(eSum, norm2(q(1:3)))
    enddo; enddo; enddo
    write(line,'(a,es10.3,a,es10.3,a)') 'worst |F - W u_g| ', fMax, ' N of a largest |F| ', eSum, ' N'
    call gate('X5', eSum > 0._R8 .and. fMax <= 1.e-12_R8*eSum, trim(line))
  end subroutine run_track


  subroutine run_absorb()
    real(R8), parameter :: MC = 2.e-6_R8, UG(3) = [10._R8, 0._R8, 0._R8], TG = 600._R8
    character(len=:), allocatable :: why
    character(len=160) :: line
    type(exchange_t) :: x
    real(R8) :: p(5), offered, absorbed(3), flow(4), errU, errE
    integer :: sweep, i, j, k
    logical :: ok

    call igloo%setup_static()
    p = absorption_packet(igloo, 1, MC, UG, TG)
    call offer_band(ok)
    call gate('A0', ok, 'the offer accepted in cells 21-40')
    do sweep = 1, 4
      if (sweep == 4) call clear_import_rates(igloo)
      call igloo%reset_state()
      call igloo%solve()
      flow(sweep) = sum(igloo%material(1)%group(1)%particle(1:igloo%material(1)%group(1)%nparticles)%mdot)
      if (sweep > 3) cycle
      absorbed(sweep) = 0._R8; offered = 0._R8; errU = 0._R8; errE = 0._R8
      do k = 1, igloo%geoblock(1)%Nz; do j = 1, igloo%geoblock(1)%Ny; do i = 21, 40
        x = exchange_rate(igloo, 1, 1, i, j, k)
        absorbed(sweep) = absorbed(sweep) + x%absorbed(1)
        if (allocated(igloo%source(1)%heldLast)) then
          if (igloo%source(1)%heldLast(1, i, j, k) > 0._R8) offered = offered + MC
        endif
        if (x%absorbed(1) > 0._R8) then
          errU = max(errU, norm2(x%absorbed(2:4)/x%absorbed(1) - p(2:4)/p(1))/norm2(UG))
          errE = max(errE, abs(x%absorbed(5)/x%absorbed(1) - p(5)/p(1))/(p(5)/p(1)))
        endif
      enddo; enddo; enddo
      if (sweep == 1) then
        write(line,'(a,es10.3,a)') 'sweep 1 absorbed ', absorbed(1), ' kg/s'
        call gate('A1', absorbed(1) == 0._R8, trim(line))
      elseif (sweep == 2) then
        write(line,'(a,es12.5,a,es12.5,a)') 'sweep 2 absorbed ', absorbed(2), ' kg/s of ', offered, &
          ' kg/s offered where drops were'
        call gate('A2a', offered > 0._R8 .and. abs(absorbed(2)/offered - 1._R8) <= 0.04_R8, trim(line))
        write(line,'(a,es10.3,a,es10.3)') 'per kg: worst velocity ', errU, ', worst energy ', errE
        call gate('A3', errU <= 1.e-10_R8 .and. errE <= 1.e-9_R8, trim(line))
        write(line,'(a,es12.5,a,es12.5,a)') 'exit flow gain ', flow(2) - flow(1), ' kg/s, absorbed ', absorbed(2), ' kg/s'
        call gate('A4', abs(flow(2) - flow(1) - absorbed(2)) <= 1.e-9_R8*absorbed(2), trim(line))
      else
        write(line,'(a,es12.5,a,es12.5,a)') 'sweep 3 absorbed ', absorbed(3), ' kg/s of ', offered, &
          ' kg/s offered where drops were'
        call gate('A2b', offered > 0._R8 .and. abs(absorbed(3)/offered - 1._R8) <= 0.005_R8, trim(line))
      endif
    enddo
    write(line,'(a,es23.15,a,es23.15)') 'exit flow ', flow(4), ' vs sweep 1 ', flow(1)
    call gate('A5', flow(4) == flow(1), trim(line))
    call add_absorption_rate(igloo, 1, 1, 61, 1, 1, p, why)
    call gate('A6a', len(why) > 0, 'cell 61 of a 60-cell block refused: '//why)
    call add_absorption_rate(igloo, 1, 1, 30, 1, 1, [0._R8, 0._R8, 0._R8, 0._R8, 0._R8], why)
    call gate('A6b', len(why) > 0, 'a zero mass refused: '//why)

  end subroutine run_absorb


  !> The offer p in cells 21-40 of block 1 (every j, k).
  subroutine offer_band(ok)
    logical, intent(out) :: ok
    character(len=:), allocatable :: why
    real(R8) :: p(5)
    integer :: i, j, k
    p = absorption_packet(igloo, 1, 2.e-6_R8, [10._R8, 0._R8, 0._R8], 600._R8)
    ok = .true.
    do k = 1, igloo%geoblock(1)%Nz; do j = 1, igloo%geoblock(1)%Ny; do i = 21, 40
      call add_absorption_rate(igloo, 1, 1, i, j, k, p, why)
      ok = ok .and. len(why) == 0
    enddo; enddo; enddo
  end subroutine offer_band


  subroutine run_births()
    real(R8), parameter :: MBIRTH(2) = [1.e-5_R8, 3.e-5_R8], DB = 20.e-6_R8
    integer,  parameter :: CELL(3,2) = reshape([10, 3, 3, 30, 2, 4], [3, 2])
    character(len=:), allocatable :: why
    character(len=160) :: line
    real(R8) :: errM, errX
    integer :: iA, iB, e, sweep, n
    logical :: ok

    call igloo%setup_static()
    iA = material_id(igloo, 'A'); iB = material_id(igloo, 'B')
    call add_birth_rate(igloo, iA, 1, 10, 3, 3, MBIRTH(1), [10._R8, 0._R8, 0._R8], 600._R8, DB, why)
    ok = len(why) > 0
    write(line,'(a)') 'A refused: '//why
    call add_birth_rate(igloo, iB, 1, 10, 3, 3, MBIRTH(1), [10._R8, 0._R8, 0._R8], 600._R8, 0._R8, why)
    ok = ok .and. len(why) > 0
    call gate('B1', ok, trim(line)//'; d = 0 refused: '//why)
    ok = .true.
    do e = 1, 2
      call add_birth_rate(igloo, iB, 1, CELL(1,e), CELL(2,e), CELL(3,e), MBIRTH(e), [5._R8, 0._R8, 0._R8], 600._R8, DB, why)
      ok = ok .and. len(why) == 0
    enddo
    call gate('B2a', ok, 'two births of B accepted')
    do sweep = 1, 3
      if (sweep == 3) call clear_import_rates(igloo)
      call igloo%reset_state()
      call igloo%solve()
      associate (gr => igloo%material(iB)%group(1))
      n = gr%nparticles
      if (sweep == 3) then
        write(line,'(a,i0)') 'parcels of B after the clear: ', n
        call gate('B5', n == 0, trim(line))
        cycle
      endif
      errM = 0._R8; errX = 0._R8
      do e = 1, min(n, 2)
        errM = max(errM, abs(gr%particle(e)%mdot - MBIRTH(e))/MBIRTH(e))
        errX = max(errX, abs(gr%particle(e)%stateVar(1) - 0.15_R8))
      enddo
      if (sweep == 1) then
        write(line,'(a,i0,a,i0,a,i0)') 'parcels of B ', n, ', IDs ', gr%particle(1)%ID, ' and ', gr%particle(min(n,2))%ID
        call gate('B2b', n == 2 .and. gr%particle(1)%ID == 1 .and. gr%particle(min(n,2))%ID == 2, trim(line))
        write(line,'(a,es10.3,a,es10.3,a)') 'worst mass-flow error ', errM, ', worst distance to the outlet ', errX, ' m'
        call gate('B3', errM <= 1.e-12_R8 .and. errX <= 1.e-9_R8, trim(line))
      else
        write(line,'(a,i0,a,es10.3)') 'second sweep: parcels of B ', n, ', worst mass-flow error ', errM
        call gate('B4', n == 2 .and. errM <= 1.e-12_R8 .and. errX <= 1.e-9_R8, trim(line))
      endif
      end associate
    enddo
  end subroutine run_births

end program test_exchange
