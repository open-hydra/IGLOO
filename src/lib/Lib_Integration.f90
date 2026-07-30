module Lib_Integration
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private
  public :: integrate


contains

  subroutine integrate(part,geoblock,gasblock,srcblock,eulblock,     &
                        hTab,cpTab,rhoTab,mupTab,sigTab,psatTab, child,addChild)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    use IGLOO_data_phases, only: obj_child
    use IGLOO_variables,  only: unitTraj,unitExit,unitScat,iprint,dtprint, &
                                nb,nspecies,toll,ord2,mesh2D,threshold,    &
                                eulerSwitch,sourceSwitch, phaseChange,     &
                                bodyAccel, srcBodyForce, axisym,           &
                                trajOn, scatOn, dNscat, sixOverPi
    use IGLOO_particles,  only: obj_particle, eps
    use IGLOO_bcBox,      only: axisymFold
    use IGLOO_data_block, only: obj_block, obj_flowblock, obj_eulerblock, obj_sourceblock
    use IGLOO_Lib_Properties
    use Lib_Equations
    use Lib_RHS, only: packAuxVars, nauxvar, packAuxState, unpackAuxState,     &
                       packEventVar, unpackEventVar,                           &
                       rhsStandard, rhsEvaporation, rhsBreakupOnly, rhsEvapBreakup, &
                       rhsAlCombustion,                                             &
                       ind_d, ind_rho, ind_m0, ind_n0, ind_sig, ind_mup, ind_m,     &
                       ind_sxi, ind_sb1, ind_sb2,                                   &
                       ind_evd, ind_evn, ind_evm0, ind_ev1, ind_ev2,                &
                       nauxstate, neventvar, nbrkst,                                &
                       mod_brkSelect, mod_propFlags, mod_model,                     &
                       mod_bp, mod_bpMethod, mod_bpScale
    use oslo,              only: Run_ODESolver
    use IGLOO_Lib_Breakup, only: nchild
    implicit none
    class(obj_particle),   intent(inout) :: part
    type(obj_block),       intent(in)    :: geoblock(nb)
    type(obj_flowblock),   intent(in)    :: gasblock(nb)
    type(obj_sourceblock), intent(inout) :: srcblock(nb)
    type(obj_eulerblock),  intent(inout) :: eulblock(nb)
    real(R8),              intent(in)    :: hTab(:),cpTab(:),rhoTab(:),     &
                                            mupTab(:),sigTab(:),psatTab(:)
    type(obj_child),       intent(out), optional :: child
    logical,               intent(out), optional :: addChild
    !> local variables
    real(R8), allocatable :: gas(:,:), gasVert(:,:), gasState(:)
    real(R8) :: vert(3,8), t1, t2, tStart, deltat, tlimit, tprint
    real(R8) :: Ein, Eout, Pin(3), Pout(3), massIn, massOut, vol
    real(R8) :: entryPos(3)   ! cell-entry position, for the closed-form body-force work (models 1,3)
    real(R8) :: wAcc          ! scatter-cloud npdot-weight accumulator (host-associated into solout)
    integer  :: neq, nsp, ng, b, i, j, k, ngVert, iu, iw, it, iter, nCross, nE, nStall
    real(R8) :: posPrev(3)    ! last outer-iter position, for the zero-progress guard
    logical  :: doLoop, IamOut, newGas, eventType, eventFlag, exitLoop, startedOut, burnedOut
    logical  :: consumed   !> droplet ended INSIDE the domain with mass still on it
    logical  :: atGasBoundary, wasBoundary   ! ord2 C2: geo consulted only at gas-boundary cells
    real(R8) :: geoHexNorms(3,2,6), geoHexCentroids(3,2,6)
    real(R8) :: gasHexNorms(3,2,6), gasHexCentroids(3,2,6)
    logical  :: geoHexDegen(2,6), gasHexDegen(2,6)
    integer,  parameter :: maxIter=500000, nMaxCell=10, nMaxStall=100
    !> accessed by ODEsystem, rhs1-4, solout via host association
    real(R8) :: y(part%neq), oldLocal(part%neq)
    real(R8) :: auxLocal(nauxvar), childState(nchild)
    logical  :: addChildLocal   !> breakupEvent needs a present target even when the caller omits addChild
    real(R8), allocatable :: stateLocal(:), oldStLocal(:)
    real(R8), allocatable :: eventLocal(:), oldEvLocal(:)
    real(R8) :: timeLocal, din, dout, deltaS(3), dir(3), taup
    integer  :: err, nDL, innerIter, jSlot
    integer,  parameter :: maxInnerIter=10, nStep=10
    real(R8), parameter :: safety=0.99, dtMin=1.e-14_R8, tauFactor=100._R8
    real(R8), parameter :: escapeDist=1.e-6_R8   ! min displacement per sliver-escape hop (startedOut)
    !> Burnout tolerance: absolute threshold on droplet MASS [kg]. A consuming droplet
    !  (evaporation/combustion) is declared fully consumed at m <= mBurnTol, stopping on a
    !  still-good state and staying out of the m->0 region where d=(6m/(pi*rho))^(1/3) has
    !  unbounded slope and the state vector goes NaN. Mass is the ODE state y(8), so the
    !  test is exact. 1e-15 kg is d ~ 1.24 um for water.
    !  Two limits pull opposite ways if this value is changed: particles injected below it
    !  are killed at injection (so sub-micron injection needs it lowered), while a case whose
    !  solver cannot reach it ends through the [WARNING] path instead (outputs stay correct).
    !  A threshold relative to the injected mass would satisfy both.
    real(R8), parameter :: mBurnTol=1.e-15_R8

    ! write(*,*) "      Processing particle n.", part%ID !> DEBUG PRINTING
    tlimit = huge(1._R8) !> steady-state by default (temporary)
    neq = part%neq       !> save local copy of neq
    atGasBoundary = .true.   !> conservative default: geo consulted until proven interior
    nsp = 1              !> nsp = nspecies(part%Iinj(1)) !> for multi-species gas
    iu  = nsp + 1
    ng  = iu + 7 ! +1 if eddy viscosity (mut) included

    !> ARRAYS ALLOCATION
    addChildLocal = .false.
    if (present(addChild)) addChild = .false.
    allocate(gasState(ng))
    if     (ord2.and.mesh2D) then; allocate(gas(ng,4)); allocate(gasVert(3,4)); ngVert=4
    elseif (     ord2      ) then; allocate(gas(ng,8)); allocate(gasVert(3,8)); ngVert=8
    else;                          allocate(gas(ng,1)); allocate(gasVert(3,1)); ngVert=1
    endif
    !> Allocate with at-least-size-1 `nauxstate`/`neventvar` 
    allocate(stateLocal(max(1,nauxstate))); allocate(oldStLocal(max(1,nauxstate)))
    allocate(eventLocal(max(1,neventvar))); allocate(oldEvLocal(max(1,neventvar)))

    !> PARCEL INITIALIZATION
    if (part%time==0._R8) then
      part%Ncell = 0
      part%gone  = .false.
      part%wasin = .false.
      part%angle = 0._R8
      part%i     = part%iInj
      part%stateVar(1:3) = part%pInj
      call part%updateCell(geoblock)

      if (part%gone) then
        write(*,'(A,I4,A)')                                                         &
                '[WARNING] Particle ',part%ID,' has been initialized out of domain!'
        return
      endif

      b = part%i(1)
      if (ord2) then; call part%findDualCell(gasblock(b)); part%igasOld = part%igas; part%xi0 = 0.5_R8
                      call setAtGasBoundary()
      else;                part%igas = part%i(2:4); endif
      call gasblock(b)%gasProperties(gas,part%igas)
      if (all(part%iInj==[0,0,0,0])) then
        part%stateVar(4:6) = part%vInj   ! DB hand-off (pin_particles can't write stateVar pre-allocation)
        if (any(part%stateVar(4:6)>=threshold)) part%stateVar(4:6) = gas(iu:iu+2,1)
        if (part%tp==0._R8) part%tp = gas(iu+3,1)
      else
        call part%initializePart(geoblock,gasblock)
      endif
      call geoblock(b)%getVertices(part%i(2:4),vert, norms=geoHexNorms, centroids=geoHexCentroids, degen=geoHexDegen)   ! geo cell (NOT the gas dual index igas)
      if (ord2) call gasblock(b)%getVertices(part%igas,gasVert,mesh2D, gasHexNorms, gasHexCentroids, gasHexDegen)
      if (part%varCp) then; part%stateVar(7) = lookupTab(hTab,part%tp); else; part%stateVar(7) = part%tp; endif
      call part%computeMass(rhoTab); part%m0 = part%m; part%npdot = part%mdot/part%m; part%npold = part%npdot

      select case(part%model)
      case(2,5); part%stateVar(8) = part%m;   part%nOde = 8
      case(3); part%stateVar(8) = part%npdot; part%nOde = 8
      case(4); part%stateVar(8) = part%m;
               part%stateVar(9) = part%npdot; part%nOde = 9
      end select
      if (eulerSwitch)    part%stateVar(part%nOde+1:part%neq) = 0._R8
      if (part%bodyAccum) part%stateVar(part%neq-1 :part%neq) = 0._R8

      if (trajOn) write(unit=unitTraj,fmt='(7F12.6,2E13.6E2,I8)') part%stateVar(1:6), part%tp, part%d, part%m, part%ID
    endif
    !> Pack auxiliary variables once (constant throughout integration)
    call packAuxVars(part, nauxvar, auxLocal)

    !> TRAJECTORY INTEGRATION LOOP
    iter = 0; nCross=0; tprint = part%time + dtprint
    nStall = 0; posPrev = part%stateVar(1:3) - 1._R8
    !> Scatter accumulator, seeded with a per-ID phase offset (golden-ratio low-discrepancy)
    !  so neighbouring streams don't drop points in lockstep => no banding. Deterministic.
    wAcc = dNscat * mod(real(part%ID,R8)*0.6180339887498949_R8, 1._R8)
    do while (part%time<tlimit.and.iter<maxIter)
      iter = iter+1
      call part%initializeCell(eulerSwitch)
      if (sourceSwitch) call part%computeSource(massIn,Pin,Ein)
      entryPos = part%stateVar(1:3)   ! cell-entry position (for the closed-form body-force work)
      i = part%i(2); j = part%i(3); k = part%i(4);
      if (ord2) then; deltat = part%computeDeltat(gasVert)
                      if (atGasBoundary) deltat = min(deltat, part%computeDeltat(vert))
      else;           deltat = part%computeDeltat(vert); endif
      !> Cap deltat by particle drag relaxation taup = rho_p*d^2/(18*mu_g).
      if (part%d > 0._R8 .and. any(gas(6,:) > toll)) taup   = part%rho * part%d**2 / (18._R8 * maxval(gas(6,:)))
      if (.not. ieee_is_finite(deltat)) deltat = tauFactor * taup

      tStart = part%time
      t2     = tStart + 20.0_R8 * deltat
      doLoop = .true.
      consumed = .false.
      innerIter = 0
      do while (doLoop)
        innerIter = innerIter + 1
        if (innerIter > maxInnerIter) then
          write(*,'(a,i0,a,i0,a)')                          &
            '[WARNING] Inner loop > ',maxInnerIter,' iter, part=',part%ID,' ==> marking gone'
          doLoop    = .false.
          part%gone = .true.
          exit
        endif
        t1     = part%time
        call ODEsystem()
      enddo
            
      if (present(addChild)) then
        if (addChild) then
          child%ipos = part%i
          child%igas = part%igas
          child%pos  = part%stateVar(1:3)
          child%temp = part%tp
          child%time = t1
        endif
      endif
      
      if (sourceSwitch) then
        call part%computeSource(massOut,Pout,Eout)
        !> Two-phase conservation on burnout: the cell source is the net parcel flux
        !  (in - out), so a droplet consumed inside the cell has zero outgoing flux and its
        !  whole remnant transfers to the gas. Consumption only (models 2/5) -- a particle
        !  that exits the domain keeps its outgoing flux.
        if (consumed) then
          massOut = 0._R8
          Pout    = 0._R8
          Eout    = 0._R8
        endif
        !> Body-force source-reaction correction: (Pin-Pout) and (Ein-Eout) absorbed the body force
        if (srcBodyForce) then
          select case(part%model)
          case(2,5)      ! J,W at neq-1,neq;
            Pin = Pin + (part%npdot * part%stateVar(part%neq-1)) * bodyAccel
            Ein = Ein + (part%npdot * part%stateVar(part%neq))
          case(4)        ! euler-on: J at neq-2 (neq-1 is F(16)=npdot)
            jSlot = part%neq - 1; if (eulerSwitch) jSlot = part%neq - 2
            Pin = Pin + part%stateVar(jSlot) * bodyAccel
            Ein = Ein + part%stateVar(part%neq)
          case default   ! models 1,3: closed form, no accumulator
            Pin = Pin + (part%mdot * part%Tstay) * bodyAccel
            Ein = Ein + (part%mdot * dot_product(bodyAccel, part%stateVar(1:3) - entryPos))
          end select
        endif
        if (ord2) then
          call computeSrcField(part, srcblock(b), part%igas(1), part%igas(2), part%igas(3), &
                               massIn, massOut, Pin, Pout, Ein, Eout)
        else
          call computeSrcField(part, srcblock(b), i, j, k, massIn, massOut, Pin, Pout, Ein, Eout)
        endif
      endif

      !> Source forcing terms. ord2=true → accumulate on gasblock (staggered)
      !  cells using part%igas; ord2=false → keep geoblock accumulation as today.
      if (eulerSwitch) then
        if (ord2) then
          call computeEulField(part, eulblock(b), part%igas(1), part%igas(2), part%igas(3))
        else
          call computeEulField(part, eulblock(b), i, j, k)
        endif
      endif

      if (ord2) then
        !> Gas crossing: advance igas on the dual mesh (gas-primary tracking).
        if (newGas) then
          call part%advanceGasCell(gasblock(b))
          call gasblock(b)%getVertices(part%igas,gasVert,mesh2D, gasHexNorms, gasHexCentroids, gasHexDegen)
          if (.not.mesh2D) then
            !> 3D: vert numbering from getVertices
            part%xi0(1) = (part%stateVar(1) - gasVert(1,1)) / (gasVert(1,5) - gasVert(1,1) + toll)
            part%xi0(2) = (part%stateVar(2) - gasVert(2,1)) / (gasVert(2,4) - gasVert(2,1) + toll)
            part%xi0(3) = (part%stateVar(3) - gasVert(3,1)) / (gasVert(3,2) - gasVert(3,1) + toll)
            part%xi0 = max(0._R8, min(1._R8, part%xi0))
          endif
          call gasblock(b)%gasProperties(gas,part%igas)
        endif
        !> Refresh boundary flag
        wasBoundary = atGasBoundary
        call setAtGasBoundary()
        if (atGasBoundary .and. .not.wasBoundary) then
          if (part%findParticle(geoblock)) then
            b = part%i(1)
            call geoblock(b)%getVertices(part%i(2:4),vert, norms=geoHexNorms, centroids=geoHexCentroids, degen=geoHexDegen)
          else
            call part%updateCell(geoblock); call handleGeoEvent()
          endif
        endif
        !> Geo domain-exit / BC: IamOut is only ever true while already at a boundary
        !  cell (solout gates it). bcDef may reflect -> handleGeoEvent re-resolves igas.
        if (IamOut) then
          call part%updateCell(geoblock); call handleGeoEvent()
        endif
        !> Stuck detection on gas-cell residency (geo index is stale in the interior).
        if (newGas .or. IamOut) then; nCross = nCross + 1; part%Ncell = 0
        else;                         part%Ncell = part%Ncell + 1; endif
      else
        !> ord1: geo-driven (unchanged).
        if (IamOut) then
          call part%updateCell(geoblock); b = part%i(1)
          call geoblock(b)%getVertices(part%i(2:4),vert, norms=geoHexNorms, centroids=geoHexCentroids, degen=geoHexDegen)
        endif
        if (newGas) then
          part%igas = part%i(2:4)
          call gasblock(b)%gasProperties(gas,part%igas)
        endif
        if (all(part%i==part%iold)) then; part%Ncell = part%Ncell+1
        elseif (IamOut) then;   nCross = nCross + 1;   part%Ncell = 0; endif
      endif

      !> Axisymmetric wedge: fold the particle back into the sector if it crossed faces 5/6.
      if (axisym) call axisymFold(part%stateVar)

      if (part%Ncell>nMaxCell) then
        write(*,'(A,I4,A,4I4,A)') '       ==> Particle ',part%ID,' stuck in cell:',part%iold
        part%gone=.true.
      endif

      !> Zero-progress guard: a sliver-pinned particle "crosses" every iter (resetting Ncell) with
      !  no net motion -> invisible to nMaxCell; catch it on displacement instead.
      if (norm2(part%stateVar(1:3)-posPrev) > eps) then; nStall = 0
      else;                                              nStall = nStall + 1; endif
      posPrev = part%stateVar(1:3)
      if (nStall > nMaxStall) then
        write(*,'(A,I4,A)') '       ==> Particle ',part%ID,' no net progress ==> marking gone'
        part%gone = .true.
      endif

      if ((dtprint>0._R8.and.part%time>=tprint).or.(mod(nCross,iprint)==0.and.part%Ncell==0)) then
        if (trajOn) write(unit=unitTraj,fmt='(7F12.6,2E13.6E2,I8)') part%stateVar(1:6), part%tp, part%d, part%m, part%ID
        tprint = tprint + dtprint
      endif
      if (part%gone) then
        part%time = 0._R8
        write(unit=unitExit,fmt='(6F12.6,2E13.6E2,I8)')                                      &
              part%stateVar(1:3), part%tp, norm2(part%stateVar(4:6)), part%angle, part%mdot, part%Af, part%ID
        return
      endif

    enddo

    !> Outer-loop maxIter silent-exit. If we ran out of iterations without the particle being marked gone,
    !  flag and write to unitExit so callers see a recognizable terminal state.
    if (iter >= maxIter .and. .not. part%gone) then
      write(*,'(A,I4,A,I0,A)') '       ==> Particle ',part%ID,' hit outer maxIter (',maxIter,'); flagging gone'
      part%gone = .true.
      part%time = 0._R8
      write(unit=unitExit,fmt='(6F12.6,2E13.6E2,I8)')                                          &
            part%stateVar(1:3), part%tp, norm2(part%stateVar(4:6)), part%angle, part%mdot, part%Af, part%ID
    endif

  contains

    !============================================================================================!
    !  ODEsystem — internal procedure, accesses integrate's locals via host association          !
    !============================================================================================!

    subroutine ODEsystem()
      implicit none
      real(R8) :: dtNew, dtSafe

      deltat = safety*deltat
      newGas = .false.; IamOut = .false.; exitLoop = .false.; eventFlag = .false.; startedOut = .false.
      burnedOut = .false.
      eventType = part%brkupEvent
      addChildLocal = .false.; childState = 0._R8
      if (present(addChild)) addChild = .false.
      y = part%oldState
      timeLocal = part%time
      oldLocal  = y
      if (allocated(stateLocal)) then
        call packAuxState(part, nauxstate, stateLocal)
        !> Segment-start init for the RT timer (told/tc in brkupState): without it a
        !  first-solout crossing rolls the timer back to a stale prior segment.
        oldStLocal = stateLocal
      endif
      if (allocated(eventLocal)) then
        call packEventVar(part, neventvar, eventLocal)
        !> Segment-start init for the oscillator (y,yDot): without it a first-solout
        !  crossing rolls deformation back and the droplet never breaks.
        oldEvLocal = eventLocal
      endif
      if (.not.ord2) gasState = gas(:,1)
      err = 0
      select case(part%model)
      case(1); nDL = 8;  call Run_ODESolver(neq, t1, t2, y, rhs1, err, deltat, solout)
      case(2); nDL = 9;  call Run_ODESolver(neq, t1, t2, y, rhs2, err, deltat, solout)
      case(3); nDL = 9;  call Run_ODESolver(neq, t1, t2, y, rhs3, err, deltat, solout)
      case(4); nDL = 10; call Run_ODESolver(neq, t1, t2, y, rhs4, err, deltat, solout)
      case(5); nDL = 9;  call Run_ODESolver(neq, t1, t2, y, rhs5, err, deltat, solout)
      end select
      if (err < 0) then
        write(*,'(a,i0,a,i0,a,es12.4,a)')                                   &
          '[IGLOO] Run_ODESolver err=',err,' part=',part%ID,                   &
          ' time=',part%time,' => marking gone'
        doLoop    = .false.
        part%gone = .true.
        return
      endif
      !> Expected end of life: the droplet evaporated/burned away on a still-good state.
      !  Physical, not an error => silent. Falls through to the normal finalize so
      !  updatePart() syncs part%d/%m/%tp; returning early would leave them stale.
      if (burnedOut) then
        part%gone = .true.
        consumed  = .true.     !> remnant transfers to the gas; see the source block
      endif

      !> Unexpected non-finite state: with the burnout guard above, a clean
      !  evaporation/combustion run never reaches the degenerate m->0 region, so a NaN here is a
      !  REAL anomaly (solver breakdown, bad gas data, ...) and must be reported. Never propagate
      !  it: computeSource would smear NaN through the whole Eulerian source (mollifier stencil)
      !  and write a NaN outloc record. Fall back to the last good accepted step.
      if (any(y/=y)) then
        write(*,'(A,I4,A)') '[WARNING] Particle ',part%ID,                                &
                            ' non-finite state ==> reverting to last good step, marking gone'
        !> Report how much droplet was left: if the absolute mass at failure is only somewhat
        !  above mBurnTol, the solver simply could not integrate down to the floor (raise it, or
        !  switch to the relative form documented at the parameter); a much larger remnant means
        !  a genuine anomaly (solver breakdown, bad gas data) worth investigating.
        if (mod_model==2 .or. mod_model==4 .or. mod_model==5)                              &
          write(*,'(A,ES10.3,A,ES10.3,A)') '          last good m/m0 =',                   &
                oldLocal(8)/part%m0, '  (burnout fires at m =', mBurnTol, ' kg)'
        part%stateVar = oldLocal
        part%time     = timeLocal
        part%oldState = oldLocal
        part%gone     = .true.
        doLoop        = .false.
        !> A consumption model losing a droplet in-domain still has to hand its remnant to the
        !  gas, exactly as a clean burnout does (see the source block); other models keep the
        !  pre-existing behaviour, since their mass does not transfer by phase change.
        if (mod_model==2 .or. mod_model==5) consumed = .true.
        !> sync d/m/tp from the REVERTED (good) state; the normal finalize is skipped here
        !  because stateLocal/eventLocal may themselves be contaminated.
        call part%updatePart(rhoTab,hTab,eulerSwitch)
        return
      endif
      part%stateVar = y
      part%time     = timeLocal
      part%oldState = oldLocal
      !> Full [t1,t2] consumed with no interrupt (solver reached XEND): finalize the segment;
      !  re-invoking on the null interval [t2,t2] hands SDIRK4 H=0 -> IDID=-1 -> spurious kill.
      !  Outer loop re-budgets deltat from the current (possibly decelerated) velocity.
      if (.not.(IamOut .or. newGas .or. eventFlag)) exitLoop = .true.
      !> startedOut escape hop: finalize at the moved position (no refinement — it would re-integrate
      !  from the sliver start and discard the motion); updateCell then relocates from there.
      if (startedOut) exitLoop = .true.
      if (exitLoop) then
        doLoop = .false.
        if (eventType) then
          if (allocated(eventLocal)) then
            part%npold = oldEvLocal(ind_evn)
            call unpackEventVar(eventLocal, part, neventvar)
            !> An event resized the droplet, but models 1/2 freeze per-droplet d and m in
            !  auxLocal at injection: resync, else the next segment recomputes d from the
            !  stale aux and reverts the breakup.
            if (ind_d > 0) auxLocal(ind_d) = part%d
            if (ind_m > 0 .and. mod_model /= 3) then
              part%m          = part%rho * part%d**3 / sixOverPi
              auxLocal(ind_m) = part%m
            endif
            !> Model 3 keeps the per-droplet count in the ODE state y(8), with d derived from
            !  it and the parcel mass-flow. Take the event (d, npdot) as truth: push npdot into
            !  y(8) and re-derive mdot to match, which conserves parcel mass for both sub-events
            !  (RT split leaves mdot alone; KH shed drops it by (dParent/dp)^3).
            !  Guarded on eventFlag: on a no-event finalize part%npdot is stale, and writing it
            !  into y(8) would freeze the continuous KH stripping.
            if (eventFlag .and. mod_model == 3 .and. part%d > 0._R8 .and. ind_m > 0) then
              part%stateVar(8) = part%npdot
              part%oldState(8) = part%npdot
              part%mdot        = part%npdot * part%rho * part%d**3 / sixOverPi
              auxLocal(ind_m)  = part%mdot
            endif
          endif
          if (present(addChild)) then
            addChild = addChildLocal
            if (addChild) then
              child%vel   = childState(1:3)
              child%diam  = childState(4)
              child%npdot = childState(5)
            endif
          endif
        endif
        if (allocated(stateLocal)) call unpackAuxState(stateLocal, part, nauxstate)
        if (eulerSwitch) then
          !> Euler moments end at nE; with a body force only W (slot neq) trails them — drop it.
          nE = neq; if (part%bodyAccum) nE = neq - 1
          part%deltaL           = y(nDL)
          part%intE(1:(nE-nDL)) = y(nDL+1:nE)
        endif
        part%time  = t1
        part%Tstay = t1 - tStart
        call part%updatePart(rhoTab,hTab,eulerSwitch)
      else
        if (allocated(eventLocal)) call unpackEventVar(oldEvLocal, part, neventvar)
        if (allocated(stateLocal)) call unpackAuxState(oldStLocal, part, nauxstate)
        deltaS = part%stateVar(1:3) - part%oldState(1:3)
        dout = norm2(deltaS)
        if (dout>eps) then
          dir  = deltaS/dout
          if     (IamOut) then; din = part%computeDs(vert,dir)
          elseif (newGas) then; din = part%computeDs(gasVert,dir)
          else;                 din = dout/nStep; endif
          !> computeDs is only reached once a crossing is already flagged (IamOut/newGas). din<=0
          !  means its Möller-Trumbore refine failed to pin the just-crossed face (isPointInsideCell
          !  vs triangle-split tolerance at sub-1e-5 proximity), NOT that there is no face. Keep
          !  converging at the interior rate so deltat & dout shrink together (~nStep cheap steps to
          !  re-cross) instead of collapsing hmax (-> NMAX) or overshooting (-> maxInnerIter).
          if (din > 0._R8) then; dtNew = max(deltat*min(1._R8/nStep,din/dout),dtMin)
          else;                  dtNew = deltat/nStep; endif
        else
          dtNew = deltat/nStep
        endif
        dtSafe = min(deltat,eps/nStep/(norm2(part%oldState(4:6))+toll))/nStep
        if (ord2) then;      deltat = part%computeDeltat(gasVert)
          if (atGasBoundary) deltat = min(deltat, part%computeDeltat(vert))
        else;                deltat = part%computeDeltat(vert); endif
        deltat = max(min(deltat,dtNew),dtMin,dtSafe)
        !> Trap a non-finite/zero deltat after the update
        if (.not. (deltat >= dtMin .and. deltat < huge(deltat))) then
          print*,'DBG[ODEsystem] bad deltat after update: part=',part%ID, &
                 ' time=',part%time,' deltat=',deltat,' dout=',dout,' din=',din,' dtNew=',dtNew
          print*,'  stateVar(1:nDL)=',part%stateVar(1:nDL)
          error stop '[ERROR] deltat update produced bad value'
        endif
      endif

    end subroutine ODEsystem

    subroutine rhs1(neq, time, Z, F)
      integer, intent(in) :: neq; real(R8), intent(in) :: time, Z(neq); real(R8), intent(out) :: F(neq)
      call rhsStandard(neq,time,Z,F, auxLocal,nauxvar, stateLocal,nauxstate, hTab,rhoTab, gas,gasVert,gasState, ng,ngVert)
    end subroutine rhs1
    subroutine rhs2(neq, time, Z, F)
      integer, intent(in) :: neq; real(R8), intent(in) :: time, Z(neq); real(R8), intent(out) :: F(neq)
      call rhsEvaporation(neq,time,Z,F, auxLocal,nauxvar, stateLocal,nauxstate, hTab,rhoTab,psatTab, gas,gasVert,gasState, ng,ngVert)
    end subroutine rhs2
    subroutine rhs3(neq, time, Z, F)
      integer, intent(in) :: neq; real(R8), intent(in) :: time, Z(neq); real(R8), intent(out) :: F(neq)
      call rhsBreakupOnly(neq,time,Z,F, auxLocal,nauxvar, stateLocal,nauxstate, hTab,rhoTab,mupTab,sigTab, gas,gasVert,gasState, ng,ngVert)
    end subroutine rhs3
    subroutine rhs4(neq, time, Z, F)
      integer, intent(in) :: neq; real(R8), intent(in) :: time, Z(neq); real(R8), intent(out) :: F(neq)
      call rhsEvapBreakup(neq,time,Z,F, auxLocal,nauxvar, stateLocal,nauxstate, hTab,rhoTab,mupTab,sigTab,psatTab, gas,gasVert,gasState, ng,ngVert)
    end subroutine rhs4
    subroutine rhs5(neq, time, Z, F)
      integer, intent(in) :: neq; real(R8), intent(in) :: time, Z(neq); real(R8), intent(out) :: F(neq)
      call rhsAlCombustion(neq,time,Z,F, auxLocal,nauxvar, stateLocal,nauxstate, hTab,rhoTab, gas,gasVert,gasState, ng,ngVert)
    end subroutine rhs5

    !****************************************************************************************************!
    !*  Solver output callback                                                                          *!
    !****************************************************************************************************!

    subroutine solout(NR,XOLD,X,Y,N,IRTRN)
      use IGLOO_RayFaceIntersection3D, only: isPointInsideCell
      use IGLOO_Lib_Breakup,           only: breakupEvent
      use Lib_Equations,               only: interphase
      use IGLOO_Lib_Properties,        only: comp_TfromTab, lookupTab
      implicit none
      integer  :: NR, N
      real(R8) :: X, XOLD, Y(N)
      integer  :: IRTRN
      real(R8) :: Vdif(3), vel, Re, Fdrag(3), Qdot, tp, acc(3), d, rho, sigma, mup, m, np
      real(R8) :: brkupState(max(nbrkst,1))
      real(R8), parameter :: oneThird=0.3333333333333333_R8, sixOverPi=1.90985931710274403_R8

      IRTRN = 1

      !> Scatter cloud: accumulate the npdot-weight of THIS accepted step (every step, incl.
      !  the one that ends on a crossing below, so the sampling stays mesh-independent).
      if (scatOn .and. dNscat>0._R8) then
        select case(mod_model)
        case(3);      np = y(8)            ! model 3: npdot is stateVar(8)
        case(4);      np = y(9)            ! model 4: npdot is stateVar(9)
        case default; np = part%npdot      ! models 1,2: npdot constant
        end select
        wAcc = wAcc + np*(x - xold)
      endif

      !> Geo containment only when the gas cell is on a boundary  (ord2==true)
      if (.not.ord2 .or. atGasBoundary) then
            IamOut = .not. isPointInsideCell(y(1:3),vert,mesh2D,geoHexNorms,geoHexCentroids,geoHexDegen,part%exitFace)
      else; IamOut = .false.
      endif
      if     ( ord2 ) then; newGas = .not.isPointInsideCell(y(1:3),gasVert,mesh2D,gasHexNorms,gasHexCentroids,gasHexDegen,part%gasExitFace)
      elseif (IamOut) then; newGas = .true.; endif

      !> Locator-lost segment (updateCell cycle-breaker fired: no cell owns the start point) =
      !  containment sliver, not a crossing: suppress containment interrupts until the particle
      !  has moved a real distance (escapeDist) off the sliver, instead of looping on zero progress.
      if (NR == 1) startedOut = (IamOut .or. newGas) .and. part%lost
      if (startedOut .and. norm2(y(1:3)-part%oldState(1:3)) <= escapeDist) then
        IamOut = .false.; newGas = .false.
      endif

      !> Burnout test, models 2 and 5 only: there y(8) falls by consumption alone, so it is an
      !  honest "how much is left" and the droplet can be declared consumed while the state is
      !  still good. Model 4 is excluded because breakup shedding also lowers y(8), so a small
      !  mass may mean "stripped into children" -- it falls back to the [WARNING] path instead.
      !  Models 1/3 never lose mass to consumption.
      if (mod_model==2 .or. mod_model==5) then
        if (y(8) <= mBurnTol) burnedOut = .true.   ! y(8) IS the droplet mass [kg]
      endif

      exitLoop = (norm2(y(1:3)-oldLocal(1:3))<eps).or.(deltat<dtMin).or.any(y/=y).or.(deltat /= deltat)

      if (eventType) then
        if (mod_propFlags(1)) then; tp = comp_TfromTab(hTab,y(7)); else; tp = y(7);                 endif
        if (mod_propFlags(2)) then; rho   = lookupTab(rhoTab, tp); else; rho   = auxLocal(ind_rho); endif
        if (mod_propFlags(3)) then; sigma = lookupTab(sigTab, tp); else; sigma = auxLocal(ind_sig); endif
        if (mod_propFlags(4)) then; mup   = lookupTab(mupTab, tp); else; mup   = auxLocal(ind_mup); endif
        select case(mod_model)
        case(3);      m = auxLocal(ind_m)/y(8)
        case(2,4,5);  m = y(8)
        case default; m = auxLocal(ind_m)
        end select
        if (ind_d > 0) then; d = auxLocal(ind_d); else; d = (sixOverPi*m/rho)**oneThird; endif

        Vdif = gasState(2:4)-y(4:6)
        vel  = norm2(Vdif)
        Re   = gasState(1)*vel*d/gasState(6)
        call interphase(gasState(:),ng,Vdif,vel,tp,d,Re,1._R8,Fdrag)
        acc  = Fdrag/m

        eventLocal(ind_evd) = d
        if (ind_sb1 > 0) brkupState(1:nbrkst) = stateLocal(ind_sb1:ind_sb2)
        !> Advance the analytic oscillator by the completed accepted-step interval (x-xold);
        !  deltat is the outer per-segment cap and is untied to elapsed time.
        call breakupEvent(eventLocal, neventvar, brkupState, nbrkst,    &
                          sigma,mup,rho,gasState(1),vel,Re,t1,acc,y(4:6),x-xold, &
                          mod_brkSelect, mod_bp,mod_bpMethod,mod_bpScale,        &
                          eventFlag,childState,addChildLocal,exitLoop)
        if (ind_sb1 > 0) stateLocal(ind_sb1:ind_sb2) = brkupState(1:nbrkst)
      endif

      if (IamOut .or. newGas .or. eventFlag .or. burnedOut) then
        IRTRN = -2724
        return
      endif

      !> Never latch a non-finite state as the "last good" one: it is the fallback the outer
      !  loop reverts to on burnout, and the scatter cloud below would emit NaN markers.
      if (any(y/=y)) return

      timeLocal  = x
      oldLocal   = y
      if (allocated(stateLocal)) oldStLocal = stateLocal
      if (allocated(eventLocal)) oldEvLocal = eventLocal

      !> Scatter cloud: this is a REAL accepted interior state (crossings returned above). Emit
      !  ONE marker per weight-quantum dNscat, carrying the remainder => distinct real states,
      !  no interpolation. Fields recovered from y exactly as the breakup branch does.
      if (scatOn .and. dNscat>0._R8 .and. wAcc >= dNscat) then
        if (mod_propFlags(1)) then; tp  = comp_TfromTab(hTab,y(7)); else; tp  = y(7);              endif
        if (mod_propFlags(2)) then; rho = lookupTab(rhoTab, tp);    else; rho = auxLocal(ind_rho); endif
        select case(mod_model)
        case(3);      m = auxLocal(ind_m)/y(8)
        case(2,4,5);  m = y(8)
        case default; m = auxLocal(ind_m)
        end select
        if (ind_d > 0) then; d = auxLocal(ind_d); else; d = (sixOverPi*m/rho)**oneThird; endif
        write(unit=unitScat,fmt='(7F12.6,2E13.6E2,I8)') y(1:3), y(4:6), tp, d, m, part%ID
        wAcc = wAcc - dNscat
      endif

    end subroutine solout

    !****************************************************************************************************!
    !*  Particle "position at boundary" handling                                                        *!
    !****************************************************************************************************!

    !> True when the gas dual cell sits on a block boundary
    subroutine setAtGasBoundary()
      atGasBoundary =  part%igas(1)==1 .or. part%igas(1)==gasblock(b)%Nx+1 .or. &
                       part%igas(2)==1 .or. part%igas(2)==gasblock(b)%Ny+1
      if (.not.mesh2D) atGasBoundary = atGasBoundary .or.                      &
                       part%igas(3)==1 .or. part%igas(3)==gasblock(b)%Nz+1
    end subroutine setAtGasBoundary

    !> After updateCell processed at geo crossing, refresh the geo position
    !  if the particle reflected (not gone), re-resolve the gas cell
    !  from the new position (gasExitFace is stale post-reflection -> reset it).
    subroutine handleGeoEvent()
      b = part%i(1)
      call geoblock(b)%getVertices(part%i(2:4),vert, norms=geoHexNorms, centroids=geoHexCentroids, degen=geoHexDegen)
      if (.not.part%gone) then
        part%gasExitFace = 0
        call part%advanceGasCell(gasblock(b))
        call gasblock(b)%getVertices(part%igas,gasVert,mesh2D, gasHexNorms, gasHexCentroids, gasHexDegen)
        call gasblock(b)%gasProperties(gas,part%igas)
        call setAtGasBoundary()
      endif
    end subroutine handleGeoEvent


    !****************************************************************************************************!
    !*  Eulerian & source field accumulators                                                            *!
    !****************************************************************************************************!

    subroutine computeSrcField(particle,sblock,ii,jj,kk,masIn,masOut,momIn,momOut,enIn,enOut)
      implicit none
      class(obj_particle),   intent(in)    :: particle
      type(obj_sourceblock), intent(inout) :: sblock
      real(R8),              intent(inout) :: masIn, momIn(3), enIn
      integer,               intent(in)    :: ii, jj, kk
      real(R8),              intent(in)    :: masOut, momOut(3), enOut
      integer :: nmat

      nmat = particle%mID
      !> Mass deposits when the material loses mass to the gas: evaporation OR combustion (model 5)
      if (phaseChange .or. particle%model==5) then
        !$OMP ATOMIC UPDATE
        sblock%sourceMass(nmat,ii,jj,kk) = sblock%sourceMass(nmat,ii,jj,kk) + (masIn-masOut)
      else
        sblock%sourceMass(nmat,ii,jj,kk) = 0.0
      endif
      !$OMP ATOMIC UPDATE
      sblock%sourceMom(1,ii,jj,kk) = sblock%sourceMom(1,ii,jj,kk) + (momIn(1)-momOut(1))
      !$OMP ATOMIC UPDATE
      sblock%sourceMom(2,ii,jj,kk) = sblock%sourceMom(2,ii,jj,kk) + (momIn(2)-momOut(2))
      !$OMP ATOMIC UPDATE
      sblock%sourceMom(3,ii,jj,kk) = sblock%sourceMom(3,ii,jj,kk) + (momIn(3)-momOut(3))
      !$OMP ATOMIC UPDATE
      sblock%sourceEn(ii,jj,kk)    = sblock%sourceEn(ii,jj,kk) + (enIn-enOut)

    end subroutine computeSrcField

    subroutine computeEulField(particle,eblock,ii,jj,kk)
      implicit none
      class(obj_particle),  intent(in)    :: particle
      type(obj_eulerblock), intent(inout) :: eblock
      integer,              intent(in)    :: ii, jj, kk
      real(R8) :: rho, factor, np, moment(3), energy

      if (particle%deltaL<toll) return
      !> Normalize by the DEPOSITION cell volume (ii,jj,kk):
      !  gas dual cell under ord2, geo cell (== i,j,k) otherwise.
      if (ord2) then; vol = gasblock(b)%cellVol(ii, jj, kk)
      else;           vol = geoblock(b)%cellVol(ii, jj, kk); endif

      factor = 1._R8/vol
      select case(part%model)
      case(4)
        np  = particle%intE(6)*factor
        rho = particle%intE(5)*factor
        factor = factor*particle%Tstay/particle%deltaL
      case(3)
        np  = particle%intE(5)            *factor
        rho = particle%mdot*particle%Tstay*factor
        factor = rho/particle%deltaL
      case(2,5)
        np  = particle%npdot*particle%Tstay*factor
        rho = particle%intE(5)             *factor
        factor = np /particle%deltaL
      case default
        np  = particle%npold*particle%Tstay*factor
        rho = particle%mdot *particle%Tstay*factor
        factor = rho/particle%deltaL
      end select
      moment = particle%intE(1:3)*factor
      energy = particle%intE( 4 )*factor !> T*rho_p (cpVariable=false) or h*rho_p (cpVariable=true)

      !> Pure += accumulators. velocity/temperature store numerators (Σ moment, Σ energy)
      !$OMP ATOMIC UPDATE
      eblock%density(ii,jj,kk)     = eblock%density(ii,jj,kk)     + rho
      !$OMP ATOMIC UPDATE
      eblock%np(ii,jj,kk)          = eblock%np(ii,jj,kk)          + np
      !$OMP ATOMIC UPDATE
      eblock%velocity(1,ii,jj,kk)  = eblock%velocity(1,ii,jj,kk)  + moment(1)
      !$OMP ATOMIC UPDATE
      eblock%velocity(2,ii,jj,kk)  = eblock%velocity(2,ii,jj,kk)  + moment(2)
      !$OMP ATOMIC UPDATE
      eblock%velocity(3,ii,jj,kk)  = eblock%velocity(3,ii,jj,kk)  + moment(3)
      !$OMP ATOMIC UPDATE
      eblock%temperature(ii,jj,kk) = eblock%temperature(ii,jj,kk) + energy

    end subroutine computeEulField

  end subroutine integrate

end module Lib_Integration
