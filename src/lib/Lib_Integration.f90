module Lib_Integration
  use, intrinsic :: iso_fortran_env, only : R8 => real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private
  public :: integrate


contains

  !> Integrate one particle from its current state to domain exit (or burnout) through the
  !  steady gas field, one cell segment at a time, depositing source and Eulerian moments.
  subroutine integrate(part,geoblock,gasblock,srcblock,eulblock,     &
                        hTab,cpTab,rhoTab,mupTab,sigTab,psatTab, shed,noShed)
    use, intrinsic :: iso_fortran_env, only : R8 => real64
    use IGLOO_data_phases, only: obj_shed, shedList
    use IGLOO_variables,  only: unitTraj,unitExit,unitScat,iprint,dtprint, &
                                nb,toll,ord2,mesh2D,threshold,            &
                                eulerSwitch,sourceSwitch, phaseChange,     &
                                bodyAccel, srcBodyForce, axisym,           &
                                nSectorFold, nMultiFold,                   &
                                trajOn, scatOn, dNscat, sixOverPi
    use IGLOO_particles,  only: obj_particle, eps
    use IGLOO_bcBox,      only: axisymFold, sectorDs
    use IGLOO_data_block, only: obj_block, obj_flowblock, obj_eulerblock, obj_sourceblock
    use IGLOO_Lib_Properties
    use Lib_Equations
    use Lib_RHS, only: packAuxVars, nauxvar, packAuxState, unpackAuxState,     &
                       packEventVar, unpackEventVar,                           &
                       rhsStandard, rhsEvaporation, rhsBreakupOnly, rhsEvapBreakup, &
                       rhsAlCombustion,                                             &
                       ind_d, ind_rho, ind_sig, ind_mup, ind_m,                     &
                       ind_sb1, ind_sb2,                                            &
                       ind_evd, ind_evn,                                            &
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
    !> This parent's shed list (intent(inout) required).
    type(shedList),        intent(inout), optional :: shed
    !> True on the last generation pass: sheds are suppressed and the parent keeps the mass.
    logical,               intent(in),    optional :: noShed
    !> local variables
    type(obj_shed) :: shedRec   !> staged child record, pushed once complete
    real(R8), allocatable :: gas(:,:), gasVert(:,:), gasState(:)
    real(R8) :: vert(3,8), t1, t2, tStart, deltat, tlimit, tprint
    real(R8) :: Ein, Eout, Pin(3), Pout(3), massIn, massOut, vol
    real(R8) :: entryPos(3)   ! cell-entry position, for the closed-form body-force work (models 1,3)
    real(R8) :: pMid(3)       ! segment mid position: azimuth of the meridian-frame source deposit
    real(R8) :: wAcc          ! scatter-cloud npdot-weight accumulator (host-associated into solout)
    !> Scatter-cloud output buffer, flushed by flushScat in one write. No initializer: it would imply SAVE.
    integer, parameter :: SCATLEN = 128, SCATCAP = 4096
    character(len=SCATLEN) :: scatBuf(SCATCAP)
    integer  :: nScat, iScat
    integer  :: neq, nsp, ng, b, i, j, k, ngVert, iu, iw, it, iter, nCross, nE, nStall
    real(R8) :: posPrev(3)    ! last outer-iter position, for the zero-progress guard
    logical  :: doLoop, IamOut, newGas, eventType, eventFlag, exitLoop, startedOut, burnedOut
    logical  :: sectorOut  !> wedge: the segment ended on a k-plane
    logical  :: foldOnly   !> sector fold with no cell crossing
    integer  :: nSect
    logical  :: consumed   !> droplet ended INSIDE the domain with mass still on it
    logical  :: atGasBoundary, wasBoundary   ! ord2: geo consulted only at gas-boundary cells
    real(R8) :: geoHexNorms(3,2,6), geoHexCentroids(3,2,6)
    real(R8) :: gasHexNorms(3,2,6), gasHexCentroids(3,2,6)
    logical  :: geoHexDegen(2,6), gasHexDegen(2,6)
    integer,  parameter :: maxIter=500000, nMaxCell=10, nMaxStall=100
    !> accessed by ODEsystem, rhs1-4, solout via host association
    real(R8) :: y(part%neq), oldLocal(part%neq)
    real(R8) :: auxLocal(nauxvar), childState(nchild)
    logical  :: addChildLocal   !> breakupEvent's shed flag
    logical  :: childDone       !> sheds suppressed for this whole call
    real(R8), allocatable :: stateLocal(:), oldStLocal(:)
    real(R8), allocatable :: eventLocal(:), oldEvLocal(:)
    !> ETAB product-velocity kick, carried out of the aborted step as a delta.
    real(R8) :: kickDV(3)
    logical  :: kickPend
    real(R8) :: timeLocal, din, dout, deltaS(3), dir(3), taup
    integer  :: err, nDL, innerIter, jSlot
    integer,  parameter :: maxInnerIter=10, nStep=10
    !> Hard cap on child parcels one parent may shed in one call.
    integer,  parameter :: maxShed=1000
    real(R8), parameter :: safety=0.99, dtMin=1.e-14_R8, tauFactor=100._R8
    real(R8), parameter :: escapeDist=1.e-6_R8   ! min displacement per sliver-escape hop (startedOut)
    !> Burnout threshold on droplet mass [kg]: a consuming droplet at m <= mBurnTol is declared consumed.
    real(R8), parameter :: mBurnTol=1.e-15_R8

    nScat = 0            !> before every exit path, so flushScat is always well-defined
    tlimit = huge(1._R8) !> steady-state by default (temporary)
    neq = part%neq       !> save local copy of neq
    atGasBoundary = .true.   !> conservative default: geo consulted until proven interior
    nsp = 1              !> nsp = nspecies(part%Iinj(1)) !> for multi-species gas
    iu  = nsp + 1
    ng  = iu + 7 ! +1 if eddy viscosity (mut) included

    !> ARRAYS ALLOCATION
    addChildLocal = .false.
    childDone     = .false.
    if (present(noShed)) childDone = noShed
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
        call flushScat()
        return
      endif

      b = part%i(1)
      if (ord2) then; call part%findDualCell(gasblock(b)); part%igasOld = part%igas; part%xi0 = 0.5_R8
                      call setAtGasBoundary()
      else;                part%igas = part%i(2:4); endif
      call gasblock(b)%gasProperties(gas,part%igas)   ! NOT hoistable: consumed just below
      if (all(part%iInj==[0,0,0,0])) then
        part%stateVar(4:6) = part%vInj   ! DB hand-off (pin_particles can't write stateVar pre-allocation)
        if (any(part%stateVar(4:6)>=threshold)) part%stateVar(4:6) = gas(iu:iu+2,1)
        if (part%tp==0._R8) part%tp = gas(iu+3,1)
      else
        call part%initializePart(geoblock,gasblock)
      endif
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
    !> Block index and cell-entry state for every entry (children enter with time /= 0).
    b = part%i(1)
    call refreshCellEntry()
    !> Pack auxiliary variables once (constant throughout integration)
    call packAuxVars(part, nauxvar, auxLocal)

    !> TRAJECTORY INTEGRATION LOOP
    iter = 0; nCross=0; tprint = part%time + dtprint
    nStall = 0; posPrev = part%stateVar(1:3) - 1._R8
    !> Scatter accumulator, seeded with a per-ID golden-ratio phase offset.
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
      !> Drag relaxation time taup; a non-finite deltat falls back to tauFactor*taup.
      if (part%d > 0._R8 .and. any(gas(6,:) > toll)) taup   = part%rho * part%d**2 / (18._R8 * maxval(gas(6,:)))
      if (.not. ieee_is_finite(deltat)) deltat = tauFactor * taup

      tStart = part%time
      t2     = tStart + 20.0_R8 * deltat
      doLoop = .true.
      consumed = .false.
      innerIter = 0
      do while (doLoop)
        innerIter = innerIter + 1
        !> Hard cap on the segment loop.
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
            
      if (sourceSwitch) then
        call part%computeSource(massOut,Pout,Eout)
        !> A droplet consumed inside the cell has zero outgoing flux.
        if (consumed) then
          massOut = 0._R8
          Pout    = 0._R8
          Eout    = 0._R8
        endif
        !> Body-force source-reaction correction.
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
        !> Wedge: rotate the momentum exchange into the meridian frame at the segment's mid azimuth.
        if (axisym) then
          pMid = 0.5_R8 * (entryPos + part%stateVar(1:3))
          Pin  = toMeridian(Pin,  pMid)
          Pout = toMeridian(Pout, pMid)
        endif
        if (ord2) then
          call computeSrcField(part, srcblock(b), part%igas(1), part%igas(2), part%igas(3), &
                               massIn, massOut, Pin, Pout, Ein, Eout)
        else
          call computeSrcField(part, srcblock(b), i, j, k, massIn, massOut, Pin, Pout, Ein, Eout)
        endif
      endif

      !> Eulerian moments, deposited on the gas dual cell (ord2) or the geo cell.
      if (eulerSwitch) then
        if (ord2) then
          call computeEulField(part, eulblock(b), part%igas(1), part%igas(2), part%igas(3))
        else
          call computeEulField(part, eulblock(b), i, j, k)
        endif
      endif

      !> Wedge sector exit: fold the azimuth back into the band after the deposits, before the cell logic.
      foldOnly = sectorOut .and. .not.(newGas .or. IamOut)
      if (sectorOut) then
        call axisymFold(part%stateVar, nSect=nSect)
        if (nSect > 0) then
          !$OMP ATOMIC
          nSectorFold = nSectorFold + 1
        else
          part%Ncell = part%Ncell + 1   ! no fold: counted as residency
        endif
        if (nSect > 1) then
          !$OMP ATOMIC
          nMultiFold = nMultiFold + 1
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
        !> Geo domain exit / BC handling.
        if (IamOut) then
          call part%updateCell(geoblock); call handleGeoEvent()
        endif
        !> Residency count on the gas cell; a pure sector fold is neither a crossing nor residency.
        if (newGas .or. IamOut) then; nCross = nCross + 1; part%Ncell = 0
        elseif (.not.foldOnly) then;  part%Ncell = part%Ncell + 1; endif
      else
        !> ord1: geo-driven.
        if (IamOut) then
          call part%updateCell(geoblock); b = part%i(1)
          call geoblock(b)%getVertices(part%i(2:4),vert, norms=geoHexNorms, centroids=geoHexCentroids, degen=geoHexDegen)
        endif
        if (newGas) then
          part%igas = part%i(2:4)
          call gasblock(b)%gasProperties(gas,part%igas)
        endif
        if (all(part%i==part%iold)) then; if (.not.foldOnly) part%Ncell = part%Ncell+1
        elseif (IamOut) then;   nCross = nCross + 1;   part%Ncell = 0; endif
      endif

      !> Wedge safety net: fold a segment that ended past a k-plane together with an x-y crossing.
      if (axisym) then
        call axisymFold(part%stateVar, nSect=nSect)
        if (nSect > 0) then
          !$OMP ATOMIC
          nSectorFold = nSectorFold + 1
        endif
        if (nSect > 1) then
          !$OMP ATOMIC
          nMultiFold = nMultiFold + 1
        endif
      endif

      if (part%Ncell>nMaxCell) then
        write(*,'(A,I4,A,4I4,A)') '       ==> Particle ',part%ID,' stuck in cell:',part%iold
        part%gone=.true.
      endif

      !> Zero-progress guard on displacement.
      if (norm2(part%stateVar(1:3)-posPrev) > eps) then; nStall = 0
      else;                                              nStall = nStall + 1; endif
      posPrev = part%stateVar(1:3)
      if (nStall > nMaxStall) then
        write(*,'(A,I4,A)') '       ==> Particle ',part%ID,' no net progress ==> marking gone'
        part%gone = .true.
      endif

      if ((dtprint>0._R8.and.part%time>=tprint).or.(mod(nCross,iprint)==0.and.part%Ncell==0.and..not.foldOnly)) then
        if (trajOn) write(unit=unitTraj,fmt='(7F12.6,2E13.6E2,I8)') part%stateVar(1:6), part%tp, part%d, part%m, part%ID
        tprint = tprint + dtprint
      endif
      if (part%gone) then
        part%time = 0._R8
        write(unit=unitExit,fmt='(6F12.6,2E13.6E2,I8)')                                      &
              part%stateVar(1:3), part%tp, norm2(part%stateVar(4:6)), part%angle, part%mdot, part%Af, part%ID
        call flushScat()
        return
      endif

    enddo

    !> Outer-loop maxIter exit: flag the particle gone and write its exit record.
    if (iter >= maxIter .and. .not. part%gone) then
      write(*,'(A,I4,A,I0,A)') '       ==> Particle ',part%ID,' hit outer maxIter (',maxIter,'); flagging gone'
      part%gone = .true.
      part%time = 0._R8
      write(unit=unitExit,fmt='(6F12.6,2E13.6E2,I8)')                                          &
            part%stateVar(1:3), part%tp, norm2(part%stateVar(4:6)), part%angle, part%mdot, part%Af, part%ID
    endif
    !> Every exit path flushes the scatter buffer.
    call flushScat()

  contains

    !> Emit the buffered scatter records in one write statement.
    subroutine flushScat()
      if (nScat <= 0) return
      write(unit=unitScat,fmt='(A)') (trim(scatBuf(iScat)), iScat=1,nScat)
      nScat = 0
    end subroutine flushScat


    !> Integrate one trajectory segment: run the ODE solver on [t1,t2] under solout interrupts,
    !  then finalize the segment (event/aux state, euler moments) or refine deltat and retry.
    subroutine ODEsystem()
      implicit none
      real(R8) :: dtNew, dtSafe

      deltat = safety*deltat
      newGas = .false.; IamOut = .false.; exitLoop = .false.; eventFlag = .false.; startedOut = .false.
      burnedOut = .false.; sectorOut = .false.
      kickPend = .false.; kickDV = 0._R8
      eventType = part%brkupEvent
      addChildLocal = .false.; childState = 0._R8
      !> childDone is reset at integrate entry only, never per segment.
      y = part%oldState
      timeLocal = part%time
      oldLocal  = y
      if (allocated(stateLocal)) then
        call packAuxState(part, nauxstate, stateLocal)
        !> Segment-start snapshot of the aux state.
        oldStLocal = stateLocal
      endif
      if (allocated(eventLocal)) then
        call packEventVar(part, neventvar, eventLocal)
        !> Segment-start snapshot of the event state.
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
      !> Burnout: the droplet was consumed on a still-good state; falls through to the normal finalize.
      if (burnedOut) then
        part%gone = .true.
        consumed  = .true.     !> remnant transfers to the gas; see the source block
      endif

      !> Non-finite state: revert to the last good accepted step and mark the particle gone.
      if (any(y/=y)) then
        write(*,'(A,I4,A)') '[WARNING] Particle ',part%ID,                                &
                            ' non-finite state ==> reverting to last good step, marking gone'
        !> Report the remaining mass fraction for mass-evolving models.
        if (mod_model==2 .or. mod_model==4 .or. mod_model==5)                              &
          write(*,'(A,ES10.3,A,ES10.3,A)') '          last good m/m0 =',                   &
                oldLocal(8)/part%m0, '  (burnout fires at m =', mBurnTol, ' kg)'
        part%stateVar = oldLocal
        part%time     = timeLocal
        part%oldState = oldLocal
        part%gone     = .true.
        doLoop        = .false.
        !> Consumption models hand the remnant to the gas as on a clean burnout.
        if (mod_model==2 .or. mod_model==5) consumed = .true.
        !> Sync d/m/tp from the reverted state.
        call part%updatePart(rhoTab,hTab,eulerSwitch)
        return
      endif
      part%stateVar = y
      part%time     = timeLocal
      part%oldState = oldLocal
      !> Apply the ETAB velocity kick to both the segment result and the resume state.
      if (kickPend) then
        part%stateVar(4:6) = part%stateVar(4:6) + kickDV
        part%oldState(4:6) = part%oldState(4:6) + kickDV
      endif
      !> Segment consumed with no interrupt: finalize, never re-invoke on the null interval.
      if (.not.(IamOut .or. newGas .or. eventFlag .or. sectorOut)) exitLoop = .true.
      !> Sliver-escape hop: finalize at the moved position without refinement.
      if (startedOut) exitLoop = .true.
      if (exitLoop) then
        doLoop = .false.
        if (eventType) then
          if (allocated(eventLocal)) then
            part%npold = oldEvLocal(ind_evn)
            call unpackEventVar(eventLocal, part, neventvar)
            !> Resync the frozen per-droplet d and m in auxLocal after an event.
            if (ind_d > 0) auxLocal(ind_d) = part%d
            if (ind_m > 0 .and. mod_model /= 3) then
              part%m          = part%rho * part%d**3 / sixOverPi
              auxLocal(ind_m) = part%m
            endif
            !> Model 3: push the event's npdot into y(8) and re-derive mdot from (d, npdot).
            if (eventFlag .and. mod_model == 3 .and. part%d > 0._R8 .and. ind_m > 0) then
              part%stateVar(8) = part%npdot
              part%oldState(8) = part%npdot
              part%mdot        = part%npdot * part%rho * part%d**3 / sixOverPi
              auxLocal(ind_m)  = part%mdot
            endif
          endif
        endif
        if (allocated(stateLocal)) call unpackAuxState(stateLocal, part, nauxstate)
        if (eulerSwitch) then
          !> Euler moments end at nE (the body-force slot W trails them).
          nE = neq; if (part%bodyAccum) nE = neq - 1
          part%deltaL           = y(nDL)
          part%intE(1:(nE-nDL)) = y(nDL+1:nE)
        endif
        part%time  = t1
        part%Tstay = t1 - tStart
        call part%updatePart(rhoTab,hTab,eulerSwitch)
        !> Capture the child record at the shed point (after updatePart, before the crossing logic)
        !  and push it onto the parent's shed list.
        if (eventType .and. present(shed)) then
          if (addChildLocal .and. .not. childDone) then
            shedRec%vel   = childState(1:3)
            shedRec%diam  = childState(4)
            shedRec%npdot = childState(5)
            shedRec%ipos  = part%i
            shedRec%igas  = part%igas
            shedRec%pos   = part%stateVar(1:3)
            shedRec%temp  = part%tp
            shedRec%time  = t1
            call shed%push(shedRec)
            !> Shed cap: suppress further sheds, the parent keeps the mass.
            if (shed%n >= maxShed) then
              childDone = .true.
              write(*,'(A,I0,A,I0,A)') '[WARNING] Particle ',part%ID,                  &
                    ' hit the per-call shed cap (',maxShed,'); further sheds suppressed'
            endif
          endif
        endif
      else
        if (allocated(eventLocal)) call unpackEventVar(oldEvLocal, part, neventvar)
        if (allocated(stateLocal)) call unpackAuxState(oldStLocal, part, nauxstate)
        deltaS = part%stateVar(1:3) - part%oldState(1:3)
        dout = norm2(deltaS)
        if (dout>eps) then
          dir  = deltaS/dout
          if     (IamOut)    then; din = part%computeDs(vert,dir)
          elseif (newGas)    then; din = part%computeDs(gasVert,dir)
          elseif (sectorOut) then; din = sectorDs(part%oldState(1:3),dir)   ! exact: the k-plane is a plane
          else;                    din = dout/nStep; endif
          !> Refine the step toward the crossed face; on a failed face pin keep the interior rate.
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

    !> ODE right-hand sides for models 1-5, forwarding integrate's host-associated state.
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


    !> Solver output callback: after each accepted step it tests cell/sector crossings,
    !  burnout and breakup events, and emits scatter-cloud samples.
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

      !> Scatter cloud: accumulate this step's npdot-weight.
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
            IamOut = .not. isPointInsideCell(y(1:3),vert,mesh2D,geoHexNorms,geoHexCentroids,geoHexDegen,part%exitFace,sectorOut)
      else; IamOut = .false.
      endif
      if     ( ord2 ) then; newGas = .not.isPointInsideCell(y(1:3),gasVert,mesh2D,gasHexNorms,gasHexCentroids,gasHexDegen,part%gasExitFace,sectorOut)
      elseif (IamOut) then; newGas = .true.; endif

      !> Locator-lost start: suppress containment interrupts until the particle has moved escapeDist.
      if (NR == 1) startedOut = (IamOut .or. newGas) .and. part%lost
      if (startedOut .and. norm2(y(1:3)-part%oldState(1:3)) <= escapeDist) then
        IamOut = .false.; newGas = .false.; sectorOut = .false.
      endif

      !> Burnout test on the droplet mass, consumption models only.
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
        !> Model 3 derives d from the current ODE state.
        if (ind_d > 0 .and. mod_model /= 3) then
          d = auxLocal(ind_d)
        else
          d = (sixOverPi*m/rho)**oneThird
        endif

        Vdif = gasState(2:4)-y(4:6)
        vel  = norm2(Vdif)
        Re   = gasState(1)*vel*d/gasState(6)
        call interphase(gasState(:),ng,Vdif,vel,tp,d,Re,1._R8,Fdrag)
        acc  = Fdrag/m

        eventLocal(ind_evd) = d
        !> Refresh npdot from the ODE state alongside d.
        if (ind_evn > 0) then
          select case (mod_model)
          case(3); eventLocal(ind_evn) = y(8)
          case(4); eventLocal(ind_evn) = y(9)
          end select
        endif
        if (ind_sb1 > 0) brkupState(1:nbrkst) = stateLocal(ind_sb1:ind_sb2)
        !> Breakup event over the accepted-step interval (x-xold).
        kickDV = Y(4:6)
        call breakupEvent(eventLocal, neventvar, brkupState, nbrkst,    &
                          sigma,mup,rho,gasState(1),vel,Re,t1,acc,y(4:6),x-xold, &
                          mod_brkSelect, mod_bp,mod_bpMethod,mod_bpScale,        &
                          eventFlag,childState,addChildLocal,exitLoop,part%rngState,childDone)
        !> Keep the ETAB kick as a delta; the abort below discards Y.
        kickDV = Y(4:6) - kickDV
        kickPend = (mod_brkSelect == 5) .and. (maxval(abs(kickDV)) > 0._R8)
        if (ind_sb1 > 0) stateLocal(ind_sb1:ind_sb2) = brkupState(1:nbrkst)
      endif

      if (IamOut .or. newGas .or. eventFlag .or. burnedOut .or. sectorOut) then
        IRTRN = -2724
        return
      endif

      !> Never latch a non-finite state as the last good one.
      if (any(y/=y)) return

      timeLocal  = x
      oldLocal   = y
      if (allocated(stateLocal)) oldStLocal = stateLocal
      if (allocated(eventLocal)) oldEvLocal = eventLocal

      !> Scatter cloud: emit one marker per weight quantum dNscat, carrying the remainder.
      if (scatOn .and. dNscat>0._R8 .and. wAcc >= dNscat) then
        if (mod_propFlags(1)) then; tp  = comp_TfromTab(hTab,y(7)); else; tp  = y(7);              endif
        if (mod_propFlags(2)) then; rho = lookupTab(rhoTab, tp);    else; rho = auxLocal(ind_rho); endif
        select case(mod_model)
        case(3);      m = auxLocal(ind_m)/y(8)
        case(2,4,5);  m = y(8)
        case default; m = auxLocal(ind_m)
        end select
        if (ind_d > 0) then; d = auxLocal(ind_d); else; d = (sixOverPi*m/rho)**oneThird; endif
        nScat = nScat + 1
        write(scatBuf(nScat),'(7F12.6,2E13.6E2,I8)') y(1:3), y(4:6), tp, d, m, part%ID
        if (nScat == SCATCAP) call flushScat()
        wAcc = wAcc - dNscat
      endif

    end subroutine solout


    !> Cell-entry geometry and gas state of the current cell.
    subroutine refreshCellEntry()
      call geoblock(b)%getVertices(part%i(2:4),vert, norms=geoHexNorms, centroids=geoHexCentroids, degen=geoHexDegen)   ! geo cell (NOT the gas dual index igas)
      if (ord2) then
        call gasblock(b)%getVertices(part%igas,gasVert,mesh2D, gasHexNorms, gasHexCentroids, gasHexDegen)
        call setAtGasBoundary()
      endif
      call gasblock(b)%gasProperties(gas,part%igas)
    end subroutine refreshCellEntry

    !> True when the gas dual cell sits on a block boundary
    subroutine setAtGasBoundary()
      atGasBoundary =  part%igas(1)==1 .or. part%igas(1)==gasblock(b)%Nx+1 .or. &
                       part%igas(2)==1 .or. part%igas(2)==gasblock(b)%Ny+1
      if (.not.mesh2D) atGasBoundary = atGasBoundary .or.                      &
                       part%igas(3)==1 .or. part%igas(3)==gasblock(b)%Nz+1
    end subroutine setAtGasBoundary

    !> After a geo crossing: refresh the geo cell and, if the particle stayed in, re-resolve the gas cell.
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



    !> Deposit the segment's net mass/momentum/energy exchange into source cell (ii,jj,kk).
    subroutine computeSrcField(particle,sblock,ii,jj,kk,masIn,masOut,momIn,momOut,enIn,enOut)
      implicit none
      class(obj_particle),   intent(in)    :: particle
      type(obj_sourceblock), intent(inout) :: sblock
      real(R8),              intent(inout) :: masIn, momIn(3), enIn
      integer,               intent(in)    :: ii, jj, kk
      real(R8),              intent(in)    :: masOut, momOut(3), enOut
      integer :: nmat

      nmat = particle%mID
      !> Mass source only for materials that lose mass to the gas (evaporation or combustion).
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

    !> Deposit the segment's Eulerian moments (density, number, momentum, energy) into cell (ii,jj,kk).
    subroutine computeEulField(particle,eblock,ii,jj,kk)
      implicit none
      class(obj_particle),  intent(in)    :: particle
      type(obj_eulerblock), intent(inout) :: eblock
      integer,              intent(in)    :: ii, jj, kk
      real(R8) :: rho, factor, np, moment(3), energy

      if (particle%deltaL<toll) return
      !> Normalize by the deposition cell volume.
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

      !> += accumulators; velocity/temperature hold numerators until finalize.
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
