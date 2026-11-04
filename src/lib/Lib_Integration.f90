module Lib_Integration
  use, intrinsic :: iso_fortran_env, only : R8 => real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private
  public :: integrate


contains

  !> Integrate one particle from its current state to domain exit (or burnout) through the
  !  steady gas field, one cell segment at a time, depositing source and Eulerian moments and, at the exit
  !  through a wall cell, its wall impact.
  subroutine integrate(part,geoblock,gasblock,srcblock,eulblock,     &
                        hTab,cpTab,rhoTab,mupTab,sigTab,psatTab,LvTab,cpvTab, slot,shed,noShed,resid,wallblock)
    use, intrinsic :: iso_fortran_env, only : R8 => real64
    use IGLOO_data_phases, only: obj_shed, shedList
    use IGLOO_variables,  only: iprint,dtprint,dsprint,dcellSet,tEnd,timeOn,snapOn,frozenOn, &
                                nb,toll,ord2,mesh2D,threshold,            &
                                eulerSwitch,sourceSwitch, phaseChange,     &
                                bodyAccel, srcBodyForce, axisym,           &
                                nSectorFold, nMultiFold, nStencilRefresh,  &
                                nInjNoCell,                                &
                                trajOn, scatOn, dNscat, sixOverPi, pi,    &
                                planeOn, planeX, failAfter
    use IGLOO_particles,  only: obj_particle, eps
    use IGLOO_bcBox,      only: axisymFold, sectorDs, grazeStandoff
    use IGLOO_data_block, only: obj_block, obj_flowblock, obj_eulerblock, obj_sourceblock, obj_wallblock
    use IGLOO_Lib_Erosion, only: isWallCode, erosionRatio
    use IGLOO_Lib_Properties
    use Lib_Equations
    use Lib_RHS, only: packAuxVars, nauxvar, packAuxState, unpackAuxState,     &
                       packEventVar, unpackEventVar,                           &
                       rhsStandard, rhsEvaporation, rhsBreakupOnly, rhsEvapBreakup, &
                       rhsAlCombustion, rhsSolidification, rhsAlCap,                &
                       ind_d, ind_rho, ind_sig, ind_mup, ind_m, ind_rs, ind_mb,     &
                       ind_sb1, ind_sb2, ind_sph,                                   &
                       ind_evd, ind_evn, ind_ev1, ind_ev2,                          &
                       nauxstate, neventvar, nbrkst,                                &
                       mod_brkSelect, mod_propFlags, mod_model, mod_isSolid,        &
                       mod_bp, mod_bpMethod, mod_bpScale
    use oslo,              only: Run_ODESolver, no_jacobian
    use IGLOO_Lib_Breakup, only: nchild
    use IGLOO_Lib_Combustion, only: imRr
    use IGLOO_Lib_Solidification, only: solidPhaseAtInjection, eventFunction, eventValue, eventCrossed, &
                                        eventPassed, solidTransition, solidDiameter, phUndercooled, &
                                        breakupStateAtSwitch
    use IGLOO_Mod_MPI,     only: mpi_size_, region_fail, region_failed
    use IGLOO_records,     only: recList, rec_emit, rec_emit_traj, rec_push, rec_break, &
                                 STR_TRAJ, STR_EXIT, STR_PLANE, STR_SNAP, STR_SCAT, &
                                 L_ROW, L_ROW_T, L_ROW_F, L_ROW_TF, L_ROW_S, L_ROW_TS, &
                                 L_EXIT, L_EXIT_T, L_EXIT_F, L_EXIT_TF, L_PLANE, L_SCAT, L_SCAT_F
    implicit none
    class(obj_particle),   intent(inout) :: part
    type(obj_block),       intent(in)    :: geoblock(nb)
    type(obj_flowblock),   intent(in)    :: gasblock(nb)
    type(obj_sourceblock), intent(inout) :: srcblock(nb)
    type(obj_eulerblock),  intent(inout) :: eulblock(nb)
    real(R8),              intent(in)    :: hTab(:),cpTab(:),rhoTab(:),     &
                                            mupTab(:),sigTab(:),psatTab(:),LvTab(:),cpvTab(:)
    !> The parcel's slot in its loop's ordered window of the particle files (its storage index)
    integer,               intent(in)    :: slot
    !> This parent's shed list (intent(inout) required).
    type(shedList),        intent(inout), optional :: shed
    !> True on the last generation pass: sheds are suppressed and the parent keeps the mass.
    logical,               intent(in),    optional :: noShed
    !> This parent's burnout-residue list (a capped material with residue = on).
    type(shedList),        intent(inout), optional :: resid
    !> [IGLOO-Erosion]: the wall-impact sums of every block.
    type(obj_wallblock),   intent(inout), optional :: wallblock(nb)
    !> local variables
    type(obj_shed) :: shedRec   !> staged child record, pushed once complete
    type(obj_shed) :: resRec    !> the burnout residue's record
    logical  :: resBorn         !> the parcel ended consumed with its cap: a residue is born
    real(R8) :: mRes, hRes, vRes
    real(R8), allocatable :: gas(:,:), gasVert(:,:), gasState(:)
    real(R8) :: vert(3,8), t1, t2, tStart, deltat, tprint
    real(R8) :: Ein, Eout, Pin(3), Pout(3), massIn, massOut, vol
    real(R8) :: entryPos(3)   ! cell-entry position, for the closed-form body-force work (models 1,3)
    real(R8) :: pMid(3)       ! segment mid position: azimuth of the meridian-frame source deposit
    real(R8) :: wAcc          ! scatter-cloud accumulator: npRef*h per step, a marker per dNscat (host-associated)
    real(R8) :: wAccOld       ! wAcc at the last latched state, restored when a cut step is integrated again
    !> The parcel's rows per particle file, written at its end by emitParcel
    type(recList) :: lTraj, lExit, lPlane, lSnap, lScat
    !> Parcel-rule rows (print-dcell unset): the last row's t and x, the pre-transport row, the last pushed state
    logical  :: rowsNew, haveRow, wrp
    real(R8) :: tL, xL(3), preRow(12), svRow(6), tRow
    integer  :: nPre, layPre
    integer  :: nWrap         !> 201 wraps so far: the number of the trajectory image the next one starts
    integer  :: neq, nsp, ng, b, i, j, k, ngVert, iu, iw, it, iter, nCross, nE, nStall
    integer  :: nRefl         ! gas order 1: segments of this cell that ended in a reflection
    real(R8) :: posPrev(3)    ! last outer-iter position, for the zero-progress guard
    integer  :: nReflW        ! reflections since the progress checkpoint
    real(R8) :: posReflW(3)   ! the progress checkpoint: the position nReflW reflections ago
    logical  :: reflSeg, refl ! this segment reflected off a boundary face; one updateCell's flag
    logical  :: doLoop, IamOut, newGas, eventType, eventFlag, exitLoop, startedOut, burnedOut
    logical  :: sectorOut  !> wedge: the segment ended on a k-plane
    logical  :: foldOnly   !> sector fold with no cell crossing
    integer  :: nSect
    logical  :: consumed   !> droplet ended INSIDE the domain with mass still on it
    real(R8) :: shedM, shedV(3)   !> KH shed of this segment: the child's birth mass flow and velocity
    real(R8) :: mdotSeg           !> model 3: the stream flow over this segment, before a shed at its end
    real(R8) :: hShed, vShed
    logical  :: atGasBoundary, wasBoundary   ! ord2: geo consulted only at gas-boundary cells
    real(R8) :: geoHexNorms(3,2,6), geoHexCentroids(3,2,6)
    real(R8) :: gasHexNorms(3,2,6), gasHexCentroids(3,2,6)
    logical  :: geoHexDegen(2,6), gasHexDegen(2,6)
    real(R8) :: geoHexTol, gasHexTol   !> the cached cells' containment tolerances
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
    !> Solidification event of the aborted step: flag, crossing fraction, threshold.
    logical  :: solidEv
    real(R8) :: sEv
    integer  :: whichEv
    real(R8) :: timeLocal, din, dout, deltaS(3), dir(3), taup, tEvt, rhoEvt
    integer  :: err, nDL, innerIter, jSlot
    integer,  parameter :: maxInnerIter=10, nStep=10
    !> Hard cap on child parcels one parent may shed in one call.
    integer,  parameter :: maxShed=1000
    real(R8), parameter :: safety=0.99, dtMin=1.e-14_R8, tauFactor=100._R8
    real(R8), parameter :: escapeDist=1.e-6_R8   ! min displacement per sliver-escape hop (startedOut)
    !> Burnout threshold on droplet mass [kg]: a consuming droplet at m <= mBurnTol is declared consumed.
    real(R8), parameter :: mBurnTol=1.e-15_R8

    !> MPI: a parcel of this loop failed; its rank aborts after the join (empty blocks keep the window whole)
    if (region_failed()) then
      call emitParcel()
      return
    endif
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
      if (.not. injCellHolds()) call part%updateCell(geoblock)

      if (part%gone) then
        write(*,'(A,I4,A)')                                                         &
                '[WARNING] Particle ',part%ID,' has been initialized out of domain!'
        call writeStartRecord()
        call emitParcel()
        return
      endif
      !> No cell holds the injection point: the parcel is not integrated (counted per sweep).
      if (any(part%i < 1)) then
        part%lost = .true.
        part%gone = .true.
        !$OMP ATOMIC
        nInjNoCell = nInjNoCell + 1
        call writeStartRecord()
        call emitParcel()
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
      !> the injection mass: the initial cap (model 8), a frozen fraction at rho-solid, else the injected diameter
      call part%injectionMass(rhoTab); part%m0 = part%m; part%npdot = part%mdot/part%m; part%npold = part%npdot
      part%n0 = part%npdot

      select case(part%model)
      case(2,5); part%stateVar(8) = part%m
      case(8); if (part%xiOx0 > 0._R8) then
                 !> the whole particle at the injected diameter: metal and its initial cap
                 part%stateVar(8) = (1._R8 - part%xiOx0)*part%m
                 part%stateVar(9) = part%xiOx0*part%m
               else
                 part%stateVar(8) = part%m
                 part%stateVar(9) = 0._R8
               endif
      case(3); part%stateVar(8) = part%npdot
      case(4); part%stateVar(8) = part%m
               part%stateVar(9) = part%npdot
      case(6,7); call solidPhaseAtInjection(part%tp, part%Tmelt, part%Tnuc, part%solidPhase, part%stateVar(8))
                 part%fSolid = part%stateVar(8)
                 if (part%model == 7) part%stateVar(9) = part%npdot
      end select
      if (eulerSwitch)    part%stateVar(part%nOde+1:part%neq) = 0._R8
      if (part%bodyAccum) part%stateVar(part%neq-1 :part%neq) = 0._R8
    endif
    !> Block index and cell-entry state for every entry (children enter with time /= 0).
    b = part%i(1)
    call refreshCellEntry()
    !> Pack auxiliary variables once (constant throughout integration)
    call packAuxVars(part, nauxvar, auxLocal)
    !> The injection row, then a breakup child's first row (its birth state, before its first step), on the
    !  parcel's own aux and cell gas.
    rowsNew = trajOn .and. .not.dcellSet
    haveRow = .false.
    nWrap = 0
    if (part%time==0._R8 .and. trajOn) call writeTrajRow()
    if (part%newborn) then
      if (trajOn) call writeTrajRow()
      part%newborn = .false.
    endif
    tL = part%time; xL = part%stateVar(1:3)

    !> TRAJECTORY INTEGRATION LOOP
    iter = 0; nCross=0; tprint = part%time + dtprint
    nStall = 0; posPrev = part%stateVar(1:3) - 1._R8
    nRefl = 0
    nReflW = 0; posReflW = part%stateVar(1:3)
    !> Scatter cloud: the marker rate npRef/dNscat, fixed at entry (a breakup child's is its parent's); the accumulator
    !  seeded with the golden-ratio phase of the parcel's own ID.
    if (part%time==0._R8 .or. part%npRef < 0._R8) part%npRef = part%npdot
    wAcc = dNscat * mod(real(part%ID,R8)*0.6180339887498949_R8, 1._R8)
    do while (part%time<tEnd.and.iter<maxIter)
      iter = iter+1
      call part%initializeCell(eulerSwitch)
      !> Set by the segment finalize or a revert to the last accepted step: a segment that ends neither way deposits nothing.
      part%deltaL = 0._R8; part%Tstay = 0._R8
      if (eulerSwitch) part%intE = 0._R8
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
      t2     = min(tStart + 20.0_R8 * deltat, tEnd)
      doLoop = .true.
      consumed = .false.
      shedM = 0._R8; shedV = 0._R8
      mdotSeg = part%mdot
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
      !> A capped parcel consumed (burnout, solver failure, non-finite state) with a cap: its residue is born.
      resBorn = .false.
      if (present(resid) .and. consumed .and. part%model == 8) resBorn = part%stateVar(9) > 0._R8
            
      if (sourceSwitch) then
        call part%computeSource(massOut,Pout,Eout)
        !> A droplet consumed inside the cell has zero outgoing flux.
        if (consumed) then
          massOut = 0._R8
          Pout    = 0._R8
          Eout    = 0._R8
        endif
        !> A KH shed hands the child's birth flux to the child, not to the gas.
        if (shedM > 0._R8) then
          hShed   = part%stateVar(7); if (.not.part%varCp) hShed = part%cp*hShed + part%hOff
          vShed   = norm2(shedV)
          massOut = massOut + shedM
          Pout    = Pout    + shedM*shedV
          Eout    = Eout    + shedM*(hShed + 0.5_R8*vShed*vShed)
        endif
        !> A burnout residue takes its birth flux, the cap at the product's heat capacity and datum, from the books.
        if (resBorn) then
          mRes    = part%npdot*part%stateVar(9)
          hRes    = part%cpRes*part%tp + part%hOffRes
          vRes    = norm2(part%stateVar(4:6))
          massOut = massOut + mRes
          Pout    = Pout    + mRes*part%stateVar(4:6)
          Eout    = Eout    + mRes*(hRes + 0.5_R8*vRes*vRes)
        endif
        !> Body-force source-reaction correction.
        if (srcBodyForce) then
          select case(part%model)
          case(2,5,8)    ! J,W at neq-1,neq;
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
        if (rowsNew) call trajRow(preRow, nPre, layPre)
        call axisymFold(part%stateVar, nSect=nSect)
        if (nSect > 0) then
          !$OMP ATOMIC
          nSectorFold = nSectorFold + 1
          if (rowsNew) call transportRows(.false.)
        else
          part%Ncell = part%Ncell + 1   ! no fold: counted as residency
        endif
        if (nSect > 1) then
          !$OMP ATOMIC
          nMultiFold = nMultiFold + 1
        endif
      endif

      reflSeg = .false.
      if (ord2) then
        !> Gas crossing: advance igas on the dual mesh (gas-primary tracking).
        if (newGas) then
          call part%advanceGasCell(gasblock(b))
          call gasblock(b)%getVertices(part%igas,gasVert,mesh2D, gasHexNorms, gasHexCentroids, gasHexDegen, gasHexTol)
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
            if (part%i(1) /= b) then
              !> the block changed: the gas stencil and state are re-sought in the new block
              !$OMP ATOMIC
              nStencilRefresh = nStencilRefresh + 1
              call handleGeoEvent()
            else
              call geoblock(b)%getVertices(part%i(2:4),vert, norms=geoHexNorms, centroids=geoHexCentroids, degen=geoHexDegen, &
                  tol=geoHexTol)
            endif
          else
            call relocate(); call handleGeoEvent()
          endif
        endif
        !> Geo domain exit / BC handling.
        if (IamOut) then
          call relocate(); call handleGeoEvent()
        endif
        !> Residency count on the gas cell; a pure sector fold is neither a crossing nor residency.
        if (newGas .or. IamOut) then; nCross = nCross + 1; part%Ncell = 0
        elseif (.not.foldOnly) then;  part%Ncell = part%Ncell + 1; endif
      else
        !> ord1: geo-driven.
        if (IamOut) then
          call relocate(); b = part%i(1)
          call geoblock(b)%getVertices(part%i(2:4),vert, norms=geoHexNorms, centroids=geoHexCentroids, degen=geoHexDegen, &
              tol=geoHexTol)
        endif
        if (newGas) then
          part%igas = part%i(2:4)
          call gasblock(b)%gasProperties(gas,part%igas)
        endif
        if (all(part%i==part%iold)) then; if (.not.foldOnly) part%Ncell = part%Ncell+1
          if (IamOut .and. .not.part%gone) nRefl = nRefl + 1   ! a reflection is not residency
        elseif (IamOut) then;   nCross = nCross + 1;   part%Ncell = 0; nRefl = 0; endif
      endif

      !> Wedge safety net: fold a segment that ended past a k-plane together with an x-y crossing.
      if (axisym) then
        if (rowsNew) call trajRow(preRow, nPre, layPre)
        call axisymFold(part%stateVar, nSect=nSect)
        if (nSect > 0) then
          !$OMP ATOMIC
          nSectorFold = nSectorFold + 1
          if (rowsNew) call transportRows(.false.)
        endif
        if (nSect > 1) then
          !$OMP ATOMIC
          nMultiFold = nMultiFold + 1
        endif
      endif

      if (part%Ncell-nRefl>nMaxCell) then
        write(*,'(A,I4,A,4I4,A)') '       ==> Particle ',part%ID,' stuck in cell:',part%iold
        part%gone=.true.
      endif

      !> A bounce that makes no progress: nMaxStall reflections within nMaxStall standoffs.
      if (reflSeg .and. .not.part%gone) then
        nReflW = nReflW + 1
        if (nReflW == nMaxStall) then
          if (norm2(part%stateVar(1:3) - posReflW) < nMaxStall*grazeStandoff) then
            write(*,'(A,I4,A,4I4,A)') '       ==> Particle ',part%ID,' stuck in cell:',part%iold
            part%gone = .true.
          endif
          nReflW = 0; posReflW = part%stateVar(1:3)
        endif
      endif

      !> Zero-progress guard on displacement.
      if (norm2(part%stateVar(1:3)-posPrev) > eps) then; nStall = 0
      else;                                              nStall = nStall + 1; endif
      posPrev = part%stateVar(1:3)
      if (nStall > nMaxStall) then
        write(*,'(A,I4,A)') '       ==> Particle ',part%ID,' no net progress ==> marking gone'
        part%gone = .true.
      endif

      !> MPI: end every parcel in flight after a failure, before the row decision writes its last row
      if (region_failed()) part%gone = .true.
      !> Rows: the parcel rule (print-ds, else print-dtime, from the last row; a last row at a stop) or print-dcell's cadence
      if (.not.dcellSet) then
        if (rowsNew) then
          if (part%gone .or. part%time>=tEnd) then
            call writeLastRow()
          elseif (dsprint>0._R8) then
            if (norm2(part%stateVar(1:3) - xL) >= dsprint) call writeRuleRow()
          elseif (dtprint>0._R8) then
            if (part%time - tL >= dtprint) call writeRuleRow()
          endif
        endif
      elseif ((dtprint>0._R8.and.part%time>=tprint).or.(mod(nCross,iprint)==0.and.part%Ncell==0.and..not.foldOnly) &
          .or.part%time>=tEnd) then
        if (trajOn) call writeTrajRow()
        tprint = tprint + dtprint
      !> A parcel that left through a boundary face ends its trajectory at the exit, whatever the cadence.
      elseif (trajOn .and. part%gone .and. part%exitCode /= 0) then
        call writeTrajRow()
      endif
      !> time-end reached inside the domain: the parcel's state at the stop.
      if (snapOn .and. part%time>=tEnd .and. .not.part%gone) then
        if (frozenOn) then
          call rec_push(lSnap, [part%stateVar(1:6), part%tp, part%d, part%m, real(part%ID, R8), part%time, part%fSolid], &
                          L_ROW_TF)
        else
          call rec_push(lSnap, [part%stateVar(1:6), part%tp, part%d, part%m, real(part%ID, R8), part%time], L_ROW_T)
        endif
      endif
      if (part%gone) then
        !> A parcel that left through a wall cell: its impact joins the cell's wall sums.
        if (present(wallblock) .and. part%exitCode /= 0) call depositErosion()
        !> The residue's record: the state its parent's exit record prints.
        if (resBorn) then
          resRec%pos    = part%stateVar(1:3)
          resRec%vel    = part%stateVar(4:6)
          resRec%temp   = part%tp
          resRec%time   = part%time
          resRec%ipos   = part%i
          resRec%igas   = part%igas
          resRec%mass   = part%stateVar(9)
          resRec%npdot  = part%npdot
          resRec%parent = part%ID
          call resid%push(resRec)
        endif
        call writeExitRow()
        part%time = 0._R8
        call emitParcel()
        return
      endif

    enddo

    !> Outer-loop maxIter exit: flag the particle gone and write its exit record.
    if (iter >= maxIter .and. .not. part%gone) then
      write(*,'(A,I4,A,I0,A)') '       ==> Particle ',part%ID,' hit outer maxIter (',maxIter,'); flagging gone'
      part%gone = .true.
      if (rowsNew) call writeLastRow()
      call writeExitRow()
      part%time = 0._R8
    endif
    !> Time end reached inside the domain: the record of where integration stopped.
    if (.not. part%gone .and. part%time >= tEnd) then
      if (rowsNew) call writeLastRow()
      call writeExitRow()
    endif
    !> Every exit path writes the parcel's rows.
    call emitParcel()

  contains

    !> One trajectory row of the parcel's current state.
    subroutine writeTrajRow()
      real(R8) :: row(12)
      integer  :: n, lay
      call trajRow(row, n, lay)
      call pushTraj(row(1:n), lay)
    end subroutine writeTrajRow


    !> The trajectory row of the parcel's current state, n values in layout lay; out-time appends the parcel time,
    !  a solidifying material the frozen fraction, a P2T material the surface temperature.
    subroutine trajRow(row, n, lay)
      real(R8), intent(out) :: row(12)
      integer,  intent(out) :: n, lay
      row(1:6) = part%stateVar(1:6); row(7) = part%tp; row(8) = part%d; row(9) = part%m; row(10) = real(part%ID, R8)
      n = 10; lay = L_ROW
      if (part%liqSelect == 1) then
        if (timeOn) then
          row(11) = part%time; row(12) = surfaceTemp(); n = 12; lay = L_ROW_TS
        else
          row(11) = surfaceTemp(); n = 11; lay = L_ROW_S
        endif
      elseif (timeOn .and. frozenOn) then
        row(11) = part%time; row(12) = part%fSolid; n = 12; lay = L_ROW_TF
      elseif (timeOn) then
        row(11) = part%time; n = 11; lay = L_ROW_T
      elseif (frozenOn) then
        row(11) = part%fSolid; n = 11; lay = L_ROW_F
      endif
    end subroutine trajRow


    !> Push one trajectory row; the parcel rule's last row compares the parcel's state with it.
    subroutine pushTraj(row, lay)
      real(R8), intent(in) :: row(:)
      integer,  intent(in) :: lay
      call rec_push(lTraj, row, lay)
      svRow = row(1:6); tRow = part%time; haveRow = .true.
    end subroutine pushTraj


    !> The parcel rule's row at a reached limit: it becomes the reference of the next one.
    subroutine writeRuleRow()
      call writeTrajRow()
      tL = part%time; xL = part%stateVar(1:3)
    end subroutine writeRuleRow


    !> The parcel rule's last row at a stop, unless the row just pushed holds this state.
    subroutine writeLastRow()
      if (haveRow) then
        if (tRow == part%time .and. all(svRow == part%stateVar(1:6))) return
      endif
      call writeTrajRow()
    end subroutine writeLastRow


    !> A wrap or fold under the parcel rule: preRow (taken before it), then the re-entry row, the next row's reference;
    !  a wrap's re-entry row starts the next image.
    subroutine transportRows(wrap)
      logical, intent(in) :: wrap
      call pushTraj(preRow(1:nPre), layPre)
      if (wrap) call rec_break(lTraj, nWrap)
      preRow(1:6) = part%stateVar(1:6)
      call pushTraj(preRow(1:nPre), layPre)
      tL = part%time; xL = part%stateVar(1:3)
    end subroutine transportRows


    !> updateCell from the outer loop; a 201 wrap starts the next trajectory image (in every mode) and writes its two rows
    !  under the parcel rule (the first: the segment end).
    subroutine relocate()
      if (rowsNew) call trajRow(preRow, nPre, layPre)
      call part%updateCell(geoblock, reflected=refl, wrapped=wrp)
      reflSeg = reflSeg .or. refl
      if (wrp) then
        nWrap = nWrap + 1
        if (rowsNew) then
          call transportRows(.true.)
        elseif (trajOn) then
          call rec_break(lTraj, nWrap)
        endif
      endif
    end subroutine relocate


    !> The surface temperature of a P2T parcel at its current state: its right-hand side on copies of the aux and gas
    !  states (the right-hand side writes both).
    function surfaceTemp() result(tsurf)
      real(R8) :: tsurf
      real(R8) :: Fs(part%neq), stS(max(1,nauxstate)), gasS(ng)
      call packAuxState(part, nauxstate, stS)
      gasS = gas(:,1)
      select case (part%model)
      case (2)
        call rhsEvaporation(neq, part%time, part%stateVar, Fs, auxLocal,nauxvar, stS,nauxstate, &
                            hTab,rhoTab,psatTab,LvTab,cpvTab, gas,gasVert,gasS, ng,ngVert, mupTabM=mupTab, TsOut=tsurf)
      case (4)
        call rhsEvapBreakup(neq, part%time, part%stateVar, Fs, auxLocal,nauxvar, stS,nauxstate, &
                            hTab,rhoTab,mupTab,sigTab,psatTab,LvTab,cpvTab, gas,gasVert,gasS, ng,ngVert, TsOut=tsurf)
      case default
        tsurf = part%tp
      end select
    end function surfaceTemp


    !> One exit record, written before the time reset; out-time appends the parcel time.
    subroutine writeExitRow(started)
      logical, intent(in), optional :: started   !> .false.: the parcel was never integrated
      real(R8) :: mdotOut
      !> Model 4: the stream's mass flow at exit; the other models: the injected or birth flow.
      mdotOut = part%mdot
      if (part%model == 4) mdotOut = part%npdot*part%m
      if (present(started)) then; if (.not.started) mdotOut = part%mdot; endif
      if (timeOn .and. frozenOn) then
        call rec_push(lExit, [part%stateVar(1:3), part%tp, norm2(part%stateVar(4:6)), part%angle, mdotOut, part%Af, &
                         real(part%ID, R8), part%time, real(part%exitCode, R8), part%d, part%fSolid], L_EXIT_TF)
      elseif (timeOn) then
        call rec_push(lExit, [part%stateVar(1:3), part%tp, norm2(part%stateVar(4:6)), part%angle, mdotOut, part%Af, &
                         real(part%ID, R8), part%time, real(part%exitCode, R8), part%d], L_EXIT_T)
      elseif (frozenOn) then
        call rec_push(lExit, [part%stateVar(1:3), part%tp, norm2(part%stateVar(4:6)), part%angle, mdotOut, part%Af, &
                         real(part%ID, R8), real(part%exitCode, R8), part%d, part%fSolid], L_EXIT_F)
      else
        call rec_push(lExit, [part%stateVar(1:3), part%tp, norm2(part%stateVar(4:6)), part%angle, mdotOut, part%Af, &
                         real(part%ID, R8), real(part%exitCode, R8), part%d], L_EXIT)
      endif
    end subroutine writeExitRow


    !> The record of a parcel never integrated: its injection point, at rest, at t = 0, code 0.
    subroutine writeStartRecord()
      part%stateVar(1:3) = part%pInj
      part%stateVar(4:6) = 0._R8
      part%angle = 0._R8; part%Af = 0._R8; part%exitCode = 0; part%exitCell = 0
      !> A boundary-face stream derives its temperature and flow at injection, which did not take place.
      if (.not. all(part%iInj == [0,0,0,0])) then; part%tp = 0._R8; part%mdot = 0._R8; endif
      call writeExitRow(started=.false.)
    end subroutine writeStartRecord


    !> One planes record per plane the step (ta,ya)-(tb,yb) crosses: forward x_a < X <= x_b, backward x_a >= X > x_b.
    subroutine emitPlanes(ta, ya, tb, yb)
      use IGLOO_variables, only: oneThird
      real(R8), intent(in) :: ta, tb, ya(:), yb(:)
      real(R8) :: s, pos(2), vel(3), z7, tp, rho, m, d, np, mdotR, f, tc, mAlP, mOxP
      integer  :: kp, idir
      do kp = 1, size(planeX)
        if (ya(1) < planeX(kp) .and. planeX(kp) <= yb(1)) then
          idir = 1
        elseif (ya(1) >= planeX(kp) .and. planeX(kp) > yb(1)) then
          idir = -1
        else
          cycle
        endif
        s   = (planeX(kp) - ya(1))/(yb(1) - ya(1))
        pos = ya(2:3) + s*(yb(2:3) - ya(2:3))
        vel = ya(4:6) + s*(yb(4:6) - ya(4:6))
        z7  = ya(7) + s*(yb(7) - ya(7))
        tc  = ta + s*(tb - ta)
        if (mod_propFlags(1)) then; tp = comp_TfromTab(hTab, z7); else; tp = z7; endif
        if (mod_propFlags(2)) then; rho = tableValue(rhoTab, Tmin, Tmax, tp); else; rho = auxLocal(ind_rho); endif
        !> The stream's number rate, droplet mass and mass flow at the crossing.
        select case(mod_model)
        case(3)
          np = ya(8) + s*(yb(8) - ya(8)); m = mdotSeg/np; mdotR = mdotSeg
        case(7)
          np = ya(9) + s*(yb(9) - ya(9)); m = mdotSeg/np; mdotR = mdotSeg
        case(2,5)
          np = part%npdot; m = ya(8) + s*(yb(8) - ya(8)); mdotR = np*m
        case(8)
          np = part%npdot; mAlP = ya(8) + s*(yb(8) - ya(8)); mOxP = ya(9) + s*(yb(9) - ya(9))
          m = mAlP + mOxP; mdotR = np*m
        case(4)
          np = ya(9) + s*(yb(9) - ya(9)); m = ya(8) + s*(yb(8) - ya(8)); mdotR = np*m
        case default
          np = part%npdot; m = auxLocal(ind_m); mdotR = part%mdot
        end select
        f = 0._R8
        if (mod_isSolid) f = ya(8) + s*(yb(8) - ya(8))
        if (mod_model == 8) then
          d = (sixOverPi*(mAlP/rho + mOxP/auxLocal(ind_mb+imRr-1)))**oneThird
        elseif (ind_d > 0 .and. mod_model /= 3 .and. mod_model /= 7 .and. ind_rs == 0) then
          d = auxLocal(ind_d)
        else
          d = (sixOverPi*m/rho)**oneThird
        endif
        if (ind_rs > 0) d = solidDiameter(m, f, rho, auxLocal(ind_rs))
        call rec_push(lPlane, [real(kp, R8), planeX(kp), pos, vel, tp, d, m, mdotR, np, f, real(part%ID, R8), &
                          real(part%parentID, R8), tc, real(idir, R8)], L_PLANE)
      enddo
    end subroutine emitPlanes


    !> The scatter markers of the step (ta,ya)-(tb,yb): one at each multiple of dNscat the accumulator crosses on it,
    !  at its fraction s of the step: the position on the cubic Hermite interpolant of the end positions and
    !  velocities, the rest linear; each weighs the droplets np/r it stands for.
    subroutine emitScatter(ta, ya, tb, yb)
      use IGLOO_variables, only: oneThird
      real(R8), intent(in) :: ta, tb, ya(:), yb(:)
      real(R8) :: s, h, lev, pos(3), vel(3), z7, tp, rho, m, d, np, f, mAlS, mOxS
      integer  :: k
      if (.not. (wAcc > wAccOld .and. wAcc < huge(wAcc))) return
      h = tb - ta
      k = 0
      lev = dNscat
      do while (wAcc >= lev)
        k = k + 1
        s   = (lev - wAccOld)/(wAcc - wAccOld)
        pos = (2._R8*s**3 - 3._R8*s**2 + 1._R8)*ya(1:3) + (s**3 - 2._R8*s**2 + s)*h*ya(4:6) &
            + (3._R8*s**2 - 2._R8*s**3)*yb(1:3) + (s**3 - s**2)*h*yb(4:6)
        vel = ya(4:6) + s*(yb(4:6) - ya(4:6))
        z7  = ya(7) + s*(yb(7) - ya(7))
        if (mod_propFlags(1)) then; tp = comp_TfromTab(hTab, z7); else; tp = z7; endif
        if (mod_propFlags(2)) then; rho = tableValue(rhoTab, Tmin, Tmax, tp); else; rho = auxLocal(ind_rho); endif
        select case(mod_model)
        case(3)
          np = ya(8) + s*(yb(8) - ya(8)); m = auxLocal(ind_m)/np
        case(7)
          np = ya(9) + s*(yb(9) - ya(9)); m = auxLocal(ind_m)/np
        case(2,5)
          np = part%npdot; m = ya(8) + s*(yb(8) - ya(8))
        case(4)
          np = ya(9) + s*(yb(9) - ya(9)); m = ya(8) + s*(yb(8) - ya(8))
        case(8)
          np = part%npdot; mAlS = ya(8) + s*(yb(8) - ya(8)); mOxS = ya(9) + s*(yb(9) - ya(9)); m = mAlS + mOxS
        case default
          np = part%npdot; m = auxLocal(ind_m)
        end select
        f = 0._R8
        if (mod_isSolid) f = ya(8) + s*(yb(8) - ya(8))
        if (ind_d > 0 .and. mod_model /= 3 .and. mod_model /= 7) then; d = auxLocal(ind_d)
        else; d = (sixOverPi*m/rho)**oneThird; endif
        if (ind_rs > 0) d = solidDiameter(m, f, rho, auxLocal(ind_rs))
        if (mod_model == 8) d = (sixOverPi*(mAlS/rho + mOxS/auxLocal(ind_mb+imRr-1)))**oneThird
        if (frozenOn) then
          call rec_push(lScat, [pos, vel, tp, d, m, real(part%ID, R8), f, dNscat*(np/part%npRef)], L_SCAT_F)
        else
          call rec_push(lScat, [pos, vel, tp, d, m, real(part%ID, R8), dNscat*(np/part%npRef)], L_SCAT)
        endif
        lev = real(k + 1, R8)*dNscat
      enddo
      wAcc = wAcc - real(k, R8)*dNscat
    end subroutine emitScatter


    !> The parcel's rows, one block per particle file.
    subroutine emitParcel()
      if (trajOn) call rec_emit_traj(lTraj, slot, part%ID, part%parentID)
      call rec_emit(lExit, STR_EXIT, slot, part%ID)
      if (planeOn) call rec_emit(lPlane, STR_PLANE, slot, part%ID)
      if (snapOn) call rec_emit(lSnap, STR_SNAP, slot, part%ID)
      if (scatOn) call rec_emit(lScat, STR_SCAT, slot, part%ID)
    end subroutine emitParcel


    !> Integrate one trajectory segment: run the ODE solver on [t1,t2] under solout interrupts,
    !  then finalize the segment (event/aux state, euler moments) or refine deltat and retry.
    subroutine ODEsystem()
      implicit none
      real(R8) :: dtNew, dtSafe
      integer  :: phaseEv, iNp

      deltat = safety*deltat
      newGas = .false.; IamOut = .false.; exitLoop = .false.; eventFlag = .false.; startedOut = .false.
      burnedOut = .false.; sectorOut = .false.
      kickPend = .false.; kickDV = 0._R8
      solidEv = .false.; sEv = 0._R8; whichEv = 0
      eventType = part%brkupEvent
      addChildLocal = .false.; childState = 0._R8
      !> childDone is reset at integrate entry only, never per segment.
      !> Models 6, 7: a threshold the segment-start state has strictly passed is applied here.
      if (mod_isSolid) call solidCatchUp()
      y = part%oldState
      timeLocal = part%time
      oldLocal  = y
      wAccOld   = wAcc
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
      case(1); nDL = 8;  call Run_ODESolver(neq, t1, t2, y, rhs1, no_jacobian, 0, err, deltat, solout)
      case(2); nDL = 9;  call Run_ODESolver(neq, t1, t2, y, rhs2, no_jacobian, 0, err, deltat, solout)
      case(3); nDL = 9;  call Run_ODESolver(neq, t1, t2, y, rhs3, no_jacobian, 0, err, deltat, solout)
      case(4); nDL = 10; call Run_ODESolver(neq, t1, t2, y, rhs4, no_jacobian, 0, err, deltat, solout)
      case(5); nDL = 9;  call Run_ODESolver(neq, t1, t2, y, rhs5, no_jacobian, 0, err, deltat, solout)
      case(6); nDL = 9;  call Run_ODESolver(neq, t1, t2, y, rhs6, no_jacobian, 0, err, deltat, solout)
      case(7); nDL = 10; call Run_ODESolver(neq, t1, t2, y, rhs6, no_jacobian, 0, err, deltat, solout)
      case(8); nDL = 10; call Run_ODESolver(neq, t1, t2, y, rhs8, no_jacobian, 0, err, deltat, solout)
      end select
      !> [IGLOO-ODE] fail-after: the cap model fails once its last accepted time reaches it.
      if (failAfter > 0._R8 .and. mod_model == 8 .and. err >= 0 .and. timeLocal >= failAfter) err = -1
      if (err < 0) then
        write(*,'(a,i0,a,i0,a,es12.4,a)')                                   &
          '[IGLOO] Run_ODESolver err=',err,' part=',part%ID,                   &
          ' time=',part%time,' => marking gone'
        doLoop    = .false.
        part%gone = .true.
        !> Consumption models hand the remnant to the gas, as on the non-finite exit below.
        if (mod_model==2 .or. mod_model==4 .or. mod_model==5 .or. mod_model==8) consumed = .true.
        !> Back to the last accepted step: the exit record and the deposits carry that progress.
        part%stateVar = oldLocal
        part%time     = timeLocal
        part%oldState = oldLocal
        call part%updatePart(rhoTab,hTab,eulerSwitch)
        part%Tstay = timeLocal - tStart
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
        if (mod_model==2 .or. mod_model==4 .or. mod_model==5 .or. mod_model==8)           &
          write(*,'(A,ES10.3,A,ES10.3,A)') '          last good m/m0 =',                   &
                oldLocal(8)/part%m0, '  (burnout fires at m =', mBurnTol, ' kg)'
        part%stateVar = oldLocal
        part%time     = timeLocal
        part%oldState = oldLocal
        part%gone     = .true.
        doLoop        = .false.
        !> Consumption models hand the remnant to the gas as on a clean burnout.
        if (mod_model==2 .or. mod_model==4 .or. mod_model==5 .or. mod_model==8) consumed = .true.
        !> Sync d/m/tp from the reverted state; the residence ends at the last good step.
        call part%updatePart(rhoTab,hTab,eulerSwitch)
        part%Tstay = timeLocal - tStart
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
        !> Solidification event alone on the step: interpolate to the crossing, then jump by absolute assignment.
        if (solidEv .and. .not.(IamOut .or. newGas .or. sectorOut)) then
          y  = oldLocal + sEv*(y - oldLocal)
          t1 = timeLocal + sEv*(t1 - timeLocal)
          wAcc = wAccOld + sEv*(wAcc - wAccOld)
        endif
        !> Planes and scatter markers the finalized step crosses, before the phase jump and the breakup event.
        if (planeOn) call emitPlanes(timeLocal, oldLocal, t1, y)
        if (scatOn .and. dNscat>0._R8) call emitScatter(timeLocal, oldLocal, t1, y)
        if (solidEv .and. .not.(IamOut .or. newGas .or. sectorOut)) then
          phaseEv = nint(stateLocal(ind_sph))
          call solidTransition(whichEv, part%cp, part%cpSol, part%Tmelt, part%hFus, y(7), y(8), phaseEv)
          stateLocal(ind_sph) = real(phaseEv, R8)
          !> Breakup state at the switch: the TAB/ETAB oscillator, the KH-RT timer.
          if ((mod_brkSelect == 4 .or. mod_brkSelect == 5) .and. allocated(eventLocal)) &
            eventLocal(ind_ev1:ind_ev2) = breakupStateAtSwitch(whichEv, mod_brkSelect, t1, eventLocal(ind_ev1:ind_ev2))
          if (ind_sb1 > 0) &
            stateLocal(ind_sb1:ind_sb2) = breakupStateAtSwitch(whichEv, mod_brkSelect, t1, stateLocal(ind_sb1:ind_sb2))
          part%stateVar = y
          part%oldState = y
        endif
        if (eventType) then
          if (allocated(eventLocal)) then
            part%npold = oldEvLocal(ind_evn)
            call unpackEventVar(eventLocal, part, neventvar)
            !> The density at the segment-end state, the one the event block derived d with.
            rhoEvt = part%rho
            if (mod_propFlags(2)) then
              if (mod_propFlags(1)) then; tEvt = comp_TfromTab(hTab,y(7)); else; tEvt = y(7); endif
              rhoEvt = tableValue(rhoTab, Tmin, Tmax, tEvt)
            endif
            !> Resync the frozen per-droplet d and m in auxLocal after an event.
            if (ind_d > 0) auxLocal(ind_d) = part%d
            !> A solidifying drop whose density varies (rho-solid, a table) changes mass only at a breakup event.
            if (ind_m > 0 .and. mod_model /= 3 .and. mod_model /= 7 .and. &
                (eventFlag .or. .not. mod_isSolid .or. (ind_rs == 0 .and. .not. mod_propFlags(2)))) then
              part%m          = rhoEvt * part%d**3 / sixOverPi
              auxLocal(ind_m) = part%m
            endif
            !> Models 3, 7: push the event's npdot into the number slot and re-derive mdot from (d, npdot).
            if (eventFlag .and. (mod_model == 3 .or. mod_model == 7) .and. part%d > 0._R8 .and. ind_m > 0) then
              iNp = 8; if (mod_model == 7) iNp = 9
              part%stateVar(iNp) = part%npdot
              part%oldState(iNp) = part%npdot
              part%mdot        = part%npdot * rhoEvt * part%d**3 / sixOverPi
              auxLocal(ind_m)  = part%mdot
            endif
            !> Models 2 and 4: the event's droplet mass (and model-4 npdot) become the ODE state.
            if (eventFlag .and. (mod_model == 2 .or. mod_model == 4) .and. part%d > 0._R8) then
              part%m           = rhoEvt * part%d**3 / sixOverPi
              part%stateVar(8) = part%m
              part%oldState(8) = part%m
              if (mod_model == 4) then
                part%stateVar(9) = part%npdot
                part%oldState(9) = part%npdot
              endif
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
            shedRec%phase = part%solidPhase
            shedRec%parent = part%ID
            shedRec%npRef = part%npRef
            call shed%push(shedRec)
            !> the child's birth mass flow, in its computeMass/computeSource operation order
            shedM = shedRec%npdot * (pi/6._R8*part%rho*shedRec%diam**3._R8)
            shedV = shedRec%vel
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
        !> The cut step is integrated again: its scatter weight counts once.
        wAcc = wAccOld
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
          !> MPI: no stop inside the parallel loop; flag it, drop the parcel, abort after the join.
          if (mpi_size_ > 1) then
            call region_fail('[ERROR] deltat update produced bad value')
            doLoop = .false.; part%gone = .true.
            return
          endif
          error stop '[ERROR] deltat update produced bad value'
        endif
      endif

    end subroutine ODEsystem

    !> Models 6, 7: apply the phase change whose threshold the segment-start state has strictly passed.
    subroutine solidCatchUp()
      implicit none
      real(R8) :: g
      integer  :: which

      call eventFunction(part%solidPhase, part%oldState(7), part%oldState(8), part%Tnuc, part%Tmelt, g, which)
      if (which == 0 .or. .not.eventPassed(g)) return
      call solidTransition(which, part%cp, part%cpSol, part%Tmelt, part%hFus, &
                           part%oldState(7), part%oldState(8), part%solidPhase)
      part%stateVar(7:8) = part%oldState(7:8)
      !> Breakup state at the switch, as the segment finalize sets it.
      part%brkupVar(1:2) = breakupStateAtSwitch(which, mod_brkSelect, part%time, part%brkupVar(1:2))
    end subroutine solidCatchUp

    !> ODE right-hand sides for models 1-8, forwarding integrate's host-associated state.
    subroutine rhs1(neq, time, Z, F)
      integer, intent(in) :: neq; real(R8), intent(in) :: time, Z(neq); real(R8), intent(out) :: F(neq)
      call rhsStandard(neq,time,Z,F, auxLocal,nauxvar, stateLocal,nauxstate, hTab,rhoTab, gas,gasVert,gasState, ng,ngVert)
    end subroutine rhs1
    subroutine rhs2(neq, time, Z, F)
      integer, intent(in) :: neq; real(R8), intent(in) :: time, Z(neq); real(R8), intent(out) :: F(neq)
      call rhsEvaporation(neq,time,Z,F, auxLocal,nauxvar, stateLocal,nauxstate, hTab,rhoTab,psatTab,LvTab,cpvTab, &
                          gas,gasVert,gasState, ng,ngVert, mupTabM=mupTab)
    end subroutine rhs2
    subroutine rhs3(neq, time, Z, F)
      integer, intent(in) :: neq; real(R8), intent(in) :: time, Z(neq); real(R8), intent(out) :: F(neq)
      call rhsBreakupOnly(neq,time,Z,F, auxLocal,nauxvar, stateLocal,nauxstate, hTab,rhoTab,mupTab,sigTab, gas,gasVert,gasState, ng,ngVert)
    end subroutine rhs3
    subroutine rhs4(neq, time, Z, F)
      integer, intent(in) :: neq; real(R8), intent(in) :: time, Z(neq); real(R8), intent(out) :: F(neq)
      call rhsEvapBreakup(neq,time,Z,F, auxLocal,nauxvar, stateLocal,nauxstate, hTab,rhoTab,mupTab,sigTab,psatTab, &
                          LvTab,cpvTab, gas,gasVert,gasState, ng,ngVert)
    end subroutine rhs4
    subroutine rhs5(neq, time, Z, F)
      integer, intent(in) :: neq; real(R8), intent(in) :: time, Z(neq); real(R8), intent(out) :: F(neq)
      call rhsAlCombustion(neq,time,Z,F, auxLocal,nauxvar, stateLocal,nauxstate, hTab,rhoTab, gas,gasVert,gasState, ng,ngVert)
    end subroutine rhs5
    subroutine rhs6(neq, time, Z, F)
      integer, intent(in) :: neq; real(R8), intent(in) :: time, Z(neq); real(R8), intent(out) :: F(neq)
      call rhsSolidification(neq,time,Z,F, auxLocal,nauxvar, stateLocal,nauxstate, rhoTab,mupTab,sigTab, &
                             gas,gasVert,gasState, ng,ngVert)
    end subroutine rhs6
    subroutine rhs8(neq, time, Z, F)
      integer, intent(in) :: neq; real(R8), intent(in) :: time, Z(neq); real(R8), intent(out) :: F(neq)
      call rhsAlCap(neq,time,Z,F, auxLocal,nauxvar, stateLocal,nauxstate, gas,gasVert,gasState, ng,ngVert)
    end subroutine rhs8


    !> Solver output callback: after each accepted step it tests cell/sector crossings,
    !  burnout and breakup events, and emits scatter-cloud samples.
    subroutine solout(NR,XOLD,X,Y,N,IRTRN)
      use IGLOO_RayFaceIntersection3D, only: isPointInsideCell
      use IGLOO_Lib_Breakup,           only: breakupEvent
      use Lib_Equations,               only: interphase
      use IGLOO_Lib_Properties,        only: comp_TfromTab, lookupTab, tableValue, Tmin, Tmax
      implicit none
      integer  :: NR, N
      real(R8) :: X, XOLD, Y(N)
      integer  :: IRTRN
      real(R8) :: Vdif(3), vel, Re, Fdrag(3), Qdot, tp, acc(3), d, rho, sigma, mup, m
      real(R8) :: brkupState(max(nbrkst,1))
      real(R8) :: gOld, gEnd
      real(R8) :: vEvt(3)   !> the velocity breakupEvent works on (Y stays the solver's)
      logical  :: brkOn     !> the breakup event acts on this step
      real(R8), parameter :: oneThird=0.3333333333333333_R8, sixOverPi=1.90985931710274403_R8

      IRTRN = 1

      !> Scatter cloud: accumulate this step's marker weight at the parcel's constant rate.
      if (scatOn .and. dNscat>0._R8) wAcc = wAcc + part%npRef*(x - xold)

      !> Geo containment only when the gas cell is on a boundary  (ord2==true)
      if (.not.ord2 .or. atGasBoundary) then
            IamOut = .not. isPointInsideCell(y(1:3),vert,mesh2D,geoHexNorms,geoHexCentroids,geoHexDegen,part%exitFace,sectorOut, &
                                            geoHexTol)
      else; IamOut = .false.
      endif
      if     ( ord2 ) then; newGas = .not.isPointInsideCell(y(1:3),gasVert,mesh2D,gasHexNorms,gasHexCentroids,gasHexDegen,part%gasExitFace, &
                                                         sectorOut,gasHexTol)
      elseif (IamOut) then; newGas = .true.; endif

      !> Locator-lost start: suppress containment interrupts until the particle has moved escapeDist.
      if (NR == 1) startedOut = (IamOut .or. newGas) .and. part%lost
      if (startedOut .and. norm2(y(1:3)-part%oldState(1:3)) <= escapeDist) then
        IamOut = .false.; newGas = .false.; sectorOut = .false.
      endif

      !> Burnout test on the droplet mass (model 8: the metal's), consumption models only.
      if (mod_model==2 .or. mod_model==4 .or. mod_model==5 .or. mod_model==8) then
        if (y(8) <= mBurnTol) burnedOut = .true.   ! y(8) IS the droplet mass [kg]
      endif

      !> Solidification: the step ends strictly past the phase's threshold from a start at or before it.
      if (mod_isSolid) then
        call eventFunction(nint(stateLocal(ind_sph)), y(7), y(8), part%Tnuc, part%Tmelt, gEnd, whichEv)
        if (whichEv > 0) then
          gOld = eventValue(whichEv, oldLocal(7), oldLocal(8), part%Tnuc, part%Tmelt)
          if (eventCrossed(gOld, gEnd)) then
            solidEv = .true.
            sEv     = gOld/(gOld - gEnd)
          endif
        endif
      endif

      exitLoop = (norm2(y(1:3)-oldLocal(1:3))<eps).or.(deltat<dtMin).or.any(y/=y).or.(deltat /= deltat)

      !> Breakup acts on the liquid only: the phase in force over the step gates the event.
      brkOn = eventType
      if (brkOn .and. mod_isSolid) brkOn = nint(stateLocal(ind_sph)) <= phUndercooled
      if (brkOn) then
        if (mod_propFlags(1)) then; tp = comp_TfromTab(hTab,y(7)); else; tp = y(7);                 endif
        if (mod_propFlags(2)) then; rho   = tableValue(rhoTab, Tmin, Tmax, tp); else; rho   = auxLocal(ind_rho); endif
        if (mod_propFlags(3)) then; sigma = tableValue(sigTab, Tmin, Tmax, tp); else; sigma = auxLocal(ind_sig); endif
        if (mod_propFlags(4)) then; mup   = tableValue(mupTab, Tmin, Tmax, tp); else; mup   = auxLocal(ind_mup); endif
        select case(mod_model)
        case(3);      m = auxLocal(ind_m)/y(8)
        case(7);      m = auxLocal(ind_m)/y(9)
        case(2,4,5);  m = y(8)
        case(8);      m = y(8) + y(9)
        case default; m = auxLocal(ind_m)
        end select
        !> Models 3, 7 and a rho-solid material derive d from the current ODE state.
        if (ind_d > 0 .and. mod_model /= 3 .and. mod_model /= 7 .and. ind_rs == 0) then
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
          case(3);   eventLocal(ind_evn) = y(8)
          case(4,7); eventLocal(ind_evn) = y(9)
          end select
        endif
        if (ind_sb1 > 0) brkupState(1:nbrkst) = stateLocal(ind_sb1:ind_sb2)
        !> Breakup event over the accepted-step interval (x-xold).
        vEvt = Y(4:6)
        call breakupEvent(eventLocal, neventvar, brkupState, nbrkst,    &
                          sigma,mup,rho,gasState(1),vel,Re,t1,acc,vEvt,x-xold, &
                          mod_brkSelect, mod_bp,mod_bpMethod,mod_bpScale,        &
                          eventFlag,childState,addChildLocal,exitLoop,part%rngState,childDone)
        !> The ETAB kick as a delta, applied once by the segment finalize.
        kickDV = vEvt - Y(4:6)
        kickPend = (mod_brkSelect == 5) .and. (maxval(abs(kickDV)) > 0._R8)
        if (ind_sb1 > 0) stateLocal(ind_sb1:ind_sb2) = brkupState(1:nbrkst)
      endif

      if (IamOut .or. newGas .or. eventFlag .or. burnedOut .or. sectorOut .or. solidEv) then
        IRTRN = -2724
        return
      endif

      !> Never latch a non-finite state as the last good one.
      if (any(y/=y)) return

      !> Planes and scatter markers the accepted step crosses, emitted once as it is latched.
      if (planeOn) call emitPlanes(timeLocal, oldLocal, x, y)
      if (scatOn .and. dNscat>0._R8) call emitScatter(timeLocal, oldLocal, x, y)

      timeLocal  = x
      oldLocal   = y
      if (allocated(stateLocal)) oldStLocal = stateLocal
      if (allocated(eventLocal)) oldEvLocal = eventLocal

      !> Latched with the state, after the emission: a refinement never writes a marker twice.
      wAccOld = wAcc

    end subroutine solout


    !> Cell-entry geometry and gas state of the current cell.
    subroutine refreshCellEntry()
      call geoblock(b)%getVertices(part%i(2:4),vert, norms=geoHexNorms, centroids=geoHexCentroids, degen=geoHexDegen, &
          tol=geoHexTol)   ! geo cell (NOT the gas dual index igas)
      if (ord2) then
        call gasblock(b)%getVertices(part%igas,gasVert,mesh2D, gasHexNorms, gasHexCentroids, gasHexDegen, gasHexTol)
        call setAtGasBoundary()
      endif
      call gasblock(b)%gasProperties(gas,part%igas)
    end subroutine refreshCellEntry

    !> A pinned parcel whose point lies in its injection cell by the test that pinned it (the 12-plane
    !  test of the cell, also on a one-layer mesh) takes that cell, with no search.
    logical function injCellHolds()
      use IGLOO_RayFaceIntersection3D, only: isPointInsideHexahedron12
      real(R8) :: v(3,8)
      injCellHolds = .false.
      if (any(part%iInj < 1)) return
      call geoblock(part%iInj(1))%getVertices(part%iInj(2:4), v)
      if (.not. isPointInsideHexahedron12(part%pInj, v)) return
      part%i = part%iInj; part%wasin = .true.; part%lost = .false.
      injCellHolds = .true.
    end function injCellHolds

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
      call geoblock(b)%getVertices(part%i(2:4),vert, norms=geoHexNorms, centroids=geoHexCentroids, degen=geoHexDegen, tol=geoHexTol)
      if (.not.part%gone) then
        part%gasExitFace = 0
        call part%advanceGasCell(gasblock(b))
        call gasblock(b)%getVertices(part%igas,gasVert,mesh2D, gasHexNorms, gasHexCentroids, gasHexDegen, gasHexTol)
        call gasblock(b)%gasProperties(gas,part%igas)
        call setAtGasBoundary()
      endif
    end subroutine handleGeoEvent


    !> The wall impact of a parcel that left through an eroding wall cell: its flow npdot*m and its eroded flow, the flow
    !  times its erosion ratio at its speed and angle from the surface at the crossing point, join the cell's sums.
    subroutine depositErosion()
      integer  :: wb, wf, wm, wn, wmat
      real(R8) :: mdotW, ermW, vmag, nmag, sinA, nn(3)

      wb = part%exitCell(1); wf = part%exitCell(2); wm = part%exitCell(3); wn = part%exitCell(4)
      if (wb < 1 .or. wf < 1) return
      if (.not. wallblock(wb)%on(wf)) return
      if (.not. isWallCode(geoblock(wb)%face(wf)%cell(wm,wn))) return
      wmat  = part%mID
      mdotW = part%npdot*part%m
      vmag  = norm2(part%vImp)
      nn    = geoblock(wb)%face(wf)%cell(wm,wn)%normal
      nmag  = norm2(nn)
      sinA  = 0._R8
      if (vmag > 0._R8 .and. nmag > 0._R8) sinA = min(1._R8, abs(dot_product(part%vImp, nn))/(vmag*nmag))
      ermW  = mdotW*erosionRatio(wmat, vmag, asin(sinA), part%d, part%rho)
      !$OMP ATOMIC UPDATE
      wallblock(wb)%face(wf)%acc(1,wmat,wm,wn) = wallblock(wb)%face(wf)%acc(1,wmat,wm,wn) + ermW
      !$OMP ATOMIC UPDATE
      wallblock(wb)%face(wf)%acc(2,wmat,wm,wn) = wallblock(wb)%face(wf)%acc(2,wmat,wm,wn) + mdotW
    end subroutine depositErosion



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
      if (phaseChange .or. particle%model==5 .or. particle%model==8) then
        !$OMP ATOMIC UPDATE
        sblock%sourceMass(nmat,ii,jj,kk) = sblock%sourceMass(nmat,ii,jj,kk) + (masIn-masOut)
      else
        !> Every rank stores the same 0 (fields pre-zeroed); a non-zero constant would be summed over the ranks.
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
      case(3,7)
        np  = particle%intE(5)            *factor
        rho = mdotSeg*particle%Tstay*factor
        factor = rho/particle%deltaL
      case(2,5,8)
        np  = particle%npdot*particle%Tstay*factor
        rho = particle%npdot*particle%intE(5)*factor   !> intE(5) = int m dt per droplet; times the number rate
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
