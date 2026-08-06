module IGLOO_module
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  use IGLOO_data_block, only: obj_block, obj_flowblock, obj_sourceblock, obj_eulerblock
  use IGLOO_data_phases, only: obj_material
  implicit none

  type, public :: obj_IGLOO
    type(obj_block),       allocatable :: geoblock(:)
    type(obj_flowblock),   allocatable :: gasblock(:)
    type(obj_sourceblock), allocatable :: source(:)
    type(obj_eulerblock),  allocatable :: euler(:,:)
    type(obj_material),    allocatable :: material(:)
    logical :: eulSwitch, srcSwitch
    ! integer            :: iprint
    ! real(R8)           :: ds, mdotMax, dtprint
  contains
    procedure, pass(self) :: setup
    procedure, pass(self) :: solve
    procedure, pass(self) :: getSourceTerms
    procedure, pass(self) :: writeout
  end type obj_IGLOO

contains

  subroutine setup(self, external_gas)
    use omp_lib
    use IGLOO_IO
    use IGLOO_IC
    use IGLOO_particles
    use IGLOO_variables
    use IGLOO_allocation
    use IGLOO_IO_INI, only: read_IGLOO_input
    use IGLOO_Lib_Statistics, only: initRandomSeed
    use Lib_ORION_data
    implicit none
    class(obj_IGLOO), intent(inout)        :: self
    type(orion_data), intent(in), optional :: external_gas
    type(orion_data)     :: own_gas
    character(len=llen)  :: gasfile
    real(8), allocatable :: pos0(:,:), vel0(:,:), mdot(:), diam(:), temp0(:)
    integer              :: m, g, fam, nthreads=1
    character(len=2)     :: method

    call print_header()

# if defined (_OPENMP)
    !$OMP PARALLEL
    nthreads = OMP_GET_NUM_THREADS()
    !$OMP END PARALLEL
# endif
    if (nthreads>1) then
      write(*,*)" OpenMP threads = ", nthreads
    endif

    call read_IGLOO_input(method,self%srcSwitch,self%eulSwitch,gasfile,pos0,vel0,temp0,mdot,diam)
    if (present(external_gas)) then
      call copyORION(external_gas,own_gas)
    else
      write(*,*)' >> Background flow field => ',trim(gasfile)
      call read_TECsolfile(gasfile,own_gas)
    endif

    call read_cdp_properties(IGLOO_phase_prefix,self%material)
    call allocate_blocks(own_gas,self%material,self%geoblock,self%gasblock, &
                         self%source,self%euler,self%srcSwitch,self%eulSwitch)
    !> 2D mesh: an out-of-plane body-acceleration component is unphysical (mesh2D is set in
    !  allocate_blocks). Zero it and refresh the gating flags. (read_general set them already.)
    if (mesh2D .and. bodyForce .and. bodyAccel(3) /= 0._R8) then
      write(*,*) ' [WARNING] 2D mesh: zeroing out-of-plane body-accel z-component'
      bodyAccel(3) = 0._R8
      bodyForce    = any(bodyAccel /= 0._R8)
      srcBodyForce = sourceSwitch .and. bodyForce
    endif
    call import_gas(own_gas,self%gasblock)
    call read_cdp_bc_file(IGLOO_phase_prefix,self%material,self%geoblock, &
                          self%gasblock,self%source,self%euler,self%srcSwitch,self%eulSwitch)
    call initRandomSeed(rng_seed)   ! before pinning: stochastic injection diameters
    fam = 0
    do m = 1, nm
      self%material(m)%mID = m
      allocate(self%material(m)%group(1:self%material(m)%ngroups))
      call self%material(m)%assign_material2group()
      do g = 1, self%material(m)%ngroups
        fam = fam + 1
        associate(gr => self%material(m)%group(g))
        gr%gID   = g
        gr%famID = fam
        call pin_particles(gr, self%geoblock, self%gasblock, method, pos0,vel0,temp0,mdot,diam, fam, gr%nparticles)
        call gr%assign_group2particle()
        end associate
      enddo
    enddo
    nfam = fam
  end subroutine setup


  subroutine solve(self)
    use Lib_Integration,  only: integrate
    use IGLOO_particles,  only: obj_particle
    use IGLOO_variables,  only: ode_word,iopt,rtol,atol, ord2, nm, nb,     &
                                eulerSwitch, sourceSwitch, unitTraj, unitScat, &
                                unitExit, IGLOO_phase_prefix, srcBodyForce,    &
                                trajOn, scatOn, dNscat, trajSample, threshold, &
                                probeOn, probeIDs
    use IGLOO_IC,         only: initialize_fields
    use IGLOO_allocation, only: allocateAccumulators
    use IGLOO_Lib_Properties, only: lookupTab   !> A23c child hand-off (enthalpy slot)
    use oslo
    use omp_lib
    implicit none
    class(obj_IGLOO), intent(inout) :: self
    integer, parameter :: maxLoop=5
    integer :: m, g, ip, iota, ch, e, start, nEnd, oldStart, oldEnd, newSize, b, fam
    integer :: loopCounter, childCounter, maxShedSeen
    ! logical, allocatable :: famDone(:)
    real(R8), allocatable :: relTol(:), absTol(:)
    !> Scatter-cloud weight-quantum (dNscat) auto-sizing scratch.
    integer  :: c, nStreams, nValid
    real(R8) :: Lref, xmin(3), xmax(3), Ndot, Vsum, Vref, tauRef, vsp
    type(obj_particle) :: ptmp   ! throwaway copy for the inject-only pre-pass

    sourceSwitch = self%srcSwitch
    eulerSwitch  = self%eulSwitch
    !> Allocate the source/euler accumulators (sized per ord2). Idempotent:
    !  re-running solve reuses the existing allocation.
    call allocateAccumulators(self%source, self%euler, sourceSwitch, eulerSwitch)
    call initialize_fields(self%source,self%euler,sourceSwitch,eulerSwitch)
    if (ord2) then
      do m = 1, nb
        call self%gasblock(m)%fillGhostGradient()
      enddo
    endif

    !> Domain length scale (bounding-box diagonal) for the scatter weight-quantum estimate.
    !  The reference transit time tauRef = Lref/Vref is finished per material below, using the
    !  resolved particle injection speed Vref (more representative than the gas mean).
    xmin =  huge(1._R8); xmax = -huge(1._R8)
    do b = 1, nb
      do c = 1, 3
        xmin(c) = min(xmin(c), minval(self%geoblock(b)%node(c,:,:,:)))
        xmax(c) = max(xmax(c), maxval(self%geoblock(b)%node(c,:,:,:)))
      enddo
    enddo
    Lref = norm2(xmax - xmin)

    do m = 1, nm
      material: associate(mat => self%material(m))
      if (trajOn) then
        open(newunit=unitTraj,file='OUTPUT/'//trim(IGLOO_phase_prefix)//'trajectories-'//trim(mat%matName)//'.dat')
        write(unitTraj,*) 'variables="X","Y","Z","U","V","W","T","d<sub>p","m<sub>p","ID"'
      endif
      open(newunit=unitExit,file='OUTPUT/'//trim(IGLOO_phase_prefix)//'outloc-'//trim(mat%matName)//'.dat')
      write(unitExit,*) 'variables="X","Y","Z","T","|u<sub>p</sub>|","<greek>a</greek>","mdot","Af","ID"'

      !> Scatter cloud: one flat point-cloud zone per material; auto-size the weight quantum
      !  dNscat (real droplets/point) so the cloud holds ~trajSample points per stream. The
      !  population estimate (Ndot*tauRef) is crude; it sets only the count, not the shape.
      dNscat = 0._R8
      if (scatOn) then
        open(newunit=unitScat,file='OUTPUT/'//trim(IGLOO_phase_prefix)//'scatter-'//trim(mat%matName)//'.dat')
        write(unitScat,*) 'variables="X","Y","Z","U","V","W","T","d<sub>p","m<sub>p","ID"'
        write(unitScat,'(A,A,A)')'Zone T="Mat ',trim(mat%matName),' scatter"'
        !> Serial inject-only pre-pass: resolve each stream's npdot (= Σ droplet rate Ndot) and
        !  injection speed (=> Vref) without the cell search. Serial => deterministic dNscat.
        Ndot = 0._R8; nStreams = 0; Vsum = 0._R8; nValid = 0
        do g = 1, mat%ngroups
          associate(gr => mat%group(g))
          call gr%setup_particleODE()
          do ip = 1, gr%nparticles
            !> Throwaway copy: resolveInjectionRate writes mdot/v/tp on it; keep the real
            !  particle pristine so the real sweep is bit-unchanged. It deliberately avoids the
            !  geometric cell search (which would pollute the threadprivate myRay).
            ptmp = gr%particle(ip)
            call ptmp%resolveInjectionRate(self%geoblock, self%gasblock)
            if (ptmp%npdot > 0._R8) Ndot = Ndot + ptmp%npdot
            vsp = norm2(ptmp%stateVar(4:6))
            if (vsp < threshold) then; Vsum = Vsum + vsp; nValid = nValid + 1; endif
          enddo
          nStreams = nStreams + gr%nparticles
          end associate
        enddo
        Vref   = 1._R8; if (nValid>0 .and. Vsum>0._R8) Vref = Vsum/real(nValid,R8)
        tauRef = Lref / max(Vref, 1.e-30_R8)           ! domain transit time at injection speed
        if (Ndot>0._R8 .and. nStreams>0) &
          dNscat = Ndot*tauRef / real(max(trajSample*nStreams,1), R8)
        if (.not. (dNscat>0._R8)) dNscat = 0._R8       ! guard NaN/neg => scatter silently skipped
        write(*,'(A,A,A,ES10.3,A,I0,A)') '     >> scatter cloud [',trim(mat%matName),       &
              ']: dNscat=',dNscat,' droplets/pt (nominal ~',trajSample,                      &
              '/stream; realized varies with residence — raise fsample-traj for a denser cloud)'
      endif
      
      write(*,*)" Compute particles dynamics for material: ",trim(mat%matName)
      if (mat%cpVariable) write(*,*) ' >> Solving enthalpy equation'
      do g = 1, mat%ngroups
        group: associate(gr => mat%group(g))
        call gr%setup_particleODE()
        if (allocated(relTol)) deallocate(relTol); allocate(relTol(gr%neq)); relTol(:) = rtol
        if (allocated(absTol)) deallocate(absTol); allocate(absTol(gr%neq)); absTol(:) = atol
        !> Exclude the body-force accumulators from the error norm (pure output quadratures): only W
        !  (neq) when euler-on, both J,W when euler-off. Mass-evolving models 2/4/5 only.
        if (srcBodyForce .and. (gr%evapSelect > 0 .or. gr%combSelect > 0)) then
          if (eulerSwitch) then       ! neq-1 is a real euler moment, keep it controlled
            relTol(gr%neq)          = 1.e30_R8; absTol(gr%neq)          = 1.e30_R8
          else                        ! J,W both appended
            relTol(gr%neq-1:gr%neq) = 1.e30_R8; absTol(gr%neq-1:gr%neq) = 1.e30_R8
          endif
        endif
        call setup_odesolver(N=gr%neq,solver=ode_word,RT=relTol,AT=absTol,iopt=iopt)

        write(*,'(A,I3,A,I6)')"     Group",g," => number of particles = ", gr%nparticles
        write(*,'(A,I2,A)') '     >> ODE system (neq=', gr%particle(1)%neq, '):'
        if (gr%evapSelect/=0) write(*,*) '       - evaporation --> ', trim(gr%evapWord)
        if (gr%combSelect/=0) write(*,*) '       - combustion  --> Beckstead d^n burn law'
        if (gr%brkupEqOde   ) write(*,*) '       - breakup     --> ', trim(gr%brkupWord)
        if (gr%evapSelect==0 .and. gr%combSelect==0 .and. .not.gr%brkupEqOde) &
          write(*,*) '       - constant particle mass and size'
        if (eulerSwitch) write(*,*) '       - eulerian field included'

        if (trajOn) write(unitTraj,'(A,A,A,I3,A)')'Zone T="Mat ',trim(mat%matName),' Group',g,'"'
        write(unitExit ,'(A,A,A,I3,A)')'Zone T="Mat ',trim(mat%matName),' Group',g,'"'

        if (gr%brkupHasChild) then
          if (allocated(gr%shed)) deallocate(gr%shed)
          allocate(gr%shed(1:gr%nparticles))
          gr%nactive = gr%nparticles
          start = 1;  nEnd = gr%nactive
          loopCounter = 0

          do while (loopCounter < maxLoop)
            loopCounter = loopCounter + 1
            !> Grow BEFORE resetting the window -- the other order indexes past the end on the
            !  pass that added children. move_alloc (not deallocate/allocate) so each list keeps
            !  the capacity it reserved on an earlier pass.
            if (size(gr%shed) < nEnd) call resizeShedArray(gr%shed, nEnd)
            !> Empty each list in the window and give it room for the common case, so nothing
            !  has to allocate inside the parallel region below.
            do ip = start, nEnd
              gr%shed(ip)%n = 0
              call gr%shed(ip)%reserve(1)
            enddo

            !$OMP PARALLEL DO SCHEDULE(DYNAMIC)
            do ip = start, nEnd
              if (probeOn) then; if (.not.any(gr%particle(ip)%ID==probeIDs)) cycle; endif
              call integrate(gr%particle(ip),self%geoblock,self%gasblock,self%source,self%euler(:,gr%famID), &
                             mat%hTab,mat%cpTab,mat%rhoTab,mat%mupTab,mat%sigTab,mat%psatTab,    &
                             gr%shed(ip), noShed=(loopCounter==maxLoop))
            enddo
            !$OMP END PARALLEL DO

            !> Total shed EVENTS, not parents that shed -- `newSize` has to size the particle
            !  array by the number of children actually created. Summed serially: the drain
            !  below walks the same ascending `ip`, so the child ordering is fixed by the
            !  traversal itself, independent of thread count and SCHEDULE(DYNAMIC).
            oldStart     = start
            childCounter = 0
            maxShedSeen  = 0
            do ip = start, nEnd
              childCounter = childCounter + gr%shed(ip)%n
              maxShedSeen  = max(maxShedSeen, gr%shed(ip)%n)
            enddo

            if (childCounter > 0) then
              write(*,*)"       Loop",loopCounter," => number of children = ", childCounter
              !> Visible tripwire: how hard the shed path is actually working. A number
              !  climbing toward the `maxShed` guard in Lib_Integration means the case is
              !  approaching runaway shedding well before the guard has to fire.
              if (maxShedSeen > 1) write(*,'(A,I0,A)')                                    &
                    "                             (max ",maxShedSeen," sheds from one parcel)"

              !> Grow particle array if needed (geometric 2x)
              oldEnd  = nEnd
              newSize = gr%nactive + childCounter
              if (newSize > size(gr%particle)) then
                call resizeParticleArray(gr%particle, max(2*size(gr%particle), newSize))
              endif

              start = gr%nactive + 1
              nEnd  = newSize
              call gr%setup_particleODE(start, nEnd)
              call gr%assign_group2particle(start, nEnd)

              !> Drain the shed lists in ascending parent index, then push order within each.
              !  This is the same traversal the old single-slot compaction performed, so with
              !  the one-shed cap on it visits exactly the same records in the same order --
              !  which is what keeps `kid%ID = oldEnd + ch` landing on the same parent.
              !  Deterministic by construction, not by luck: no sort, no thread-number
              !  dependence, unaffected by SCHEDULE(DYNAMIC).
              !
              !> A23c: a child inherits time > 0, so `integrate` SKIPS its `part%time==0`
              !  initialization block entirely. Everything that block would have established
              !  must therefore be set here, or the child enters the solver with an
              !  uninitialized cell index and mass state (which segfaulted on the first
              !  getVertices). This is the child hand-off that was never written -- the path
              !  had no trigger before A23a/A23b, so it had never executed for any model.
              ch = 0
              do ip = oldStart, oldEnd
              do e  = 1, gr%shed(ip)%n
                ch   = ch + 1
                iota = oldEnd + ch
                associate(kid => gr%particle(iota), src => gr%shed(ip)%item(e))
                kid%ID            = iota
                kid%stateVar(1:3) = src%pos
                kid%stateVar(4:6) = src%vel
                kid%tp            = src%temp
                kid%d             = src%diam
                kid%npdot         = src%npdot
                kid%time          = src%time
                !> cell/gas locators inherited from the parent at the shed point
                kid%i       = src%ipos
                kid%igas    = src%igas
                kid%igasOld = src%igas
                kid%xi0     = 0.5_R8
                !> fresh trajectory bookkeeping
                kid%Ncell = 0; kid%gone = .false.; kid%wasin = .false.; kid%angle = 0._R8
                kid%lost  = .false.
                !> parcel mass state: m from (rho, d), then the parcel flow mdot = npdot*m
                call kid%computeMass(mat%rhoTab)
                kid%m0   = kid%m
                kid%mdot = kid%npdot * kid%m
                !> ODE state: enthalpy/temperature slot, then the model-3 count slot
                if (kid%varCp) then
                  kid%stateVar(7) = lookupTab(mat%hTab, kid%tp)
                else
                  kid%stateVar(7) = kid%tp
                endif
                select case (kid%model)
                case(2,5); kid%stateVar(8) = kid%m
                case(3);   kid%stateVar(8) = kid%npdot
                case(4);   kid%stateVar(8) = kid%m; kid%stateVar(9) = kid%npdot
                end select
                if (eulerSwitch)    kid%stateVar(kid%nOde+1:kid%neq) = 0._R8
                if (kid%bodyAccum)  kid%stateVar(kid%neq-1 :kid%neq) = 0._R8
                kid%npold   = kid%npdot
                kid%oldState = kid%stateVar
                end associate
              enddo
              enddo
              gr%nactive = newSize
            else
              exit
            endif
          enddo

          !> Trim excess capacity + sync nparticles
          if (size(gr%particle) > gr%nactive) then
            call resizeParticleArray(gr%particle, gr%nactive)
          endif
          gr%nparticles = gr%nactive
        else
          !$OMP PARALLEL DO SCHEDULE(DYNAMIC)
          do ip = 1,gr%nparticles
            if (probeOn) then; if (.not.any(gr%particle(ip)%ID==probeIDs)) cycle; endif
            call integrate(gr%particle(ip),self%geoblock,self%gasblock,self%source,self%euler(:,gr%famID), &
                            mat%hTab,mat%cpTab,mat%rhoTab,mat%mupTab,mat%sigTab,mat%psatTab)
          enddo
          !$OMP END PARALLEL DO
        endif
        end associate group
      enddo
      end associate material
      if (trajOn) close(unitTraj)
      if (scatOn) close(unitScat)
      close(unitExit)
    enddo

    !> End-of-solve finalization.
    !  - obj_eulerblock%finalize normalizes +-accumulated numerators into
    !    weighted averages (and inverts h→T for cpVariable groups), then if
    !    ord2=true reduces gasblock-shape arrays to geoblock-shape via the
    !    sub-octant volume weighting.
    !  - obj_sourceblock%finalize is a no-op when ord2=false; when ord2=true
    !    it performs the same reduction on the source-mass/mom/en arrays.
    if (sourceSwitch) then
      do b = 1, nb
        call self%source(b)%finalize(self%geoblock(b))
      enddo
    endif
    if (eulerSwitch) then
      ! allocate(famDone(size(self%euler, 2)))
      ! famDone = .false.
      do m = 1, nm
        do g = 1, self%material(m)%ngroups
          fam = self%material(m)%group(g)%famID
          ! if (.not. famDone(fam)) then
          do b = 1, nb
            call self%euler(b, fam)%finalize(self%geoblock(b), &
                                              self%material(m)%hTab, &
                                              self%material(m)%cpVariable)
          enddo
            ! famDone(fam) = .true.
          ! endif
        enddo
      enddo
      ! deallocate(famDone)
    endif

    write(*,*)" Stop condition : All particles out of domain!"

  end subroutine solve


  pure function getSourceTerms(self,gas,vel,Tp,rhop,np,mID) result(FdragQdot)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    use IGLOO_Lib_Properties, only: lookupTab
    use IGLOO_variables,      only: pi, dragSelect, heatSelect
    use IGLOO_Lib_Drag
    use IGLOO_Lib_Heat
    implicit none
    class(obj_IGLOO), intent(in) :: self
    integer,          intent(in) :: mID
    real(R8),         intent(in) :: gas(9), vel(3), Tp, np, rhop
    real(R8) :: FdragQdot(4), vdiff(3), Fdrag(3)
    real(R8) :: rhog, Tg, gamma, mug, kg, Rg
    real(R8) :: Re, Ma, Tr, Pr, Nu, Cd, slip, Qdot, d, rhom
    real(R8), parameter :: oneThird =0.33333333333333333_R8, &
                           sixOverPi=1.90985931710274403_R8, &
                           piOver8  =0.39269908169872415_R8
    
    if (self%material(mID)%rhoVariable) then
      rhom = lookupTab(self%material(mID)%rhoTab,Tp)
    else
      rhom = self%material(mID)%rho
    endif
    rhog  = gas(1)
    vdiff = gas(2:4)-vel
    Tg = gas(5); gamma = gas(6); Rg = gas(7); mug = gas(8); kg = gas(9)
    slip  = norm2(vdiff)
    d  = (sixOverPi*rhop/(rhom*np+1e-20))**oneThird
    Re = rhog*slip*d/mug
    Ma = slip/sqrt(gamma*Rg*Tg)
    Tr = Tp/Tg
    Pr = mug*gamma*Rg/((gamma-1)*kg)
    Cd = drag(Re,Ma,gamma,Tr,dragSelect)
    Nu = heat(Re,Pr,Ma,heatSelect)
    Fdrag = piOver8*Cd*d*Re*mug*vdiff
    Qdot  = Nu*kg*pi*d*(Tg-Tp)
    Qdot  = Qdot + Fdrag(1)*vel(1)+Fdrag(2)*vel(2)+Fdrag(3)*vel(3)
    FdragQdot = [Fdrag, Qdot]*np

  end function getSourceTerms


  subroutine writeout(self)
    use IGLOO_IO
    implicit none
    class(obj_IGLOO), intent(inout) :: self

    call write_outfield(self%material,self%geoblock,self%source,self%euler,self%srcSwitch,self%eulSwitch)

  end subroutine writeout


  subroutine print_header()
    write(*,*)
    write(*,*) ' ============================================================================= '
    write(*,*) '|                  ///    /////    //       /////    /////                    |'
    write(*,*) '|                  ///   //        //      //  //   //  //                    |'
    write(*,*) '|                  ///   //  ///   //      //  //   //  //                    |'
    write(*,*) '|                  ///   //   //   //      //  //   //  //                    |'
    write(*,*) '|                  ///   //////    //////  /////    /////                     |'
    write(*,*) '|-----------------------------------------------------------------------------|'
    write(*,*) '|            Integration of a General Lagrangian One-Way ODE set              |'
    write(*,*) ' ============================================================================= '
    write(*,*)
  end subroutine print_header


  subroutine resizeParticleArray(arr, newCapacity)
    use IGLOO_particles, only: obj_particle
    implicit none
    type(obj_particle), allocatable, intent(inout) :: arr(:)
    integer, intent(in) :: newCapacity
    type(obj_particle), allocatable :: tmp(:)
    integer :: n
    n = min(size(arr), newCapacity)
    allocate(tmp(newCapacity))
    tmp(1:n) = arr(1:n)
    call move_alloc(tmp, arr)
  end subroutine resizeParticleArray


  !> Same idiom for the per-parent shed lists. Intrinsic assignment deep-copies each list's
  !> allocatable `item` component, so previously-reserved capacity survives the growth.
  subroutine resizeShedArray(arr, newCapacity)
    use IGLOO_data_phases, only: shedList
    implicit none
    type(shedList), allocatable, intent(inout) :: arr(:)
    integer, intent(in) :: newCapacity
    type(shedList), allocatable :: tmp(:)
    integer :: n
    n = min(size(arr), newCapacity)
    allocate(tmp(newCapacity))
    tmp(1:n) = arr(1:n)
    call move_alloc(tmp, arr)
  end subroutine resizeShedArray


end module IGLOO_module
