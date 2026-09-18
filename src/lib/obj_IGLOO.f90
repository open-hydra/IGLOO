!> Top-level IGLOO object: geometry, gas, feedback fields and materials, with the
!  setup / reset_state / solve / writeout entry points.
module IGLOO_module
  use, intrinsic :: iso_fortran_env, only : R8 => real64
  use IGLOO_data_block, only: obj_block, obj_flowblock, obj_sourceblock, obj_eulerblock
  use IGLOO_data_phases, only: obj_material
  implicit none

  type, public :: obj_IGLOO
    type(obj_block),       allocatable :: geoblock(:)
    type(obj_flowblock),   allocatable :: gasblock(:)
    !> Valid only between a solve and the next reset_state (which deallocates them).
    type(obj_sourceblock), allocatable :: source(:)
    type(obj_eulerblock),  allocatable :: euler(:,:)
    type(obj_material),    allocatable :: material(:)
    logical :: eulSwitch, srcSwitch
    !> Sweep bookkeeping: setup_static ran; reset_state ran since the last solve.
    logical :: staticDone    = .false.
    logical :: stateIsFresh  = .false.
    !> Output generation: -1 before the first reset_state, 0 on the first sweep (untagged files).
    integer :: sweep = -1
  contains
    procedure, pass(self) :: setup
    procedure, pass(self) :: setup_static
    procedure, pass(self) :: reset_state
    procedure, pass(self) :: sweepTag
    procedure, pass(self) :: solve
    procedure, pass(self) :: getSourceTerms
    procedure, pass(self) :: writeout
  end type obj_IGLOO

contains

  !> Single-sweep entry point: setup_static followed by reset_state.
  subroutine setup(self, external_gas)
    use Lib_ORION_data
    implicit none
    class(obj_IGLOO), intent(inout)        :: self
    type(orion_data), intent(in), optional :: external_gas

    call self%setup_static(external_gas)
    call self%reset_state()

  end subroutine setup


  !> Once-only setup: input parsing, gas import (`external_gas` or the Tecplot file), block and
  !  geometry allocation, BC tagging and particle pinning; replicated identically on every rank.
  subroutine setup_static(self, external_gas)
    use omp_lib
    use IGLOO_IO
    use IGLOO_IC
    use IGLOO_particles
    use IGLOO_variables
    use IGLOO_allocation
    use IGLOO_IO_INI, only: read_IGLOO_input
    use IGLOO_Lib_Statistics, only: initRandomSeed
    use IGLOO_Mod_MPI,        only: mpi_is_root, mpi_size_
    use Lib_ORION_data
    implicit none
    class(obj_IGLOO), intent(inout)        :: self
    type(orion_data), intent(in), optional :: external_gas
    type(orion_data)     :: own_gas
    character(len=llen)  :: gasfile
    real(8), allocatable :: pos0(:,:), vel0(:,:), mdot(:), diam(:), temp0(:)
    integer              :: m, g, fam, nthreads
    character(len=2)     :: method

    if (self%staticDone) then
      write(*,*) ' [IGLOO] setup_static already completed -- repeat call ignored'
      write(*,*) '         (mesh, gas import and pinning are static; use reset_state per sweep)'
      return
    endif

    !> Informational prints are root-gated; contract violations print on every rank.
    if (mpi_is_root) call print_header()

    nthreads = 1
# if defined (_OPENMP)
    nthreads = OMP_GET_MAX_THREADS()   !> outside the region: no shared write to race on
# endif
    if (nthreads>1 .and. mpi_is_root) then
      write(*,*)" OpenMP threads = ", nthreads
    endif
    !> Rank count, printed only when actually decomposed.
    if (mpi_size_>1 .and. mpi_is_root) then
      write(*,'(A,I0)')" MPI ranks = ", mpi_size_
    endif

    call read_IGLOO_input(method,self%srcSwitch,self%eulSwitch,gasfile,pos0,vel0,temp0,mdot,diam)
    if (present(external_gas)) then
      call copyORION(external_gas,own_gas)
    else
      if (mpi_is_root) write(*,*)' >> Background flow field => ',trim(gasfile)
      call read_TECsolfile(gasfile,own_gas)
    endif

    call read_cdp_properties(IGLOO_phase_prefix,self%material)
    call allocate_blocks(own_gas,self%material,self%geoblock,self%gasblock, &
                         self%source,self%euler,self%srcSwitch,self%eulSwitch)
    !> 2D mesh: drop the out-of-plane body-acceleration component.
    if (mesh2D .and. bodyForce .and. bodyAccel(3) /= 0._R8) then
      write(*,*) ' [WARNING] 2D mesh: zeroing out-of-plane body-accel z-component'
      bodyAccel(3) = 0._R8
      bodyForce    = any(bodyAccel /= 0._R8)
      srcBodyForce = sourceSwitch .and. bodyForce
    endif
    call import_gas(own_gas,self%gasblock)
    !> 2.5D notice: DB injection with an out-of-plane velocity.
    if (mesh2D .and. mpi_is_root) then
      if (allocated(vel0)) then
        if (any(vel0(:,3) /= 0._R8 .and. abs(vel0(:,3)) < threshold)) &
          write(*,'(A)') '  >> DB injection carries wp /= 0 (2.5D)'
      endif
    endif
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
        !> Pin-time census, never overwritten.
        gr%nInjected = gr%nparticles
        end associate
      enddo
    enddo
    nfam = fam
    self%staticDone = .true.

  end subroutine setup_static


  !> Per-sweep state restore: puts the pinned population back to its injection state and
  !  re-imports the gas when `external_gas` is given. Must run before every `solve`.
  subroutine reset_state(self, external_gas)
    use IGLOO_variables,      only: nm, nb
    use IGLOO_allocation,     only: import_gas
    use IGLOO_Lib_Statistics, only: rngSeedFor
    use Lib_ORION_data
    implicit none
    class(obj_IGLOO), intent(inout)        :: self
    type(orion_data), intent(in), optional :: external_gas
    type(orion_data) :: own_gas
    integer :: m, g, ip, b, fam

    !> Next output generation (shared by the solve/writeout pair).
    self%sweep = self%sweep + 1

    !> Re-import the gas field values (the mesh is static); cell%mdotGas is not refreshed.
    if (present(external_gas)) then
      call copyORION(external_gas, own_gas)
      call import_gas(own_gas, self%gasblock)
    endif

    !> Drop last sweep's source/euler accumulators (all of them); solve re-allocates them at the right shape.
    do b = 1, nb
      if (allocated(self%source(b)%sourceMass)) deallocate(self%source(b)%sourceMass)
      if (allocated(self%source(b)%sourceMom))  deallocate(self%source(b)%sourceMom)
      if (allocated(self%source(b)%sourceEn))   deallocate(self%source(b)%sourceEn)
      do fam = 1, size(self%euler, 2)
        if (allocated(self%euler(b,fam)%density))     deallocate(self%euler(b,fam)%density)
        if (allocated(self%euler(b,fam)%velocity))    deallocate(self%euler(b,fam)%velocity)
        if (allocated(self%euler(b,fam)%temperature)) deallocate(self%euler(b,fam)%temperature)
        if (allocated(self%euler(b,fam)%np))          deallocate(self%euler(b,fam)%np)
      enddo
    enddo

    do m = 1, nm
      do g = 1, self%material(m)%ngroups
        associate(gr => self%material(m)%group(g))
        !> Counts: drop last sweep's children and restore the pin-time census.
        gr%nparticles = gr%nInjected
        gr%nactive    = gr%nInjected
        if (size(gr%particle) /= gr%nInjected) &
          call resizeParticleArray(gr%particle, gr%nInjected)
        !> Shed lists are re-allocated by solve.
        if (allocated(gr%shed)) deallocate(gr%shed)

        !> Restore each particle's injection state (same stochastic diameter draw every sweep).
        do ip = 1, gr%nInjected
          associate(part => gr%particle(ip))
          part%ID = ip
          !> Re-seed the particle's RNG stream from (rng_seed, famID, ID).
          part%rngState = rngSeedFor(gr%famID, ip)
          part%d  = part%dInj
          !> Assigned/DB streams only: BC streams re-derive tp and mdot at injection.
          if (all(part%iInj == [0,0,0,0])) then
            part%tp   = part%tpInj
            part%mdot = part%mdotInj
          endif
          !> time == 0 is integrate's "needs injection init" flag.
          part%time        = 0._R8
          part%lost        = .false.
          part%exitFace    = 0
          part%gasExitFace = 0
          part%Af          = 0._R8
          part%gone  = .false.
          part%wasin = .false.
          part%Ncell = 0
          part%angle = 0._R8
          part%i     = part%iInj
          part%iold  = [0,0,0,0]
          end associate
        enddo

        !> ODE sizing first, then group properties.
        call gr%setup_particleODE()
        call gr%assign_group2particle()
        end associate
      enddo
    enddo

    self%stateIsFresh = .true.

  end subroutine reset_state


  !> Filename tag of the current output generation: empty on sweep 0, '-sweep<N>' afterwards.
  function sweepTag(self) result(tag)
    implicit none
    class(obj_IGLOO), intent(in)  :: self
    character(len=:), allocatable :: tag
    character(len=16) :: buf

    if (self%sweep <= 0) then
      tag = ''
    else
      write(buf,'(A,I0)') '-sweep', self%sweep
      tag = trim(buf)
    endif

  end function sweepTag


  !> Integrate every particle of every material through the gas field: opens the per-material
  !  output files, sizes the scatter quantum, runs the OpenMP particle loop (with the child
  !  generation loop for shedding models), then reduces and finalizes the source/euler fields.
  subroutine solve(self)
    use Lib_Integration,  only: integrate
    use IGLOO_particles,  only: obj_particle
    use IGLOO_variables,  only: ode_word,iopt,rtol,atol, ord2, nm, nb,     &
                                eulerSwitch, sourceSwitch, unitTraj, unitScat, &
                                unitExit, IGLOO_phase_prefix, srcBodyForce,    &
                                trajOn, scatOn, dNscat, trajSample, threshold, &
                                probeOn, probeIDs, axisym, nSectorFold, nMultiFold
    use IGLOO_IC,         only: initialize_fields
    use IGLOO_allocation, only: allocateAccumulators
    use IGLOO_Lib_Properties, only: lookupTab   !> child hand-off (enthalpy slot)
    use IGLOO_Lib_Statistics, only: rngSeedFor
    use IGLOO_Mod_MPI,        only: mpi_size_, mpi_abort_all, owns_particle,           &
                                    mpi_allreduce_sum_i4_array, mpi_is_root, rank_suffix, &
                                    reduce_accumulators
    use oslo
    use omp_lib
    implicit none
    class(obj_IGLOO), intent(inout) :: self
    integer, parameter :: maxLoop=5
    integer :: m, g, ip, iota, ch, e, start, nEnd, oldStart, oldEnd, newSize, b, fam
    integer :: loopCounter, maxShedSeen
    !> Child-ID band of the current generation. No initializer here: it would imply SAVE.
    integer :: curBase, curCount, newBase, childLocal, childGlobal, pidx, idBase
    integer :: foldStat(2)
    integer, allocatable :: nShedAll(:), shedOff(:)
    ! logical, allocatable :: famDone(:)
    real(R8), allocatable :: relTol(:), absTol(:)
    !> Scatter-cloud weight-quantum (dNscat) auto-sizing scratch.
    integer  :: c, nStreams, nValid
    real(R8) :: Lref, xmin(3), xmax(3), Ndot, Vsum, Vref, tauRef, vsp
    type(obj_particle) :: ptmp   ! throwaway copy for the inject-only pre-pass

    !> Refuse to run on a stale particle state.
    if (.not. self%stateIsFresh) then
      write(*,'(A)') ' [IGLOO::solve] particle state is not fresh.'
      write(*,'(A)') '   Call reset_state() before each solve() (setup() does it for the first).'
      error stop 1
    endif
    self%stateIsFresh = .false.
    nSectorFold = 0; nMultiFold = 0

    sourceSwitch = self%srcSwitch
    eulerSwitch  = self%eulSwitch
    !> Allocate and zero the source/euler accumulators.
    call allocateAccumulators(self%source, self%euler, sourceSwitch, eulerSwitch)
    call initialize_fields(self%source,self%euler,sourceSwitch,eulerSwitch)
    if (ord2) then
      do m = 1, nb
        call self%gasblock(m)%fillGhostGradient()
      enddo
    endif

    !> Domain length scale (bounding-box diagonal) for the scatter weight-quantum estimate.
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
      !> Per-material output files: <kind>-<material><sweeptag><ranksuffix>.dat (rank suffix empty at one rank).
      if (trajOn) then
        open(newunit=unitTraj,file='OUTPUT/'//trim(IGLOO_phase_prefix)//'trajectories-'//trim(mat%matName)//self%sweepTag()//trim(rank_suffix())//'.dat')
        write(unitTraj,*) 'variables="X","Y","Z","U","V","W","T","d<sub>p","m<sub>p","ID"'
      endif
      open(newunit=unitExit,file='OUTPUT/'//trim(IGLOO_phase_prefix)//'outloc-'//trim(mat%matName)//self%sweepTag()//trim(rank_suffix())//'.dat')
      write(unitExit,*) 'variables="X","Y","Z","T","|u<sub>p</sub>|","<greek>a</greek>","mdot","Af","ID"'

      !> Scatter cloud: one zone per material; auto-size the weight quantum dNscat to ~trajSample points per stream.
      dNscat = 0._R8
      if (scatOn) then
        open(newunit=unitScat,file='OUTPUT/'//trim(IGLOO_phase_prefix)//'scatter-'//trim(mat%matName)//self%sweepTag()//trim(rank_suffix())//'.dat')
        write(unitScat,*) 'variables="X","Y","Z","U","V","W","T","d<sub>p","m<sub>p","ID"'
        write(unitScat,'(A,A,A)')'Zone T="Mat ',trim(mat%matName),' scatter"'
        !> Serial inject-only pre-pass over ALL particles (replicated on every rank so dNscat is
        !  identical everywhere): resolves each stream's npdot and injection speed.
        Ndot = 0._R8; nStreams = 0; Vsum = 0._R8; nValid = 0
        do g = 1, mat%ngroups
          associate(gr => mat%group(g))
          call gr%setup_particleODE()
          do ip = 1, gr%nparticles
            !> Throwaway copy keeps the real particle pristine.
            ptmp = gr%particle(ip)
            !> DB/assigned streams: seed the velocity as integrate does.
            if (all(ptmp%iInj == [0,0,0,0])) ptmp%stateVar(4:6) = ptmp%vInj
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
        if (mpi_is_root) &
          write(*,'(A,A,A,ES10.3,A,I0,A)') '     >> scatter cloud [',trim(mat%matName),     &
              ']: dNscat=',dNscat,' droplets/pt (nominal ~',trajSample,                      &
              '/stream; realized varies with residence — raise fsample-traj for a denser cloud)'
      endif

      if (mpi_is_root) then
        write(*,*)" Compute particles dynamics for material: ",trim(mat%matName)
        if (mat%cpVariable) write(*,*) ' >> Solving enthalpy equation'
      endif
      do g = 1, mat%ngroups
        group: associate(gr => mat%group(g))
        call gr%setup_particleODE()
        if (allocated(relTol)) deallocate(relTol); allocate(relTol(gr%neq)); relTol(:) = rtol
        if (allocated(absTol)) deallocate(absTol); allocate(absTol(gr%neq)); absTol(:) = atol
        !> Exclude the body-force accumulators (output quadratures) from the error norm.
        if (srcBodyForce .and. (gr%evapSelect > 0 .or. gr%combSelect > 0)) then
          if (eulerSwitch) then       ! neq-1 is a real euler moment, keep it controlled
            relTol(gr%neq)          = 1.e30_R8; absTol(gr%neq)          = 1.e30_R8
          else                        ! J,W both appended
            relTol(gr%neq-1:gr%neq) = 1.e30_R8; absTol(gr%neq-1:gr%neq) = 1.e30_R8
          endif
        endif
        call setup_odesolver(N=gr%neq,solver=ode_word,RT=relTol,AT=absTol,iopt=iopt)

        if (mpi_is_root) then
          write(*,'(A,I3,A,I6)')"     Group",g," => number of particles = ", gr%nparticles
          write(*,'(A,I2,A)') '     >> ODE system (neq=', gr%particle(1)%neq, '):'
          if (gr%evapSelect/=0) write(*,*) '       - evaporation --> ', trim(gr%evapWord)
          if (gr%combSelect/=0) write(*,*) '       - combustion  --> Beckstead d^n burn law'
          if (gr%brkupEqOde   ) write(*,*) '       - breakup     --> ', trim(gr%brkupWord)
          if (gr%evapSelect==0 .and. gr%combSelect==0 .and. .not.gr%brkupEqOde) &
            write(*,*) '       - constant particle mass and size'
          if (eulerSwitch) write(*,*) '       - eulerian field included'
        endif

        if (trajOn) write(unitTraj,'(A,A,A,I3,A)')'Zone T="Mat ',trim(mat%matName),' Group',g,'"'
        write(unitExit ,'(A,A,A,I3,A)')'Zone T="Mat ',trim(mat%matName),' Group',g,'"'

        if (gr%brkupHasChild) then
          if (allocated(gr%shed)) deallocate(gr%shed)
          allocate(gr%shed(1:gr%nparticles))
          gr%nactive = gr%nparticles
          start = 1;  nEnd = gr%nactive
          loopCounter = 0
          !> Generation 0 = the pinned parents, IDs 1..nInjected (the ID space is per group).
          curBase  = 0
          curCount = gr%nInjected   !> the immutable pin-time census, never the live counter

          do while (loopCounter < maxLoop)
            loopCounter = loopCounter + 1
            !> Grow the shed-list array before resetting the window.
            if (size(gr%shed) < nEnd) call resizeShedArray(gr%shed, nEnd)
            !> Empty each list in the window, with room reserved outside the parallel region.
            do ip = start, nEnd
              gr%shed(ip)%n = 0
              call gr%shed(ip)%reserve(1)
            enddo

            !$OMP PARALLEL DO SCHEDULE(DYNAMIC)
            do ip = start, nEnd
              if (probeOn) then; if (.not.any(gr%particle(ip)%ID==probeIDs)) cycle; endif
              !> MPI ownership filter applies to the parent generation only.
              if (loopCounter == 1) then
                if (.not. owns_particle(ip)) cycle
              endif
              call integrate(gr%particle(ip),self%geoblock,self%gasblock,self%source,self%euler(:,gr%famID), &
                             mat%hTab,mat%cpTab,mat%rhoTab,mat%mupTab,mat%sigTab,mat%psatTab,    &
                             gr%shed(ip), noShed=(loopCounter==maxLoop))
            enddo
            !$OMP END PARALLEL DO

            !> Shed census of the current generation, indexed by position in its global ID band.
            oldStart = start
            oldEnd   = nEnd
            if (allocated(nShedAll)) deallocate(nShedAll, shedOff)
            allocate(nShedAll(curCount), shedOff(curCount))   !> curCount is globally agreed
            nShedAll   = 0
            childLocal = 0
            do ip = oldStart, oldEnd
              pidx = gr%particle(ip)%ID - curBase
              if (pidx < 1 .or. pidx > curCount) &
                call mpi_abort_all('shed census: parcel outside its generation ID band')
              nShedAll(pidx) = gr%shed(ip)%n
              childLocal     = childLocal + gr%shed(ip)%n
            enddo
            !> Pass-1 invariants: generation-0 IDs are the storage index; non-owned parents shed nothing.
            if (loopCounter == 1) then
              do ip = oldStart, oldEnd
                if (gr%particle(ip)%ID /= ip) &
                  call mpi_abort_all('generation-0 IDs are not the 1..nInjected storage index')
                if (.not.owns_particle(ip) .and. gr%shed(ip)%n /= 0) &
                  call mpi_abort_all('non-owned parent shed a child (ownership gate leak)')
              enddo
            endif

            call mpi_allreduce_sum_i4_array(nShedAll, curCount)   !> the ONLY in-loop collective

            !> Exclusive prefix over the global census: total children and each parent's offset.
            childGlobal = 0
            maxShedSeen = 0
            do e = 1, curCount
              shedOff(e)  = childGlobal
              childGlobal = childGlobal + nShedAll(e)
              maxShedSeen = max(maxShedSeen, nShedAll(e))
            enddo
            newBase = curBase + curCount
            if (childGlobal > huge(1) - newBase) &
              call mpi_abort_all('child ID band overflowed int32')

            if (childGlobal > 0) then
              !> Reports the global child count, once.
              if (mpi_is_root) then
                write(*,*)"       Loop",loopCounter," => number of children = ", childGlobal
                if (maxShedSeen > 1) write(*,'(A,I0,A)')                                  &
                    "                             (max ",maxShedSeen," sheds from one parcel)"
              endif

              !> Grow the particle array (geometric 2x) for this rank's children.
              newSize = gr%nactive + childLocal
              if (newSize > size(gr%particle)) then
                call resizeParticleArray(gr%particle, max(2*size(gr%particle), newSize))
              endif

              start = gr%nactive + 1
              nEnd  = newSize
              call gr%setup_particleODE(start, nEnd)
              call gr%assign_group2particle(start, nEnd)

              !> Drain the shed lists in ascending parent index and hand each child its full
              !  initial state (integrate skips its injection-init block for time > 0).
              ch = 0
              do ip = oldStart, oldEnd
              idBase = newBase + shedOff(gr%particle(ip)%ID - curBase)
              do e  = 1, gr%shed(ip)%n
                ch   = ch + 1
                iota = oldEnd + ch
                associate(kid => gr%particle(iota), src => gr%shed(ip)%item(e))
                !> Global ID from the census band; iota is the local storage slot.
                kid%ID            = idBase + e
                if (mpi_size_ == 1 .and. kid%ID /= iota) &
                  call mpi_abort_all('child-ID census broke serial equivalence')
                !> Own RNG stream, seeded from the global ID.
                kid%rngState      = rngSeedFor(gr%famID, kid%ID)
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
              !> Advance the band to the new generation.
              curBase  = newBase
              curCount = childGlobal
            else
              exit
            endif
          enddo

          !> Trim excess capacity and sync nparticles (both rank-local after a solve).
          if (size(gr%particle) > gr%nactive) then
            call resizeParticleArray(gr%particle, gr%nactive)
          endif
          gr%nparticles = gr%nactive
        else
          !$OMP PARALLEL DO SCHEDULE(DYNAMIC)
          do ip = 1,gr%nparticles
            if (probeOn) then; if (.not.any(gr%particle(ip)%ID==probeIDs)) cycle; endif
            !> MPI ownership filter.
            if (.not. owns_particle(ip)) cycle
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

    !> Allreduce the partial grid sums across ranks before finalize, which is nonlinear.
    call reduce_accumulators(self%source, self%euler, sourceSwitch, eulerSwitch)

    !> Finalize the source/euler fields: normalize the moments, reduce to geoblock shape under ord2.
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

    !> Wedge fold counters.
    if (axisym) then
      foldStat = [nSectorFold, nMultiFold]
      call mpi_allreduce_sum_i4_array(foldStat, 2)
      if (mpi_is_root .and. sum(foldStat) > 0) write(*,'(A,I0,A,I0,A)') '     - wedge sector folds: ', &
                                                 foldStat(1), ' (multi-sector: ', foldStat(2), ')'
    endif
    if (mpi_is_root) write(*,*)" Stop condition : All particles out of domain!"

  end subroutine solve


  !> Drag force and heat rate exerted on the gas by np particles of material mID at the given state.
  pure function getSourceTerms(self,gas,vel,Tp,rhop,np,mID) result(FdragQdot)
    use, intrinsic :: iso_fortran_env, only : R8 => real64
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


  !> Write the grid fields (root only) and merge this sweep's per-rank particle files.
  subroutine writeout(self)
    use IGLOO_IO
    use IGLOO_Mod_MPI, only: mpi_is_root, mpi_barrier_env
    implicit none
    class(obj_IGLOO), intent(inout) :: self

    !> Every rank has closed its shards before root merges them.
    call mpi_barrier_env()

    !> Root only from here on.
    if (.not. mpi_is_root) return

    call write_outfield(self%material,self%geoblock,self%source,self%euler,self%srcSwitch, &
                        self%eulSwitch, tag=self%sweepTag())

    !> Collapse this sweep's per-rank .dat shards into the serial layout (no-op at one rank).
    call merge_rank_particle_files(self%material, tag=self%sweepTag())

  end subroutine writeout


  !> Print the IGLOO banner.
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


  !> Resize a particle array to newCapacity, preserving the leading entries.
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


  !> Resize a shed-list array to newCapacity, preserving the leading entries and their capacity.
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
