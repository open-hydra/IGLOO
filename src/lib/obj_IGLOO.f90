module IGLOO_module
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  use IGLOO_data_block, only: obj_block, obj_flowblock, obj_sourceblock, obj_eulerblock
  use IGLOO_data_phases, only: obj_material
  implicit none

  type, public :: obj_IGLOO
    type(obj_block),       allocatable :: geoblock(:)
    type(obj_flowblock),   allocatable :: gasblock(:)
    !> ⚠ Valid only BETWEEN a solve and the next reset_state: reset_state deallocates both
    !> (they are re-allocated at the shape the coming sweep needs). Read them before it.
    type(obj_sourceblock), allocatable :: source(:)
    type(obj_eulerblock),  allocatable :: euler(:,:)
    type(obj_material),    allocatable :: material(:)
    logical :: eulSwitch, srcSwitch
    !> Repeatability bookkeeping. `setup_static` is once-only (mesh, gas import and pinning
    !  are static); `reset_state` restores the per-sweep particle state and must run before
    !  every `solve`. A second `solve` without it silently integrates finished particles.
    logical :: staticDone    = .false.
    logical :: stateIsFresh  = .false.
    !> Output-file generation, counting reset_state calls from 0. Deliberately starts at -1 so
    !  the FIRST reset_state lands on 0, whose tag is empty: a single-sweep standalone run
    !  keeps byte-identical filenames, and the 49 oracles that glob fixed names still match.
    integer :: sweep = -1
    ! integer            :: iprint
    ! real(R8)           :: ds, mdotMax, dtprint
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

  !> Back-compat entry point: the documented hydra hook and all e2e cases call this.
  !  A driver that integrates repeatedly should call setup_static once and reset_state
  !  before each solve instead (see src/app/IGLOO.f90 for the single-sweep shape).
  subroutine setup(self, external_gas)
    use Lib_ORION_data
    implicit none
    class(obj_IGLOO), intent(inout)        :: self
    type(orion_data), intent(in), optional :: external_gas

    call self%setup_static(external_gas)
    call self%reset_state()

  end subroutine setup


  !> Once-only: input parsing, gas field, block/geometry allocation, BC tagging and
  !  particle pinning. Everything here is independent of how many sweeps follow.
  subroutine setup_static(self, external_gas)
    use omp_lib
    use IGLOO_IO
    use IGLOO_IC
    use IGLOO_particles
    use IGLOO_variables
    use IGLOO_allocation
    use IGLOO_IO_INI, only: read_IGLOO_input
    use IGLOO_Lib_Statistics, only: initRandomSeed
    use IGLOO_Mod_MPI,        only: mpi_is_root
    use Lib_ORION_data
    implicit none
    class(obj_IGLOO), intent(inout)        :: self
    type(orion_data), intent(in), optional :: external_gas
    type(orion_data)     :: own_gas
    character(len=llen)  :: gasfile
    real(8), allocatable :: pos0(:,:), vel0(:,:), mdot(:), diam(:), temp0(:)
    integer              :: m, g, fam, nthreads
    character(len=2)     :: method

    !> Not re-runnable: allocate_blocks declares its outputs intent(inout), allocatable and
    !  allocates unconditionally, so a second call aborts on a double allocate. Say so
    !  instead of letting the runtime produce that message.
    if (self%staticDone) then
      write(*,*) ' [IGLOO] setup_static already completed -- repeat call ignored'
      write(*,*) '         (mesh, gas import and pinning are static; use reset_state per sweep)'
      return
    endif

    !> Root-gated banners. The "already completed" message above is NOT gated: it is a contract
    !  violation, and a non-root rank swallowing one would strand the others at the next
    !  collective. Same rule for solve's stateIsFresh failure.
    if (mpi_is_root) call print_header()

    nthreads = 1
# if defined (_OPENMP)
    nthreads = OMP_GET_MAX_THREADS()   !> outside the region: no shared write to race on
# endif
    if (nthreads>1 .and. mpi_is_root) then
      write(*,*)" OpenMP threads = ", nthreads
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
        !> Immutable pin-time census. `solve` overwrites nparticles with nactive (children
        !  folded in), so without this the count is unrecoverable and a second sweep
        !  re-injects the previous sweep's children as if they were originals.
        gr%nInjected = gr%nparticles
        end associate
      enddo
    enddo
    nfam = fam
    self%staticDone = .true.

  end subroutine setup_static


  !> Per-sweep state restore: puts the pinned population back the way pinning left it, so
  !  the next `solve` integrates the same particles from injection again. Must run before
  !  every `solve`; `solve` enforces that through `stateIsFresh`.
  !
  !  `assign_group2particle` lives here rather than in `setup_static` because it is the
  !  per-sweep half of pinning: it is `pure`, writes only `self%particle(i)%*`, and among
  !  other things re-zeroes `brkupVar` for originals (whose only other zeroing is the
  !  child-range call in `solve`, which never reaches an original).
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

    !> Next output generation. Counting here rather than in solve keeps the tag stable across
    !  the solve/writeout pair, which write different files for the SAME sweep.
    self%sweep = self%sweep + 1

    !> Refresh the background field when the parent hands one in. Without this an embedding
    !  integrates a frozen gas forever, which defeats the point of being re-runnable.
    !  The mesh is static, so only the field values are re-imported.
    !  ⚠ cell%mdotGas is NOT refreshed: its seeder early-returns when non-zero
    !  (obj_block.f90:1186), so bcdef-401 streams keep the mass flow they were seeded with.
    !  Under an evolving gas the injected mass flow will not track it. Physics change with
    !  its own gate; recorded, not silently inherited.
    if (present(external_gas)) then
      call copyORION(external_gas, own_gas)
      call import_gas(own_gas, self%gasblock)
    endif

    !> F7 -- the accumulators must be dropped, not reused. Under ord2 they are allocated at
    !  gasblock shape (1..Nx+1) but `finalize` move_allocs geoblock-shaped arrays (1..Nx) over
    !  them at end of solve (obj_block.f90:237-239, :318-321). allocateAccumulators guards on
    !  allocation STATUS, not shape, so a second solve would keep the geoblock-shaped arrays
    !  and then deposit through part%igas -- up to Nx+1 -- past the end, via !$OMP ATOMIC
    !  UPDATE. That is silent heap corruption in a release build, and it was measurable:
    !  vie-plait's source.tec moved 12 decades above its run-to-run noise floor.
    !
    !  Deallocate the COMPLETE set, unconditionally. allocateAccumulators guards on one array
    !  each (sourceMass, density) but allocateSRC/allocateEUL allocate all seven, so freeing
    !  only the guard array makes the next allocate fatal. Not ord2-gated either: it costs
    !  nothing on sweep 1 (nothing is allocated yet, so every guard below is false) and
    !  solve's allocateAccumulators + initialize_fields rebuild and zero them anyway.
    !
    !  Caller contract: self%source and self%euler are valid only BETWEEN a solve and the
    !  next reset_state.
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
        !> solve re-allocates the shed lists per group when the model has children.
        if (allocated(gr%shed)) deallocate(gr%shed)

        !> Fan the pin-time capture back out to the live fields the last sweep consumed.
        !  Restoring d = dInj deliberately re-injects the SAME stochastic draw every sweep:
        !  frozen realization vs resampled ensemble is a physics choice and belongs behind
        !  an explicit input key, not a refactor side effect.
        do ip = 1, gr%nInjected
          associate(part => gr%particle(ip))
          part%ID = ip
          !> Re-seed this particle's RNG stream. Deterministic in (rng_seed, famID, ID), so every
          !  sweep replays the same draws -- which is what makes an RNG-consuming model
          !  repeatable at all. famID comes from the group (set in setup_static), not from
          !  part%famID, which assign_group2particle only fills after this loop.
          part%rngState = rngSeedFor(gr%famID, ip)
          part%d  = part%dInj
          !> Assigned/DB streams only. BC streams re-derive tp and mdot from live gas in
          !  initializePart, and their *Inj fields are 0 -- restoring those would inject
          !  a zero-temperature, zero-mass-flow particle.
          if (all(part%iInj == [0,0,0,0])) then
            part%tp   = part%tpInj
            part%mdot = part%mdotInj
          endif
          !> time == 0 IS the "needs injection init" flag integrate keys off
          !  (Lib_Integration.f90:119), and every exit path zeroes it -- which is what makes
          !  a finished particle re-injectable at all.
          part%time        = 0._R8
          part%lost        = .false.
          part%exitFace    = 0
          part%gasExitFace = 0
          part%Af          = 0._R8
          !> Redundant with integrate's injection-init block, but reset_state should read as
          !  a statement of the fresh state, not as a transcript of integrate's internals.
          part%gone  = .false.
          part%wasin = .false.
          part%Ncell = 0
          part%angle = 0._R8
          part%i     = part%iInj
          !> `iold` has no default initializer and only the DB pin path ever wrote it, so
          !  bc-pinned particles carry stale/undefined values here. Latent today: at
          !  injection findParticle succeeds (injViable guarantees it) and updateCell
          !  returns before reaching the `retry` seed that reads iold. Defined now.
          part%iold  = [0,0,0,0]
          end associate
        enddo

        !> Same order as solve's child hand-off: ODE sizing first, then properties.
        call gr%setup_particleODE()
        call gr%assign_group2particle()
        end associate
      enddo
    enddo

    self%stateIsFresh = .true.

  end subroutine reset_state


  !> Filename tag for the current output generation: empty on the first sweep (so standalone
  !  runs are byte-identical to before this existed), '-sweep<N>' afterwards. Applied to every
  !  output file INCLUDING the .tec ones -- source.tec is what hydra consumes, and it clobbers
  !  just as readily as the trajectory dumps.
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
    !> Child-ID band. curBase/curCount describe the GLOBAL generation held in the current
    !  window: its IDs are (curBase, curBase+curCount]. childLocal sizes LOCAL storage,
    !  childGlobal fixes the next band and the rank-uniform termination test. Deliberately
    !  uninitialized -- an initializer here would make them implicit-SAVE (work package A).
    integer :: curBase, curCount, newBase, childLocal, childGlobal, pidx, idBase
    integer, allocatable :: nShedAll(:), shedOff(:)
    ! logical, allocatable :: famDone(:)
    real(R8), allocatable :: relTol(:), absTol(:)
    !> Scatter-cloud weight-quantum (dNscat) auto-sizing scratch.
    integer  :: c, nStreams, nValid
    real(R8) :: Lref, xmin(3), xmax(3), Ndot, Vsum, Vref, tauRef, vsp
    type(obj_particle) :: ptmp   ! throwaway copy for the inject-only pre-pass

    !> Repeatability contract. Every exit path zeroes part%time, so a finished particle is
    !  indistinguishable from one that never started: calling solve twice without
    !  reset_state does not fail, it silently integrates a corrupted population (stale d,
    !  last sweep's children counted as originals, accumulators of the wrong shape). Fail
    !  loudly instead of returning plausible garbage.
    if (.not. self%stateIsFresh) then
      write(*,'(A)') ' [IGLOO::solve] particle state is not fresh.'
      write(*,'(A)') '   Call reset_state() before each solve() (setup() does it for the first).'
      error stop 1
    endif
    self%stateIsFresh = .false.

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
      !> Per-rank shards. rank_suffix() composes AFTER the sweep tag, giving
      !  `trajectories-A-sweep1.rank2.dat`: `<kind>-<material><sweeptag>` stays the logical file
      !  identity and `.rank<r>` is a pure shard marker, so Phase 4 merges per sweep by globbing
      !  `<logical>.rank*.dat`. Reversing the order would make that glob straddle sweeps.
      !  Empty at one rank, so serial filenames are byte-identical. Each rank writes its own
      !  `variables=` and Zone headers, so every shard stays independently Tecplot-loadable.
      if (trajOn) then
        open(newunit=unitTraj,file='OUTPUT/'//trim(IGLOO_phase_prefix)//'trajectories-'//trim(mat%matName)//self%sweepTag()//trim(rank_suffix())//'.dat')
        write(unitTraj,*) 'variables="X","Y","Z","U","V","W","T","d<sub>p","m<sub>p","ID"'
      endif
      open(newunit=unitExit,file='OUTPUT/'//trim(IGLOO_phase_prefix)//'outloc-'//trim(mat%matName)//self%sweepTag()//trim(rank_suffix())//'.dat')
      write(unitExit,*) 'variables="X","Y","Z","T","|u<sub>p</sub>|","<greek>a</greek>","mdot","Af","ID"'

      !> Scatter cloud: one flat point-cloud zone per material; auto-size the weight quantum
      !  dNscat (real droplets/point) so the cloud holds ~trajSample points per stream. The
      !  population estimate (Ndot*tauRef) is crude; it sets only the count, not the shape.
      dNscat = 0._R8
      if (scatOn) then
        open(newunit=unitScat,file='OUTPUT/'//trim(IGLOO_phase_prefix)//'scatter-'//trim(mat%matName)//self%sweepTag()//trim(rank_suffix())//'.dat')
        write(unitScat,*) 'variables="X","Y","Z","U","V","W","T","d<sub>p","m<sub>p","ID"'
        write(unitScat,'(A,A,A)')'Zone T="Mat ',trim(mat%matName),' scatter"'
        !> Serial inject-only pre-pass: resolve each stream's npdot (= Σ droplet rate Ndot) and
        !  injection speed (=> Vref) without the cell search. Serial => deterministic dNscat.
        !  ⚠ MPI: this stays a FULL REPLICATED SWEEP -- do NOT add an owns_particle guard. Every
        !  rank computing it over ALL particles is what keeps dNscat identical everywhere; an
        !  ownership guard would turn Ndot/Vsum/nValid/nStreams into partial sums and give each
        !  rank a different scatter quantum. Pinning and initRandomSeed stay replicated for the
        !  same reason.
        Ndot = 0._R8; nStreams = 0; Vsum = 0._R8; nValid = 0
        do g = 1, mat%ngroups
          associate(gr => mat%group(g))
          call gr%setup_particleODE()
          do ip = 1, gr%nparticles
            !> Throwaway copy: resolveInjectionRate writes mdot/v/tp on it; keep the real
            !  particle pristine so the real sweep is bit-unchanged. It deliberately avoids the
            !  geometric cell search (which would pollute the threadprivate myRay).
            ptmp = gr%particle(ip)
            !> DB/assigned streams: resolveInjectionRate skips initializePart, so stateVar is
            !  allocated-but-unwritten here. Seed it as integrate does (Lib_Integration:140).
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
        !> Root-gated: every rank computes the same value, so N copies would only duplicate the
        !  line. Informational prints follow this rule throughout; WARNINGS and contract
        !  violations stay on every rank (a silent non-root failure hangs the next collective).
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
          !> Generation 0 = the pinned parents, IDs 1..nInjected (reset_state assigns the storage
          !  index, and so do all three pin sites). curBase=0 makes newBase = n0 on pass 1, which
          !  is byte-for-byte the `oldEnd` that `kid%ID = oldEnd + ch` used before the census.
          !  Per GROUP: the ID space is per group, so this cannot be hoisted out.
          curBase  = 0
          curCount = gr%nInjected   !> the immutable pin-time census, never the live counter

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
              !> Ownership -- FIRST PASS ONLY. Pass 1 iterates the replicated parents, of which
              !  each rank integrates its stripe; later passes iterate [start:nEnd] = the
              !  children this rank created, which it must integrate unconditionally. Identically
              !  true at one rank, so the serial path is unchanged.
              if (loopCounter == 1) then
                if (.not. owns_particle(ip)) cycle
              endif
              call integrate(gr%particle(ip),self%geoblock,self%gasblock,self%source,self%euler(:,gr%famID), &
                             mat%hTab,mat%cpTab,mat%rhoTab,mat%mupTab,mat%sigTab,mat%psatTab,    &
                             gr%shed(ip), noShed=(loopCounter==maxLoop))
            enddo
            !$OMP END PARALLEL DO

            !> Shed census -- total shed EVENTS, not parents that shed, because `newSize` has to
            !  size the particle array by the number of children actually created. One entry per
            !  parcel of the CURRENT generation, indexed by its position in that generation's
            !  global ID band, so MPI_SUM reconstructs the exact serial per-parent vector: every
            !  parcel is written by exactly ONE rank (pass 1 -- a non-owned parent is cycled and
            !  keeps the n=0 stored just above, which every rank writes unconditionally; passes
            !  >=2 -- the window holds only locally-created children, each on its creator).
            !  Child ordering is fixed by the ascending `ip` traversal the drain below repeats,
            !  so it is independent of thread count and SCHEDULE(DYNAMIC).
            oldStart = start
            oldEnd   = nEnd    !> hoisted: needed on every rank every pass, childLocal==0 included
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
            !> Both premises of the pass-1 authority argument, asserted rather than assumed.
            if (loopCounter == 1) then
              do ip = oldStart, oldEnd
                if (gr%particle(ip)%ID /= ip) &
                  call mpi_abort_all('generation-0 IDs are not the 1..nInjected storage index')
                if (.not.owns_particle(ip) .and. gr%shed(ip)%n /= 0) &
                  call mpi_abort_all('non-owned parent shed a child (ownership gate leak)')
              enddo
            endif

            call mpi_allreduce_sum_i4_array(nShedAll, curCount)   !> the ONLY in-loop collective

            !> Exclusive prefix over the GLOBAL vector: childGlobal fixes the next band and the
            !  rank-uniform termination test, shedOff places each parent's children inside it.
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
              !> Root-gated, and childGlobal (not childLocal) is load-bearing: khrt-e2e's
              !  check.py and check_threads.py parse this line out of run_out.txt and sum it, so
              !  it has to report the same total serial does, exactly once.
              if (mpi_is_root) then
                write(*,*)"       Loop",loopCounter," => number of children = ", childGlobal
                !> Visible tripwire: how hard the shed path is actually working. A number
                !  climbing toward the `maxShed` guard in Lib_Integration means the case is
                !  approaching runaway shedding well before the guard has to fire.
                if (maxShedSeen > 1) write(*,'(A,I0,A)')                                  &
                    "                             (max ",maxShedSeen," sheds from one parcel)"
              endif

              !> Grow particle array if needed (geometric 2x). Sized on childLocal, not
              !  childGlobal: only this rank's children live in this rank's array. `oldEnd` was
              !  hoisted above the census -- do not re-derive it here, `nEnd` moves just below.
              newSize = gr%nactive + childLocal
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
              idBase = newBase + shedOff(gr%particle(ip)%ID - curBase)
              do e  = 1, gr%shed(ip)%n
                ch   = ch + 1
                iota = oldEnd + ch
                associate(kid => gr%particle(iota), src => gr%shed(ip)%item(e))
                !> ID comes from the GLOBAL census band; `iota` stays the LOCAL storage slot. Keep
                !  the two axes apart -- at size 1 they coincide and the tripwire says so, but
                !  part%ID also seeds the scatter sampling (Lib_Integration:182), so silently
                !  re-baselining it would move which scatter points get emitted.
                kid%ID            = idBase + e
                if (mpi_size_ == 1 .and. kid%ID /= iota) &
                  call mpi_abort_all('child-ID census broke serial equivalence')
                !> Own RNG stream. Without this a child keeps the default 0 -- deterministic, but
                !  IDENTICAL for every child, which the repeatability gate cannot see. Seeded from
                !  the global ID, so the stream is rank-invariant as well as thread-invariant.
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
              !> Advance the band: the children just created are the next generation.
              curBase  = newBase
              curCount = childGlobal
            else
              exit
            endif
          enddo

          !> Trim excess capacity + sync nparticles.
          !  ⚠ MPI: both counters are now RANK-LOCAL after a solve -- each rank holds only its own
          !  stripe plus the children it created. Nothing downstream in standalone mode consumes
          !  the particle arrays, and reset_state restores both from gr%nInjected before every
          !  sweep, so the divergence cannot survive into the next one. Drive ownership and
          !  curCount from gr%nInjected, never from these.
          if (size(gr%particle) > gr%nactive) then
            call resizeParticleArray(gr%particle, gr%nactive)
          endif
          gr%nparticles = gr%nactive
        else
          !$OMP PARALLEL DO SCHEDULE(DYNAMIC)
          do ip = 1,gr%nparticles
            if (probeOn) then; if (.not.any(gr%particle(ip)%ID==probeIDs)) cycle; endif
            !> Ownership: no children here, so every particle is a replicated parent and the
            !  stripe is unconditional. Identically true at one rank.
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

    !> MPI: merge the per-rank partial grid sums BEFORE finalize, because finalize is NONLINEAR --
    !  it divides the +=-accumulated numerators by density (Favre average), so reducing afterwards
    !  would average averages. ALLREDUCE, not reduce-to-root: the finalize below then runs on
    !  identical raw sums on every rank, which keeps the post-solve object state bit-identical
    !  everywhere and is what makes the hydra embedding straightforward. The ord2
    !  gasblock->geoblock reduction inside finalize is linear and runs once on already-reduced
    !  arrays. Must sit HERE and not in writeout: reset_state deallocates all seven accumulators
    !  every sweep, so they are only guaranteed allocated between allocateAccumulators and the next
    !  reset_state -- shapes are therefore read fresh per sweep, never cached across sweeps.
    !  INVARIANT: with the Phase-2a census this is the SECOND and last collective in solve().
    call reduce_accumulators(self%source, self%euler, sourceSwitch, eulerSwitch)

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

    if (mpi_is_root) write(*,*)" Stop condition : All particles out of domain!"

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
    use IGLOO_Mod_MPI, only: mpi_is_root, mpi_barrier_env
    implicit none
    class(obj_IGLOO), intent(inout) :: self

    !> Every rank must have closed its own .dat shards before root reads or merges them (Phase 4),
    !  so the barrier stays even though the write below is root-only. No-op at one rank.
    call mpi_barrier_env()

    !> ROOT ONLY. write_outfield writes exclusively grid .tec files to one fixed name per
    !  (material, sweep) -- no rank shard -- so every rank calling it is N writers on one path.
    !  Since Phase 3 the values are correct on every rank (solve allreduces the accumulators before
    !  finalize), so root-only is the permanent answer here, not scaffolding.
    if (.not. mpi_is_root) return

    call write_outfield(self%material,self%geoblock,self%source,self%euler,self%srcSwitch, &
                        self%eulSwitch, tag=self%sweepTag())

    !> Choke-point 3: collapse THIS sweep's per-rank .dat shards into the serial layout, so the
    !  oracles and any parent post-processing see one file per (kind, material, sweep). No-op at one
    !  rank. Safe here because solve closed every shard on every rank before returning; the barrier
    !  above makes that ordering explicit instead of incidental.
    call merge_rank_particle_files(self%material, tag=self%sweepTag())

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
