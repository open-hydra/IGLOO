module IGLOO_IC
  implicit none

contains

subroutine pin_particles(group, block, gas, method, pos0,vel,temp,mdot,diam,fam,npart)
  use IGLOO_variables,    only: ds, mdotMax, dsSwitch
  use IGLOO_RayFaceIntersection3D
  use IGLOO_data_phases, only: obj_group
  use IGLOO_data_block, only: obj_block, obj_flowblock
  implicit none
  type(obj_block),     intent(inout) :: block(:)   !> inout: pin_particles_bc_ds writes cell%nInj
  type(obj_flowblock), intent(in)    :: gas(:)
  type(obj_group),     intent(inout) :: group
  character(len=2), intent(in)    :: method
  real(8), intent(in)    :: mdot(:), diam(:), temp(:)
  real(8), intent(in)    :: pos0(:,:), vel(:,:)
  integer, intent(in)    :: fam
  integer, intent(out)   :: npart
  integer                :: np

  if     (method=='FB' .and. ds==0.0 .and. mdotMax==0.0 .and. .not.dsSwitch) then
    call pin_particles_bc_center(group,block,fam,npart)
  elseif (method=='FB') then
    call pin_particles_bc_ds(group,block,gas,fam,npart)
  else
    npart = size(pos0,1)
    if (npart==0) then
      error stop ( '[ERROR] NO PARTICLES INJECTED!' )
    endif
    allocate(group%particle(1:npart))
      do np = 1, npart
        associate(part => group%particle(np))
        part%ID = np
        part%pInj = pos0(np,:)
        part%iInj = [0,0,0,0]
        part%iold = [0,0,0,0]
        part%mdot = mdot(np)
        part%d    = diam(np)
        part%tp   = temp(np)
        part%vInj = vel(np,:)   ! stateVar not yet allocated here; handed off in integrate's init
        end associate
      enddo
  endif

contains

  subroutine pin_particles_bc_center(group, block, fam, npart)
    use IGLOO_variables, only: fsample
    use IGLOO_Lib_Statistics, only: sampleDiameter
    implicit none
    type(obj_block),  intent(in)    :: block(:)
    type(obj_group),  intent(inout) :: group
    integer, intent(in)  :: fam
    integer, intent(out) :: npart
    integer  :: nn=0, p=0
    integer  :: nb, b, f, m, n, i, j, k
    real(8)  :: eps=1.0e-8_R8

    write(*,*)
    write(*,'(A,I3,A)')'  Placing particles at the bounday every', fsample ,' cells...'

    nb = size(block)
    npart = 0
    do b = 1, nb; do f = 1, 6
      associate(face => block(b)%face(f))
      do m = 1, face%Nm; do n = 1, face%Nn
          if (isInj(face%cell(m,n)%bcdef)) then
              nn = nn+1
            if (mod(nn,fsample)==0) npart = npart+1
          endif
      enddo; enddo
      end associate
    enddo; enddo
    allocate(group%particle(1:npart))
    associate(particle => group%particle)
    nn = 0
    do b = 1, nb; do f = 1, 6
      associate(face => block(b)%face(f))
      do m = 1, face%Nm; do n = 1, face%Nn
        if (isInj(face%cell(m,n)%bcdef)) then
            nn = nn+1
          if (mod(nn,fsample)==0) then
            p = p+1
            call block(b)%fmn2ijk(f,m,n,i,j,k)
            particle(p)%ID = p
            particle(p)%pInj = face%cell(m,n)%center - eps*face%cell(m,n)%normal
            particle(p)%iInj = [b, i, j, k]
            particle(p)%fInj = f
            particle(p)%Ninj = 1
            particle(p)%d = sampleDiameter( 2._R8*face%cell(m,n)%properties(fam,6), &  ! dp = 2*rp
                                                  face%cell(m,n)%properties(fam,7), &  ! sigmap
                                            nint( face%cell(m,n)%properties(fam,8) ) ) ! law code
          endif
        endif
      enddo; enddo
      end associate
    enddo; enddo
    end associate
  end subroutine pin_particles_bc_center

end subroutine pin_particles


!> Pin one-or-more injection particles per boundary cell on inflow faces.
!  Spacing is set per cell by compute_ds (per-cell ds col 9, global ds, or mdotMax).
!  2D faces (Nz==1) use the in-plane circle/line sweep; 3D faces use an
!  advancing-front (BFS) hexagonal packing seeded from each injection cell.
subroutine pin_particles_bc_ds(group, block, gas, fam, npart)
  use, intrinsic :: iso_fortran_env, only : R8 => real64
  use IGLOO_variables, only: ds, mdotMax, nb, pi
  use IGLOO_RayFaceIntersection3D, only: isPointInsideHexahedron12
  use IGLOO_VectorModule, only: rotateVector
  use IGLOO_Lib_Statistics, only: sampleDiameter
  use IGLOO_data_phases, only: obj_group
  use IGLOO_data_block,  only: obj_block, obj_flowblock
  implicit none
  type(obj_block),     intent(inout) :: block(:)
  type(obj_flowblock), intent(in)    :: gas(:)
  type(obj_group),     intent(inout) :: group
  integer,             intent(in)    :: fam
  integer,             intent(out)   :: npart
  integer  :: alloc=1e7
  integer  :: p, b, f, m, n, i, j, k, Nm, Nn
  integer  :: mInj, nInj, head, m0, n0, mc, nc, m2, n2, ic, jc, kc, i2, j2, k2, kd
  integer,      allocatable :: ind(:,:)
  real(R8),     allocatable :: Inj(:,:)
  logical,      allocatable :: singleCell(:,:)
  real(R8),     allocatable :: Rcell(:,:), Dsize(:,:)
  logical  :: startSlot, injected, skip, retry, found, single, overflow
  real(R8) :: cA, cB, xInj, yInj, const, ratio, R, dsR, dcell, dir, orient
  real(R8) :: eps, toll, delta, array(3), pCand(3), start(3), vertices(3,8)
  real(R8) :: pt(3), Pc(3), dvec(3), tng(3), nrm(3), n2v(3), Rc, verts2(3,8)

  write(*,*)
  if     ((ds/=0._R8) .and. (mdotMax==0._R8)) then
    write(*,'(A,F9.5,A)')'  Placing particles at the bounday every', ds*100 ,' cm...'
  elseif ((ds==0._R8) .and. (mdotMax/=0._R8)) then
    write(*,'(A,F9.5,A)')'  Placing particles such that maximum mdot <', mdotMax*1000 ,' g/s...'
  else
    write(*,'(A,F9.5,A)')'  Placing particles such that maximum mdot <', mdotMax*1000 ,' g/s'
    write(*,'(A,F9.5,A)')'  - maximum distance between injection points is', ds*100 ,' cm...'
  endif

  npart    = 0
  toll     = 1.0e-10_R8
  eps      = 1.0e-8_R8
  overflow = .false.
  allocate(ind(7,alloc))
  allocate(Inj(3,alloc))

  if (block(1)%Nz==1) then
    !> -------------------------------------------------------------------------
    !> 2D injection algorithm: in-plane circle/line sweep along m
    !> -------------------------------------------------------------------------
    write(*,*)' >> Using 2D injection algorithm'
    write(*,*) ''
    blockCycle: do b = 1, nb; do f = 1, 6
      pinFace: associate(face => block(b)%face(f))
      Nm = face%Nm; Nn = face%Nn
      allocate(singleCell(Nm,Nn))

      !> Pass 0 — reset tally, classify single cells, seed their one (viable) centred
      !  point. Single cells are always gated out of the sweep; a degenerate single
      !  cell (collapsing slivers, center not locatable) injects nothing (injViable).
      do n = 1, Nn; do m = 1, Nm
        face%cell(m,n)%nInj = 0
        singleCell(m,n)     = .false.
        if (.not.isInj(face%cell(m,n)%bcdef)) cycle
        call block(b)%fmn2ijk(f,m,n,i,j,k)
        call block(b)%getVertices([i,j,k],vertices)
        call compute_ds(block(b),gas(b),f,m,n,i,j,k,vertices,group,ds,mdotMax,dsR,single,dcell)
        singleCell(m,n) = single
        if (single) then
          pCand = face%cell(m,n)%center - eps*face%cell(m,n)%normal
          if (injViable(pCand,vertices,dcell)) then
            call place(b,i,j,k,f,m,n, pCand)
            if (overflow) exit blockCycle
          endif
        endif
      enddo; enddo

      !> Sweep (multi-particle cells only; single cells gated out by singleCell)
      do n = 1, Nn
        mInj = 0
        do m = 1, Nm
          startSlot = .true.
          pinCell: associate(cell => face%cell)
          if (isInj(cell(m,n)%bcdef) .and. (m >= mInj) .and. .not.singleCell(m,n)) then
            mInj = m
            nInj = n
            call block(b)%fmn2ijk(f,mInj,nInj,i,j,k)
            select case (f)
            case (1,3,5)
              start = block(b)%node(:,i-1,j-1,k-1)
            case (2)
              start = block(b)%node(:, i ,j-1,k-1)
            case (4)
              start = block(b)%node(:,i-1, j ,k-1)
            case (6)
              start = block(b)%node(:,i-1,j-1, k )
            end select
            call block(b)%getVertices([i,j,k],vertices)
            call compute_ds(block(b),gas(b),f,mInj,nInj,i,j,k,vertices,group,ds,mdotMax,dsR,single,dcell)
            R = 0.5_R8*dsR
            do while ((mInj <= Nm) .and. isInj(cell(mInj,n)%bcdef) .and. .not.singleCell(mInj,n))
              injected = .false.
              skip  = .false.
              retry = .true.
              do while ((.not.injected) .and. isInj(cell(mInj,n)%bcdef) .and. .not.singleCell(mInj,n))
                associate(cc => cell(mInj,nInj)%center, &
                          cn => cell(mInj,nInj)%normal)

                !> Choose direction (increasing/decreasing x coordinate)
                if (startSlot) then
                  array  = [0.0_R8,1.0_R8,0.0_R8]
                  orient = dot_product(cn, array)
                  dir = sign(1._R8, orient)
                endif

                !> Compute Intersection (circumference-line)
                if (abs(cn(1)) < toll) then
                  xInj = start(1) + dir*sqrt(R*R - (cc(2)-start(2))**2)
                  yInj = cc(2)
                elseif (abs(cn(2)) < toll) then
                  xInj = cc(1)
                  yInj = start(2) + dir*sqrt(R*R - (cc(1)-start(1))**2)
                else
                  ratio = cn(1)/cn(2)
                  const = cc(1)*ratio + cc(2) - start(2)
                  cA = 1 + ratio*ratio                      !> A
                  cB = const*ratio + start(1)               !> -B/2
                  delta = cA*R*R - ((cc(2)-start(2))+ratio*(cc(1)-start(1)))**2
                  xInj  = (cB + dir*sqrt(delta))/cA
                  yInj  = ratio*(cc(1)-xInj) + cc(2)
                endif
                if (isnan(xInj) .or. isnan(yInj)) then
                  error stop ( '[ERROR] Computing injection points: NaN' )
                endif
                array = [xInj, yInj, 0.0_R8]
                pCand = array - eps*cn

                !> Injection point inside cell verification
                if (isPointInsideHexahedron12(pCand,vertices)) then
                  call place(b,i,j,k,f,mInj,nInj,pCand)
                  if (overflow) exit
                  startSlot = .false.
                  injected  = .true.
                  R = dsR
                elseif (skip.and.retry) then
                  dir   = -dir
                  retry = .false.
                  if (startSlot) startSlot = .false.
                else
                  skip  = .true.
                  retry = .true.
                  mInj  = mInj + 1
                  if (mInj>Nm) exit
                  call block(b)%fmn2ijk(f,mInj,nInj,i,j,k)
                  call block(b)%getVertices([i,j,k],vertices)
                  if (mdotMax/=0._R8) then
                    select case (f)
                    case (1,3,5)
                      start = block(b)%node(:,i-1,j-1,k-1)
                    case (2)
                      start = block(b)%node(:, i ,j-1,k-1)
                    case (4)
                      start = block(b)%node(:,i-1, j ,k-1)
                    case (6)
                      start = block(b)%node(:,i-1,j-1, k )
                    end select
                    call compute_ds(block(b),gas(b),f,mInj,nInj,i,j,k,vertices,group,ds,mdotMax,dsR,single,dcell)
                    R = 0.5_R8*dsR
                  endif
                endif
                end associate
              enddo
              if (injected) start = array
              if (overflow .or. mInj>Nm) exit
            enddo
          endif
          if (overflow) exit blockCycle
          end associate pinCell
        enddo
      enddo

      !> Pass 2 — coverage guarantee: every (non-degenerate) injection cell the
      !  sweep missed gets one viable centred particle.
      do n = 1, Nn; do m = 1, Nm
        if (.not.isInj(face%cell(m,n)%bcdef)) cycle
        if (singleCell(m,n) .or. face%cell(m,n)%nInj>0) cycle
        call block(b)%fmn2ijk(f,m,n,i,j,k)
        call block(b)%getVertices([i,j,k],vertices)
        call compute_ds(block(b),gas(b),f,m,n,i,j,k,vertices,group,ds,mdotMax,dsR,single,dcell)
        pCand = face%cell(m,n)%center - eps*face%cell(m,n)%normal
        if (injViable(pCand,vertices,dcell)) then
          call place(b,i,j,k,f,m,n, pCand)
          if (overflow) exit blockCycle
        endif
      enddo; enddo

      deallocate(singleCell)
      end associate pinFace
    enddo; enddo blockCycle
    if (allocated(singleCell)) deallocate(singleCell)

  else
    !> -------------------------------------------------------------------------
    !> 3D injection algorithm: advancing-front (BFS) hexagonal packing.
    !  Each placed point spawns 6 hex-neighbour candidates one spacing R away,
    !  obtained by rotating the cell's local tangent about its normal (Rodrigues).
    !  check_distance enforces the minimum spacing; locate_inj_cell finds which
    !  injection cell a candidate falls in. Single cells get one centred point
    !  (Pass 0); uncovered cells are back-filled (Pass 2).
    !> -------------------------------------------------------------------------
    write(*,*)' >> Using 3D injection algorithm'
    write(*,*) ''
    blockCycle3D: do b = 1, nb; do f = 1, 6
      associate(face => block(b)%face(f))
      Nm = face%Nm; Nn = face%Nn
      allocate(singleCell(Nm,Nn), Rcell(Nm,Nn), Dsize(Nm,Nn))

      !> Pass 0 — reset tally, cache per-cell spacing & size, seed viable single cells.
      do n = 1, Nn; do m = 1, Nm
        face%cell(m,n)%nInj = 0
        singleCell(m,n)     = .false.
        Rcell(m,n)          = 0._R8
        Dsize(m,n)          = 0._R8
        if (.not.isInj(face%cell(m,n)%bcdef)) cycle
        call block(b)%fmn2ijk(f,m,n,i,j,k)
        call block(b)%getVertices([i,j,k],vertices)
        call compute_ds(block(b),gas(b),f,m,n,i,j,k,vertices,group,ds,mdotMax,dsR,single,dcell)
        Rcell(m,n)      = dsR
        Dsize(m,n)      = dcell
        singleCell(m,n) = single
        if (single) then
          pCand = face%cell(m,n)%center - eps*face%cell(m,n)%normal
          if (injViable(pCand,vertices,dcell)) then
            call place(b,i,j,k,f,m,n, pCand)
            if (overflow) exit blockCycle3D
          endif
        endif
      enddo; enddo

      !> Pass 1 — advancing front over the remaining (non-single) injection cells.
      !  The append arrays Inj/ind double as the BFS queue: head walks the tail.
      head = npart
      do n0 = 1, Nn; do m0 = 1, Nm
        if (.not.isInj(face%cell(m0,n0)%bcdef)) cycle
        if (singleCell(m0,n0) .or. face%cell(m0,n0)%nInj>0) cycle
        call block(b)%fmn2ijk(f,m0,n0,i,j,k)
        call block(b)%getVertices([i,j,k],vertices)
        pCand = face%cell(m0,n0)%center - eps*face%cell(m0,n0)%normal
        if (.not.injViable(pCand,vertices,Dsize(m0,n0))) cycle
        call place(b,i,j,k,f,m0,n0, pCand)
        if (overflow) exit blockCycle3D
        do while (head < npart)
          head = head + 1
          pt = Inj(:,head)
          mc = ind(6,head); nc = ind(7,head)
          call block(b)%fmn2ijk(f,mc,nc,ic,jc,kc)
          nrm = face%cell(mc,nc)%normal
          tng = cellTangent(b,f,ic,jc,kc,nrm)
          Rc  = Rcell(mc,nc)
          do kd = 0, 5
            dvec = rotateVector(tng, nrm, real(kd,R8)*pi/3._R8)
            Pc   = pt + Rc*dvec
            call locate_inj_cell(b,f,Nm,Nn,mc,nc,Pc, m2,n2,verts2,found)
            if (.not.found) cycle
            if (singleCell(m2,n2)) cycle
            !> project the candidate onto the target cell plane, push just inside
            n2v = face%cell(m2,n2)%normal
            Pc  = Pc - dot_product(Pc-face%cell(m2,n2)%center, n2v)*n2v - eps*n2v
            if (.not.isPointInsideHexahedron12(Pc,verts2)) cycle
            if (.not.check_distance(npart,Rc,Inj,Pc)) cycle
            call block(b)%fmn2ijk(f,m2,n2,i2,j2,k2)
            call place(b,i2,j2,k2,f,m2,n2, Pc)
            if (overflow) exit
          enddo
          if (overflow) exit
        enddo
        if (overflow) exit blockCycle3D
      enddo; enddo

      !> Pass 2 — coverage guarantee: any (non-degenerate) injection cell the front
      !  never reached gets one viable centred particle.
      do n = 1, Nn; do m = 1, Nm
        if (.not.isInj(face%cell(m,n)%bcdef)) cycle
        if (singleCell(m,n) .or. face%cell(m,n)%nInj>0) cycle
        call block(b)%fmn2ijk(f,m,n,i,j,k)
        call block(b)%getVertices([i,j,k],vertices)
        pCand = face%cell(m,n)%center - eps*face%cell(m,n)%normal
        if (.not.injViable(pCand,vertices,Dsize(m,n))) cycle
        call place(b,i,j,k,f,m,n, pCand)
        if (overflow) exit blockCycle3D
      enddo; enddo

      deallocate(singleCell, Rcell, Dsize)
      end associate
    enddo; enddo blockCycle3D
    if (allocated(singleCell)) deallocate(singleCell)
    if (allocated(Rcell))      deallocate(Rcell)
    if (allocated(Dsize))      deallocate(Dsize)
  endif

  !> --- Final assembly: one obj_particle per pinned point; Ninj = cell tally ---
  allocate(group%particle(1:npart))
  associate(particle => group%particle)
  do p = 1, npart
    particle(p)%ID   = p
    particle(p)%pInj = Inj(:,p)
    particle(p)%iInj = ind(1:4,p)
    particle(p)%fInj = ind( 5 ,p)
    b = ind(1,p); f = ind(5,p); m = ind(6,p); n = ind(7,p)
    particle(p)%Ninj = block(b)%face(f)%cell(m,n)%nInj
    particle(p)%d = sampleDiameter( 2._R8*block(b)%face(f)%cell(m,n)%properties(fam,6), &  ! dp = 2*rp
                                          block(b)%face(f)%cell(m,n)%properties(fam,7), &  ! sigmap
                                    nint( block(b)%face(f)%cell(m,n)%properties(fam,8) ) ) ! law code
  enddo
  end associate

contains

  !> Viability guard for a centred (single / coverage) placement. Rejects a cell
  !  that is degenerate for injection so the integrator never receives a particle
  !  it cannot locate or that sits in a collapsing sliver:
  !    (1) tangential size >= dsDegen (configurable floor; 0 disables the check), and
  !    (2) the candidate point passes the *exact* per-cell acceptance test the
  !        integrator's findParticle uses (axis-aligned bbox + z-range + point-in-
  !        cell). Matching it guarantees findParticle re-locates the point in its
  !        own cell instead of falling into the boundary-crossing path (which can
  !        overshoot a domain-edge cell to index Ny+1 and crash getVertices).
  logical function injViable(point,verts,dcell)
    use IGLOO_variables,             only: mesh2D, dsDegen
    use IGLOO_RayFaceIntersection3D, only: isPointInsideCell
    implicit none
    real(R8), intent(in) :: point(3), verts(3,8), dcell
    real(R8) :: cmin(3), cmax(3)
    integer  :: d

    injViable = .false.
    if (dcell < dsDegen) return
    do d = 1, 3
      cmin(d) = minval(verts(d,:)); cmax(d) = maxval(verts(d,:))
    enddo
    if (point(1) < cmin(1) .or. point(1) > cmax(1)) return
    if (point(2) < cmin(2) .or. point(2) > cmax(2)) return
    if (point(3) < cmin(3) .or. point(3) > cmax(3)) return
    injViable = isPointInsideCell(point, verts, mesh2D)
  end function injViable

  !> Append one pinned point to Inj/ind and bump the owning cell's nInj tally.
  !  Single shared overflow guard for every placement site.
  subroutine place(b,i,j,k,f,m,n,point)
    implicit none
    integer,  intent(in) :: b,i,j,k,f,m,n
    real(R8), intent(in) :: point(3)

    if (npart >= alloc) then
      if (.not.overflow) write(*,*) '[WARNING] Exceeded maximum particle number (', alloc,')'
      overflow = .true.
      return
    endif
    npart = npart + 1
    Inj(:,npart) = point
    ind(:,npart) = [b,i,j,k,f,m,n]
    block(b)%face(f)%cell(m,n)%nInj = block(b)%face(f)%cell(m,n)%nInj + 1
  end subroutine place

  !> Find the injection cell of face (b,f) that contains point P. Tests the
  !  current cell and its 8 neighbours first, then a full face scan.
  subroutine locate_inj_cell(b,f,Nm,Nn,mc,nc,P, m2,n2,verts2,found)
    use IGLOO_RayFaceIntersection3D, only: isPointInsideHexahedron12
    implicit none
    integer,  intent(in)  :: b,f,Nm,Nn,mc,nc
    real(R8), intent(in)  :: P(3)
    integer,  intent(out) :: m2,n2
    real(R8), intent(out) :: verts2(3,8)
    logical,  intent(out) :: found
    integer  :: ii,jj,kk, mm, nq, dm, dn
    real(R8) :: vv(3,8), nr(3)

    found = .false.; m2 = 0; n2 = 0; verts2 = 0._R8
    !> neighbourhood first (current cell included at dm=dn=0)
    do dn = -1, 1; do dm = -1, 1
      mm = mc+dm; nq = nc+dn
      if (mm<1 .or. mm>Nm .or. nq<1 .or. nq>Nn) cycle
      if (.not.isInj(block(b)%face(f)%cell(mm,nq)%bcdef)) cycle
      call block(b)%fmn2ijk(f,mm,nq,ii,jj,kk)
      call block(b)%getVertices([ii,jj,kk],vv)
      nr = block(b)%face(f)%cell(mm,nq)%normal
      if (isPointInsideHexahedron12(P - eps*nr, vv)) then
        m2 = mm; n2 = nq; verts2 = vv; found = .true.; return
      endif
    enddo; enddo
    !> full face scan
    do nq = 1, Nn; do mm = 1, Nm
      if (.not.isInj(block(b)%face(f)%cell(mm,nq)%bcdef)) cycle
      call block(b)%fmn2ijk(f,mm,nq,ii,jj,kk)
      call block(b)%getVertices([ii,jj,kk],vv)
      nr = block(b)%face(f)%cell(mm,nq)%normal
      if (isPointInsideHexahedron12(P - eps*nr, vv)) then
        m2 = mm; n2 = nq; verts2 = vv; found = .true.; return
      endif
    enddo; enddo
  end subroutine locate_inj_cell

  !> Local in-plane tangent of cell (i,j,k) on face f: first m-edge of the cell,
  !  projected onto the cell plane and normalised. (Same edge picks as the legacy
  !  fixed-frame dSin, but evaluated per cell.)
  function cellTangent(b,f,i,j,k,normal) result(t)
    implicit none
    integer,  intent(in) :: b,f,i,j,k
    real(R8), intent(in) :: normal(3)
    real(R8) :: t(3), tn

    select case (f)
    case (1); t = block(b)%node(:,i-1, j ,k-1) - block(b)%node(:,i-1,j-1,k-1)
    case (2); t = block(b)%node(:, i , j ,k-1) - block(b)%node(:, i ,j-1,k-1)
    case (3); t = block(b)%node(:, i ,j-1,k-1) - block(b)%node(:,i-1,j-1,k-1)
    case (4); t = block(b)%node(:, i , j ,k-1) - block(b)%node(:,i-1, j ,k-1)
    case (5); t = block(b)%node(:, i ,j-1,k-1) - block(b)%node(:,i-1,j-1,k-1)
    case (6); t = block(b)%node(:, i ,j-1, k ) - block(b)%node(:,i-1,j-1, k )
    end select
    t  = t - dot_product(t,normal)*normal
    tn = norm2(t)
    if (tn > toll) t = t/tn
  end function cellTangent

end subroutine pin_particles_bc_ds


logical function check_distance(npart,R,Inj,newInj)
  use, intrinsic :: iso_fortran_env, only : R8 => real64
  implicit none
  real(R8), intent(in) :: Inj(3,npart), newInj(3)
  real(R8), intent(in) :: R
  integer,  intent(in) :: npart
  integer :: i

  check_distance = .true.
  if (npart>=1) then
    do i = 1, npart
      if (norm2(newInj-Inj(:,i))<(R/2)) then
        check_distance = .false.
        return
      endif
    enddo
  endif

end function check_distance

subroutine compute_ds(block,gas,f,m,n,i,j,k,vertices,group,ds,mdotMax,R,single,dcell)
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  use IGLOO_data_phases, only: obj_group
  use IGLOO_data_block,     only: obj_block, obj_flowblock
  implicit none
  class(obj_block),    intent(in)  :: block
  type(obj_flowblock), intent(in)  :: gas
  type(obj_group),     intent(in)  :: group
  real(R8),            intent(in)  :: vertices(3,8)
  real(8),             intent(in)  :: ds, mdotMax
  integer,             intent(in)  :: f, m, n, i, j, k
  real(8),             intent(out) :: R
  logical,             intent(out) :: single   !> .true. => exactly one centred particle in this cell
  real(R8),            intent(out) :: dcell     !> cell tangential (mid-edge) size [m]
  real(R8) :: dl(3), dsCell, dsEff
  integer  :: NpMin

  single = .false.

  ! Cell tangential size (mid-edge along the m-direction): needed both for the
  ! mdotMax spacing and to detect ds >= cell (=> a single particle).
  select case (f)
  case (1); dl = 0.5_R8*(block%node(:,i-1, j ,k-1) - block%node(:,i-1,j-1,k-1) + &
                         block%node(:,i-1, j , k ) - block%node(:,i-1,j-1, k ))
  case (2); dl = 0.5_R8*(block%node(:, i , j ,k-1) - block%node(:, i ,j-1,k-1) + &
                         block%node(:, i , j , k ) - block%node(:, i ,j-1, k ))
  case (3); dl = 0.5_R8*(block%node(:, i ,j-1,k-1) - block%node(:,i-1,j-1,k-1) + &
                         block%node(:, i ,j-1, k ) - block%node(:,i-1,j-1, k ))
  case (4); dl = 0.5_R8*(block%node(:, i , j ,k-1) - block%node(:,i-1, j ,k-1) + &
                         block%node(:, i , j , k ) - block%node(:,i-1, j , k ))
  case (5); dl = 0.5_R8*(block%node(:,i-1, j ,k-1) - block%node(:,i-1,j-1,k-1) + &
                         block%node(:, i , j ,k-1) - block%node(:, i ,j-1,k-1))
  case (6); dl = 0.5_R8*(block%node(:,i-1, j , k ) - block%node(:,i-1,j-1, k ) + &
                         block%node(:, i , j , k ) - block%node(:, i ,j-1, k ))
  end select
  dcell = norm2(dl)

  ! Per-cell ds (bc.txt col 9, meters) overrides the global ds; <=0 falls back to it.
  dsCell = block%face(f)%cell(m,n)%properties(group%famID, 9)
  if (dsCell > 0._R8) then; dsEff = dsCell; else; dsEff = ds; endif

  if (mdotMax==0._R8) then
    if (dsEff <= 0._R8 .or. dsEff >= dcell) then
      single = .true.    ! ds >= cell (or unset): exactly one centred particle
      R = dcell
    else
      R = dsEff
    endif
    return
  endif

  associate(cell => block%face(f)%cell(m,n), &
            prop => block%face(f)%cell(m,n)%properties(group%famID,1:7))
  select case (cell%bcdef)
  case (401)
    ! mdot_p = krho/(1-krhoTot)/Ninj · (mdotGas + mdotPart)
    ! NpMin (~ Ninj when mdot_p hits mdotMax) inverts that relation.
    NpMin = ceiling(prop(1) / (1._R8 - cell%krhoTot) / mdotMax * &
                    (abs(cell%mdotGas) + cell%mdotPart))
  case (402, 403)
    ! prop(1) is gp [kg/(s·m²)]; mdot_p = gp*area/Ninj => NpMin = gp*area/mdotMax
    NpMin = ceiling(prop(1) * cell%area / mdotMax)
  case default
    write(*,'(A,I0)') ' [compute_ds] unsupported cell%bcdef=', cell%bcdef
    error stop 1
  end select
  end associate
  if (NpMin <= 1) then
    single = .true.    ! mdotMax already satisfied by one particle
    R = dcell
    return
  endif
  R = dcell/NpMin
  if (dsEff > 0._R8) R = min(dsEff,R)

end subroutine compute_ds

logical function isInj(bcdef)
  implicit none
  integer, intent(in) :: bcdef

  ! Cells that inject mass: ATLAS 401-420 (inlet/outlet) and 501-502 (SRM grain).
  select case (bcdef)
  case (401:420, 501:502)
    isInj = .true.
  case default
    isInj = .false.
  end select

end function isInj

subroutine initialize_fields(sourceblock,eulerblock,srcSwitch,eulSwitch)
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  use IGLOO_variables,  only: nb
  use IGLOO_data_block, only: obj_sourceblock, obj_eulerblock
  implicit none
  logical,               intent(in)    :: srcSwitch, eulSwitch
  type(obj_sourceblock), intent(inout) :: sourceblock(nb)
  type(obj_eulerblock),  intent(inout) :: eulerblock(:,:)
  integer :: b, fam

  do b = 1, nb
    if (srcSwitch) then
      sourceblock(b)%sourceMass(:,:,:,:) = 0.0_R8
      sourceblock(b)%sourceMom(:,:,:,:)  = 0.0_R8
      sourceblock(b)%sourceEn(:,:,:)     = 0.0_R8
    endif
    if (eulSwitch) then
      do fam = 1, size(eulerblock(b,:))
        eulerblock(b,fam)%density(:,:,:)     = 0.0_R8
        eulerblock(b,fam)%velocity(:,:,:,:)  = 0.0_R8
        eulerblock(b,fam)%temperature(:,:,:) = 0.0_R8
        eulerblock(b,fam)%np(:,:,:)          = 0.0_R8
      enddo
    endif
  enddo

end subroutine initialize_fields

end module IGLOO_IC
