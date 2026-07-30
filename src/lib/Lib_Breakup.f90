module IGLOO_Lib_Breakup
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    private
    public :: assign_breakup, assign_scaleFactor
    public :: breakupOde
    public :: breakupEvent

    !> Breakup parameters — variable-size array, set by assign_breakup per model
    integer,  public :: nbp = 0
    real(R8), allocatable, public :: bp(:)
    integer,  public :: bpMethod = 1          !> TAB sub-method selector
    real(R8), public :: bpScale  = 1._R8      !> TAB precomputed scale factor

    !> childState layout (fixed): [vel(3), diam, npdot]
    integer(I4), parameter, public :: nchild = 5

contains

    subroutine assign_scaleFactor()
        implicit none
        select case (bpMethod)
        case (2); bpScale = gamma(1._R8 + 2._R8/bp(4)) / gamma(1._R8 + 3._R8/bp(4))
        case default; bpScale = 1._R8
        end select
    end subroutine assign_scaleFactor

    subroutine assign_breakup(breakup_word,brkupSelect,brkupEqOde,brkupEvent,brkupHasChild)
        implicit none
        character(len=*), intent(in)  :: breakup_word
        integer,          intent(out) :: brkupSelect
        logical,          intent(out) :: brkupEqOde
        logical,          intent(out) :: brkupEvent
        logical,          intent(out) :: brkupHasChild

        brkupEqOde    = .false.
        brkupEvent    = .false.
        brkupHasChild = .false.
        select case (breakup_word)
        case('NoBreakup')
            brkupSelect = 0
        case('Pilch-Erdman')
            brkupSelect = 1
            brkupEqOde  = .true.
            ! nbp = 2   !> [Cd, B]
        case('Reitz-Diawakar')
            brkupSelect = 2
            brkupEqOde  = .true.
            ! nbp = 4   !> [WeBag, Cb, Cstrip, Cs]
        case('Reitz-KHRT')
            brkupSelect   = 3
            brkupEqOde    = .true.
            brkupEvent    = .true.
            brkupHasChild = .true.
        case('TAB')
            brkupSelect = 4
            brkupEvent  = .true.
            ! nbp = 4   !> [Comega, Cmu, WeCrit, nSpread]
        case('ETAB')
            brkupSelect = 5
            brkupEvent  = .true.
            ! nbp = 6   !> [k1, k2, WeCrit, WeTrans, Comega, Cmu]
        case default
            write(*,*)
            write(*,*)
            write(*,*) "Wrong breakup model input ---> "//trim(breakup_word)
            write(*,*) "Choose one of the following :"
            write(*,*) "- Pilch-Erdman "
            write(*,*) "- Reitz-Diawakar "
            write(*,*) "- Reitz-KHRT (Kevin-Helmholtz & Rayleigh-Taylor)"
            write(*,*) "- TAB        (Taylor Analogy Breakup)"
            write(*,*) "- ETAB       (Enhanced TAB)"
            write(*,*)
            stop
        end select

        ! if (allocated(bp)) deallocate(bp)
        ! if (nbp > 0) then; allocate(bp(nbp)); bp = 0._R8; endif

    end subroutine assign_breakup

!**********************************************************************************************************************!
!********************** MODELS WITH CONTINUOUS PARCEL MODIFICATION - FUNCTIONS IN THE ODE SYSTEM **********************!
!**********************************************************************************************************************!
    function breakupOde(dp,sigma,mup,rhop,npdot,rhog,vel,Re,time,acc,vp,var,brkupSelect,bp) result(npdotDot)
        use IGLOO_variables, only: toll
        implicit none
        real(R8), intent(in)    :: dp, sigma, mup, rhop, npdot, rhog, vel, Re, time, acc(3), vp(3)
        real(R8), intent(inout) :: var(2)
        integer,  intent(in)    :: brkupSelect
        real(R8), intent(in)    :: bp(:)
        real(R8) :: npdotDot

        if (vel<toll) then; npdotDot = 0._R8; return; endif
        select case (brkupSelect)
        case (1); npdotDot = PilchErdman(dp,sigma,mup,rhop,npdot,rhog,vel, bp)
        case (2); npdotDot = ReitzDiwakar(dp,sigma,rhop,npdot,rhog,vel,Re, bp)
        case (3); npdotDot = ReitzKHRT(dp,sigma,mup,rhop,npdot,rhog,vel,time,acc,vp,var(1),var(2), bp)
        end select

    end function breakupOde


    function PilchErdman(dp,sigma,mup,rhop,npdot,rhog,vel, bp) result(npdotDot)
        use IGLOO_variables, only: toll
        implicit none
        real(R8), intent(in) :: dp, sigma, mup, rhop, npdot, rhog, vel
        real(R8), intent(in) :: bp(:)   !> [Cd, B]
        real(R8) :: npdotDot
        real(R8) :: We, Wec, Oh, TBT, VdV, dStable, tau

        npdotDot = 0._R8
        We  = dp * rhog * vel**2 / sigma
        Oh  = mup / sqrt( rhop * dp * sigma )
        Wec = 12 * (1._R8 + 1.077_R8 * Oh**1.6_R8)

        TBT = huge(1._R8)
        if (We <= Wec) return

        if     (We > 2670._R8) then
            TBT = 5.5_R8
        elseif (We >  351._R8) then
            TBT = 0.766_R8*(We - 12._R8)**0.25_R8
        elseif (We >   45._R8) then
            TBT =  14.1_R8*(We - 12._R8)**(-0.25_R8)   !> A16 fix: PE87 45<We<=351 is NEGATIVE (Fig.7 continuous)
        elseif (We >   18._R8) then
            TBT =  2.45_R8*(We - 12._R8)**0.25_R8
        elseif (We >   12._R8) then
            TBT =    6._R8*(We - 12._R8)**(-0.25_R8)
        else
            return
        endif

        if (Oh>0.1_R8 .and. We<228._R8) TBT = 4.5_R8*(1._R8 + 1.2_R8*Oh**1.64_R8)
        VdV = sqrt(rhog/rhop) * (0.75_R8*bp(1)*TBT + 3._R8*bp(2)*TBT**2)
        dStable = Wec*sigma/(rhog*vel*vel*(1._R8-VdV)**2+toll)

        if (dp < dStable) return

        tau = TBT*dp/(vel*sqrt(rhog/rhop))
        npdotDot = - 3*npdot/dp * (dStable-dp)/tau

    end function PilchErdman

    function ReitzDiwakar(dp,sigma,rhop,npdot,rhog,vel,Re, bp) result(npdotDot)
        implicit none
        real(R8), intent(in) :: dp, sigma, rhop, npdot, rhog, vel, Re
        real(R8), intent(in) :: bp(:)   !> [WeBag, Cb, Cstrip, Cs]
        real(R8) :: npdotDot, radius, We, tau, dStable

        npdotDot = 0._R8
        radius = 0.5_R8*dp
        We  =  radius * rhog * vel**2 / sigma

        if (We > bp(1)) then
            if (We > bp(3)*sqrt(Re)) then
                tau     = bp(4)*0.5_R8*dp*sqrt(rhop/rhog)/vel
                dStable = (2._R8*bp(3)*sigma)**2*Re/(dp*rhog**2*vel**4)
            else
                tau     = bp(2)*dp*sqrt(0.0625_R8*rhop*dp/sigma)
                dStable = 2._R8*sigma*bp(1)/(rhog*vel**2)
            endif
            npdotDot = - 3*npdot/dp * (dStable-dp)/tau
        endif

    end function ReitzDiwakar

    function ReitzKHRT(dp,sigma,mup,rhop,npdot,rhog,vel,time,acc,vp,told,tc, bp) result(npdotDot)
        use IGLOO_variables, only: pi, toll
        implicit none
        real(R8), intent(in)    :: dp, sigma, mup, rhop, npdot, rhog, vel, time, acc(3), vp(3)
        real(R8), intent(inout) :: told, tc
        real(R8), intent(in)    :: bp(:)   !> [B0, B1, Ctau, CRT, mShedLim, WeLimit]
        real(R8) :: npdotDot
        real(R8) :: radius, WeGas, WePart, Oh, Tay, omegaKH, lambdaKH, tauKH, dStable
        real(R8) :: gt, force, omegaRT, lambdaRT, tauRT
        logical  :: noRTbreakup

        npdotDot = 0._R8
        radius  = 0.5_R8*dp
        WeGas  = radius * rhog * vel**2 / sigma
        WePart = WeGas/rhog * rhop
        !> Z = sqrt(We_l)/Re_l = mu_l/sqrt(rho_l sigma r); T = Z sqrt(We_g) (Reitz-87 Eqs 4-5)
        Oh  = mup/sqrt(rhop*sigma*radius)
        Tay = Oh*sqrt(WeGas)
        if (Oh>1 .or. rhog/rhop>0.1_R8) write(*,*) ' [WARNING] ReitzKHRT hypothesis falling'

        omegaKH  = (0.34_R8 + 0.38_R8*WeGas**1.5_R8)/((1._R8 + Oh)*(1._R8 + 1.4_R8*Tay**0.6_R8)) &
                    *sqrt(sigma/(rhop*radius**3))   ! Reitz 1987 KH growth rate [1/s]
        lambdaKH = radius*9.02_R8*(1._R8 + 0.45_R8*sqrt(Oh))*(1._R8 + 0.4_R8*Tay**0.7_R8)        &
                    /(1._R8 + 0.87_R8*WeGas**1.67_R8)**0.6_R8
        tauKH    = 3.726_R8*bp(2)*radius/(lambdaKH*omegaKH)
        dStable  = 2._R8*bp(1)*lambdaKH

        gt    = (acc(1)*vp(1)+acc(2)*vp(2)+acc(3)*vp(3))/norm2(vp)
        force = abs(gt*(rhog - rhop))
        omegaRT  = sqrt(2._R8*force**1.5_R8/(3._R8*sqrt(3._R8*sigma)*(rhop+rhog)))
        lambdaRT = 2._R8*pi*bp(4)/(sqrt(force/(3._R8*sigma))+toll)
        tauRT    = bp(3)/(omegaRT+toll)

        if ((tc > 0._R8) .or. (lambdaRT < dp)) then; tc = tc + time - told; told = time; endif

        noRTbreakup = .not.((tc > tauRT).and.(lambdaRT < dp))
        if (noRTbreakup.and.(dStable < dp).and.(WeGas > bp(6))) then
            npdotDot = - 3*npdot/dp * (dStable-dp)/tauKH
        endif

    end function ReitzKHRT

!**********************************************************************************************************************!
!********************** MODELS WITH PARCEL MODIFICATION "EVENTS" - SOUBROUTINE IN THE ODE SOLOUT **********************!
!**********************************************************************************************************************!
    subroutine breakupEvent(eventvar, neventvar, brkupState, nbrkst, &
                             sigma,mup,rhop,rhog,vel,Re,time,acc,vp,dt, &
                             brkupSelect,bp,bpMethod,bpScale, &
                             event,childState,addChild,exitLoop)
        use IGLOO_variables, only: toll
        implicit none
        integer,  intent(in)    :: neventvar, nbrkst
        real(R8), intent(inout) :: eventvar(neventvar)
        real(R8), intent(inout) :: brkupState(max(nbrkst,1))
        real(R8), intent(in)    :: sigma, mup, rhop, rhog, vel, Re, time, dt, acc(3)
        real(R8), intent(inout) :: vp(3)   !> ETAB kicks the parcel velocity (vperp) at the apply
        integer,  intent(in)    :: brkupSelect
        real(R8), intent(in)    :: bp(:)
        integer,  intent(in)    :: bpMethod
        real(R8), intent(in)    :: bpScale
        real(R8), intent(out)   :: childState(nchild)
        logical,  intent(out)   :: event, addChild
        logical,  intent(in)    :: exitLoop

        childState = 0._R8
        if (vel<toll) then; event = .false.; addChild = .false.; return; endif
        !> eventvar(1) = dp, eventvar(2) = npdot (always, for all models)
        select case (brkupSelect)
        case (3)   !> KHRT: eventvar=[dp,npdot,m0,KHidx], brkupState=[told,tc]
            call ReitzKHRTevent(eventvar(1),sigma,mup,rhop,eventvar(2),rhog,vel, &
                                time,acc,vp,eventvar(3),brkupState(1),brkupState(2),      &
                                eventvar(4), bp, event,childState,addChild,exitLoop)
        case (4)   !> TAB: eventvar=[dp,npdot,y,yDot]
            addChild = .false.
            call TABmodel(eventvar(1),sigma,mup,rhop,eventvar(2),rhog,vel,dt, &
                          eventvar(3),eventvar(4), bp,bpMethod,bpScale, event,exitLoop)
        case (5)   !> ETAB: eventvar=[dp,npdot,y,yDot]; vp kicked by vperp at the apply
            addChild = .false.
            call ETABmodel(eventvar(1),sigma,mup,rhop,eventvar(2),rhog,vel,Re,dt, &
                           eventvar(3),eventvar(4),vp, bp, event,exitLoop)
        end select

    end subroutine breakupEvent
    

    subroutine ReitzKHRTevent(dp,sigma,mup,rhop,npdot,rhog,vel,time,acc,vp,m0,told,tc,KHindex, &
                               bp, event,childState,addChild,exitLoop)
        use IGLOO_variables, only: pi, toll
        implicit none
        real(R8), intent(in)    :: sigma, mup, rhop, rhog, vel, time, acc(3), vp(3)
        real(R8), intent(inout) :: dp, npdot, m0, told, tc, KHindex
        real(R8), intent(in)    :: bp(:)   !> [B0, B1, Ctau, CRT, mShedLim, WeLimit]
        real(R8), intent(out)   :: childState(nchild)  !> [vel(3), diam, npdot]
        logical,  intent(out)   :: event, addChild
        logical,  intent(in)    :: exitLoop
        real(R8) :: radius, WeGas, WePart, Oh, Tay, omegaKH, lambdaKH, tauKH, dStable!, dOld
        real(R8) :: gt, force, omegaRT, lambdaRT, tauRT, nDrops, mShed!, lengthScale
        real(R8) :: ae3, be3, ce3, de3, qe3, pe3, d3, ue3, ve3, dParent, mc, nChild

        event    = .false.
        addChild = .false.
        radius   = 0.5_R8*dp
        !> We definitions based on the RADIUS
        WeGas  = radius * rhog * vel**2 / sigma
        WePart = WeGas/rhog * rhop
        !> Z/T per Reitz-87 Eqs 4-5
        Oh  = mup/sqrt(rhop*sigma*radius)
        Tay = Oh*sqrt(WeGas)
        if (Oh>1 .or. rhog/rhop>0.1_R8) write(*,*) ' [WARNING] ReitzKHRT hypothesis falling'

        omegaKH  = (0.34_R8 + 0.38_R8*WeGas**1.5_R8)/((1._R8 + Oh)*(1._R8 + 1.4_R8*Tay**0.6_R8)) &
                    *sqrt(sigma/(rhop*radius**3))   ! Reitz 1987 KH growth rate [1/s]
        lambdaKH = radius*9.02_R8*(1._R8 + 0.45_R8*sqrt(Oh))*(1._R8 + 0.4_R8*Tay**0.7_R8)        &
                    /(1._R8 + 0.87_R8*WeGas**1.67_R8)**0.6_R8
        tauKH    = 3.726_R8*bp(2)*radius/(lambdaKH*omegaKH)
        dStable  = 2._R8*bp(1)*lambdaKH

        gt    = (acc(1)*vp(1)+acc(2)*vp(2)+acc(3)*vp(3))/norm2(vp)
        force = abs(gt*(rhog - rhop))
        omegaRT  = sqrt(2._R8*force**1.5_R8/(3._R8*sqrt(3._R8*sigma)*(rhop+rhog)))
        lambdaRT = 2._R8*pi*bp(4)/(sqrt(force/(3._R8*sigma))+toll)
        tauRT    = bp(3)/(omegaRT+toll)

        if ((tc > 0._R8) .or. (lambdaRT < dp)) then; tc = tc + time - told; told = time; endif

        if ((tc > tauRT) .and.(lambdaRT < dp)) then
            tc = - 0._R8
            event  = .true.
            if (exitLoop) then
                !> children of diameter 2r_c = lambdaRT (Beale-Reitz Eq 11) => cubic count by mass conservation
                nDrops = (dp/lambdaRT)**3
                dp     = (dp**3/nDrops)**(1._R8/3._R8)
                npdot  = nDrops*npdot
                m0     = pi/6._R8*rhop*dp**3
            endif
        elseif (dStable < dp) then
            if (WeGas > bp(6)) then
                mShed = m0 - pi/6._R8*rhop*dp**3
                if (mShed/m0 > bp(5)) then
                    ae3 = 1._R8
                    be3 = -dStable
                    ce3 = 0._R8
                    de3 = dp*dp*(dStable - dp)
                    qe3 = (be3/(3._R8*ae3))**3 - be3*ce3/(6._R8*ae3*ae3) + de3/(2._R8*ae3);
                    pe3 = (3._R8*ae3*ce3 - be3*be3)/(9._R8*ae3*ae3);
                    d3  = qe3*qe3 + pe3*pe3*pe3;

                    if (d3 >= 0._R8) then
                        d3  = sqrt(d3)
                        ue3 = (-qe3 + d3)**(1._R8/3._R8)
                        ve3 = (-qe3 - d3)**(1._R8/3._R8)
                        dParent = ue3 + ve3 - be3/3._R8
                        mc  = npdot*(dp**3-dParent**3)
                        nChild = mc/dStable**3
                        if (nChild >= npdot) then
                            event       = .true.
                            if (exitLoop) then
                                addChild       = .true.
                                childState(1:3) = vp   ! child velocity = parent (no normal component in KHRT sheds)
                                childState(4)   = dStable
                                childState(5)   = nChild
                                dp = dParent
                                m0 = pi/6._R8*rhop*dp**3
                            endif
                        endif
                    endif
                endif
            endif
        elseif (KHindex < 0.5_R8) then
            ! lengthScale = min(lambdaKH, 2._R8*pi*vel/omegaKH)
            ! event = .true.
            ! if (exitLoop) then
            !     dOld  = dp
            !     dp    = (1.5_R8*dp*dp*lengthScale)**(1._R8/3._R8)
            !     npdot = npdot*(dOld/dp)**3
            !     m0    = pi/6._R8*rhop*dp**3
            !     KHindex = 1._R8
            ! endif
        endif

    end subroutine ReitzKHRTevent

    subroutine TABmodel(dp,sigma,mup,rhop,npdot,rhog,vel,dt,y,yDot, &
                         bp,method,scaleFactor, event,exitLoop)
        use IGLOO_variables,      only: pi, toll
        use IGLOO_Lib_Statistics, only: ChiSquare, RosinRammler, RonsinRammlerMoments
        implicit none
        real(R8), intent(in)    :: sigma, mup, rhop, rhog, vel, dt
        real(R8), intent(inout) :: dp, npdot, y, yDot
        real(R8), intent(in)    :: bp(:)        !> [Comega, Cmu, WeCrit, nSpread]
        integer,  intent(in)    :: method       !> TAB sub-method selector
        real(R8), intent(in)    :: scaleFactor   !> precomputed scale factor
        logical,  intent(out)   :: event
        logical,  intent(in)    :: exitLoop
        real(R8) :: twoPi, radius, rdt, omega, We, WeCr, y1, y2, a, phi, quad, tb
        real(R8) :: coste, theta, k, rs, rMin, rNew, rSize, nSpread, m3, m2, r32, dOld
        integer  :: n, iter
        integer, parameter :: maxIter=20

        event  = .false.
        twoPi  = 2._R8*pi
        radius = 0.5_R8*dp

        rdt   = 0.5_R8*bp(2)*mup/(rhop*radius**2)         !> inverse of characteristic dumping time
        omega = bp(1)*sigma/(rhop*radius**3) - rdt*rdt    !> oscillation frequency squared

        if (omega>0._R8) then
            omega = sqrt(omega)
            We    = radius * vel**2 / sigma * rhog
            WeCr  = We/bp(3)

            y1 = y - WeCr
            y2 = yDot/omega
            a  = sqrt(y1*y1+y2*y2)

            if (a+WeCr > 1._R8) then
                phi  = acos(min(max(y1/a,-1._R8), 1._R8))
                quad = -y2/a
                if (quad < 0._R8) phi = twoPi - phi
                tb = 0._R8

                if (abs(y) < 1._R8) then
                    coste = 1._R8
                    ! if (((WeCr-a)<-1._R8).and.(yDot<0._R8)) coste = -1._R8
                    theta = acos((coste-WeCr)/a)
                    if (theta<phi) then
                        if ((twoPi-theta)>=phi) theta = - theta
                        theta = twoPi+theta
                    endif
                    tb = (theta-phi)/omega
                endif

                if (dt>tb) then
                    event = .true.
                    if (exitLoop) then
                        yDot = - a * omega * sin(omega*tb + phi)
                        k    = (1._R8 + (4._R8/3._R8)+ rhop*radius**3/(8*sigma)*yDot**2)
                        rs   = radius/k
                        rNew = 0._R8
                        rMin = 0.01_R8*rs !> minimum value threshold
                        if     (method==1) then
                            !> Chi-Square distribution with n degrees of freedom from the original paper
                            !> [WARNING: number of degrees of freedom n has been chosen "arbitrarily"]
                            n = nint(k)
                            rSize = rs/real(n + 4)
                            do while (rNew<rMin .or. rNew>=radius)
                                rNew = ChiSquare(n)*rSize
                            enddo
                        elseif (method==2) then
                            !> Rosin-Rammler distribution with size parameter rs and spread parameter n=3.5
                            !> Most used probability density distribution (DEFAULT CHOICE)
                            nSpread = bp(4)
                            rSize   = rs*scaleFactor
                            do iter = 1, maxIter
                                ! truncated moments
                                m3 = RonsinRammlerMoments(3, rSize, nSpread, rMin, radius)
                                m2 = RonsinRammlerMoments(2, rSize, nSpread, rMin, radius)
                                ! resulting D32 at the current rSize
                                if (m2 > toll) then; r32 = m3/m2
                                else;                r32 = rs 
                                endif
                                ! convergence test
                                if (abs(r32-rs)/rs< 0.0001_R8) exit
                                ! fixed-point update of rSize
                                rSize = rSize * (rs / r32)
                            enddo
                            do while (rNew<rMin .or. rNew>=radius)
                                rNew = RosinRammler(rSize,nSpread)
                            enddo
                        elseif (method==3) then
                            !> LogNormal/Normal distributions (to be implemented)
                        endif
                        dOld  = dp
                        dp    = 2._R8*rNew
                        npdot = npdot*(dOld/dp)**3
                        y     = 0._R8
                        yDot  = 0._R8
                    endif
                else
                    call YupdateTAB(omega,dt,y1,rdt,WeCr,y,yDot)
                endif
            else
                call YupdateTAB(omega,dt,y1,rdt,WeCr,y,yDot)
            endif
        else
            y    = 0._R8
            yDot = 0._R8
        endif            

    end subroutine TABmodel

    subroutine ETABmodel(dp,sigma,mup,rhop,npdot,rhog,vel,Re,dt,y,yDot,vp, &
                          bp, event,exitLoop)
        use IGLOO_variables, only: pi, toll
        implicit none
        real(R8), intent(in)    :: sigma, mup, rhop, rhog, vel, Re, dt
        real(R8), intent(inout) :: dp, npdot, y, yDot, vp(3)
        real(R8), intent(in)    :: bp(:)   !> [k1, k2, WeCrit, WeTrans, Comega, Cmu]
        logical,  intent(out)   :: event
        logical,  intent(in)    :: exitLoop
        real(R8) :: AWe, twoPi, radius, rdt, omega, We, WeCr, y1, y2, a, phi, quad, tb
        real(R8) :: theta, rNew, sqrtWe, Kbr, dOld
        real(R8) :: Cd0, Asq, vperp, omega0, pHat(3), e1(3), e2(3), psi, vpn

        AWe = (bp(2)/bp(1)*sqrt(bp(4))-1._R8)/bp(4)**4
        event  = .false.
        twoPi  = 2._R8*pi
        radius = 0.5_R8*dp

        rdt   = 0.5_R8*bp(6)*mup/(rhop*radius**2)        !> inverse of characteristic dumping time
        omega = bp(5)*sigma/(rhop*radius**3) - rdt*rdt   !> oscillation frequency squared

        if (omega>0._R8) then
            omega = sqrt(omega)
            We    = radius * vel**2 / sigma * rhog
            WeCr  = We/bp(3)
            
            y1 = y - WeCr
            y2 = yDot/omega
            a  = sqrt(y1*y1+y2*y2)

            if (a+WeCr > 1._R8) then
                phi  = acos(min(max(y1/a,-1._R8), 1._R8))
                quad = -y2/a
                if (quad < 0._R8) phi = twoPi - phi
                tb = 0._R8

                if (abs(y) < 1._R8) then
                    theta = acos((1._R8-WeCr)/a)
                    if (theta<phi) then
                        if ((twoPi-theta)>=phi) theta = - theta
                        theta = twoPi+theta
                    endif
                    tb = (theta-phi)/omega
                    if (dt>tb) then
                        y    = 1._R8
                        yDot = - a * omega * sin(omega*tb + phi)
                    endif
                endif

                if (dt>tb) then
                    sqrtWe = AWe*We**4+1._R8
                    Kbr = bp(1)*omega*sqrtWe
                    !> stripping branch carries k2 (Tanner Eq.6); the AWe smoother reaches k2*sqrt(WeTrans)
                    !  at WeTrans by design, so k2 here keeps Kbr continuous for k1/=k2
                    if (We > bp(4)) Kbr = bp(2)*omega*sqrt(We)
                    rNew = radius*exp(-Kbr/omega*acos(min(max(1._R8-1._R8/WeCr,-1._R8),1._R8)))
                    if (rNew<radius) then
                        event = .true.
                        if (exitLoop) then
                            !> Product normal velocity vperp = A*xdot, xdot = a*ydot/2 (Tanner-97
                            !  Eqs 8-10): A^2 = 3(1 - a/rSMR + 5*Cd*We/72)*omega0^2/ydot^2 with the
                            !  undamped omega0^2 = Comega*sigma/(rho_l a^3) -- ydot cancels in A*xdot.
                            !  Cd = incompressible sphere estimate (energy-balance coefficient only);
                            !  direction: random azimuth in the plane normal to the parent path.
                            if (Re > 1000._R8) then; Cd0 = 0.43_R8
                            else;                    Cd0 = 24._R8/(Re+toll)*(1._R8+0.15_R8*Re**0.687_R8); endif
                            omega0 = sqrt(bp(5)*sigma/(rhop*radius**3))
                            Asq    = 3._R8*(1._R8 - radius/rNew + 5._R8*Cd0*We/72._R8)
                            vpn    = norm2(vp)
                            if (Asq > 0._R8 .and. vpn > toll) then
                                vperp = 0.5_R8*radius*omega0*sqrt(Asq)
                                pHat  = vp/vpn
                                !> basis normal to the path: Gram-Schmidt vs the least-aligned axis
                                e1 = 0._R8; e1(minloc(abs(pHat),1)) = 1._R8
                                e1 = e1 - dot_product(e1,pHat)*pHat
                                e1 = e1/norm2(e1)
                                e2 = [pHat(2)*e1(3)-pHat(3)*e1(2), &
                                      pHat(3)*e1(1)-pHat(1)*e1(3), &
                                      pHat(1)*e1(2)-pHat(2)*e1(1)]
                                call random_number(psi); psi = twoPi*psi
                                vp = vp + vperp*(cos(psi)*e1 + sin(psi)*e2)
                            endif
                            dOld  = dp
                            dp    = 2._R8*rNew
                            npdot = npdot*(dOld/dp)**3
                            y     = 0._R8
                            yDot  = 0._R8
                        endif
                    endif
                else
                    call YupdateTAB(omega,dt,y1,rdt,WeCr,y,yDot)
                endif
            else
                call YupdateTAB(omega,dt,y1,rdt,WeCr,y,yDot)
            endif
        else
            y    = 0._R8
            yDot = 0._R8
        endif

    end subroutine ETABmodel


    ! function/subroutine SHFmodel(...) result(npdotDot)
    ! end function/subroutine SHFmodel

!**********************************************************************************************************************!
!****************************************** FUNCTION USED BY BREAKUP MODELS  ******************************************!
!**********************************************************************************************************************!

    subroutine YupdateTAB(omega,dt,y1,rdt,WeCr,y,yDot)
        implicit none
        real(R8), intent(in)    :: omega, dt, y1, rdt, WeCr
        real(R8), intent(inout) :: y, yDot  
        real(R8) :: tempVar, sint, cost

        tempVar = omega*dt
        sint = sin(tempVar)
        cost = cos(tempVar)

        y    =     WeCr     + exp(-dt*rdt)*(       y1*cost + (yDot+y1*rdt)*sint/omega)
        yDot = (WeCr-y)*rdt + exp(-dt*rdt)*(-omega*y1*sint + (yDot+y1*rdt)*cost)

    endsubroutine YupdateTAB


end module IGLOO_Lib_Breakup