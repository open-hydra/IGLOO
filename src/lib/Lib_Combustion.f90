!> Beckstead d^n burn law for metal particles; the oxide cap's exposed-surface factor and stretch.
module IGLOO_Lib_Combustion
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    use IGLOO_variables, only : pi
    implicit none
    private
    public :: becksteadRate
    public :: buildLobeTable, lobeF, capStretch, capStretchX0
    public :: buildXiLaw, xiGallier

    !> Metal parameter slots; the order is mirrored by Lib_RHS ind_mb../nmetal.
    integer(I4), parameter, public :: nmp = 13
    integer(I4), parameter, public :: imKb=1, imN=2, imXe=3, imBeta=4, imXi=5
    integer(I4), parameter, public :: imTig=6, imHf=7, imCps=8, imQc=9
    integer(I4), parameter, public :: imRr=10, imCr=11, imCs=12   !> the cap: product density, heat capacity, stretch
    integer(I4), parameter, public :: imXl=13                     !> the cap's retained fraction: 0 xi-cap, 1 Gallier21

    !> X_eff exponent of the Beckstead (2005) correlation.
    real(R8), parameter, public :: aXeff = 1.0_R8

    !> Oxide formed per kg of aluminium burned, M(Al2O3)/(2 M(Al)).
    real(R8), parameter, public :: nuOx = (2._R8*26.9815385_R8 + 3._R8*15.999_R8)/(2._R8*26.9815385_R8)

    !> The lobe temperature of the exposed-surface factor [K] (Gallier et al. 2021, 1 atm).
    real(R8), parameter, public :: TLobe = 2570._R8
    integer(I4), parameter, public :: NLOBE = 513     !> lobe table nodes, uniform in the lobe angle
    integer(I4), parameter, public :: NSTRETCH = 4096 !> midpoint cells of the stretch
    !> The lobe table in increasing f_ox: f_ox, (1 - f_ox)^(1/3) and F; set once by buildLobeTable.
    real(R8) :: lobeFr(NLOBE) = 0._R8, lobeUr(NLOBE) = 0._R8, lobeFv(NLOBE) = 0._R8
    logical  :: lobeBuilt = .false.

    !> Gallier et al. 2021's deposition law, Figs. 11 and 7 (tests/combustion/burn-box-cap-gallier/reference/): h the
    !  retained fraction against f_ox at 70 um, g the regression line through f_dep [%] against d_Al [um], at 70 um 1.
    integer(I4), parameter :: NXH = 6, NXG = 6
    real(R8), parameter :: xhF(NXH) = [0._R8, 0.1_R8, 0.3_R8, 0.5_R8, 0.7_R8, 0.9_R8]
    real(R8), parameter :: xhV(NXH) = [0.069522_R8, 0.101342_R8, 0.150355_R8, 0.192828_R8, 0.280447_R8, 0.520769_R8]
    real(R8), parameter :: xgD(NXG) = [20._R8, 50._R8, 70._R8, 100._R8, 220._R8, 400._R8]
    real(R8), parameter :: xgP(NXG) = [7.8595_R8, 7.3922_R8, 6.9633_R8, 6.1459_R8, 3.2639_R8, 1.6087_R8]
    !> g = 1 + xgS (d - 70), d in um held in [20, 400]; set once by buildXiLaw
    real(R8) :: xgS = 0._R8
    logical  :: xiLawBuilt = .false.

contains

    !> Beckstead d^n burn law: d^n(t) = d0^n - Keff*t => mdot = -(rho*pi*Keff/(2n))*d^(3-n),
    !  Keff = K-burn*X-eff; inert (mdot = 0) below T-ign. Ref: Beckstead, CESW 41(5):533 (2005).
    pure subroutine becksteadRate(Tp, dp, rhop, mp, mdot, qrel)
        implicit none
        real(R8), intent(in)  :: Tp, dp, rhop
        real(R8), intent(in)  :: mp(nmp)
        real(R8), intent(out) :: mdot, qrel

        mdot = 0._R8
        qrel = 0._R8
        if (Tp < mp(imTig)) return   ! ignition gate
        if (dp <= 0._R8) return      ! burnout/degenerate guard
        mdot = -(rhop * pi * mp(imKb) * mp(imXe)**aXeff / (2._R8 * mp(imN))) &
               * dp**(3._R8 - mp(imN))
        qrel = -mp(imBeta) * mp(imQc) * mdot   ! heat release to the particle [W]

    end subroutine becksteadRate

    !> The exposed-surface factor F(f_ox) of a drop carrying its oxide as one lobe, at TLobe with Gallier's
    !  property laws: NLOBE nodes uniform in the lobe angle alpha, stored in increasing f_ox; once per run.
    !  Ref: Gallier et al., Combust. Flame 228:142 (2021), eqs. 9-14.
    subroutine buildLobeTable()
        implicit none
        real(R8) :: rAl, rOx, sAl, sOx, sAO, th1, th2, delta, alo, ahi, a, g, h
        integer  :: k, j

        if (lobeBuilt) return
        rAl = 3226._R8 - 0.604_R8*TLobe
        rOx = 3056._R8 - 0.97_R8*(TLobe - 2327._R8)
        sAl = 1.267_R8 - 2.6e-4_R8*(TLobe - 933._R8)
        sOx = 0.66_R8 - 6.0e-5_R8*(TLobe - 2327._R8)
        sAO = 0.687_R8 - 1.6e-4_R8*(TLobe - 933._R8)
        th1 = acos((sOx**2 - sAO**2 - sAl**2)/(2._R8*sAO*sAl))
        th2 = acos((sAl**2 - sAO**2 - sOx**2)/(2._R8*sAO*sOx))
        delta = th1 + th2 - pi
        alo = delta - 0.5_R8*pi        ! the pole: all oxide
        ahi = 0.5_R8*pi                ! no oxide
        do k = 0, NLOBE - 1
            j = NLOBE - k
            if (k == 0) then
                a = alo + 1.e-12_R8*(ahi - alo)
            elseif (k == NLOBE - 1) then
                a = ahi
            else
                a = alo + (ahi - alo)*real(k, R8)/real(NLOBE - 1, R8)
            endif
            g = lobeG(a, th1)
            if (k == 0) then
                lobeFr(j) = 1._R8
            elseif (k == NLOBE - 1) then
                lobeFr(j) = 0._R8
            else
                h = lobeH(a, th1, delta)
                lobeFr(j) = rOx*h/(rAl*g + rOx*h)
            endif
            if (k == NLOBE - 1) then
                lobeFv(j) = 1._R8
            else
                lobeFv(j) = 2._R8*(1._R8 + sin(a))/(2._R8*g)**(2._R8/3._R8)
            endif
            lobeUr(j) = (1._R8 - lobeFr(j))**(1._R8/3._R8)
        enddo
        lobeBuilt = .true.

    end subroutine buildLobeTable

    !> 3 V_Al/(pi R^3) of the lobe at alpha (Gallier eq. 12).
    pure function lobeG(a, th1) result(g)
        implicit none
        real(R8), intent(in) :: a, th1
        real(R8) :: g, b

        b = th1 - a
        g = (1._R8 + sin(a))**2*(2._R8 - sin(a)) - cos(a)**3/cos(b)**3*(2._R8 + sin(b))*(1._R8 - sin(b))**2

    end function lobeG

    !> 3 V_ox/(pi R^3) of the lobe at alpha (Gallier eq. 13).
    pure function lobeH(a, th1, delta) result(h)
        implicit none
        real(R8), intent(in) :: a, th1, delta
        real(R8) :: h, b, e

        b = th1 - a
        e = a - delta
        h = cos(a)**3/cos(e)**3*(1._R8 - sin(e))**2*(2._R8 + sin(e)) &
          + cos(a)**3/cos(b)**3*(2._R8 + sin(b))*(1._R8 - sin(b))**2

    end function lobeH

    !> F(f_ox) from the lobe table: the cell by bisection on f_ox, the weight linear in u = (1 - f_ox)^(1/3).
    pure function lobeF(fox) result(Fx)
        implicit none
        real(R8), intent(in) :: fox
        real(R8) :: Fx, fc, u, w
        integer  :: lo, hi, mid

        fc = min(max(fox, 0._R8), 1._R8)
        lo = 1; hi = NLOBE
        do while (hi - lo > 1)
            mid = (lo + hi)/2
            if (lobeFr(mid) <= fc) then; lo = mid; else; hi = mid; endif
        enddo
        u  = (1._R8 - fc)**(1._R8/3._R8)
        w  = (u - lobeUr(lo))/(lobeUr(lo+1) - lobeUr(lo))
        Fx = (1._R8 - w)*lobeFv(lo) + w*lobeFv(lo+1)

    end function lobeF

    !> The stretch c(xi, n) = int_0^1 ds/F(f_ox(s)), s = (d_Al/d0)^n, which keeps Beckstead's burn time d0^n/K_eff
    !  under F; midpoint rule on NSTRETCH cells over lobeF.
    pure function capStretch(xi, n) result(c)
        implicit none
        real(R8), intent(in) :: xi, n
        real(R8) :: c, s, x, fox
        integer  :: k

        c = 0._R8
        do k = 0, NSTRETCH - 1
            s   = (real(k, R8) + 0.5_R8)/real(NSTRETCH, R8)
            x   = s**(3._R8/n)
            fox = xi*nuOx*(1._R8 - x)/(x + xi*nuOx*(1._R8 - x))
            c   = c + 1._R8/lobeF(fox)
        enddo
        c = c/real(NSTRETCH, R8)

    end function capStretch

    !> The stretch with an initial cap: c = (d_Al0/d_p0)^n int_0^1 ds/F(f_ox(s)), so t_b = d_p0^n/K_eff.
    pure function capStretchX0(xi, n, x0, rhoAl, rhoRes) result(c)
        implicit none
        real(R8), intent(in) :: xi, n, x0, rhoAl, rhoRes
        real(R8) :: c, s, x, fox, r0
        integer  :: k

        r0 = x0/(1._R8 - x0)
        c  = 0._R8
        do k = 0, NSTRETCH - 1
            s   = (real(k, R8) + 0.5_R8)/real(NSTRETCH, R8)
            x   = s**(3._R8/n)
            fox = (r0 + xi*nuOx*(1._R8 - x))/(x + r0 + xi*nuOx*(1._R8 - x))
            c   = c + 1._R8/lobeF(fox)
        enddo
        c = c/real(NSTRETCH, R8)
        c = c*((1._R8 - x0)/rhoAl/((1._R8 - x0)/rhoAl + x0/rhoRes))**(n/3._R8)

    end function capStretchX0

    !> g's slope from the least-squares line f = a + b d through Fig. 7's points: g = f(d)/f(70) = 1 + b (d - 70)/f(70).
    subroutine buildXiLaw()
        implicit none
        real(R8) :: sx, sy, sxx, sxy, a, b, n

        if (xiLawBuilt) return
        n   = real(NXG, R8)
        sx  = sum(xgD); sy = sum(xgP); sxx = sum(xgD*xgD); sxy = sum(xgD*xgP)
        b   = (n*sxy - sx*sy)/(n*sxx - sx*sx)
        a   = (sy - b*sx)/n
        xgS = b/(a + 70._R8*b)
        xiLawBuilt = .true.

    end subroutine buildXiLaw

    !> The retained fraction xi = g(d_Al) h(f_ox) in [0, 1]: g at the metal diameter dAl [m], held at its 20 and 400 um
    !  values outside; h linear between Fig. 11's points, held past f_ox = 0.9.
    pure function xiGallier(dAl, fox) result(xi)
        implicit none
        real(R8), intent(in) :: dAl, fox
        real(R8) :: xi, d, g, h, w
        integer  :: k

        d = min(max(dAl*1.e6_R8, xgD(1)), xgD(NXG))
        g = 1._R8 + xgS*(d - 70._R8)
        if (fox <= xhF(1)) then
            h = xhV(1)
        elseif (fox >= xhF(NXH)) then
            h = xhV(NXH)
        else
            k = 1
            do while (fox > xhF(k+1))
                k = k + 1
            enddo
            w = (fox - xhF(k))/(xhF(k+1) - xhF(k))
            h = xhV(k) + w*(xhV(k+1) - xhV(k))
        endif
        xi = min(max(g*h, 0._R8), 1._R8)

    end function xiGallier

end module IGLOO_Lib_Combustion
