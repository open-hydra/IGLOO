!> Solidification of a molten droplet with supercooling and recalescence (ODE model 6).
module IGLOO_Lib_Solidification
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    private
    public :: solidPhaseAtInjection, nucleationJump, plateauRate, hSolid
    public :: eventValue, eventFunction, solidTransition

    !> Phases: liquid, undercooled liquid, freezing plateau at T-melt, solid.
    integer(I4), parameter, public :: phLiquid=0, phUndercooled=1, phPlateau=2, phSolid=3

contains

    !> Phase and frozen fraction at injection: solid at or below T-nuc, undercooled below T-melt, else liquid.
    pure subroutine solidPhaseAtInjection(Tp, Tmelt, Tnuc, phase, f)
        implicit none
        real(R8), intent(in)  :: Tp, Tmelt, Tnuc
        integer,  intent(out) :: phase
        real(R8), intent(out) :: f

        if (Tp <= Tnuc) then
            phase = phSolid;       f = 1._R8
        elseif (Tp < Tmelt) then
            phase = phUndercooled; f = 0._R8
        else
            phase = phLiquid;      f = 0._R8
        endif

    end subroutine solidPhaseAtInjection

    !> Adiabatic recalescence from Tev: back to T-melt with f0 = cpl*(Tmelt-Tev)/hFus when f0 < 1,
    !  else frozen whole at the temperature that conserves hSolid.
    pure subroutine nucleationJump(Tev, cpl, cpSol, Tmelt, hFus, Tafter, fAfter, phaseAfter)
        implicit none
        real(R8), intent(in)  :: Tev, cpl, cpSol, Tmelt, hFus
        real(R8), intent(out) :: Tafter, fAfter
        integer,  intent(out) :: phaseAfter
        real(R8) :: f0

        f0 = cpl*(Tmelt - Tev)/hFus
        if (f0 < 1._R8) then
            Tafter = Tmelt;                                   fAfter = f0;    phaseAfter = phPlateau
        else
            Tafter = Tmelt - (cpl*(Tmelt - Tev) - hFus)/cpSol; fAfter = 1._R8; phaseAfter = phSolid
        endif

    end subroutine nucleationJump

    !> Frozen-fraction rate on the plateau: m*hFus*df/dt = -QdotW (QdotW [W] > 0 heats the droplet).
    pure function plateauRate(QdotW, m, hFus) result(dfdt)
        implicit none
        real(R8), intent(in) :: QdotW, m, hFus
        real(R8) :: dfdt

        dfdt = -QdotW/(m*hFus)

    end function plateauRate

    !> Specific enthalpy with the latent heat of fusion; hOff is the table datum.
    pure function hSolid(T, f, phase, cpl, cpSol, hFus, Tmelt, hOff) result(h)
        implicit none
        real(R8), intent(in) :: T, f, cpl, cpSol, hFus, Tmelt, hOff
        integer,  intent(in) :: phase
        real(R8) :: h

        select case (phase)
        case (phPlateau); h = cpl*Tmelt + hOff - f*hFus
        case (phSolid);   h = cpl*Tmelt + hOff - hFus - cpSol*(Tmelt - T)
        case default;     h = cpl*T + hOff
        end select

    end function hSolid

    !> Event function of threshold `which` (1 nucleation, 2 fully frozen, 3 re-melted); positive before it.
    pure function eventValue(which, T, f, Tnuc) result(g)
        implicit none
        integer,  intent(in) :: which
        real(R8), intent(in) :: T, f, Tnuc
        real(R8) :: g

        select case (which)
        case (1);     g = T - Tnuc
        case (2);     g = 1._R8 - f
        case (3);     g = f
        case default; g = 1._R8
        end select

    end function eventValue

    !> Threshold that ends the current phase (0: none) and its event function at (T, f).
    pure subroutine eventFunction(phase, T, f, Tnuc, g, which)
        implicit none
        integer,  intent(in)  :: phase
        real(R8), intent(in)  :: T, f, Tnuc
        real(R8), intent(out) :: g
        integer,  intent(out) :: which

        select case (phase)
        case (phLiquid, phUndercooled); which = 1
        case (phPlateau)
            if (1._R8 - f <= f) then; which = 2; else; which = 3; endif
        case default;                   which = 0
        end select
        g = eventValue(which, T, f, Tnuc)

    end subroutine eventFunction

    !> Phase change `which` applied to (T, f, phase), conserving hSolid.
    pure subroutine solidTransition(which, cpl, cpSol, Tmelt, hFus, T, f, phase)
        implicit none
        integer,  intent(in)    :: which
        real(R8), intent(in)    :: cpl, cpSol, Tmelt, hFus
        real(R8), intent(inout) :: T, f
        integer,  intent(inout) :: phase
        real(R8) :: Tev

        select case (which)
        case (1)
            Tev = T
            call nucleationJump(Tev, cpl, cpSol, Tmelt, hFus, T, f, phase)
        case (2)
            T = Tmelt - (f - 1._R8)*hFus/cpSol; f = 1._R8; phase = phSolid
        case (3)
            T = Tmelt - f*hFus/cpl;             f = 0._R8; phase = phLiquid
        end select

    end subroutine solidTransition

end module IGLOO_Lib_Solidification
