module IGLOO_Lib_Combustion
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    use IGLOO_variables, only : pi
    implicit none
    private
    public :: becksteadRate

    !> Metal parameter array: contiguous layout (order mirrored by Lib_RHS ind_mb../nmetal)
    integer(I4), parameter, public :: nmp = 11
    integer(I4), parameter, public :: imKb=1, imN=2, imXe=3, imBeta=4, imXi=5
    integer(I4), parameter, public :: imTig=6, imTm=7, imHf=8, imTn=9, imCps=10, imQc=11

    !> Beckstead (2005) final correlation: t_b ~ 1/X_eff (linear exponent)
    real(R8), parameter, public :: aXeff = 1.0_R8

contains

    !===================================================================!
    !  Beckstead d^n burn law                                           !
    !  d^n(t) = d0^n - Keff*t  =>  mdot = -(rho*pi*Keff/(2n))*d^(3-n)   !
    !  Keff = K-burn * X-eff;  inert (mdot=0 exactly) below T-ign       !
    !  Ref: Beckstead, Combust. Explos. Shock Waves 41(5):533 (2005)    !
    !===================================================================!
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

end module IGLOO_Lib_Combustion
