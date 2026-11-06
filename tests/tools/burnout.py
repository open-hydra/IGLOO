"""The burnout rule as the oracles state it: a consuming parcel (models 2, 4, 5, 8) is consumed at the first accepted step
whose end has its mass slot at or below m_b = min(1e-6 m0, 1e-20 kg), m0 the drop's reference mass (its injection mass,
a child's birth mass, a KH-RT drop's mass at its last event; model 8 the whole particle against the metal)."""
import math

EPS_BURN = 1.0e-6        # the relative bound, a fraction of m0
M_BURN_ABS = 1.0e-20     # the absolute bound [kg]


def burn_mass(m0):
    """m_b [kg] for a reference mass m0 [kg]."""
    return min(EPS_BURN * m0, M_BURN_ABS)


def burn_diam(m0, rho):
    """The diameter [m] of a sphere of density rho [kg/m^3] holding m_b."""
    return (6.0 * burn_mass(m0) / (math.pi * rho)) ** (1.0 / 3.0)
