#!/usr/bin/env python3
"""The boiling-zgr-* oracles: one water drop at Re = 0 (kv = 1) in the uniform box gas, integrated independently of
IGLOO with scipy's solve_ivp.

Above the boiling temperature at the local pressure, T_B(p) = 1/(1/T_boil - ln(p/p_atm) R_u/(L_v M_v)), the
Zuo-Gomes-Rutland branch (boiling = ZGR) with the vapour leaving at T_B:
    dT = T - T_B,  dT_f = max(dT, 0.5 K)
    alpha_s = 760 dT_f^0.26 (dT_f < 5), 27 dT_f^2.33 (< 25), 13800 dT_f^0.39   [W m-2 K-1]  (LiquidEvaporationBoil.C)
    G_f = alpha_s dT_f pi d^2/L_v;  A = c_p,g (T_g - T_B)/L_v;  T_f = T_B + (T_g - T_B)/3, k_f = k_g (T_f/T_g)^0.7
    B = pi k_f d Nu/c_p,g;  G (1 + G_f/G) = B ln(1 + A (1 + G_f/G))  (G = 0 for A <= 0);  S = G + G_f
    Q_s = S c_p,g (T_g - T_B)/(exp(S/B) - 1);  m c dT/dt = Q_s - S (L_v - c dT);  dm/dt = -S
Below T_B, the classical model (CEM, Sh = Nu = 2, Le = 1, Y_inf = 0, the clamp's x_s cap) with the vapour leaving at the
mean temperature: m c dT/dt = pi d k_g Nu (T_g - T) + mdot L_v, mdot = -pi d (k_g/c_p,g) Sh ln(1 + B_M).
d = (6 m/(pi rho_l))^(1/3).  Gas: rho 1.2, R 287, gamma 1.4, k 0.026 (tools/make_box_case.py); liquid rho 958, c 4186.
"""
import math
import sys

try:
    from scipy.integrate import solve_ivp
    from scipy.optimize import brentq
except ImportError:            # a missing dependency is a failure, never a skip
    print("[FAIL] boiling-zgr oracle needs scipy (scipy.integrate.solve_ivp, scipy.optimize.brentq)")
    sys.exit(1)

RU, PATM, XSCAP = 8314.46, 101325.0, 1e-12
RHOG, RG, GAM, KG = 1.2, 287.0, 1.4, 0.026
CPG = GAM*RG/(GAM - 1.0)
MG = RU/RG
RHOL, CL, LV, MV = 958.0, 4186.0, 2.26e6, 18.0
NU = 2.0
MBURN = 1e-15                  # IGLOO's burnout mass [kg]
EPS = 2.220446049250313e-16


def t_boil(p, tboil):
    return 1.0/(1.0/tboil - math.log(p/PATM)/(LV*MV/RU))


def alpha_s(dtf):
    if dtf < 5.0:
        return 760.0*dtf**0.26
    if dtf < 25.0:
        return 27.0*dtf**2.33
    return 13800.0*dtf**0.39


def diam(m):
    return (6.0*max(m, 1e-300)/(math.pi*RHOL))**(1.0/3.0)


def heat_fed(gf, a, b):
    if a <= 0.0:
        return 0.0
    h = lambda lx: gf*(1.0 + math.exp(lx))/math.exp(lx) - b*math.log1p(a*(1.0 + math.exp(lx)))
    return gf/math.exp(brentq(h, -60.0, 60.0, xtol=1e-15, rtol=4*EPS, maxiter=500))


def zgr_rate(T, m, tb, tg):
    """m c dT/dt and the efflux S of the boiling branch."""
    d = diam(m)
    dt = T - tb
    dtf = max(dt, 0.5)
    gf = alpha_s(dtf)*dtf*math.pi*d*d/LV
    a = CPG*(tg - tb)/LV
    kf = KG*((tb + (tg - tb)/3.0)/tg)**0.7
    b = math.pi*kf*d*NU/CPG
    s = heat_fed(gf, a, b) + gf
    qs = s*CPG*(tg - tb)/math.expm1(s/b)
    return qs - s*(LV - CL*dt), s


def cem_rate(T, m, tg, tboil):
    """m c dT/dt and mdot of the classical model below T_B."""
    d = diam(m)
    p = RHOG*RG*tg
    psat = PATM*math.exp(-(LV*MV/RU)*(1.0/T - 1.0/tboil))
    xs = 1.0 - XSCAP if psat >= p else min(psat/p, 1.0 - XSCAP)
    ys = xs*MV/(xs*MV + (1.0 - xs)*MG)
    mdot = -math.pi*d*(KG/CPG)*NU*math.log1p(ys/(1.0 - ys))
    return math.pi*d*KG*NU*(tg - T) + mdot*LV, mdot


class Drop:
    """The drop injected at T0, d0, carried at u = u_g = U; .state(t) -> (T, m), .slope(t) -> (dT/dt, dm/dt)."""

    def __init__(self, T0, d0, U, tg, tboil, tend):
        self.tg, self.tboil, self.U = tg, tboil, U
        self.tb = t_boil(RHOG*RG*tg, tboil)
        self.m0 = RHOL*math.pi/6.0*d0**3
        self.tc = None                       # the time T crosses T_B (none if it does not before tend)
        self.tburn = None                    # the burnout time (m = 1e-15 kg)

        def rhs_zgr(t, y):
            e, s = zgr_rate(y[0], y[1]*self.m0, self.tb, tg)
            return [e/(y[1]*self.m0*CL), -s/self.m0]

        def rhs_cem(t, y):
            e, md = cem_rate(y[0], max(y[1], 1e-300)*self.m0, tg, tboil)
            return [e/(max(y[1], 1e-300)*self.m0*CL), md/self.m0]

        self.rhs_zgr, self.rhs_cem = rhs_zgr, rhs_cem
        cross = lambda t, y: y[0] - self.tb
        cross.terminal, cross.direction = True, -1
        s1 = solve_ivp(rhs_zgr, (0.0, tend), [T0, 1.0], method='LSODA', rtol=1e-12, atol=[1e-12, 1e-15],
                       dense_output=True, events=cross)
        self.s1, self.s2 = s1, None
        if len(s1.t_events[0]):
            self.tc = s1.t_events[0][0]
            burn = lambda t, y: y[1]*self.m0 - MBURN
            burn.terminal = True
            s2 = solve_ivp(rhs_cem, (self.tc, tend), list(s1.y_events[0][0]), method='Radau', rtol=1e-10,
                           atol=[1e-8, 1e-14], dense_output=True, events=burn)
            self.s2 = s2
            if len(s2.t_events[0]):
                self.tburn = s2.t_events[0][0]

    def state(self, t):
        if self.tc is None or t <= self.tc:
            T, mu = self.s1.sol(t)
        else:
            T, mu = self.s2.sol(t)
        return T, mu*self.m0

    def slope(self, t):
        T, m = self.state(t)
        f = self.rhs_zgr if (self.tc is None or t <= self.tc) else self.rhs_cem
        dT, dmu = f(t, [T, m/self.m0])
        return dT, dmu*self.m0


def e13_half_ulp(v):
    """Half the last digit of an E13.6 print of v (0.dddddd E+xx)."""
    if v == 0.0:
        return 0.5e-6*1e-300
    return 0.5e-6*10.0**(math.floor(math.log10(abs(v))) + 1)


def log_count(path, text):
    try:
        return sum(1 for line in open(path) if text in line)
    except FileNotFoundError:
        return -1
