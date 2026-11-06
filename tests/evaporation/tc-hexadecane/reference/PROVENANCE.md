# TC2012 Fig. 11 digitized reference

Digitized (WebPlotDigitizer) from **Tonini, S.; Cossali, G. E., "An analytical model of
liquid drop evaporation in gaseous environment," Int. J. Thermal Sciences 57 (2012) 45-53**,
DOI 10.1016/j.ijthermalsci.2012.01.017, **Figure 11** — non-dimensional drop size and
temperature profiles for an n-hexadecane drop (R_d0=10 um, T_d0=300 K, T_inf=600 K,
vapour-free air), the **present-model (TC)** curves.

- `tc2012_fig11_d2.csv`: (tau, R^2/R_d0^2), tau = t*Dv/R_d0^2 (non-dimensional time).
- `tc2012_fig11_T.csv`:  (tau, T/T_d0).

Non-gating overlay only (verify.py). The tight gate (check.py) is IGLOO-vs-the-TC kernel.

## Case property provenance
IGLOO's TC kernel is TC2012's **present model** (eq. 16) -- proven: its small-rate limit is
`m_hat=(P_vs-chi)/[(T~s+1)/2]`, eq. 16's stated form (not Stefan-Fuchs eq. 2b). The case
properties are sourced as follows (not fitted to the overlay):

- **p_sat(T)** -- **TC2012's OWN Table 1** (paper's own data) gives n-hexadecane's saturation
  curve as `(T, Pvs/PN)` anchors
  `[(473.86,0.1),(510.11,0.3),(529.58,0.5),(543.5,0.7),(554.52,0.9),(560,1.0)]*Patm`. The `Psat`
  column of `INPUT/properties.dat` is that curve on the table's 1-K nodes: `ln p` linear in `1/T`
  between consecutive anchors, the end segments continued beyond them, written `%.6e`
  (`tests/tools/make_psat_table.py --anchors 473.86 0.1 510.11 0.3 529.58 0.5 543.5 0.7 554.52 0.9
  560 1.0 --merge INPUT/properties.dat`). Exact at every anchor; 0 (double underflow) on the rows below 10 K.
  The segments' effective latent heats (`Lv = slope*Ru/Mv`):

  | segment [K] | `Lv*Mv/Ru` [K] | `Lv` [kJ/kg] |
  |---|---|---|
  | 473.86-510.11 | 7325.7 | 269.0 |
  | 510.11-529.58 | 7087.6 | 260.2 |
  | 529.58-543.5 | 6957.3 | 255.4 |
  | 543.5-554.52 | 6873.1 | 252.4 |
  | 554.52-560 | 5970.4 | 219.2 (short segment; Table-1 rounding, kept: it fixes psat(560)=Patm) |

  A single Clausius-Clapeyron fit through the anchors (anchored at boiling) has `Lv=2.58e5` and
  reproduces them to <=2.5 %; at 490 K it agrees with the table to 0.1 %.
- **Lv = 2.2695e5 J/kg** -- Table 1's `dHv(Tb) = 226.95 kJ/kg`, the latent heat at the boiling
  point: the energy sink. It does not set the saturation curve (the column does).
- **cp_l(T)** (the `Cp` column of `INPUT/properties.dat`, varying; the `Enthalpy` column its trapezoid
  sum from h(1 K) = cp(1 K) * 1 K, a relative datum with hOff = 0) -- liquid n-hexadecane at 1 atm from
  the Helmholtz equation of state of **R. Romeo, E. W. Lemmon, "Equations of State for n-Hexadecane
  and n-Docosane", Int. J. Thermophys. 43 (2022), doi:10.1007/s10765-022-03059-0** (open NIST
  manuscript, https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=934414): the ideal part, eqs. (9)
  and (10) with Table 1, the residual part, eq. (12) with Table 3, Tc = 722.1 K, rho_c = 1.000 mol/dm^3,
  M = 226.441 g/mol. `tests/tools/hexadecane_eos.py --merge INPUT/properties.dat` writes it on the
  integer kelvins of 292..560 K (the triple point is 291.33 K, TC2012's boiling point 560 K) and holds
  the end values outside: 2216 J/kg/K at 300 K, 2538 at 400 K, 2739 at 454 K, 2874 at 490 K, 3135 at
  560 K. The equation's stated uncertainty in cp is 0.25 %; the cp data it was fitted to reach 454 K
  (Petit and Ter Minassian, 298-454 K, its Table 6), so the column above 454 K, the drop's plateau
  included, is the equation extrapolated. Checks of the implementation: its liquid density at
  298.15 K, 770.13 kg/m^3; its analytic derivatives against finite differences, 7e-7.
  TC2012 takes its properties constant (p. 47: "Fluid properties (Dv, mu, k, c) are assumed to be
  constant") and among the liquid's names only the density as temperature dependent (p. 52); the
  column is the liquid, not TC2012's model. Against the constant 2800 J/kg/K (NIST/Chemeo values at
  298 and 400 K extrapolated to 490 K) the replica's plateau stays 493.69 K, its heat-up fraction goes
  from 0.600 to 0.591 (the figure's 0.637) and its lifetime is 3.0 % shorter: measured, not tuned.
- **rho_l(T)** -- Table-1 boiling anchor 569.9 @ 560 K, linear `767-0.758*(T-300)`. The equation of
  state above gives 633.0 kg/m^3 at 490 K against this law's 623.0 (1.6 %); the column stays on
  TC2012's anchors.

## Gas state and gas-film properties
TC2012 gives no gas property values (p. 47 above; Table 1 lists the fuels only) and runs its cases at
the atmospheric pressure (p. 48, p. 50, Fig. 12(b)); the Fig. 11 caption names no pressure, read as
1 atm. The properties are those of air taken at the film temperature T_r = T_s + (T_inf - T_s)/3
(Hubbard, Denny & Mills 1975, the rule IGLOO's CEM-B uses) at the plateau, T_s = 493.69 K:
T_r = 529.1 K. Only the temperature half of the 1/3 rule is applied: the film's composition, the
vapour's share of the air-vapour mixture, is not.

- **rho_inf = 0.58841 kg/m^3** -- p_atm/(R T_inf) = 101325/(287*600) (solfile `Roi`): IGLOO's
  pressure is rho R T.
- **k_g = 0.04256 W/m/K, c_pg = 1035.8 J/kg/K** -- air at 1 atm, Incropera & DeWitt, Fundamentals of
  Heat and Mass Transfer, Table A.4: k = 40.7 and 43.9 mW/m/K, c_p = 1.030 and 1.040 kJ/kg/K at 500 and
  550 K, linear between (solfile `KL`, and `GAM` = c_pg/(c_pg - R) = 1.38328).
- **D_v = 1.2415e-5 m^2/s** -- n-hexadecane in air at T_r and 1 atm by the Fuller-Schettler-Giddings
  correlation, D = 1e-7 T^1.75 sqrt(1/M_A + 1/M_B)/(p_atm (V_A^(1/3) + V_B^(1/3))^2), M_A = 226.45,
  M_B = 28.97, the atomic diffusion volumes C 15.9, H 2.31 (V_A = 16*15.9 + 34*2.31 = 332.94) and air
  V_B = 19.7 (Poling, Prausnitz & O'Connell 2001, eq. 11-4.4 and Table 11-1; the form with p in bar
  and the constant 0.00143 gives 0.2 % less). Its error for organic vapours reaches about 15 % (Tang
  et al., Atmos. Chem. Phys. 14, 9233, 2014), about 0.8 in Le and 4 K in the plateau.
- **Le = 5.625** -- k_g/(rho_inf c_pg D_v): the product rho_inf D_v of TC2012's eq. 1
  (`mdot = 4 pi R_d D_v rho_inf m_hat`), which IGLOO's TC kernel forms from Le (D_v = k_g/(rho_g c_pg Le)).
- **cpv = 2582 J/kg/K** -- n-hexadecane vapour at T_r, Burcat's C16H34 NASA-7 fit (note
  "HexadecanS5/01", `utils/ATLAS/database/thermo/burcat.yaml` of hydra's ATLAS): 584.7 J/mol/K.
- mu_g stays 1.8e-5 Pa s: at Re = 0 it enters no rate.

With these the run's plateau is 493.69 K against the digitized 493.72 K, and the 0-D replica of
the kernel (`zero_d.py`) reproduces it to the printed digit and the heating path to 3e-4 K (see
`INFO.md`): the replica checks the run against its own kernel and inputs, not the inputs against
TC2012. The figure's time scale is not reproduced: the run lives tau = t D_v/R_d0^2 = 512.9 against
the figure's 212.8.
