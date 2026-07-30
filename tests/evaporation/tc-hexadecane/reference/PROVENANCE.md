# TC2012 Fig. 11 digitized reference

Digitized (WebPlotDigitizer) from **Tonini, S.; Cossali, G. E., "An analytical model of
liquid drop evaporation in gaseous environment," Int. J. Thermal Sciences 57 (2012) 45-53**,
DOI 10.1016/j.ijthermalsci.2012.01.017, **Figure 11** — non-dimensional drop size and
temperature profiles for an n-hexadecane drop (R_d0=10 um, T_d0=300 K, T_inf=600 K,
vapour-free air), the **present-model (TC)** curves.

- `tc2012_fig11_d2.csv`: (tau, R^2/R_d0^2), tau = t*Dv/R_d0^2 (non-dimensional time).
- `tc2012_fig11_T.csv`:  (tau, T/T_d0).

Non-gating overlay only (verify.py). The tight gate (check.py) is IGLOO-vs-the-TC kernel.

## Case property provenance (Tier-1 reconstruction, 2026-07-28)
IGLOO's TC kernel is TC2012's **present model** (eq. 16) -- proven: its small-rate limit is
`m_hat=(P_vs-chi)/[(T~s+1)/2]`, eq. 16's stated form (not Stefan-Fuchs eq. 2b). The case
properties are sourced as follows (not fitted to the overlay):

- **psat / Lv = 2.58e5 J/kg** -- from **TC2012's OWN Table 1** (paper's own data): Table 1
  gives n-hexadecane's saturation curve as `(T, Pvs/PN)` anchors
  `[(473.86,0.1),(510.11,0.3),(529.58,0.5),(543.5,0.7),(554.52,0.9),(560,1.0)]*Patm`. A
  Clausius-Clapeyron fit through them (anchored at boiling) has effective `Lv=2.58e5` and
  reproduces the anchors to <=2.5%. This ~= Watson `Lv` at the ~490 K wet-bulb. The prior
  `Lv=2.9e5` made CC-psat run 20-24% BELOW Table 1. (Table-1 boiling `dHv(Tb)=226.95 kJ/kg`
  is the low-T end; the psat-slope value 258 is the operating-range effective latent heat.)
- **cp_l = 2800 J/kg/K** (constant, `INPUT/properties.dat` Cp column) -- **NIST WebBook /
  Chemeo** n-hexadecane LIQUID isobaric heat capacity: `Cp(298 K)=499-500 J/mol*K = 2205
  J/kg*K` (= the prior 2200, i.e. the 298 K value), `Csat(400 K)=572 J/mol*K = 2527`,
  linear-extrapolated to the ~490 K operating point `~= 635 J/mol*K = 2804 J/kg*K -> 2800`.
  TC2012 uses CONSTANT gas-film properties (Dv/mu/k/c; only rho_l is T-dependent), so a
  constant cp_l at the operating T is the faithful choice.
  Sources: NIST Chemistry WebBook (webbook.nist.gov, hexadecane C544763);
  Chemeo (chemeo.com/cid/30-657-9/Hexadecane); Tc=722-723 K.
- **rho_l(T)** -- Table-1 boiling anchor 569.9 @ 560 K, linear `767-0.758*(T-300)`.
- **Le = 2.5** -- n-hexadecane vapour (Dv-justified); kept physical.

Tier-1 reproduces the heating shape (heat-frac 0.536 vs Fig.11 0.540) and plateau to 0.7%
(490.2 vs ~493.6 K). Tier 2 (a psat curve DECOUPLED from the sink Lv, via the variable-psat
capability) closes the residual -- see `plan-bucket/tc-hexadecane-tier2-decoupled-psat.md`.
