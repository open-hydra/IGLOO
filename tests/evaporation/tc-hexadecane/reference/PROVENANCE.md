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
  point: the energy sink. It no longer sets the saturation curve (the column does).
- **cp_l = 2800 J/kg/K** (constant, `INPUT/properties.dat` Cp column) -- **NIST WebBook /
  Chemeo** n-hexadecane LIQUID isobaric heat capacity: `Cp(298 K)=499-500 J/mol*K = 2205
  J/kg*K`, `Csat(400 K)=572 J/mol*K = 2527`, linear-extrapolated to the ~490 K operating point
  `~= 635 J/mol*K = 2804 J/kg*K -> 2800`. Both sources list liquid Cp only up to 313 K besides
  that 400 K value, so no anchor supports a higher value at 490 K; 2900 would bring the 0-D
  heat-frac closer to the figure (0.629 against 0.618, figure 0.637) and is not used for that
  reason alone. TC2012 uses CONSTANT gas-film properties (Dv/mu/k/c; only rho_l is
  T-dependent), so a constant cp_l at the operating T is the faithful choice.
  Sources: NIST Chemistry WebBook (webbook.nist.gov, hexadecane C544763);
  Chemeo (chemeo.com/cid/30-657-9/Hexadecane); Tc=722-723 K.
- **rho_l(T)** -- Table-1 boiling anchor 569.9 @ 560 K, linear `767-0.758*(T-300)`.
- **Le = 2.5** -- n-hexadecane vapour (Dv-justified); kept physical.

With these the run's plateau is 492.6 K against the digitized ~493.7 K, and the 0-D replica of
the kernel (`zero_d.py`) reproduces it to the printed digit (see `INFO.md`).
