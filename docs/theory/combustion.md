# Metal Combustion

Aluminum combustion is the first model of the metal track: a material with
`[GPB-PhaseX] combustion = Beckstead` in `input.ini` is routed to the dedicated
ODE family `model = 5` (`src/lib/Lib_RHS.f90::rhsAlCombustion`), whose closure
lives in `src/lib/Lib_Combustion.f90`.  The state layout is identical to the
evaporation family (model 2): position/velocity in Z(1:6), temperature (or
enthalpy) at Z(7), the particle **Al mass** at Z(8).  Combustion is mutually
exclusive per material with evaporation (a loud warning drops evaporation) and
with breakup (hard error; the coupling is deferred).

---

## Beckstead d^n burn law

Beckstead's correlation of ~400 aluminum burn-time data points gives
$t_b \propto d^n$ with $n \approx 1.8$ (range 1.5–1.8) — deliberately *below*
the diffusion-controlled d²-law exponent, which is why a fuel-evaporation
closure must not be reused here.  IGLOO posits the implied regression

$$
d^n(t) = d_0^n - K_\mathrm{eff}\,t
$$

and differentiates it into an instantaneous mass rate (with
$m = \rho_p \pi d^3/6$):

$$
\dot{m} = -\frac{\rho_p\,\pi\,K_\mathrm{eff}}{2\,n}\; d^{\,3-n}
$$

so the analytic burn time $t_b = d_0^n/K_\mathrm{eff}$ is exact by
construction.  For $n = 1.8$ the rate scales as $d^{1.2}$.

### How X-eff enters K

IGLOO's gas field carries no species, so the local oxidizer composition is
replaced by a **constant effective oxidizer mole fraction** per material
(input `X-eff`), following Beckstead's effectiveness ranking
(O₂ ≈ 2× H₂O ≈ 5× CO₂):

$$
X_\mathrm{eff} = C_{\mathrm{O}_2} + 0.6\,C_{\mathrm{H}_2\mathrm{O}} + 0.22\,C_{\mathrm{CO}_2},
\qquad
K_\mathrm{eff} = K_\mathrm{burn}\; X_\mathrm{eff}^{\,1.0}
$$

The effective oxidizer enters *linearly* ($X_\mathrm{eff}^{1.0}$) — Beckstead's
final correlation ([Beck05], Eq. p.541: $t_b X_\mathrm{eff}\,p^{0.1}T_0^{0.2} =
0.00735\,D^{1.8}$); the 0.9 sometimes quoted is Belyaev's *earlier* $D^{1.5}$ fit
(exposed as the named parameter `aXeff` in `Lib_Combustion`).  The correlation's *weak* pressure and initial
temperature dependence ($p^{0.1}\,T_0^{0.2}$) is **not** applied in-code — it
is folded into the user-supplied `K-burn` (units m^`n-burn`/s, defined at
`X-eff = 1`).  Per-cell oxidizer composition is deferred until the ORION gas
import grows species fields.

### Ignition gate

Below the ignition temperature `T-ign` (default 2350 K — the Al₂O₃-shell
melting anchor) the particle is **inert**: $\dot m = 0$ exactly (bitwise), and
only drag and convective heat act.  The gate is $T_p < T_\mathrm{ign}$, so a
particle at exactly `T-ign` burns.

### Heat release

Combustion *releases* heat to the particle — the sign is opposite to the
evaporation latent-heat sink hard-coded in the model-2 $F(7)$.  The energy
equation of `rhsAlCombustion` is

$$
\frac{dT_p}{dt} \;(\text{or } \tfrac{dh_p}{dt})
= \frac{\dot{Q}_\mathrm{conv} + \beta_\mathrm{part}\, q_\mathrm{comb}\, |\dot m| \cdot c_p\text{-factor}}{m}
$$

where $\dot Q_\mathrm{conv}$ is the standard interphase convective term (which
already carries the `cpFactor` unit convention: 1 for the enthalpy state,
$1/c_p$ for the temperature state — see [Heat Transfer](heat.md)), and the
release term is scaled by the same factor so both live in consistent units.

!!! warning "beta-part sensitivity (plan risk R5)"
    The heat-partition fraction $\beta_\mathrm{part}$ (share of the reaction
    enthalpy deposited on the condensed phase rather than in the detached
    flame) is **weakly constrained in the literature**.  It is therefore an
    explicit input with default **0** (no heat-back unless configured) rather
    than a buried constant, and `q-comb` (heat of combustion per unit Al mass,
    J/kg; ≈31 MJ/kg for Al) also defaults to 0.  Results in hot-particle
    regimes can be sensitive to $\beta_\mathrm{part}\,q_\mathrm{comb}$;
    sweep it before trusting particle-temperature predictions.

### Robustness

Near burnout an implicit-solver trial step can drive $m \le 0$, making
$d = (6m/\rho\pi)^{1/3}$ NaN; like models 3/4, any non-finite RHS entry is
scrubbed to the `1e30` penalty so SDIRK4 rejects the trial instead of
poisoning the integration.  Burnout itself (terminating the parent and
spawning the inert alumina residual) is phase M2.

---

## Inputs

All in the material's `[GPB-PhaseX]` section (see the
[registry](../user/registry.md)):

| Key | Meaning | Default |
| :--- | :--- | :---: |
| `combustion` | `Beckstead` switches the material to the metal track | absent = off |
| `K-burn` | burn-rate coefficient $K$ at $X_\mathrm{eff}=1$ (m^`n-burn`/s) | required > 0 |
| `n-burn` | diameter exponent $n$ | 1.8 |
| `X-eff` | effective oxidizer mole fraction (0, 1] | 1.0 |
| `T-ign` | ignition temperature (K) | 2350 |
| `beta-part` | heat-partition fraction to the particle | 0 |
| `q-comb` | heat of combustion (J/kg) | 0 |

Euler feedback and the mass source follow the model-2 rules: the instantaneous
stream rate $\dot n_p\,m$ deposits the burned Al mass on the gas
(see [Eulerian Feedback](eulerian-feedback.md)).

---

## V&V

Unit family `tests/combustion/` (closed-form complex-step oracle, d^n
regression identity, bitwise ignition gate, informational Beckstead
correlation points) and e2e box case `tests/combustion/burn-box/` (exact
closed-form $d^n(x)$, independent RK4 of the heat-release energy balance,
mass-telescoping audit).  See `tests/combustion/INFO.md`.

Reference: Beckstead, M. W., "Correlating Aluminum Burning Times,"
*Combust. Explos. Shock Waves* 41(5):533–546, 2005.
