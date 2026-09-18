# Breakup

Five breakup models are available, selected by `[IGLOO-Models] breakup` in
`input.ini` and assigned by `assign_breakup` in `src/lib/Lib_Breakup.f90`.  They split
into two mechanistic classes:

- **Continuous-rate (ODE)**: add a $\dot{n}_p$ equation to the ODE system; active models set
  `brkupEqOde = .true.` → model 3 or 4.  `breakupOde` is evaluated inside the RHS.
- **Event-based (solout callback)**: modify parcel diameter/count at a discrete breakup
  event detected in the `solout` callback.  These set `brkupEvent = .true.` and run under
  model 1 or 2; they do **not** set `brkupEqOde`.

| Model | Keyword | `brkupEqOde` | `brkupEvent` | Spawns children |
| :--- | :--- | :---: | :---: | :---: |
| Pilch-Erdman | `Pilch-Erdman` | yes | no | no |
| Reitz-Diwakar | `Reitz-Diawakar` | yes | no | no |
| Reitz-KHRT | `Reitz-KHRT` | yes | yes | yes |
| TAB | `TAB` | no | yes | no |
| ETAB | `ETAB` | no | yes | no |

Model constants are read from `[IGLOO-Models]` under the names listed per model below
(defaults in the [registry](../user/registry.md)); the same set applies to every material.
The liquid properties they need — surface tension $\sigma$ and viscosity $\mu_p$ — come
from `INPUT/properties.dat` or the `[IGLOO-Properties] sigma` / `mu` overrides.

---

## Continuous-rate models

### Pilch-Erdman

```
breakup = Pilch-Erdman
```

Parameters: `Cd`, `B` (calibration constants).

Computes the dimensionless breakup time $T_{BT}$ from a five-branch $\mathrm{We}$–$\mathrm{Oh}$
map, then the stable diameter via

$$
V_{\Delta V} = \sqrt{\rho_g/\rho_p}\,(0.75\,C_d\,T_{BT} + 3\,B\,T_{BT}^2), \qquad
d_\mathrm{stable} = \frac{W_{e,c}\,\sigma}{\rho_g\,v_\mathrm{slip}^2\,(1 - V_{\Delta V})^2 + \varepsilon}
$$

The breakup-time scale $\tau = T_{BT}\,d / (v_\mathrm{slip}\,\sqrt{\rho_g/\rho_p})$, and the
parcel-number rate

$$
\dot{n}_p = -\frac{3\,n_p}{d}\,\frac{d_\mathrm{stable} - d}{\tau}
$$

Only fires when $\mathrm{We} > \mathrm{We}_c = 12(1 + 1.077\,\mathrm{Oh}^{1.6})$ and
$d > d_\mathrm{stable}$.

---

### Reitz-Diwakar

```
breakup = Reitz-Diawakar
```

Parameters: `WeBag`, `Cb`, `Cstrip`, `Cs`.

Distinguishes bag breakup ($\mathrm{We} > We_\mathrm{bag}$) from stripping breakup
($\mathrm{We} > C_\mathrm{strip}\sqrt{\mathrm{Re}}$):

| Regime | $\tau$ | $d_\mathrm{stable}$ |
| :--- | :--- | :--- |
| Stripping | $C_s\,d\,\sqrt{\rho_p/\rho_g}/(2 v)$ | $(2\,C_\mathrm{strip}\,\sigma)^2\,\mathrm{Re}/(\rho_g^2\,v^4\,d)$ |
| Bag | $C_b\,d\,\sqrt{\rho_p\,d/(16\sigma)}$ | $2\,\sigma\,\mathrm{We}_\mathrm{bag}/(\rho_g\,v^2)$ |

The parcel rate uses the same form as Pilch-Erdman:
$\dot{n}_p = -3n_p(d_\mathrm{stable} - d)/(d\,\tau)$.  The default constants are the
calibrated values of Reitz & Diwakar (1987): $C_s = 20$ (curve-fit) and $C_b = \pi$.

---

### Reitz-KHRT

```
breakup = Reitz-KHRT
```

Parameters: `B0`, `B1`, `Ctau`, `CRT`, `mShedLim`, `WeLimit`.

Kelvin-Helmholtz + Rayleigh-Taylor model.  With $r = d/2$, the gas Weber number
$\mathrm{We}_g$, $\mathrm{Oh} = \mu_p/\sqrt{\rho_p\,\sigma\,r}$ and
$\mathrm{Ta} = \mathrm{Oh}\sqrt{\mathrm{We}_g}$ (Reitz 1987, eqs. 4–5), the KH growth
rate, wavelength and timescale are

$$
\omega_\mathrm{KH} = \frac{0.34 + 0.38\,\mathrm{We}_g^{1.5}}{(1+\mathrm{Oh})(1+1.4\,\mathrm{Ta}^{0.6})} \sqrt{\frac{\sigma}{\rho_p\,r^3}}, \qquad
\lambda_\mathrm{KH} = 9.02\,r\,\frac{(1 + 0.45\sqrt{\mathrm{Oh}})(1 + 0.4\,\mathrm{Ta}^{0.7})}{(1 + 0.87\,\mathrm{We}_g^{1.67})^{0.6}}
$$

$$
\tau_\mathrm{KH} = \frac{3.726\,B_1\,r}{\lambda_\mathrm{KH}\,\omega_\mathrm{KH}}, \qquad
d_\mathrm{stable} = 2\,B_0\,\lambda_\mathrm{KH}
$$

RT growth rate, wavelength, and timescale:

$$
\omega_\mathrm{RT} = \sqrt{\frac{2\,|g_t(\rho_g - \rho_p)|^{3/2}}{3\sqrt{3\sigma}\,(\rho_p + \rho_g)}}, \qquad
\lambda_\mathrm{RT} = 2\pi\,C_{RT}\,\sqrt{\frac{3\sigma}{|g_t(\rho_g - \rho_p)|}}, \qquad
\tau_\mathrm{RT} = \frac{C_\tau}{\omega_\mathrm{RT}}
$$

where $g_t = (\mathbf{a}\cdot\mathbf{v}_p)/|\mathbf{v}_p|$ is the acceleration projected
along the particle velocity, and $\mathbf{a}$ is the drag acceleration.

The continuous ODE contribution is the KH stripping rate, active while no RT breakup is
pending, $d > d_\mathrm{stable}$ and $\mathrm{We}_g >$ `WeLimit`:

$$
\dot{n}_p = -\frac{3\,n_p}{d}\,\frac{d_\mathrm{stable} - d}{\tau_\mathrm{KH}}
$$

Two discrete events ride on top of it in the `solout` callback:

- **RT breakup** — once the RT clock exceeds $\tau_\mathrm{RT}$ with
  $\lambda_\mathrm{RT} < d$, the parent shatters into drops of diameter
  $\lambda_\mathrm{RT}$ (Beale & Reitz 1999, eq. 11): $n_p$ is multiplied by
  $(d/\lambda_\mathrm{RT})^3$ and $d$ reduced accordingly, conserving parcel mass.  No
  child parcel is created.
- **KH shed** — the mass stripped since the last event accumulates on the parent; once it
  exceeds the fraction `mShedLim` of the per-drop mass and the product count would be at
  least the parent count (Reitz 1987 product-parcel rule), a child parcel of diameter
  $d_\mathrm{stable}$ carrying the stripped mass is appended (`brkupHasChild = .true.`,
  `childState = [v_p, d_\mathrm{stable}, n_\mathrm{child}]`, velocity equal to the
  parent's) and the parent's count is restored.

---

## Event-based models

### TAB — Taylor Analogy Breakup

```
breakup = TAB
```

Parameters: `Comega`, `Cmu`, `WeCrit`, `n` (Rosin–Rammler spread), `method`.

Models the droplet as a forced, damped oscillator for the non-dimensional deformation $y$
(O'Rourke & Amsden 1987), integrated in closed form over each accepted ODE step:

$$
y(t) = y_\mathrm{eq} + e^{-t/t_d}
  \!\left[(y_0 - y_\mathrm{eq})\cos\omega t
  + \frac{(y_0 - y_\mathrm{eq})/t_d + \dot{y}_0}{\omega}\sin\omega t\right]
$$

with the paper's constants $C_k = 8$ (`Comega`), $C_d = 5$ (`Cmu`), $C_f = 1/3$,
$C_b = 1/2$ and $K = 10/3$.  Breakup is detected when $|y| > 1$ within the current
step.  At breakup the child drop radius is sampled from either a chi-square distribution
(`method = 1`, the original paper) or a Rosin-Rammler distribution (`method = 2`, the
default) whose scale parameter is iterated to match the $D_{32}$ implied by energy
conservation; the parcel count is rescaled to conserve mass.  Any other `method` value
is refused at setup.

---

### ETAB — Enhanced TAB

```
breakup = ETAB
```

Parameters: `k1`, `k2`, `WeCrit`, `WeTrans`, `Comega`, `Cmu`.

ETAB (Tanner 1997) keeps the TAB oscillator for the breakup *time* and replaces the
stochastic product sizing with a deterministic exponential cascade driven by a
Weber-dependent breakup rate $K_\mathrm{br}$:

$$
\frac{K_\mathrm{br}}{\omega} = \begin{cases}
k_1\,(A_{We}\,\mathrm{We}^4 + 1) & \mathrm{We} \le \mathrm{We}_\mathrm{trans} \quad\text{(bag)} \\[4pt]
k_2\,\sqrt{\mathrm{We}} & \mathrm{We} > \mathrm{We}_\mathrm{trans} \quad\text{(stripping)}
\end{cases}
$$

with $A_{We} = (k_2\sqrt{\mathrm{We}_\mathrm{trans}}/k_1 - 1)/\mathrm{We}_\mathrm{trans}^4$
chosen so that $K_\mathrm{br}$ is continuous at the transition.  The product radius is

$$
r_\mathrm{child} = r\,\exp\!\left(-\frac{K_\mathrm{br}}{\omega}\,\arccos\!\left(1 - \frac{1}{\mathrm{We}_\mathrm{cr}}\right)\right)
$$

and the product drops receive a velocity kick normal to the path,
$v_\perp = A\,\dot x$ (Tanner 1997, eqs. 8–10), at a random azimuth.  The parcel count
is rescaled to conserve mass; no child parcel is created.

---

## Child-particle spawning

Reitz-KHRT is the only model that creates new particles (`brkupHasChild = .true.`).
`integrate` in `src/lib/Lib_Integration.f90` receives the optional `shed` list and `noShed`
flag; each parent owns a growable `shedList` (`gr%shed(:)`).  `obj_IGLOO%solve` (outer
`maxLoop`) drains the lists into fully initialised `obj_particle` objects and appends them
to the next integration pass.  Every parcel — parents and children — draws its
stochastic samples from its own RNG stream seeded from `(seed, family, ID)`, so results
do not depend on thread scheduling or MPI rank count.

---

## V&V

Routine-level comparisons of all five models with their source papers, and the
end-to-end Weber-sweep cases (`tab-e2e`, `etab-e2e`, `pilch-erdman-e2e`,
`reitz-diwakar-e2e`, `khrt-e2e`): [../vv/literature.md](../vv/literature.md),
[../vv/e2e.md](../vv/e2e.md).
