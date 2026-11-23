# Breakup

Five breakup models are available, chosen per material ([Materials & Models](../user/models.md#breakup)). They work in
two ways:

- **At a continuous rate**: Pilch–Erdman, Reitz–Diwakar and the KH stripping of Reitz KH-RT add an equation for the
  number rate $\dot n_p$: the parcel's drops shrink and multiply at a constant stream mass (systems 3, 4 and 7).
- **By events**: TAB, ETAB and the RT breakup and KH shedding of Reitz KH-RT change the parcel at a discrete instant,
  between two integration steps. TAB and ETAB add no equation, so their material stays in system 1, 2 or 6.

| Model | Continuous rate | Events | Creates new parcels |
| :--- | :---: | :---: | :---: |
| Pilch–Erdman | yes | no | no |
| Reitz–Diwakar | yes | no | no |
| Reitz KH-RT | yes | yes | yes |
| TAB | no | yes | no |
| ETAB | no | yes | no |

An event sets a new diameter and scales the number rate so that the stream mass is kept; on an evaporating drop it
also resets the drop's mass. On a solidifying material breakup acts on the liquid only, and stops at nucleation
([Solidification](solidification.md#breakup-on-the-liquid)).

The models need the liquid's surface tension $\sigma$ and viscosity $\mu_p$, constants or tables in temperature read
at the drop's temperature ([Material Properties](../user/properties.md)). Each model's constants have published
defaults, which the input may change.

---

## Continuous-rate models

### Pilch-Erdman

Constants: $C_d$ and $B$.

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

### Reitz–Diwakar

Constants: $\mathrm{We}_\mathrm{bag}$, $C_b$, $C_\mathrm{strip}$, $C_s$.

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

### Reitz KH-RT

Constants: $B_0$, $B_1$, $C_\tau$, $C_{RT}$, the shed limit and the Weber limit.

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
pending, $d > d_\mathrm{stable}$ and $\mathrm{We}_g$ is above the Weber limit (6):

$$
\dot{n}_p = -\frac{3\,n_p}{d}\,\frac{d_\mathrm{stable} - d}{\tau_\mathrm{KH}}
$$

Two discrete events ride on top of it:

- **RT breakup** — once the RT clock exceeds $\tau_\mathrm{RT}$ with
  $\lambda_\mathrm{RT} < d$, the parent shatters into drops of diameter
  $\lambda_\mathrm{RT}$ (Beale & Reitz 1999, eq. 11): $n_p$ is multiplied by
  $(d/\lambda_\mathrm{RT})^3$ and $d$ reduced accordingly, conserving parcel mass.  No
  child parcel is created.
- **KH shed** — the mass stripped since the last event accumulates on the parent; once it
  exceeds the shed limit (3 % by default) of the per-drop mass and the product count would be at
  least the parent count (Reitz 1987 product-parcel rule), a child parcel of diameter
  $d_\mathrm{stable}$, with the parent's velocity, carries the stripped mass away, and the parent's count is
  restored.  Without evaporation the stream
  mass is constant and the stripped mass is the per-drop deficit since the last event.  With
  evaporation the drop mass $m$ also falls by evaporation, which never changes
  $\dot n$: with $n_0$ the number rate right after the last event, the stripped mass flow is
  $(\dot n - n_0)\,m$, the shed fires when it exceeds the shed limit times $n_0 m_0$ and
  $n_\mathrm{child} = (\dot n - n_0)(d/d_\mathrm{stable})^3 \ge n_0$, and the parent returns
  to $n_0$, so evaporated mass never enters a child.  The stripped mass leaves with the child and
  never enters the gas source: the parent's outgoing flux in the shed segment includes the
  child's birth flux, which the child deposits along its own path (see
  [Eulerian feedback](eulerian-feedback.md)).

The model is valid for $\mathrm{Oh} \le 1$ and $\rho_g/\rho_p \le 0.1$. IGLOO runs it outside that range too, and
the run log counts the evaluations that left it.

---

## Event-based models

### TAB — Taylor Analogy Breakup

Constants: $C_k$, $C_d$, the critical Weber number, the product law and its spread $n$.

Models the droplet as a forced, damped oscillator for the non-dimensional deformation $y$
(O'Rourke & Amsden 1987), integrated in closed form over each accepted ODE step:

$$
y(t) = y_\mathrm{eq} + e^{-t/t_d}
  \!\left[(y_0 - y_\mathrm{eq})\cos\omega t
  + \frac{(y_0 - y_\mathrm{eq})/t_d + \dot{y}_0}{\omega}\sin\omega t\right]
$$

with the paper's constants $C_k = 8$, $C_d = 5$, $C_f = 1/3$,
$C_b = 1/2$ and $K = 10/3$.  Breakup is detected when $|y| > 1$ within the current
step.  At breakup the product radius is drawn from the mass distribution of a chi-square
(the exponential law of O'Rourke and Amsden's KIVA-II implementation) or a
Rosin–Rammler (the default) size law whose Sauter radius is the $r_s$ of energy
conservation.  Both are drawn inside $[0.01\,r_s, r)$: the Rosin–Rammler scale is iterated so
that the truncated law keeps $r_s$; the chi-square law has $r_s$ untruncated and keeps the
Sauter radius of the window, at least $0.946\,r_s$ (at $r/r_s = 7/3$, its smallest value).  The
parcel keeps its mass, so the drops it represents follow the size law.

---

### ETAB — Enhanced TAB

Constants: $k_1$, $k_2$, the critical and the transition Weber numbers, $C_k$, $C_d$.

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
$v_\perp = A\,\dot x$ (Tanner 1997, eqs. 8–10), at a random azimuth, added to the parcel velocity once, at the event.  The parcel count is rescaled to conserve mass; no child parcel is created.

---

## New parcels

Reitz KH-RT is the only model that creates parcels. A child starts where its parent shed it and is integrated after
the parents, in the same sweep, until no new child appears. Every parcel, parent or child, draws its random numbers
from its own stream, seeded from the run's seed, its family and its number, so the results do not depend on how the
parcels are shared among threads or processes.

---

## V&V

Every model against its source paper, the Weber sweeps end to end, the TAB and ETAB products beside Hsiang and
Faeth's correlations, and the KH-RT parent drop beside their core drop: [Breakup V&V](../vv/breakup.md).
