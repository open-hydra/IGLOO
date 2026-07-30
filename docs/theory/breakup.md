# Breakup

Five breakup models are available, selected by `[IGLOO-Models] breakup-model` in
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
| Reitz-Diawakar | `Reitz-Diawakar` | yes | no | no |
| Reitz-KHRT | `Reitz-KHRT` | yes | yes | yes |
| TAB | `TAB` | no | yes | no |
| ETAB | `ETAB` | no | yes | no |

---

## Continuous-rate models

### Pilch-Erdman

```
breakup-model = Pilch-Erdman
```

Parameters: `bp(1)=Cd`, `bp(2)=B` (calibration constants, 2 total).

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

### Reitz-Diawakar

```
breakup-model = Reitz-Diawakar
```

Parameters: `bp(1)=WeBag`, `bp(2)=Cb`, `bp(3)=Cstrip`, `bp(4)=Cs` (4 total).

Distinguishes bag breakup ($\mathrm{We} > We_\mathrm{bag}$) from stripping breakup
($\mathrm{We} > C_\mathrm{strip}\sqrt{\mathrm{Re}}$):

| Regime | $\tau$ | $d_\mathrm{stable}$ |
| :--- | :--- | :--- |
| Stripping | $C_s\,d\,\sqrt{\rho_p/\rho_g}/(2 v)$ | $(2\,C_\mathrm{strip}\,\sigma)^2\,\mathrm{Re}/(\rho_g^2\,v^4\,d)$ |
| Bag | $C_b\,d\,\sqrt{\rho_p\,d/(4\sigma)}$ | $2\,\sigma\,\mathrm{We}_\mathrm{bag}/(\rho_g\,v^2)$ |

The parcel rate uses the same form as Pilch-Erdman:
$\dot{n}_p = -3n_p(d_\mathrm{stable} - d)/(d\,\tau)$.

---

### Reitz-KHRT

```
breakup-model = Reitz-KHRT
```

Parameters: `bp(1)=B0`, `bp(2)=B1`, `bp(3)=Ctau`, `bp(4)=CRT`, `bp(5)=mShedLim`,
`bp(6)=WeLimit` (6 total).

Kelvin-Helmholtz + Rayleigh-Taylor model.  Computes KH growth rate $\omega_\mathrm{KH}$,
wavelength $\lambda_\mathrm{KH}$, and timescale $\tau_\mathrm{KH}$:

$$
\omega_\mathrm{KH} = \frac{0.34 + 0.38\,\mathrm{We}_g^{1.5}}{(1+\mathrm{Oh})(1+1.4\,\mathrm{Ta}^{0.6})} \cdot \sqrt{\frac{\sigma}{\rho_p\,r^3}}
$$

!!! success "Fixed (A2, 2026-07-07)"
    `sqrt` and grouping corrected at both `Lib_Breakup.f90` sites (lines 192–193 and
    280–281).  Previously coded as $\sigma r^3/\rho_p$ — missing the `sqrt` and $r^3$
    on the wrong side — giving units m⁶/s² instead of s⁻¹ and making `tauKH`
    dimensionally garbage.  Compare with TAB/ETAB (lines 373/478), which always used
    the correct grouping.  Gated by `test_breakup_khrt` (code ≡ Reitz-1987 chain to
    1e-12 on the pure-KH path; RT quantities and event path deferred).

RT growth rate, wavelength, and timescale:

$$
\omega_\mathrm{RT} = \sqrt{\frac{2\,|g_t(\rho_g - \rho_p)|^{3/2}}{3\sqrt{3\sigma(\rho_p + \rho_g)}}}, \qquad
\lambda_\mathrm{RT} = 2\pi\,C_{RT}\,\sqrt{\frac{3\sigma}{|g_t(\rho_g - \rho_p)|} + \varepsilon}, \qquad
\tau_\mathrm{RT} = \frac{C_\tau}{\omega_\mathrm{RT} + \varepsilon}
$$

where $g_t = (\mathbf{a}\cdot\mathbf{v}_p)/|\mathbf{v}_p|$ is the acceleration projected
along the particle velocity, and $\mathbf{a}$ is the drag acceleration.

`brkupHasChild = .true.`: the KH stripping event appends a child parcel with
`childState = [v_p, d_\mathrm{stable}, n_\mathrm{child}]`.  The RT event directly updates
$d$ and $n_p$ on the parent (no child).

---

## Event-based models

### TAB — Taylor Analogy Breakup

```
breakup-model = TAB
```

Parameters: `bp(1)=Comega`, `bp(2)=Cmu`, `bp(3)=WeCrit`, `bp(4)=nSpread` (4 total).

Models the droplet as a forced, damped oscillator for the non-dimensional deformation $y$.
Breakup is detected when $|y| > 1$ within the current timestep.  At breakup the child
drop radius is sampled from either a chi-square distribution (method=1) or a
Rosin-Rammler distribution (method=2, default), with scale parameter iterated to match
the $D_{32}$ implied by energy conservation.

Constants $C_k=8, C_d=5, C_f=1/3, C_b=1/2$ (from O'Rourke and Amsden 1987) are embedded
in `bp`.  All values verified against the original paper (no bug).

---

### ETAB — Enhanced TAB

```
breakup-model = ETAB
```

Parameters: `bp(1)=k1`, `bp(2)=k2`, `bp(3)=WeCrit`, `bp(4)=WeTrans`, `bp(5)=Comega`,
`bp(6)=Cmu` (6 total).

ETAB replaces the child-size distribution with an exponential drop-size law driven by the
oscillation amplitude and a Weber-dependent breakup rate $K_\mathrm{br}$:

$$
K_\mathrm{br} = k_1\,\omega\,\sqrt{\mathrm{We}} \qquad (\mathrm{We} > \mathrm{We}_\mathrm{trans})
$$

otherwise $K_\mathrm{br} = k_1\,\omega\,(\mathrm{AWe}\cdot\mathrm{We}^4 + 1)^{1/2}$ in
the low-We limit.

---

## Child-particle spawning

Reitz-KHRT is the only model that creates new particles (`brkupHasChild = .true.`).
`integrate` in `src/lib/Lib_Integration.f90` receives `child` and `addChild` optional
arguments.  `obj_IGLOO%solve` (outer `maxLoop`) collects spawned children into
`gr%child(:)`, compacts them into fully initialised `obj_particle` objects, and appends
them to the next integration pass.

---

## V&V

Reitz-KHRT, TAB, and ETAB literature comparisons: [../vv/literature.md](../vv/literature.md).
