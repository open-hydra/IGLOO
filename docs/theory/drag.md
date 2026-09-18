# Drag

The drag force on each particle is computed in `Lib_Equations.f90::interphase` as

$$
\mathbf{F}_\mathrm{drag} = \frac{\pi}{8}\,C_d\,d\,\mathrm{Re}\,\mu_g\,(\mathbf{v}_g - \mathbf{v}_p)
$$

where $d$ is the particle diameter, $\mathrm{Re} = \rho_g\,|\mathbf{v}_g - \mathbf{v}_p|\,d/\mu_g$
is the particle Reynolds number, and $C_d(\mathrm{Re}, \mathrm{Ma}, \gamma, T_r)$ is selected
by `assign_drag` in `Lib_Drag.f90`.  The keyword string from `[IGLOO-Models] drag` in
`input.ini` maps directly to one of the 13 correlations below.  The compressible
correlations also receive the slip Mach number $\mathrm{Ma} = |\mathbf{v}_g - \mathbf{v}_p| /
\sqrt{\gamma R_g T_g}$, the gas $\gamma$, and the temperature ratio $T_r = T_p / T_g$.

Every correlation carries the sentinel $\varepsilon = 10^{-20}$ (a per-function `toll` parameter)
in each denominator that could vanish at $\mathrm{Re} = 0$ or $\mathrm{Ma} = 0$; it is
omitted from the formulas below where it only guards a division.

---

## Incompressible correlations

### Stokes

```
drag = Stokes
```

$$
C_d = \frac{24}{\mathrm{Re} + \varepsilon}
$$

Pure viscous limit.  With $\varepsilon = 10^{-20}$ the drag is linear in the slip velocity
to machine precision for any $\mathrm{Re}$, which is what makes the closed-form
verification cases exact.

---

### Schlichting

```
drag = Schlichting
```

$$
C_d = \frac{24}{\mathrm{Re} + \varepsilon}\left(1 + \frac{3\,\mathrm{Re}}{16}\right)
$$

One-term Oseen correction; exact for $\mathrm{Re}\ll 1$.

---

### Schiller-Naumann

```
drag = Schiller-Naumann
```

$$
C_d = \frac{24}{\mathrm{Re} + \varepsilon}\left(1 + 0.15\,\mathrm{Re}^{0.687}\right)
$$

Valid up to $\mathrm{Re} \approx 800$.

---

### Morsi-Alexander

```
drag = Morsi-Alexander
```

$$
C_d = a_1 + \frac{a_2}{\mathrm{Re} + \varepsilon} + \frac{a_3}{\mathrm{Re}^2 + \varepsilon}
$$

Eight Reynolds-number ranges; the $a_i$ coefficients as coded in `Lib_Drag.f90`:

| Re range | $a_1$ | $a_2$ | $a_3$ |
| :--- | ---: | ---: | ---: |
| $[0,\,0.1]$ | 0 | 24 | 0 |
| $(0.1,\,1]$ | 3.69 | 22.73 | 0.0903 |
| $(1,\,10]$ | 1.222 | 29.1667 | −3.8889 |
| $(10,\,100]$ | 0.6167 | 46.5 | −116.67 |
| $(100,\,1000]$ | 0.3644 | 98.33 | −2778 |
| $(1000,\,5000]$ | 0.357 | 148.62 | −4.75×10⁴ |
| $(5000,\,10000]$ | 0.46 | −490.546 | 5.787×10⁵ |
| $(10000,\,\infty)$ | 0.5191 | −1662.5 | 5.4167×10⁶ |

---

### Wen-Yu

```
drag = Wen-Yu
```

$$
C_d = \begin{cases}
\dfrac{24}{\mathrm{Re}+\varepsilon}(1 + 0.15\,\mathrm{Re}^{0.687}) & \mathrm{Re} \le 1000 \\[6pt]
0.43 & \mathrm{Re} > 1000
\end{cases}
$$

Schiller–Naumann below $\mathrm{Re} = 1000$ and the constant plateau $0.43$ above it
(Shimada 2006, eq. 16).  The low branch evaluates to $\approx 0.4383$ at
$\mathrm{Re} = 1000$, so a small step remains at the handoff — a property of the published
form.  Wen–Yu is also the incompressible base $C_{d,0}$ of the Crowe, Hermsen and
Carlson–Hoglund correlations below.

---

### Newton

```
drag = Newton
```

$$
C_d = 0.45
$$

Constant, valid in the Newton regime ($\mathrm{Re} \in [1000,\,2\times10^5]$).

---

### Chang

```
drag = Chang
```

$$
C_d = \frac{24}{\mathrm{Re}+\varepsilon}\left(1 + 0.15\,\mathrm{Re}^{0.687}\right) + \frac{0.42}{1 + 42500\,(\mathrm{Re}+\varepsilon)^{-1.16}}
$$

Continuous across all $\mathrm{Re}$; algebraically identical to Clift–Gauvin below
($24\times0.0175 = 0.42$), differing only in where $\varepsilon$ enters.

---

### Clift-Gauvin

```
drag = Clift-Gauvin
```

$$
C_d = \frac{24}{\mathrm{Re}+\varepsilon}\left(1 + 0.15\,\mathrm{Re}^{0.687} + \frac{0.0175\,\mathrm{Re}}{1 + 4.25\times10^{4}\,\mathrm{Re}^{-1.16}}\right)
$$

Continuous across all $\mathrm{Re}$ ranges.

---

### Putnam

```
drag = Putnam
```

$$
C_d = \begin{cases}
\dfrac{24}{\mathrm{Re}+\varepsilon}\left(1 + \dfrac{\mathrm{Re}^{2/3}}{6}\right) & \mathrm{Re} < 1000 \\[6pt]
0.4392 & \mathrm{Re} \ge 1000
\end{cases}
$$

The plateau $0.4392$ is the published value (Shimada 2006, eq. 17); the low branch gives
$24\,(1 + 100/6)/1000 = 0.424$ at $\mathrm{Re} = 1000$, so the published form carries a
$\approx 3.5\,\%$ step at the handoff.

---

## Compressible correlations

### Henderson

```
drag = Henderson
```

Henderson (1976); the branch bodies below are the forms as coded (from the Shimada 2006
catalog, eqs. 20–22).  With $s = \mathrm{Ma}\sqrt{\gamma/2}$, the subsonic branch
($\mathrm{Ma} \le 1$) is

$$
C_{d,1} = \frac{24}{\mathrm{Re} + s\left[4.33 + \dfrac{3.65 - 1.53\,T_r}{1 + 0.353\,T_r}\,e^{-0.247\,\mathrm{Re}/s}\right]}
+ e^{-0.5\,\mathrm{Ma}/\sqrt{\mathrm{Re}}}\left[\frac{4.5 + 0.38\,(0.03\,\mathrm{Re} + 0.48\sqrt{\mathrm{Re}})}{1 + 0.03\,\mathrm{Re} + 0.48\sqrt{\mathrm{Re}}} + 0.1\,\mathrm{Ma}^2 + 0.2\,\mathrm{Ma}^8\right]
+ 0.6\,s\left(1 - e^{-\mathrm{Ma}/\mathrm{Re}}\right)
$$

and the supersonic branch ($\mathrm{Ma} \ge 1.75$), as coded, is

$$
C_{d,2} = \frac{0.9 + \dfrac{0.34}{\mathrm{Ma}^2} + 1.86\sqrt{\dfrac{\mathrm{Ma}}{\mathrm{Re}}}\left(2 + \dfrac{2}{s^2} + \dfrac{1.058\sqrt{T_r}}{s}\right) - \dfrac{1}{s^4}}{1 + 1.86\sqrt{\mathrm{Ma}/\mathrm{Re}}}
$$

In the transonic interval $1 < \mathrm{Ma} < 1.75$ the two branches are bridged linearly,

$$
C_d = C_{d,1}(\mathrm{Ma}{=}1) + \tfrac{4}{3}\,(\mathrm{Ma}-1)\,\bigl[C_{d,2}(\mathrm{Ma}{=}1.75) - C_{d,1}(\mathrm{Ma}{=}1)\bigr]
$$

where $4/3 = 1/(1.75 - 1)$ makes the bridge continuous at both ends.

---

### Crowe

```
drag = Crowe
```

Shimada 2006, eq. 23.  With $C_{d,0}$ the Wen–Yu value,

$$
C_d = 2 + (C_{d,0} - 2)\,\exp\!\left(-\frac{3.07\sqrt{\gamma}\,\mathrm{Ma}\,g(\mathrm{Re})}{\mathrm{Re}}\right)
+ \frac{h(\mathrm{Ma}, T_r)}{\mathrm{Ma}\sqrt{\gamma}}\,\exp\!\left(-\frac{\mathrm{Re}}{2\,\mathrm{Ma}}\right)
$$

$$
g(\mathrm{Re}) = 10^{\,1.25\,[1 + \tanh(0.77\log_{10}\mathrm{Re} - 1.92)]}, \qquad
h(\mathrm{Ma}, T_r) = 2.3 + 1.7\sqrt{T_r} - 2.3\tanh(1.17\log_{10}\mathrm{Ma})
$$

The first exponential recovers $C_{d,0}$ in the continuum limit; the second is the
rarefaction (free-molecular) term, which vanishes for $\mathrm{Re}/\mathrm{Ma} \gg 1$.

---

### Hermsen

```
drag = Hermsen
```

Shimada 2006, eq. 24.  Same structure as Crowe with

$$
g(\mathrm{Re}) = \frac{1 + \mathrm{Re}\,(12.278 + 0.548\,\mathrm{Re})}{1 + 11.278\,\mathrm{Re}}, \qquad
h(\mathrm{Ma}, T_r) = \frac{5.6}{1 + \mathrm{Ma}} + 1.7\sqrt{T_r}
$$

---

### Carlson-Hoglund

```
drag = Carlson-Hoglund
```

With $C_{d,0}$ the Wen–Yu value,

$$
C_d = C_{d,0}\;\frac{1 + \exp\!\left(-\dfrac{0.427}{\mathrm{Ma}^{4.63}} - \dfrac{3}{\mathrm{Re}^{0.88}}\right)}
{1 + \dfrac{\mathrm{Ma}}{\mathrm{Re}}\left(3.82 + 1.28\,e^{-1.25\,\mathrm{Re}/\mathrm{Ma}}\right)}
$$

Compressible correlation valid across the full Mach range.

---

## Selection summary

| Keyword | Compressibility | Re range | Notes |
| :--- | :---: | :--- | :--- |
| `Stokes` | — | $\mathrm{Re}\to 0$ | theoretical limit; exactly linear drag |
| `Schlichting` | — | $\mathrm{Re}\ll 1$ | Oseen correction |
| `Schiller-Naumann` | — | $\mathrm{Re}<800$ | common default |
| `Morsi-Alexander` | — | all | 8-range piecewise |
| `Wen-Yu` | — | all | plateau 0.43 above Re = 1000 |
| `Newton` | — | $10^3$–$2\times10^5$ | constant 0.45 |
| `Chang` | — | all | continuous; ≡ Clift-Gauvin |
| `Clift-Gauvin` | — | all | continuous |
| `Putnam` | — | all | plateau 0.4392 above Re = 1000 |
| `Henderson` | yes | all | three Mach regimes, linear transonic bridge |
| `Crowe` | yes | all | continuum + rarefaction terms on Wen–Yu |
| `Hermsen` | yes | all | continuum + rarefaction terms on Wen–Yu |
| `Carlson-Hoglund` | yes | all | Mach correction on Wen–Yu |

---

## V&V

Routine-level comparison of every correlation against its source expression, and the
end-to-end Stokes relaxation case: [Literature tests](../vv/literature.md),
[Stokes drag](../vv/drag-stokes.md).
