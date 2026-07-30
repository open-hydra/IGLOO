# Drag

The drag force on each particle is computed in `Lib_Equations.f90::interphase` as

$$
\mathbf{F}_\mathrm{drag} = \frac{\pi}{8}\,C_d\,d\,\mathrm{Re}\,\mu_g\,(\mathbf{v}_g - \mathbf{v}_p)
$$

where $d$ is the particle diameter, $\mathrm{Re} = \rho_g\,|\mathbf{v}_g - \mathbf{v}_p|\,d/\mu_g$
is the particle Reynolds number, and $C_d(\mathrm{Re}, \mathrm{Ma}, \gamma, T_r)$ is selected
by `assign_drag` in `Lib_Drag.f90`.  The keyword string from `[IGLOO-Models] drag-model` in
`input.ini` maps directly to one of the 13 correlations below.

---

## Correlations

### Stokes

```
drag-model = Stokes
```

$$
C_d = \frac{24}{\mathrm{Re} + \varepsilon}
$$

Pure viscous limit.  $\varepsilon = 10^{-20}$ (the module-level `toll`) prevents division by zero
at $\mathrm{Re}=0$; this sentinel appears in every correlation below.

---

### Schlichting

```
drag-model = Schlichting
```

$$
C_d = \frac{24}{\mathrm{Re} + \varepsilon}\left(1 + \frac{3\,\mathrm{Re}}{16}\right)
$$

One-term Oseen correction; exact for $\mathrm{Re}\ll 1$.

---

### Schiller-Naumann

```
drag-model = Schiller-Naumann
```

$$
C_d = \frac{24}{\mathrm{Re} + \varepsilon}\left(1 + 0.15\,\mathrm{Re}^{0.687}\right)
$$

Valid up to $\mathrm{Re} \approx 800$.

---

### Morsi-Alexander

```
drag-model = Morsi-Alexander
```

$$
C_d = a_1 + \frac{a_2}{\mathrm{Re} + \varepsilon} + \frac{a_3}{\mathrm{Re}^2 + \varepsilon}
$$

Eight Reynolds-number ranges; the $a_i$ coefficients are as coded in `Lib_Drag.f90`:

| Re range | $a_1$ | $a_2$ | $a_3$ |
| :--- | ---: | ---: | ---: |
| $[0,\,0.1)$ | 0 | 24 | 0 |
| $[0.1,\,1)$ | 3.69 | 22.73 | 0.0903 |
| $[1,\,10)$ | 1.222 | 29.1667 | −3.8889 |
| $[10,\,100)$ | 0.6167 | 46.5 | −116.67 |
| $[100,\,1000)$ | 0.3644 | 98.33 | −2778 |
| $[1000,\,5000)$ | 0.357 | 148.62 | −4.75×10⁴ |
| $[5000,\,10000)$ | 0.46 | −490.546 | 57874.52 |
| $[10000,\,\infty)$ | 0.5191 | −1662.5 | 5.4167×10⁵ |

---

### Wen-Yu

```
drag-model = Wen-Yu
```

$$
C_d = \begin{cases}
\dfrac{24}{\mathrm{Re}+\varepsilon}(1 + 0.15\,\mathrm{Re}^{0.687}) & \mathrm{Re} \le 1000 \\[6pt]
0.44 & \mathrm{Re} > 1000
\end{cases}
$$

Fixed (2026-07-07, gated XD5): high-Re constant corrected to the published Rowe plateau
`0.44` (was `0.43`).  The Schiller–Naumann low branch evaluates to ≈ 0.4383 at
Re = 1000, so a small inherent step remains at the handoff — this is a property of the
published form, not a continuity bug.

---

### Newton

```
drag-model = Newton
```

$$
C_d = 0.45
$$

Constant, valid in the Newton regime ($\mathrm{Re} \in [1000,\,2\times10^5]$).

---

### Chang

```
drag-model = Chang
```

$$
C_d = \begin{cases}
\dfrac{24}{\mathrm{Re}+\varepsilon}(1 + 0.0175\,\mathrm{Re}) & \mathrm{Re} < 1000 \\[6pt]
0.43 & \mathrm{Re} \ge 1000
\end{cases}
$$

Algebraically identical to Clift-Gauvin at the coded coefficient values
($24\times0.0175=0.42$, breakpoint constant $0.43$).

---

### Clift-Gauvin

```
drag-model = Clift-Gauvin
```

$$
C_d = \frac{24}{\mathrm{Re}+\varepsilon}\left(1 + 0.15\,\mathrm{Re}^{0.687}\right) + \frac{0.42}{1 + 42500\,\mathrm{Re}^{-1.16}}
$$

Continuous across all Re ranges.

---

### Putnam

```
drag-model = Putnam
```

$$
C_d = \begin{cases}
\dfrac{24}{\mathrm{Re}+\varepsilon}\left(1 + \dfrac{\mathrm{Re}^{2/3}}{6}\right) & \mathrm{Re} < 1000 \\[6pt]
0.424 & \mathrm{Re} \ge 1000
\end{cases}
$$

!!! success "Fixed (A5, 2026-07-07)"
    High-Re constant corrected from `0.4392` to `0.424`, matching the low-Re branch at
    Re = 1000 ($24\cdot(1+100/6)/1000 = 0.424$ exactly).  Gated by XD3 (jump < 1e-6).

---

### Henderson

```
drag-model = Henderson
```

Compressibility-aware correlation; distinguishes subsonic ($\mathrm{Ma}<1$), transonic
($1 \le \mathrm{Ma} < 1.75$), and supersonic ($\mathrm{Ma} \ge 1.75$) regimes.

!!! success "Fixed (A6, 2026-07-07)"
    Transonic blend corrected to $C_{d,1} + \tfrac{4}{3}\,(\mathrm{Ma}-1)\,(C_{d,2}-C_{d,1})$.
    The interval $\mathrm{Ma}\in[1,\,1.75]$ has width 0.75, so the slope must be $1/0.75=4/3$;
    the old coefficient `0.75` left a 44 % step in $C_d$ at Ma = 1.75.  Gated by XD4
    (jump < 1e-6).

---

### Crowe

```
drag-model = Crowe
```

!!! success "Fixed (A1, 2026-07-07)"
    Rarefaction exponential moved from denominator to numerator: the corrected form is
    `hfun / (Ma*sqrt(G)) * exp(-Re/(2*Ma))`.  Previously the exp was in the denominator,
    underflowing to `toll=1e-20` for Re/Ma ≳ 90 and giving $C_d \approx 10^{20}$.
    Gated by XD1 (Cd finite ≈ 0.44 at Re = 1e3, Ma = 2).

---

### Hermsen

```
drag-model = Hermsen
```

!!! success "Fixed (A1, 2026-07-07)"
    Same denominator-exp defect as Crowe, fixed identically at `Lib_Drag.f90:287`.
    Gated by XD2 (Cd finite).

---

### Carlson-Hoglund

```
drag-model = Carlson-Hoglund
```

Compressible correlation valid across the full Mach range.

---

## Selection summary

| Keyword | Compressibility | Re range | Notes |
| :--- | :---: | :--- | :--- |
| `Stokes` | — | $\mathrm{Re}\to 0$ | theoretical limit |
| `Schlichting` | — | $\mathrm{Re}\ll 1$ | Oseen correction |
| `Schiller-Naumann` | — | $\mathrm{Re}<800$ | common default |
| `Morsi-Alexander` | — | all | 8-range piecewise |
| `Wen-Yu` | — | all | plateau at Re=1000 |
| `Newton` | — | $10^3$–$2\times10^5$ | constant |
| `Chang` | — | all | ≡ Clift-Gauvin (coded) |
| `Clift-Gauvin` | — | all | continuous |
| `Putnam` | — | all | |
| `Henderson` | yes | all | |
| `Crowe` | yes | all | |
| `Hermsen` | yes | all | |
| `Carlson-Hoglund` | yes | all | |
