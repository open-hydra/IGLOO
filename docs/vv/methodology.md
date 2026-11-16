# Methodology

How the end-to-end cases are built and why their gates can be trusted. The unit tests follow the
same rules on a single routine.

---

## Four rules

1. **Independent references.** No reference value comes from the code under test.
2. **Tolerances from theory.** Every tolerance is derived (output precision, truncation bounds,
   statistical margins) before the test is run.
3. **Production is never patched to pass a test.** When a test disagrees with the code, the
   disagreement is a defect to investigate; the oracle wins until an *external* reference (the
   paper, a closed form, a second implementation) proves it wrong.
4. **Every gate is falsifiable.** A new gate is shown to fail on a deliberately broken build (a
   wrong constant, a reverted rotation, a fold applied to position only) before it is accepted.

---

## The box case

Most end-to-end cases run on one shared box.

| Parameter | Value |
|---|---|
| Domain | $0.15 \times 0.05 \times 0.05$ m |
| Cells | $60 \times 5 \times 5$, $\Delta x = 2.5$ mm |
| Inlet | the $x = 0$ face, one parcel stream per inlet cell (25 streams) |
| Other faces | `bcdef = 400`; the parcels leave through $x = 0.15$ m |
| Gas | uniform, $u_g = 10$ m/s along $x$, $T_g = 600$ K, $\mu = 1.8\times10^{-5}$ Pa s, $k_g = 0.026$ W/(m K), unless stated |
| Particles | $d_0 = 11.89$ µm, $\rho_p = 2950$ kg/m³ unless stated |

Two inlet ratios isolate one physics at a time: an injection velocity equal to the gas velocity
removes the slip (drag is a no-op, $\mathrm{Re} = 0$, $\mathrm{Nu} = 2$), and an injection
temperature equal to the gas temperature removes the heat exchange.

Other geometries:

- a non-uniform compressive gas on the box (`vie-plait`);
- a one-cell annular wedge with solid-body swirl (the `swirl-wedge` family);
- a one-cell planar square in solid-body rotation (`vortex-cloud`);
- a real MOSE nozzle solution on an axisymmetric wedge (`db-2daxi`, `axis-200`, `wedge-fold`);
- the box cut into two blocks, joined by a connection or by overlapping chimera faces.

---

## The oracle

Each case has its own Python oracle, `check.py`. It builds the prediction from the case inputs
alone (particle properties, gas constants, body force); the run supplies only the points to test.
Where time would enter, the closed form is written in position — $x(v)$, $T(x)$ — and anchored at
each parcel's first row, so the ODE solver's step choice cannot affect the comparison.

The **kernel along the run** cases are the one exception: their rate law is integrated along the
drop temperature measured in the run, because the coupled equations have no closed form. The oracle
still shares no code with IGLOO.

A case passes only when the solver exits 0, every number in its output files is finite, and the
oracle exits 0. The solver's standard error is not read.

---

## Tolerances

The particle files print with format F12.6, so every coordinate, velocity and temperature carries a
half-unit truncation $\delta = 0.5\times10^{-6}$. For a residual $r = q_\mathrm{meas} -
q_\mathrm{pred}(x_\mathrm{meas})$ the tolerance carries that truncation through the closed form:

$$
\mathrm{tol}(x) = \delta\,\bigl(2 + \text{sensitivity terms}\bigr) + \varepsilon_\mathrm{int},
$$

the sensitivity terms being $|\partial r/\partial q_i|$ at the anchor and test points, and
$\varepsilon_\mathrm{int} = 10^{-8}$ a floor for the integrator (at `rtol = atol = 1e-11` the
SDIRK4 error sits three orders below the print, so the floor never binds).

Two guards keep a pass meaningful:

- **The resolved window.** A row is tested only where the print can still tell the prediction from
  its limit — e.g. $v \le 0.95\,u_g$ for a velocity relaxation, $|T - T_g| > 1$ K for a temperature
  relaxation. Near saturation the sensitivity diverges and a print digit maps to any error.
- **Non-vacuity.** At least 20 parcels must be tested in that window, so a transient too short to
  sample cannot pass by default.

Where a case carries a competing model (VLE against Langmuir–Knudsen, CEM against Tonini–Cossali),
the margin between the two is itself gated, so the comparison cannot pass vacuously.

---

## Figures

Each figure is drawn by the case's `verify.py` from the last run and copied into the documentation by
the full test run. A figure is a picture of a gate, never a gate itself: some draw a reference
(experiment, a correlation outside its fit) that is plotted but not gated, and each case card says
so when it does.
