---
title: IGLOO
hide:
  - navigation
  - toc
---

<div class="igloo-hero" markdown>

<div class="igloo-hero__text" markdown>

<p class="igloo-hero__eyebrow">HYDRA CFD SUITE</p>

<h1 class="igloo-hero__title">IGLOO</h1>

<p class="igloo-hero__lead">Open-source Lagrangian particle solver. Each particle's ODE state is integrated through a
steady gas field, from injection to the exit of the domain.</p>

<p class="igloo-hero__cmd"><span>$</span> ./install.sh build --compilers=gnu --use-openmp</p>

[Get started](getting-started/index.md){ .md-button .md-button--primary }
[Verification](vv/index.md){ .md-button .igloo-hero__ghost }

<p class="igloo-hero__name"><strong>I</strong>ntegration of a <strong>G</strong>eneral <strong>L</strong>agrangian
<strong>O</strong>ne-way <strong>O</strong>DE set</p>

</div>

<!-- The hero figure is docs/images/home-hero.png, drawn by extra.css as a background: without the file the
     hero shows only its own background. -->
<div class="igloo-hero__shot" aria-hidden="true"><div class="igloo-hero__figure"></div></div>

</div>

<p class="igloo-kicker">The solver</p>

## Particles, one trajectory at a time { .igloo-heading }

<p class="igloo-sub" markdown>IGLOO follows every parcel through the gas, from its injection to the exit of the
domain. Start with the [Quick Start](getting-started/quick-start.md): a verification case runs in a few
minutes.</p>

<div class="igloo-tools" markdown>

[:material-weather-windy: <b>Drag and heat</b> <i>Incompressible and compressible drag laws, Nusselt correlations and grey-body radiation.</i>](theory/drag.md){ .igloo-tool }
[:material-water-percent: <b>Evaporation</b> <i>d²-law, Abramzon-Sirignano and more, with an equilibrium or a Langmuir-Knudsen interface.</i>](theory/evaporation.md){ .igloo-tool }
[:material-fire: <b>Breakup and combustion</b> <i>Five secondary-breakup models with child parcels; aluminium burn with ignition and burnout.</i>](theory/breakup.md){ .igloo-tool }
[:material-grid: <b>Injection and tracking</b> <i>Patch injection with size distributions; cell tracking through multi-block 2-D, axisymmetric and 3-D grids.</i>](theory/injection.md){ .igloo-tool }
[:material-swap-horizontal: <b>Feedback to the gas</b> <i>Eulerian and source fields deposited on the gas mesh, smoothed by a mollifier.</i>](theory/eulerian-feedback.md){ .igloo-tool }
[:material-chip: <b>Parallel runs</b> <i>OpenMP and MPI over particles; the answer does not depend on the number of ranks.</i>](development/parallelization.md){ .igloo-tool }

</div>

<p class="igloo-kicker"></p>

## Hydra CFD Suite { .igloo-heading }

<p class="igloo-sub">Hydra couples solvers of gases, real fluids, solids, particles and radiation in one
simulation. IGLOO is its Lagrangian particle solver.</p>

<div class="hydra" markdown>

<div class="hydra__case" markdown>

<p class="hydra__label">Your case</p>

<p class="hydra__file" markdown>:material-grid: structured mesh</p>

<p class="hydra__file" markdown>:material-database: [databases](https://open-hydra.github.io/ATLAS/databases/)</p>

</div>

<div class="hydra__atlas" markdown>

<p class="hydra__name">ATLAS</p>

<p class="hydra__tools" markdown>[GPB](https://open-hydra.github.io/ATLAS/user-guide/gpb/) [BCB](https://open-hydra.github.io/ATLAS/user-guide/bcb/) [ICB](https://open-hydra.github.io/ATLAS/user-guide/icb/) [STB](https://open-hydra.github.io/ATLAS/user-guide/stb/) [MDB](https://open-hydra.github.io/ATLAS/user-guide/mdb/) [KAnT](https://open-hydra.github.io/ATLAS/user-guide/kant/)</p>

</div>

<div class="hydra__link">
<span class="hydra__link-line"></span>
<span class="hydra__link-bottom">input files for each run</span>
</div>

<div class="hydra__suite" markdown>

<p class="hydra__label" markdown>[Hydra](https://github.com/open-hydra) · the coupled simulation</p>

<div class="hydra__solvers" markdown>

[<b>MOSE</b> compressible flow <i>ideal-gas</i>](https://github.com/open-hydra/MOSE){ .hydra__solver .hydra__solver--gas }
[<b>GROOT</b> radiative transfer <i>no ATLAS phase</i>](https://github.com/open-hydra/GROOT){ .hydra__solver .hydra__solver--radiation }
[<b>ARES</b> real-fluid flow <i>real-fluid</i>](https://github.com/open-hydra/ARES){ .hydra__solver .hydra__solver--fluid }
[<b>FUSS</b> heat conduction in solids <i>solid</i>](https://github.com/open-hydra/FUSS){ .hydra__solver .hydra__solver--solid }
[<b>IGLOO</b> Lagrangian particles <i>dispersed</i>](https://github.com/open-hydra/IGLOO){ .hydra__solver .hydra__solver--particle .is-here }
[<b>ICE</b> Eulerian particles <i>dispersed</i>](https://github.com/open-hydra/ICE){ .hydra__solver .hydra__solver--particle }

</div>

</div>

<div class="hydra__libs" markdown>

<p class="hydra__label">Built on</p>

[<b>ORION</b> Tecplot, VTK and PLOT3D I/O](https://MarcoGrossi92.github.io/ORION/){ .hydra__lib }
[<b>OSLO</b> ODE integrators](https://github.com/MarcoGrossi92/OSLO){ .hydra__lib }
[<b>FiNeR</b> INI file parser](https://github.com/szaghi/FiNeR){ .hydra__lib }

</div>

</div>

<div class="igloo-links" markdown>

[Quick Start](getting-started/quick-start.md)
[User Guide](user/index.md)
[Input Parameters](user/registry.md)
[Theory](theory/index.md)
[Verification](vv/index.md)
[Developer Guide](development/index.md)

</div>
