# Theoretical Guide

IGLOO follows particles and droplets through a steady gas. Each computational particle is a **parcel**: it stands
for $\dot n_p$ real drops per second, all alike. A parcel carries its position, velocity, temperature and, when they
change, its mass and its number rate. These are integrated in time from the inlet to the domain exit, through the
gas the parcel meets.

The physics a parcel feels depends on its material: drag and heat exchange always, then evaporation, breakup,
burning or freezing when the material asks for them. Each combination is one system of ordinary differential
equations.

## The eight systems

| System | Material | What changes besides position, velocity, temperature |
| :--- | :--- | :--- |
| 1. inert | no phase change, no continuous breakup | nothing |
| 2. evaporating | evaporates | the drop's mass |
| 3. breaking | breaks up at a continuous rate | the number rate |
| 4. evaporating and breaking | both | the drop's mass and the number rate |
| 5. burning | burns (Beckstead or direct) | the drop's mass |
| 6. freezing | solidifies | the frozen fraction |
| 7. freezing and breaking | solidifies, and breaks up while liquid | the frozen fraction and the number rate |
| 8. burning and carrying | burns and keeps part of its product | the metal's mass and the carried mass |

The numbers are only labels: the V&V pages use them. The equations of each system are in
[Governing equations](governing-equations.md). Breakup by events (TAB, ETAB) adds no equation: it changes the
parcel between two steps, so a material that breaks this way stays in system 1, 2 or 6.

The input that selects each model: [Materials & Models](../user/models.md).

---

## Pages in this guide

<div class="grid cards" markdown>

-   :material-math-integral: **[Governing equations](governing-equations.md)**

    The parcel's state, the equations of each system, the gas seen by a parcel.

-   :material-arrow-expand-all: **[Drag](drag.md)**

    The 13 drag-coefficient correlations.

-   :material-thermometer: **[Heat transfer](heat.md)**

    Nusselt-number correlations and the grey-body radiative exchange.

-   :material-water-percent: **[Evaporation](evaporation.md)**

    d²-law, CEM, CEM-B, Abramzon–Sirignano and Tonini–Cossali; the Langmuir–Knudsen interface, Stefan blowing,
    the temperature profile inside a drop and the boiling of a superheated drop.

-   :material-fire: **[Metal combustion](combustion.md)**

    Beckstead's $d^n$ law, the direct $d^2$ law, and a drop that carries what it forms.

-   :material-snowflake: **[Solidification](solidification.md)**

    Undercooling, recalescence, the freezing plateau, melting.

-   :material-explosion: **[Breakup](breakup.md)**

    Pilch–Erdman, Reitz–Diwakar, Reitz KH-RT, TAB and ETAB.

-   :material-wall: **[Wall erosion](erosion.md)**

    Oka et al. (2005) and Neilson–Gilchrist (1968).

-   :material-spray: **[Injection](injection.md)**

    Placing parcels on an inlet, and drawing their diameters.

-   :material-transfer: **[Eulerian feedback](eulerian-feedback.md)**

    What the parcels give the gas, and the averaged fields of the particles.

-   :material-blur: **[Mollification](mollification.md)**

    Smoothing the deposited fields without losing what they hold.

-   :material-clock-fast: **[Time integration](time-integration.md)**

    Advancing a parcel cell by cell, and how a step ends.

-   :material-cube-outline: **[Geometry](geometry.md)**

    Which cell holds a parcel, and where it leaves it.

</div>
