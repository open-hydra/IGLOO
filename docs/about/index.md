# About

IGLOO (**I**ntegration of a **G**eneral **L**agrangian **O**ne-way **O**DE set)
is an open-source Fortran solver for Lagrangian particle dynamics in a steady
background gas field. It integrates each particle's ODE state — position,
velocity, temperature, diameter — from injection to exit of the domain in a single
pass. The gas field is either read from a Tecplot file or passed in memory by a
parent solver. IGLOO is the dispersed-phase sub-solver of the
[Hydra](https://github.com/open-hydra/IGLOO) CFD ecosystem and is licensed under
the [GNU General Public License v3](license.md).

---

- [Acknowledgements](acknowledgements.md) — authors, maintainers, and dependency credits
- [License](license.md) — GPL-3.0 full text and summary
