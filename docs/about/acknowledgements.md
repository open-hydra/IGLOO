# Acknowledgements

---

## Authors

IGLOO was created and is maintained by:

| Name | Role |
|------|------|
| Giacomo Passarani | Original author, present maintainer |
| Marco Grossi | Original author, present maintainer |
| Paolo Zolla | Original author |

See [`AUTHORS.md`](https://github.com/open-hydra/IGLOO/blob/master/AUTHORS.md) at
the repository root for the authoritative list.

---

## Dependencies

IGLOO builds on three companion libraries, included as Git submodules in standalone
mode and shared with the hydra suite when built as a submodule:

| Library | Authors / maintainers | Role |
|---------|-----------------------|------|
| [ORION](https://github.com/MarcoGrossi92/ORION) | Marco Grossi | Multi-format I/O (Tecplot ASCII/binary) |
| [OSlo](https://github.com/MarcoGrossi92/OSlo) | Marco Grossi | ODE solver library (DOPRI5, H-SDIRK4) |
| [FiNeR](https://github.com/szaghi/FiNeR) | Stefano Zaghi | INI configuration file parser |

Optional dependency: **TecIO** (Tecplot binary output), pulled in transitively by
ORION when `--use-tecio` is set.

Documentation is built with [MkDocs Material](https://squidfunk.github.io/mkdocs-material/)
by Martin Donath.
