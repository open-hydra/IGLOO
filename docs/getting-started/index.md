# Getting Started

Welcome to the IGLOO getting-started guide. This section covers everything needed to build the solver and run your first end-to-end case.

!!! info "IGLOO within Hydra"
    IGLOO is the Lagrangian particle solver of the **Hydra** CFD ecosystem. The pre-processor **ATLAS** (distributed separately) produces the boundary-condition and material input files that IGLOO requires. See the [Overview](../overview.md) for the full ecosystem picture.

<div class="grid cards" markdown>

-   :material-download:{ .lg .middle } __Installation__

    ---

    Build IGLOO from source — standalone or as a hydra submodule

    [:octicons-arrow-right-24: Install IGLOO](installation.md)

-   :material-rocket-launch:{ .lg .middle } __Quick Start__

    ---

    Run the Stokes-drag verification case end to end in minutes

    [:octicons-arrow-right-24: Quick start tutorial](quick-start.md)

</div>

## Scope of This Section

1. **Installation** — prerequisites, standalone and hydra-submodule build modes, incremental rebuild.
2. **Quick Start** — run one oracle-checked end-to-end case (`drag-stokes`) and confirm the solver integrates particles correctly.
