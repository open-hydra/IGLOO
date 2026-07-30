# Quick Start: Stokes Drag Relaxation

This guide walks through running the **`drag-stokes`** end-to-end case — a uniform-gas Stokes-drag relaxation with a closed-form oracle. It is the simplest IGLOO case: one-dimensional axial motion, no evaporation or breakup, verified against an exact analytic solution.

!!! note "Prerequisites"
    - IGLOO built and available at `bin/IGLOO` (see [Installation](installation.md))
    - Python 3 (for the oracle check)

---

## Physics

Particles are injected at velocity $v_0 = 0.1\,u_g$ into a uniform gas moving at $u_g = 10\ \mathrm{m/s}$. Stokes drag forces $C_d = 24/Re$, so the drag law reduces to an exact exponential relaxation:

$$v(t) = u_g + (v_0 - u_g)\,e^{-t/\tau}, \qquad \tau = \frac{\rho_p d_p^2}{18\,\mu}$$

Eliminating time gives a geometry- and time-independent $x(v)$ relation that `check.py` verifies per particle.

| Parameter | Value |
|-----------|-------|
| Gas velocity $u_g$ | $10\ \mathrm{m/s}$ (+x) |
| Gas temperature | $600\ \mathrm{K}$ |
| Gas viscosity $\mu$ | $1.8\times10^{-5}\ \mathrm{Pa\cdot s}$ |
| Particle density $\rho_p$ | $2950\ \mathrm{kg/m^3}$ |
| Particle diameter $d_p$ | $11.89\ \mathrm{\mu m}$ (Dirac distribution) |
| Injection velocity fraction $k_V$ | $0.1$ |
| Stokes relaxation time $\tau$ | $\approx 1.29\ \mathrm{ms}$ |

---

## Case Directory

The case lives at `tests/standard/drag-stokes/`:

```
tests/standard/drag-stokes/
├── input.ini          # solver configuration
├── check.py           # independent oracle (verifies closed-form x(v))
├── INPUT/
│   ├── solfile.tec    # uniform-gas box mesh (60×5×5 cells, Tecplot)
│   ├── bc.txt         # boundary condition table (ATLAS-generated)
│   ├── phase.txt      # material definition (one material: "A")
│   └── properties.dat # thermodynamic property tables
└── OUTPUT/            # created at run time
```

The gas field (`INPUT/solfile.tec`) is a 3D axis-aligned box ($x\in[0,0.15]$, $y\in[0,0.05]$, $z\in[0,0.05]$ m, 60×5×5 cells) with spatially uniform properties. Particles inject from the $x=0$ inlet face.

---

## Running the Case

### Single case

```bash
cd tests/standard/drag-stokes

# create output directory
mkdir -p OUTPUT

# run the solver (set stack and OMP threads first)
ulimit -s unlimited
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=5
../../../bin/IGLOO
```

Output files appear in `OUTPUT/` after the run completes.

### Full e2e suite

Run all oracle-checked cases at once with the test runner:

```bash
./tests/test.sh e2e
```

`test.sh` builds `build/verif/` (refreshing `bin/IGLOO`) and drives the cases through CTest. Each e2e case runs `bin/IGLOO` in its own directory with `ulimit -s unlimited`, `KMP_STACKSIZE=100M`, and `OMP_NUM_THREADS=5`, then runs its oracle. Pass/fail is determined by the oracle exit code, not by whether stderr is empty.

To clean output files and the verification build:

```bash
./tests/test.sh clean
```

---

## Understanding the Output

The solver prints a header, then per-material and per-group progress:

```
 =============================================================================
|                  ///    /////    //       /////    /////                    |
...
 =============================================================================

  OpenMP threads =  5
 >> Background flow field =>  INPUT/solfile.tec
 >> Output fields: gas coupling source & equivalent eulerian
 >> Field mollification ON: 3 binomial pass(es)
 Compute particles dynamics for material: A
     Group  1 => number of particles =   25
     >> ODE system (neq= 7):
       - constant particle mass and size
 Stop condition : All particles out of domain!
```

Output files written to `OUTPUT/`:

| File | Content |
|------|---------|
| `trajectories-A.dat` | Per-cell-crossing state: X Y Z U V W T d_p m_p ID (Tecplot ASCII) |
| `outloc-A.dat` | Exit location, speed, impact angle, area, ID per particle |
| `scatter-A.dat` | Number-density scatter cloud (npdot-weighted point cloud) |
| `source.tec` | Gas-coupling source terms on the mesh |

---

## Verifying the Result

```bash
cd tests/standard/drag-stokes
python3 -B check.py
```

The oracle:

1. Reads `OUTPUT/trajectories-A.dat`.
2. For each particle, anchors the closed-form $x(v)$ relation at the first recorded point.
3. Verifies every subsequent point within the well-conditioned relaxation window ($v\in[v_0,\,0.95\,u_g]$) against a theoretical truncation tolerance derived from the F12.6 output format.
4. Requires at least 20 independently verifiable particles.

A passing run prints:

```
[PASS] 25 particles match the Stokes closed form within theoretical tolerance.
```

!!! warning "Pass criterion"
    The e2e harness uses the oracle (`check.py` exit code), not the legacy "stderr empty" criterion. A run that crashes mid-integration will fail the oracle because `OUTPUT/trajectories-A.dat` will be absent or empty.

---

## Next Steps

- **[User Guide](../user/using.md)** — full workflow, case directory contract, hydra embedding.
- **[Input File Reference](../user/input.md)** — all `input.ini` keys.
- **[V&V Suite](../vv/index.md)** — additional physics cases (heat relaxation, body force, evaporation).
