# Input File

IGLOO is configured through a single **INI-format** file called `input.ini`, located in the case root directory. The file is parsed by FiNeR and organized into named sections.

The sections fall into two categories:

- **IGLOO-native sections** — read directly by the solver at run time.
- **ATLAS-generated sections** — consumed by ATLAS to produce `INPUT/bc.txt`, `INPUT/phase.txt`, and `INPUT/properties.dat`. IGLOO ignores these sections at run time; they are present in the same file as a matter of convention.

For the exhaustive list of every IGLOO-native key, its type, default, and allowed values, see the auto-generated **[Parameter Registry](registry.md)**.

---

## File Structure

```ini
[SECTION-NAME]
parameter = value
```

Parameters not specified take their default values. Unknown sections are silently ignored by FiNeR.

!!! warning "Parameter names are case-sensitive"
    An incorrectly spelled key is silently ignored and the default is used. Check the [Parameter Registry](registry.md) for exact spelling.

---

## IGLOO-Native Sections

| Section | Description |
|---------|-------------|
| `[IGLOO-General]` | Gas file, output modes, trajectory sampling, mollification, body acceleration |
| `[IGLOO-Models]` | Drag law, heat correlation, evaporation model, breakup model |
| `[IGLOO-Properties]` | Constant scalar properties for evaporation/breakup (per material) |
| `[IGLOO-BC]` | Injection spacing (`ds`) or explicit particle coordinates |
| `[IGLOO-ODE]` | ODE solver selection, tolerances, step count |

## ATLAS-Generated Sections (not parsed by IGLOO)

| Section | Role |
|---------|------|
| `[GPB-Phase*]` | Per-material condensed-phase properties (type, density, heat capacity) for ATLAS BC building |
| `[BCB-Block*]` | Maps each block face to a named patch (e.g. `face1 = in`) |
| `[<patch-name>]` | Per-patch BC definition (type, `krho`, `dp`, etc.) used by ATLAS to write `bc.txt` |

These are kept in the case file because they document the **injected phase** — its
properties and the patch each face carries. Mesh-generation sections (`[GRIB-*]`) are
**not** kept: IGLOO reads its geometry and gas from the Tecplot solution named by
`[IGLOO-General] gas-file`, which is another solver's *output*, so there is nothing left
for a mesh generator to specify.

---

## Worked Example

The `tests/standard/drag-stokes/input.ini` case runs a uniform-gas Stokes relaxation — the simplest possible IGLOO configuration:

```ini
; uniform-gas Stokes drag: drag=Stokes forces the exact exponential relaxation.
; gas-order=1 reads cell-center values; out-file=S writes source only.

[IGLOO-General]
gas-file    = INPUT/solfile.tec
gas-order   = 1
print-dcell = 1
out-file    = S

[IGLOO-ODE]
ode-solver   = H-sdirk4
relative-tol = 1e-11
absolute-tol = 1e-11

[IGLOO-Models]
drag = Stokes
heat = Ranz-Marshall

[IGLOO-BC]
ds = 10.0

[GPB-Phase1]           ; ATLAS section — sets rho/cp for ATLAS bc.txt generation
type = condensed-dispersed
rho  = 2950
cp   = 1250

[BCB-Block1]           ; ATLAS section — maps face labels to patch names
face1 = in
face2 = out
face3 = wall
face4 = wall
face5 = wall
face6 = wall

[out]                  ; ATLAS patch definition
type = outlet

[in]                   ; ATLAS patch definition — krho, dp go into bc.txt col 1/col 6
type = inlet
krho = 0.34
dp   = 11.89e-6

[wall]
type = wall
q    = 0.0
```

---

## Section Reference

### `[IGLOO-General]`

Controls the gas field source, output modes, and miscellaneous run-time options.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `gas-file` | string | — | Path to the Tecplot gas field (ignored when `external_gas` is provided by hydra) |
| `gas-order` | integer | `2` | Gas interpolation order: `1` = cell-center value, `2` = second-order reconstruction |
| `out-file` | string | both | Output fields: `E` = Eulerian only, `S` = source only; any other value = both |
| `fsample-traj` | integer | `100` | Trajectory sampling: record every nth particle-state (scatter and trajectory output) |
| `print-dcell` | integer | `1` | Console print frequency in cell crossings |
| `print-dtime` | real | `−1` | Console print frequency in seconds; `−1` disables time-based printing |
| `mdot-max` | real | `0` | Maximum mass flow rate per particle [g/s]; used to auto-size injection spacing |
| `out-traj` | string | `on` | Enable trajectory output (`off` to disable) |
| `out-scatter` | string | `on` | Enable scatter-cloud output (`off` to disable) |
| `seed` | integer | `42` | RNG seed for stochastic diameter sampling |
| `mollify` | string | `on` | Enable field mollification smoother (`off` to disable) |
| `mollify-passes` | integer | auto | Binomial smoother passes; default suppresses ≤4-cell deposition noise |
| `body-accel` | real(3) | `0 0 0` | Uniform body acceleration [m/s²]: `gx gy gz`; absent or all-zero = no-op |

!!! warning "on/off switches are strings, not Fortran logicals"
    FiNeR's `get(logical)` only accepts `T` or `F`. The `out-traj`, `out-scatter`, and `mollify` keys are parsed as strings; accepted off-tokens are `off`, `false`, `no`, `0`, `F`, `f` (any case). Anything else is treated as on.

### `[IGLOO-Models]`

Selects physical model closures.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `drag` | string | — | Drag law: `Stokes`, `Morsi-Alexander`, `Crowe`, `Hermsen`, `Henderson`, `Putnam` |
| `heat` | string | — | Nusselt correlation: `Ranz-Marshall`, `Kavanau-Drake` |
| `evaporation` | string | absent | Evaporation model; when present, enables phase change |
| `breakup` | string | absent | Breakup model; when present, enables secondary breakup |

Breakup model–specific tuning constants (`B0`, `B1`, `Cs`, etc.) are also read from `[IGLOO-Models]` when a breakup model is active; see the [Parameter Registry](registry.md) for per-model keys.

### `[IGLOO-Properties]`

Constant evaporation and breakup properties, one value per material in the order they appear in `INPUT/phase.txt`. Required only when the corresponding model is active and the property is not already in `INPUT/properties.dat`.

| Key | Description |
|-----|-------------|
| `psat` | Saturation pressure [Pa] |
| `Mv` | Vapor molecular weight [kg/mol] |
| `Lv` | Latent heat of vaporization [J/kg] |
| `Tboil` | Boiling temperature [K] |
| `cpv` | Vapor specific heat [J/(kg·K)] |
| `Le` | Lewis number |
| `Yinf` | Far-field vapor mass fraction |
| `sigma` | Surface tension [N/m] (breakup) |
| `mu` | Dynamic viscosity [Pa·s] (breakup) |

### `[IGLOO-BC]`

Controls injection particle placement.

**Boundary-patch injection** (default, `method = FB`): particles are placed on all faces tagged as inlet in `INPUT/bc.txt`. Controlled by:

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `ds` | real | `0` | Global injection spacing [cm]. When `> 0`, uses the advancing-front (2D sweep or 3D BFS hexagonal packing) algorithm. When `= 0`, one particle per injection cell center. |
| `ds-degen` | real | `0` | Degeneracy floor [cm]: skip injection in boundary cells whose tangential size is below this threshold. |
| `fsample` | integer | `1` | Injection cell subsampling factor. |

**Assigned-position injection** (`method = DB`): provide explicit coordinates instead of boundary-face scanning. Requires `x`, `y`, `z` (or any combination), plus `mdot` and `diam`.

| Key | Description |
|-----|-------------|
| `x`, `y`, `z` | Injection point coordinates (scalar or array) [m] |
| `mdot` | Mass flow rate per particle [g/s] (scalar or array) |
| `diam` | Particle diameter [m] (scalar or array) |
| `temp0` | Initial temperature [K] (optional, defaults to gas temperature) |
| `up`, `vp`, `wp` | Initial velocity components [m/s] (optional) |

### `[IGLOO-ODE]`

Selects and tunes the ODE integrator.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `ode-solver` | string | `H-sdirk4` | Integrator: `H-dopri5` (explicit DOPRI5) or `H-sdirk4` (implicit H-SDIRK4) |
| `relative-tol` | real | `1e-10` | Relative ODE tolerance |
| `absolute-tol` | real | `1e-10` | Absolute ODE tolerance |
| `max-steps-ode` | integer | `100000` | Maximum ODE steps per cell crossing |

`H-sdirk4` is recommended for stiff cases (evaporation, small Stokes number). `H-dopri5` is faster for non-stiff drag-only cases.

---

## Gotchas

!!! warning "`on`/`off` switches must be strings, not Fortran logicals"
    FiNeR's built-in `get(logical)` accepts only `T` or `F`. IGLOO boolean-like keys (`out-traj`, `out-scatter`, `mollify`) are parsed as character strings. Writing `out-traj = .true.` will be misinterpreted; use `out-traj = on` or simply omit the key (default is on).

!!! warning "`[GPB-Phase1] rho` is ignored at run time"
    Particle density at run time comes from `INPUT/properties.dat`, not from the `rho` key in `[GPB-Phase1]`. The `[GPB-Phase*]` sections are ATLAS input; IGLOO uses the property tables that ATLAS writes to `INPUT/properties.dat`.
