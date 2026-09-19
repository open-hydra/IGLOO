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

IGLOO reads only the `[IGLOO-*]` sections above. The sections below are the *preprocessor's* input:
whatever IGLOO needs from them reaches it through the files ATLAS writes (`phase.txt`,
`properties.dat`, `bc.txt`). The per-material model keys (`evaporation`, `combustion`, `alpha-e`,
...) are `[GPB-Phase*]` input too: ATLAS GPB writes them as `key=value` tokens after `<name> <groups>`
on the material line of the phase file, and IGLOO reads them there (an unknown key stops the run).

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
| `phase` | string | `''` | ATLAS name of the condensed phase to read (`[GPB-Phase*] name`): the files become `INPUT/<phase>-{phase.txt,properties.dat,bc.txt}` and every output gets the `<phase>-` prefix. Absent keeps the prefix in force (unnamed `INPUT/phase.txt` standalone; the parent app's assignment under hydra-MI2). Must not contain `-` |
| `gas-order` | integer | `2` | Gas interpolation order: `1` = cell-center value, `2` = second-order reconstruction |
| `out-file` | string | `E+S` | Output fields: `E` = Eulerian only, `S` = source only, `E+S` (also `S+E`, `ES`, `SE`, `ALL`, `both`, case-insensitive; or absent) = both; any other token is refused |
| `fsample-traj` | integer | `100` | Scatter-cloud density: nominal points per injection stream (sets the droplets-per-point quantum); trajectory rows follow `print-dcell`/`print-dtime`, not this key |
| `print-dcell` | integer | `1` | Console print frequency in cell crossings |
| `print-dtime` | real | `−1` | Console print frequency in seconds; `−1` disables time-based printing |
| `mdot-max` | real | `0` | Maximum mass flow rate per particle [g/s]; used to auto-size injection spacing |
| `out-traj` | string | `on` | Enable trajectory output (`off` to disable) |
| `out-scatter` | string | `on` | Enable scatter-cloud output (`off` to disable) |
| `seed` | integer | `42` | RNG seed for stochastic diameter sampling |
| `mollify` | string | `on` | Enable field mollification smoother (`off` to disable) |
| `mollify-passes` | integer | `8` | Binomial smoother passes; `0` = off |
| `body-accel` | real(3) | `0 0 0` | Uniform body acceleration [m/s²]: `gx gy gz`; absent or all-zero = no-op |

!!! warning "on/off switches are strings, not Fortran logicals"
    FiNeR's `get(logical)` only accepts `T` or `F`. The `out-traj`, `out-scatter`, and `mollify` keys are parsed as strings; accepted off-tokens are `off`, `false`, `no`, `0`, `F`, `f` (any case). Anything else is treated as on.

### `[IGLOO-Models]`

Selects physical model closures.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `drag` | string | — | Drag law: one of the 13 correlations in [Drag](../theory/drag.md) (`Stokes`, `Schiller-Naumann`, `Morsi-Alexander`, `Henderson`, `Crowe`, ...) |
| `heat` | string | — | Nusselt correlation: `Ranz-Marshall`, `Kavanau-Drake`, `JAXA1`–`JAXA4` |
| `evaporation` | string | absent | Evaporation model (`d2-law`, `CEM`, `CEM-B`, `ASM`, `TC`); when present, enables phase change |
| `interface` | string | `VLE` | Surface state: `VLE` equilibrium or `LK` Langmuir–Knudsen non-equilibrium |
| `blowing` | string | `none` | Stefan-blowing reduction of the convective heat: `none` or `LK` |
| `breakup` | string | absent | Breakup model (`Pilch-Erdman`, `Reitz-Diawakar`, `Reitz-KHRT`, `TAB`, `ETAB`); when present, enables secondary breakup |

Breakup model–specific tuning constants (`B0`, `B1`, `Cs`, etc.) are also read from `[IGLOO-Models]` when a breakup model is active; see the [Parameter Registry](registry.md) for per-model keys. Unknown model tokens, and the reserved `evaporation = LEB`, are refused at setup.

### `[IGLOO-Properties]`

Constant evaporation and breakup properties, one value per material in the order they appear in `INPUT/phase.txt`. Required only when the corresponding model is active and the property is not already in `INPUT/properties.dat`.

| Key | Description |
|-----|-------------|
| `psat` | Saturation pressure [Pa] |
| `Mv` | Vapor molar mass [kg/kmol] |
| `Lv` | Latent heat of vaporization [J/kg] |
| `Tboil` | Boiling temperature [K] |
| `cpv` | Vapor specific heat [J/(kg·K)] |
| `Le` | Lewis number |
| `Yinf` | Far-field vapor mass fraction |
| `sigma` | Surface tension [N/m] (breakup) |
| `mu` | Dynamic viscosity [Pa·s] (breakup) |

### `[IGLOO-BC]`

Controls injection particle placement.

**Boundary-patch injection** (method `FB`, selected when no `x`/`y`/`z` is given): particles are placed on all faces tagged as inlet in `INPUT/bc.txt`. Controlled by:

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `ds` | real | `0` | Global injection spacing [cm]. When `> 0`, uses the advancing-front (2D sweep or 3D BFS hexagonal packing) algorithm. When `= 0`, one particle per injection cell center. |
| `ds-degen` | real | `0` | Degeneracy floor [cm]: skip injection in boundary cells whose tangential size is below this threshold. |
| `fsample` | integer | `1` | Injection cell subsampling factor. |

**Assigned-position injection** (method `DB`, selected by the presence of `x`/`y`/`z`): provide explicit coordinates instead of boundary-face scanning. Requires `x`, `y`, `z` (or any combination), plus `mdot` and `diam`.

| Key | Description |
|-----|-------------|
| `x`, `y`, `z` | Injection point coordinates (scalar or array) [m] |
| `mdot` | Mass flow rate per particle [kg/s] (scalar or array) |
| `diam` | Particle diameter [m] (scalar or array) |
| `temp0` | Initial temperature [K] (optional, defaults to gas temperature) |
| `up`, `vp`, `wp` | Initial velocity components [m/s] (optional) |

### `[IGLOO-ODE]`

Selects and tunes the ODE integrator.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `ode-solver` | string | `H-sdirk4` | Integrator: `H-sdirk4` (implicit SDIRK4) or `H-dopri5` (explicit Dormand–Prince 5(4)); any other token is refused |
| `relative-tol` | real | `1e-10` | Relative ODE tolerance |
| `absolute-tol` | real | `1e-10` | Absolute ODE tolerance |
| `max-steps-ode` | integer | `100000` | Maximum internal steps per integrator call (`iopt(1)`) |

`H-sdirk4` is recommended for stiff cases (evaporation, small Stokes number). `H-dopri5` is faster for non-stiff drag-only cases.

---

## Gotchas

!!! warning "`on`/`off` switches must be strings, not Fortran logicals"
    FiNeR's built-in `get(logical)` accepts only `T` or `F`. IGLOO boolean-like keys (`out-traj`, `out-scatter`, `mollify`) are parsed as character strings. Writing `out-traj = .true.` will be misinterpreted; use `out-traj = on` or simply omit the key (default is on).

!!! warning "`[GPB-Phase1] rho` is ignored at run time"
    Particle density at run time comes from `INPUT/properties.dat`, not from the `rho` key in `[GPB-Phase1]`. The `[GPB-Phase*]` sections are ATLAS input; IGLOO uses the property tables that ATLAS writes to `INPUT/properties.dat`.

!!! warning "`properties.dat` zones bind to materials by ORDER"
    Zone *i* of `INPUT/properties.dat` is the material on line *i* of `INPUT/phase.txt`; the zone title (`ZONE T="..."`) is not read. ATLAS writes both files in material order, so a generated case is always consistent — the trap is a hand-edited or regenerated `phase.txt` paired with a stale `properties.dat`, which is a silent density/cp swap. IGLOO refuses a file whose zone count differs from the material count or whose zones do not share one temperature table.

!!! note "The enthalpy column carries a datum"
    Column 4 of `properties.dat` is named `Enthalpy` (relative: `cp·T`, or the SP-database integral from `Tmin`) or `Enthalpy_abs` (absolute: formation enthalpy included — thermo tables, or a fixed `cp` with `[GPB-Phase*] h0`). For a constant-`cp` material IGLOO integrates the temperature and takes the table's level once, `hOff = h(Tmin) − cp·Tmin` (0 for a relative table), and the gas coupling source credits the transferred mass at `cp·T + hOff`; a variable-`cp` material integrates the enthalpy state read from the table, which carries its datum by itself. The datum only matters when mass is transferred to a gas solver whose energy is absolute (hydra-MI2): a relative table then leaves the coupler to correct the level. A relative-tagged constant-`cp` column that is not `cp·T` is refused at setup. The datum is printed per material at setup (`enthalpy datum absolute|relative (hOff = ...)`).
