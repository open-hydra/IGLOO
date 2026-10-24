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
| `[IGLOO-Properties]` | Constant scalar properties for evaporation/breakup and the radiation emissivity (per material) |
| `[IGLOO-BC]` | Injection spacing (`ds`) or explicit particle coordinates |
| `[IGLOO-ODE]` | ODE solver selection, tolerances, step count |

## ATLAS-Generated Sections (not parsed by IGLOO)

IGLOO reads only the `[IGLOO-*]` sections above. The sections below are the *preprocessor's* input:
whatever IGLOO needs from them reaches it through the files ATLAS writes (`phase.txt`,
`properties.dat`, `bc.txt`). The per-material model keys (`evaporation`, `combustion`, `solidification`,
`breakup`, `alpha-e`, ...), the solidification constants (`h-fus`, `cp-solid`, `T-melt`, `T-nuc`, `rho-solid`; see
[Solidification](../theory/solidification.md#inputs)) and the breakup constants (`<tag>-<key>`, as `KHRT-B1`; table under
[`[IGLOO-Models]`](#igloo-models)) are `[GPB-Phase*]` input too: ATLAS GPB writes them as `key=value` tokens after
`<name> <groups>` on the material line of the phase file, and IGLOO reads them there (an unknown key stops the run).
The `<groups>` summed over the material lines is the family count, and therefore the number of table
copies `bc.txt` may carry (one copy feeds every family; see [Several families](boundary-conditions.md#several-families)).

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
| `print-dcell` | integer | `1` | Trajectory row frequency in cell crossings |
| `print-dtime` | real | `−1` | Trajectory row interval [s]; `−1` disables time-based rows |
| `time-end` | real | `−1` | Parcel end time [s]: every parcel still in the domain stops there, gets a last trajectory row at that time, and its state goes to `snapshot-<mat>.dat`; `≤ 0` or absent = off (integrate until the parcel leaves the domain) |
| `out-time` | string | `off` | `on` appends the parcel time [s] as the last column of every trajectory and exit row |
| `plane-x` | real list | — | Planes normal to x at the listed positions [m], space-separated (`plane-x = 0.001 0.18`): `planes-<mat>.dat` holds one record per crossing of each plane ([Plane crossings](using.md#plane-crossings)); absent = no file. A token that is not a finite number stops the run |
| `mdot-max` | real | `0` | Maximum mass flow rate per particle [g/s]; used to auto-size injection spacing |
| `out-traj` | string | `on` | Enable trajectory output (`off` to disable) |
| `out-scatter` | string | `on` | Enable scatter-cloud output (`off` to disable) |
| `seed` | integer | `42` | RNG seed for stochastic diameter sampling |
| `mollify` | string | `on` | Enable field mollification smoother (`off` to disable) |
| `mollify-passes` | integer | `8` | Binomial smoother passes; `0` = off |
| `body-accel` | real(3) | `0 0 0` | Uniform body acceleration [m/s²]: `gx gy gz`; absent or all-zero = no-op |
| `timing` | string | `off` | `on` prints `[IGLOO-timing]` lines after each output write: wall time per phase, parcel count and integration time per thread and per rank, the rank imbalance, the allocated array sizes (see [Parallelization](../development/parallelization.md#per-phase-timers)) |

!!! warning "on/off switches are strings, not Fortran logicals"
    FiNeR's `get(logical)` only accepts `T` or `F`. The `out-traj`, `out-scatter`, `out-time`, `timing` and `mollify` keys are parsed as strings; accepted off-tokens are `off`, `false`, `no`, `0`, `F`, `f` (any case). Anything else is treated as on.

### `[IGLOO-Models]`

Selects physical model closures.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `drag` | string | — | Drag law: one of the 13 correlations in [Drag](../theory/drag.md) (`Stokes`, `Schiller-Naumann`, `Morsi-Alexander`, `Henderson`, `Crowe`, ...), or `NoDrag` (C_d = 0, no momentum exchange) |
| `heat` | string | — | Nusselt correlation: `Ranz-Marshall`, `Kavanau-Drake`, `JAXA1`–`JAXA4`, or `NoHeat` (Nu = 0, no convective heat exchange) |
| `evaporation` | string | absent | Evaporation model (`d2-law`, `CEM`, `CEM-B`, `ASM`, `TC`); when present, enables phase change |
| `interface` | string | `VLE` | Surface state: `VLE` equilibrium or `LK` Langmuir–Knudsen non-equilibrium |
| `blowing` | string | `none` | Stefan-blowing reduction of the convective heat: `none` or `LK` |
| `liquid-conduction` | string | `ITC` | Heat conduction inside an evaporating droplet: `ITC` (a uniform temperature) or `P2T` (a parabolic profile; the material line carries `k-liq` and, on a material that does not break, `mu-liq`); see [Evaporation](../theory/evaporation.md#heat-conduction-inside-the-droplet) |
| `boiling` | string | `clamp` | An evaporating droplet at or above its boiling point at the gas pressure: `clamp` (the gas-side model with its surface fraction capped) or `ZGR` (the Zuo–Gomes–Rutland boiling branch); see [Evaporation](../theory/evaporation.md#boiling-of-a-superheated-droplet) |
| `breakup` | string | absent | Breakup model (`Pilch-Erdman`, `Reitz-Diawakar`, `Reitz-KHRT`, `TAB`, `ETAB`); when present, enables secondary breakup. The default of every material: a material's `breakup=<model>` token on its phase-file line overrides it |

Each breaking material reads the constants of its own model: its `<tag>-<key>` token on the phase-file line, else
the `[IGLOO-Models]` key of that name, else the default (see the [Parameter Registry](registry.md)). TAB and ETAB
each take their own `Comega`, `Cmu` and `WeCrit` tokens:

| `breakup=` | constant tokens |
|---|---|
| `Pilch-Erdman` | `PE-Cd`, `PE-B` |
| `Reitz-Diawakar` | `RD-WeBag`, `RD-Cb`, `RD-Cstrip`, `RD-Cs` |
| `Reitz-KHRT` | `KHRT-B0`, `KHRT-B1`, `KHRT-Ctau`, `KHRT-CRT`, `KHRT-mShedLim`, `KHRT-WeLimit` |
| `TAB` | `TAB-Comega`, `TAB-Cmu`, `TAB-WeCrit`, `TAB-method` (1 or 2), `TAB-n` |
| `ETAB` | `ETAB-k1`, `ETAB-k2`, `ETAB-WeCrit`, `ETAB-WeTrans`, `ETAB-Comega`, `ETAB-Cmu` |

A token of a model the material does not run is accepted and unused; a material that carries a constant token
logs its resolved set at setup. Unknown model tokens are refused at setup.

### `[IGLOO-Properties]`

Constant evaporation and breakup properties and the emissivity of the radiative exchange, one value per material in the order they appear in `INPUT/phase.txt`. The evaporation properties are required only when evaporation is active and the property is not already in `INPUT/properties.dat`; `sigma` and `mu` only when a material breaks and its zone of `properties.dat` has no `Sigma` / `Viscosity` column (the entry of a material that does not break is unused).

| Key | Description |
|-----|-------------|
| `psat` | Legacy, accepted and ignored — the saturation pressure is computed from Clausius-Clapeyron (`Lv`, `Mv`, `boiling-temperature`) or read from the `Psat` column of `properties.dat` |
| `Mv` | Vapor molar mass [kg/kmol], or the `Mv` column of `properties.dat`; 0 for a material that has the column |
| `Lv` | Latent heat of vaporization [J/kg], or the `Lv` column of `properties.dat`; 0 for a material that has the column |
| `boiling-temperature` (alias `Tboil`) | Boiling temperature [K]; give one name or the other, not both |
| `cpv` | Vapor specific heat [J/(kg·K)], or the `Cpv` column of `properties.dat`; 0 for a material that has the column |
| `Le` | Lewis number |
| `Yinf` | Far-field vapor mass fraction |
| `sigma` | Surface tension [N/m] of a breaking material, or the `Sigma` column of `properties.dat`; 0 for a material that has the column |
| `mu` | Liquid dynamic viscosity [Pa·s] of a breaking material, or the `Viscosity` column of `properties.dat`; 0 for a material that has the column |
| `emissivity` | Surface emissivity of the grey-body exchange with the gas temperature, $\varepsilon\,\sigma\,\pi d^2 (T_g^4 - T_p^4)$ (see [Radiation](../theory/heat.md#radiation)); default 0 = off; refused below 0, cut to 1.0 with a warning above 1 |

IGLOO reads `emissivity` from `[IGLOO-Properties]` only. ICE models the same exchange through `[ICE-Physics] emissivity`: a hydra case that runs both condensed solvers sets the two keys together, with the same value for a material both phases carry; ICE's key alone leaves IGLOO's exchange off.

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

## Property tables (`INPUT/properties.dat`)

One zone per material, in the order of the material lines of `INPUT/phase.txt`; every zone on
the same temperatures, one row per integer kelvin (`T = Tmin, Tmin+1, ..., Tmax`). The first line
holding `VARIABLES` names the columns:

```text
VARIABLES = "Temperature", "Cp", "Density", "Enthalpy", "Psat", "Lv", "Mv", "Cpv", "Sigma", "Viscosity"
```

| column | unit | |
|---|---|---|
| `Temperature` | K | first; the nodes |
| `Cp` | J/(kg K) | required |
| `Density` | kg/m³ | required; the liquid's for a solidifying material (the solid's is its `rho-solid` token, when given) |
| `Enthalpy` or `Enthalpy_abs` | J/kg | required, exactly one: `Enthalpy` is relative (`cp·T`, or the integral of cp), `Enthalpy_abs` includes the formation enthalpy (the datum, see below) |
| `Psat` | Pa, absolute | optional: the saturation pressure of an evaporating material |
| `Lv` | J/kg | optional, with `Psat`: the latent heat of vaporization |
| `Mv` | kg/kmol | optional, with `Psat`: the vapour molar mass, the same on every row |
| `Cpv` | J/(kg K) | optional, with `Psat`: the vapour specific heat |
| `Sigma` | N/m | optional: the surface tension of a breaking material |
| `Viscosity` | Pa s | optional: the liquid dynamic viscosity of a breaking material |

After `Temperature` the columns may come in any order; the lowercase names are accepted too
(`PSAT`, `LV`, `MV`, `CPV`, `SIGMA`, `VISCOSITY` as well). A column whose property is constant in T is read as a constant; one that varies
becomes a table, linear between the nodes. Outside the table the density holds its end values and
the enthalpy continues its first and last segments, both as in ICE. With a varying `Cp` the particle
state is the enthalpy, and its temperature is read back through the same table, so the table may
start at any temperature (see [Governing equations](../theory/governing-equations.md)).

**`Psat`.** An evaporating material takes its saturation pressure from the column, linear between
the nodes and at its end values outside the table; a column of zeros keeps Clausius-Clapeyron
(`Lv`, `Mv`, `boiling-temperature` of `[IGLOO-Properties]`), and a material that does not evaporate
ignores it. With the column, `Lv` is the latent-heat sink only. Setup prints
`p_sat tabulated from properties.dat (Tmin..Tmax K), psat(boiling-temperature)/Patm = ...` per
material that uses it. ATLAS GPB writes the column from a liquid/vapour pair of its thermo database
(`[GPB-Phase*] psat-vapour`), zeros for a material without a pair.

**`Lv`, `Mv`, `Cpv`.** ATLAS GPB writes them beside `Psat` from the same pair (the latent heat with
the enthalpy difference that `Psat`'s Clausius-Clapeyron slope uses), zeros for a material without
a pair. An evaporating material takes `Lv` and `Cpv` at the droplet temperature in every right-hand
side (linear between the nodes, the end values outside, as `Psat`) — the latent-heat sink, the
d²-law transfer number and the vapour heat capacity of ASM and TC — and `Mv` as its constant. Each
column replaces its `[IGLOO-Properties]` key: a material whose column is non-zero must give 0 (or
nothing) for that key, and a material without the column needs the key. For a tabulated `Lv` the
material's scalar latent heat is `Lv(298.15 K)` (the enthalpy datum's temperature). Setup prints
`vapour data tabulated from properties.dat: Lv, Mv, Cpv (Lv(298.15 K) = ..., Lv(boiling-temperature) = ...)`
per material that uses them; a material that does not evaporate ignores them.

**`Sigma`, `Viscosity`.** A breaking material takes its surface tension and its liquid viscosity
from these columns at the droplet temperature, in every right-hand side and at every breakup event
(linear between the nodes, the end values outside the table, as `Psat`); a column of zeros is
absent, and a material that does not break ignores them. Each column replaces its
`[IGLOO-Properties]` key: a material whose column is non-zero must give 0 (or nothing) for that
key, and a breaking material without the column needs the key > 0. Setup prints
`sigma tabulated from properties.dat (Tmin..Tmax K)` and `mu tabulated from properties.dat (Tmin..Tmax K)`
per material that uses them.

A table is refused at setup, with the reason and the expected header, when: there is no
`VARIABLES` line, the first column is not `Temperature`, `Cp`, `Density` or an enthalpy column is
missing, a column is named twice, both `Enthalpy` and `Enthalpy_abs` appear, a data row does not
hold a number for every named column, text follows the last row, a zone's `I=` differs from the
rows it holds, the rows are not on consecutive integer kelvins (to 1e-6 K), a value is not
finite, a density or cp is not positive, the enthalpy does not increase with T or disagrees with
cp (a constant cp: `h = cp·T + hOff` on every row; a varying cp: the trapezoid sum), or a relative
`Enthalpy` of a constant cp has an offset. The `Psat` column of an evaporating material is refused
for a value that is not finite, a negative pressure, a pressure that decreases with T (equal
neighbours, such as low-T values rounded to 0, are accepted), a constant pressure,
`boiling-temperature` outside `[Tmin, Tmax−1]`, or `psat(boiling-temperature)` more than a factor 2
away from one atmosphere (a unit slip or a wrong liquid/vapour pair; a database curve sits within
a few percent). The `Lv`, `Mv`, `Cpv` columns of an evaporating material are refused for a value
that is not finite, a column zero on some rows and not on others, a negative value, an `Mv` that
differs between rows, any of them without a valid `Psat` column, or an `Lv` whose Clausius-Clapeyron
slope at `boiling-temperature` is more than a factor 2 from the `Psat` column's
(`ln(p(i0+1)/p(i0))` against `½(Lv(i0)+Lv(i0+1))·Mv/Ru·(1/i0 − 1/(i0+1))`, `i0 = int(boiling-temperature)`,
`Mv` from its column or else from the key: `Lv` in kJ/kg is a factor 1000); a non-zero key beside a
non-zero column of the same quantity is refused as well. The `Sigma` and `Viscosity` columns of a
breaking material are refused for a value that is not finite, a column zero on some rows and not on
others, or a value that is not positive; a non-zero key beside a non-zero column is refused as well,
and so is a breaking material with neither a column nor a key > 0.

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
