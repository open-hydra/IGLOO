# Material Properties

A material's properties come from two places:

- **`INPUT/properties.dat`**, one table per material in temperature, written by ATLAS. Heat capacity, density and
  enthalpy always come from it, and any other property may.
- **`[IGLOO-Properties]`** in `input.ini`, one constant per material, for the properties the table does not give.

Give each property in one place only: a non-zero column beside a non-zero key for the same property stops the run.

---

## `[IGLOO-Properties]`

Each key takes one value per material, in the order of the phase file. A material that has the property's column
gives 0.

| Key | Unit | Needed by | Column that replaces it |
| :--- | :---: | :--- | :--- |
| `Lv` | J/kg | evaporation; `combustion=direct` | `Lv` |
| `Mv` | kg/kmol | evaporation | `Mv` |
| `cpv` | J/(kg K) | evaporation (ASM, TC) | `Cpv` |
| `boiling-temperature` (or `Tboil`) | K | evaporation; `direct` without `T-burn` | — |
| `boiling-pressure` | Pa | the pressure of `boiling-temperature`; default 101325 | — |
| `Le` | – | evaporation: the Lewis number | — |
| `Yinf` | – | evaporation: the far-field vapour mass fraction | — |
| `sigma` | N/m | breakup: surface tension | `Sigma` |
| `mu` | Pa s | breakup: liquid viscosity | `Viscosity` |
| `emissivity` | – | radiation: 0 (default) is off, above 1 is cut to 1 with a warning | `Emissivity` |
| `psat` | — | accepted and ignored: the saturation pressure is Clausius–Clapeyron's, or the `Psat` column | `Psat` |

**Saturation pressure.** Without a `Psat` column the saturation pressure follows Clausius–Clapeyron through
(`boiling-temperature`, `boiling-pressure`), using `Lv` and `Mv`. A boiling point given at another pressure than one
atmosphere is moved once at setup to the boiling point at one atmosphere on the same line, and the log prints both
([Evaporation](../theory/evaporation.md#saturation-pressure)). Setup stops the run when `Lv` or `Mv` is not positive,
when the line has no boiling point at one atmosphere, or when `boiling-pressure` is not positive and finite.

**Emissivity.** IGLOO reads its own key only. ICE models the same exchange through `[ICE-Physics] emissivity`: a hydra
case that runs both condensed solvers sets both keys, with the same value for a material both carry.

---

## `properties.dat`

One zone per material, in the order of the phase file. Every zone covers the same temperatures, one row per
integer kelvin from $T_\mathrm{min}$ to $T_\mathrm{max}$. The first line with `VARIABLES` names the columns:

```text
VARIABLES = "Temperature", "Cp", "Density", "Enthalpy", "Psat", "Lv", "Mv", "Cpv", "Sigma", "Viscosity", "Conductivity", "Emissivity"
```

| Column | Unit | |
| :--- | :---: | :--- |
| `Temperature` | K | first; the nodes |
| `Cp` | J/(kg K) | required |
| `Density` | kg/m³ | required; the liquid's for a solidifying material |
| `Enthalpy` or `Enthalpy_abs` | J/kg | required, exactly one ([below](#the-enthalpy-datum)) |
| `Psat` | Pa | evaporation: the saturation pressure |
| `Lv`, `Mv`, `Cpv` | J/kg, kg/kmol, J/(kg K) | evaporation, with `Psat`: latent heat, vapour molar mass (constant), vapour heat capacity. `Lv` alone serves `combustion=direct` |
| `Sigma`, `Viscosity` | N/m, Pa s | breakup: surface tension and liquid viscosity |
| `Conductivity` | W/(m K) | `liquid-conduction = P2T`: the liquid conductivity |
| `Viscosity` | Pa s | also `P2T` on a material that does not break |
| `Emissivity` | – | radiation, any material |

After `Temperature` the columns may come in any order, in upper or lower case. A column that is constant in
temperature is a constant; one that varies is a table, linear between the nodes. Outside the table the density and
every optional column hold their end values, and the enthalpy continues its first and last segments. A column of
zeros counts as absent.

Where each column is read:

- `Psat`, `Lv`, `Cpv`, `Sigma`: at the drop's temperature, in every evaluation.
- `Conductivity`: at the drop's mean temperature.
- `Viscosity` under `P2T`: at the surface temperature (the breakup still uses the mean temperature).
- `Emissivity`: at the temperature of each radiative term, the particle's or, under `P2T`, the surface's.
- `Lv` for `direct`: once, at the burning temperature.

ATLAS writes `Psat`, `Lv`, `Mv` and `Cpv` from a liquid/vapour pair of its database (`[GPB-Phase*] psat-vapour`), so
that the latent heat is the slope of the saturation curve; zeros for a material without a pair. For a tabulated `Lv`
the material's single latent-heat value is `Lv` at 298.15 K.

With a varying `Cp` the parcel's state is its enthalpy, read back to a temperature through the same table, so the
table may start at any temperature ([Governing equations](../theory/governing-equations.md#the-property-tables)). A
solidifying material keeps the temperature as its state and uses the table's chord as its liquid heat capacity.

The setup log prints, per material, which columns it uses (`p_sat tabulated from properties.dat ...`,
`vapour data tabulated ...`) and warns about a non-zero column that no model of the material reads.

### The enthalpy datum

The enthalpy column carries a datum:

- `Enthalpy` is **relative**: $c_p T$, or the integral of $c_p$ from $T_\mathrm{min}$;
- `Enthalpy_abs` is **absolute**: the formation enthalpy is included (a thermo table, or a constant $c_p$ with
  `[GPB-Phase*] h0`).

For a constant $c_p$ IGLOO keeps the offset $h_\mathrm{off} = h(T_\mathrm{min}) - c_p T_\mathrm{min}$ (0 for a relative
table) and credits the gas with mass at $c_p T + h_\mathrm{off}$. The datum matters only when mass goes to a gas solver
whose energy is absolute (hydra): a relative table then leaves the coupler to fix the level. A carrying material
needs an absolute table, since its parcel holds two materials, each on its own datum
([Energy](../theory/combustion.md#energy-the-absolute-datum-and-the-streams)). The log prints each material's datum.

### Checks at setup

??? note "What stops the run"
    **Any table**

    - no `VARIABLES` line, or a first column other than `Temperature`;
    - `Cp`, `Density` or an enthalpy column missing; a column named twice; both `Enthalpy` and `Enthalpy_abs`;
    - a row without a number for every column, text after the last row, or a zone whose `I=` differs from its rows;
    - rows not on consecutive integer kelvins (to $10^{-6}$ K);
    - a value that is not finite; a density or $c_p$ that is not positive;
    - an enthalpy that does not rise with $T$ or disagrees with $c_p$ (a constant $c_p$: $h = c_p T + h_\mathrm{off}$
      on every row; a varying one: the trapezoid sum); a relative table of constant $c_p$ with an offset;
    - a zone count that differs from the material count, or zones on different temperatures.

    **`Psat`, of an evaporating material**

    - a negative pressure, one that falls with $T$ (equal neighbours are accepted), or a constant one;
    - `boiling-temperature` outside $[T_\mathrm{min}, T_\mathrm{max} - 1]$;
    - $p_\mathrm{sat}$(`boiling-temperature`) more than a factor 2 from `boiling-pressure`: a unit slip or a wrong
      liquid/vapour pair (a database curve sits within a few percent).

    **`Lv`, `Mv`, `Cpv`, of an evaporating material**

    - a column zero on some rows only, a negative value, or an `Mv` that changes between rows;
    - any of them without a valid `Psat` column;
    - an `Lv` whose Clausius–Clapeyron slope at `boiling-temperature` is more than a factor 2 from the `Psat`
      column's (`Lv` in kJ/kg is a factor 1000 off).

    **`Sigma`, `Viscosity`, `Conductivity`, `Emissivity`**

    - a column zero on some rows only, or a value that is not positive;
    - a breaking material with neither the column nor the key > 0;
    - a `P2T` material with neither `Conductivity` nor `k-liq` > 0, or, when it does not break, neither
      `Viscosity` nor `mu-liq` > 0.

    **Any column with its key**: a non-zero key or line word (`k-liq`, `mu-liq`, `emissivity`, ...) beside a
    non-zero column of the same property.
