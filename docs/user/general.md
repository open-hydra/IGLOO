# Run Settings

`[IGLOO-General]` sets where the gas comes from, what IGLOO writes, and the limits of a run. The keys are grouped
below by what they do; the [Parameter Registry](registry.md) lists them all in one table.

---

## The gas and the phase

| Key | Default | Meaning |
| :--- | :---: | :--- |
| `gas-file` | — | the Tecplot file of the gas field ([below](#the-gas-file)); not read under hydra, which hands the gas over in memory |
| `gas-order` | `2` | how the gas is sampled at a parcel: `1` the value of the cell's centre, `2` a second-order reconstruction ([Gas-state interpolation](../theory/governing-equations.md#gas-state-interpolation)) |
| `phase` | — | the ATLAS name of the condensed phase (`[GPB-Phase*] name`). The input files become `INPUT/<phase>-phase.txt`, `-properties.dat`, `-bc.txt`, and every output takes the prefix `<phase>-`. It must not contain `-`. Absent: the unnamed files standalone, the name hydra assigns under hydra |
| `body-accel` | `0 0 0` | a uniform body acceleration $g_x\ g_y\ g_z$ [m/s²], gravity for instance |

### The gas file

The gas file is an ASCII Tecplot file of ordered zones, one per block, in BLOCK packing:

- a `VARIABLES` line, whose names may run on over the following lines up to the first `ZONE`;
- per block, a `ZONE` record with `I`, `J` and `K` (it may span several lines; the first line starting with a number
  starts the values);
- the node coordinates, then the variables, all at the nodes or all at the cell centres
  (`VARLOCATION=([...]=CELLCENTERED)`), the same in every zone.

Values may be separated by blanks, tabs, commas or line ends, and a line may be of any length.

A file with POINT packing, CR line ends, a value longer than 100 characters or a null value (`/`, `r*`) is read by
ORION's reader instead, which cuts every line at 1000 characters; the log then says
`[IGLOO] gas file read by ORION's reader:` and why. On a file whose lines are under 1000 characters the two
readers give the same field, bit for bit. Setting `IGLOO_TEC_FAST=0` in the environment sends every gas file to
ORION's reader.

---

## What IGLOO writes

The output files are described in [Running IGLOO](using.md#output-files).

| Key | Default | Meaning |
| :--- | :---: | :--- |
| `out-file` | `E+S` | which fields: `E` the Eulerian fields of the particles, `S` the source fields of the gas, `E+S` both (also `S+E`, `ES`, `SE`, `ALL`, `both`, any case). Another word stops the run |
| `out-traj` | `on` | the trajectory files |
| `out-scatter` | `on` | the scatter cloud |
| `out-time` | `off` | `on` adds the parcel's time [s] as the last column of every trajectory and exit row |
| `plane-x` | — | planes normal to $x$ at the listed positions [m], e.g. `plane-x = 0.001 0.18`: one record per crossing in `planes-<mat>.dat` ([Plane crossings](using.md#plane-crossings)). An empty value or a word that is not a number stops the run |
| `timing` | `off` | `on` prints the wall time of each phase, the parcel counts and the load balance after each write ([Parallelization](../development/parallelization.md#per-phase-timers)) |

### Trajectory rows

A trajectory row is written at the end of an integration segment. Three keys choose which segment ends get one:

| Key | Default | A row when |
| :--- | :---: | :--- |
| `print-dcell` | off | every `print-dcell`-th cell crossing, at `time-end` and at a boundary exit; it disables `print-ds` |
| `print-ds` | off | the parcel is this far [m] from its last row; ignored with `print-dcell` |
| `print-dtime` | off | this much time [s] has passed since the last row; ignored with `print-ds` |

- With `print-dcell` set, `print-dtime` still adds rows: at the first segment end on or after a target time
  that starts one `print-dtime` after the parcel's start and moves on by `print-dtime` at every row written.
- Rows fall on segment ends, which are cell faces: a `print-ds` equal to the cell size is decided by each step's
  overshoot ([Trajectories](using.md#trajectories)).
- A value of 0 or less turns `print-ds` and `print-dtime` off. A `print-dcell` below 1, or any of the three that is
  not a number, stops the run.

### Scatter cloud

| Key | Default | Meaning |
| :--- | :---: | :--- |
| `fsample-traj` | `100` | the number of markers per injection stream. A parcel writes $\dot n_\mathrm{ref}/\Delta N$ markers per second, $\dot n_\mathrm{ref}$ its number rate at entry and $\Delta N$ the quantum this key sets |
| `scatter-floor` | `0.1` | the floor $F \ge 0$ that keeps sparse streams visible: a parcel whose number rate is below $F$ times its material's mean writes markers at that floor (at least $F$ × `fsample-traj` per domain transit), with weights scaled down so the cloud stays unbiased. `0` turns it off. A negative value, a list or a word stops the run |

A breakup child keeps its parent's rate.

---

## Run limits

| Key | Default | Meaning |
| :--- | :---: | :--- |
| `time-end` | off | every parcel still in the domain stops at this time [s], with a last trajectory row; its state goes to `snapshot-<mat>.dat` |
| `max-reflections` | `10000` | a parcel is dropped after this many reflections at boundary faces; below 1 keeps the value in force |
| `graze-standoff` | `0.01` | after a grazing impact on a symmetry face, the distance the parcel is set back from it, in parcel diameters ([below](#reflections-on-symmetry-faces)); 0 or less is ignored with a notice |
| `corner-clear` | `on` | after a reflection, set the parcel clear of every symmetry face of its cell, not only the one it hit |
| `seed` | `42` | the seed of the random draws (injection diameters, TAB and ETAB breakup) |
| `mdot-max` | `0` | the largest mass flow of one parcel [g/s]; the inlets place more parcels to respect it ([Placing the parcels](boundary-conditions.md#placing-the-parcels-igloo-bc)); 0 is off |

`time-end`, `max-reflections`, `graze-standoff` and `corner-clear` are read again at every sweep, so a coupled run
can change them between sweeps (`time-end` only when set at start-up).

### Reflections on symmetry faces

An impact whose normal speed is below 2 % of the speed is grazing: the parcel is set back
`graze-standoff` diameters along the face normal, at least $10^{-9}$ m and at most 1 % of the cell's size along it.
This is the height a parcel pressed onto the plane falls back from. A smaller value resolves the contact more finely
at more reflections per second, about as the inverse square root: an 11.89 µm parcel bouncing along 0.146 m of plane
under gravity reflects 9 100 times at the default (the `sym-plane-bounce` case). A steeper impact is mirrored and
restarts $10^{-12}$ m inside.

In the corner of two symmetry planes, `corner-clear = on` moves the parcel along the face it hit until it is clear
of the other one too.

---

## Smoothing the fields

| Key | Default | Meaning |
| :--- | :---: | :--- |
| `mollify` | `on` | smooth the source and Eulerian fields before writing them ([Field Mollification](../theory/mollification.md)) |
| `mollify-passes` | `8` | the number of smoothing passes; the filter width is about $\sqrt{n/2}$ cells. `0` turns smoothing off |
