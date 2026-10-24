# refusals — setup-time refusals, gated (coverage audit G20)

**Purpose.** IGLOO refuses several inputs with an `error stop` at a single choke point in
`read_cdp_properties` (`src/lib/IO.f90`) — models that are parsed and threaded through but whose
physics is absent (ledger O3/O4), a misleading model combination (evaporation plan F1 guard), and
a `properties.dat` that does not match `phase.txt` (ledger O24 residual, both checks added
2026-09-16). None of these refusals had a test: a regression that turned a refusal into a silent
default — the failure mode the choke point exists to prevent — would have passed the suite.

**Harness.** `igloo_refusal_case(name dir "expected" labels)` in `tests/CMakeLists.txt` runs the
solver and hands its exit code to `tools/check_refusal.py`, which requires ALL of: exit **exactly
the compiler's `error stop` code** — probed at configure time by `tools/stop_code_probe.f90` and
handed over as `IGLOO_ERROR_STOP_RC`: 128 under ifx/ifort (exactly why `igloo_e2e_case`, whose
rc ≥ 128 test means "died on signal", cannot host a refusal), 1 under gfortran; 0 = a plain `stop`
or no refusal, 129–255 = a signal death after the message; the expected text in **`run_err.txt`** — the `error stop` payload, i.e. the string of
the stop that ended the run, never an `[ERROR]` line merely printed to stdout before some other
death; and nothing set up past the choke point — no "Placing particles" (injection), "Compute
particles dynamics" or "Stop condition" line, no `OUTPUT/*.dat`. **Non-vacuity proven**: the same
tool on drag-stokes's own output fails on every count (exit 0, no payload, injected, integrated,
output written).

| case | fixture (all else = `drag-stokes`) | expected `error stop` payload |
|---|---|---|
| `refuse-p2t` | `INPUT/` symlinked; `[IGLOO-Models] liquid-conduction = P2T` on a material that does not evaporate | `liquid-conduction=P2T requires an evaporating material` |
| `refuse-p2t-kliq` | `evaporation/d2law` (`INPUT/` symlinked: no `k-liq` on its material line) + `liquid-conduction = P2T` | `liquid-conduction=P2T requires k-liq > 0` |
| `refuse-p2t-muliq` | as `refuse-p2t-kliq` with its own `INPUT/phase.txt` (`A 1 k-liq=0.6`), the other three files symlinked | `liquid-conduction=P2T requires mu-liq > 0` |
| `refuse-p2t-mu-two` | `breakup/khrt-evap-frozen` (Reitz-KHRT, `[IGLOO-Properties] mu = 1.0e-3`) + `liquid-conduction = P2T`, its own `INPUT/phase.txt` (`A 1 k-liq=0.6 mu-liq=1.0e-3`), `khrt-e2e`'s other three files symlinked: a breaking material has one liquid viscosity, the breakup one | `liquid-conduction=P2T with breakup takes the breakup viscosity` |
| `refuse-p2t-mu-zero` | as `refuse-p2t-mu-two` with `A 1 k-liq=0.6` and `mu = 0.0` | `requires mu > 0` (`IGLOO: Breakup requires mu > 0: …`) |
| `refuse-zgr` | `INPUT/` symlinked; `[IGLOO-Models] boiling = ZGR` on a material that does not evaporate | `boiling=ZGR requires an evaporating material` |
| `refuse-lk-d2law` | `INPUT/` symlinked; `evaporation = d2-law` + `interface = LK` (F1 guard fires before the `[IGLOO-Properties]` requirements) | `interface=LK needs evaporation in` |
| `refuse-properties-zones` | own `INPUT/`: drag-stokes's `bc.txt`/`solfile.tec`, two-mat's `phase.txt` (A, B), the ONE-zone `common/properties.dat`; `input.ini` also carries two-mat's `[GPB-Phase2]` | `zone count /= number of materials` |
| `refuse-properties-range` | as above with zone B cut to 3000 rows | `zones do not share one temperature range` |
| `refuse-properties-two-enthalpy` | own 10-row `INPUT/properties.dat` whose header names `"Enthalpy"` and `"Enthalpy_abs"` | `the enthalpy datum is ambiguous` |
| `refuse-properties-bad-row` | own 10-row `properties.dat`, row 3's enthalpy typed `375O.000000` (a letter O) | `unreadable rows: line 7 does not hold 4 numbers` |
| `refuse-properties-psat` | tc-box (`bc.txt`/`solfile.tec` symlinked), its table on 1..1000 K plus a `"Psat"` column = its Clausius-Clapeyron line (`tools/make_psat_table.py --anchors 500 0.018150089922984093 600 1.0`, the first ratio `exp(-Lv Mv/Ru (1/500 - 1/600))`) with the 549 K and 550 K values swapped | `Psat column: a pressure that decreases with T` |
| `refuse-pure2d-mesh` | own `INPUT/solfile.tec`: `planar-slab`'s gas file reduced to its first node plane, `K=1` (`python3 make_pure2d.py`, run from this directory); `planar-slab`'s `bc.txt` and settings | `IGLOO: a pure 2D mesh (one node plane) is not supported: give IGLOO the mesh of one cell layer` — refused at mesh import (`allocation.f90`): ORION reads the one-plane zone as a block of `Nk = 0` cells |
| `refuse-pure2d-bc` | own `INPUT/bc.txt`: `planar-slab`'s records of faces 1–4 with 5-integer headers `b i j f code`, as ATLAS writes a pure-2D mesh (`make_pure2d.py`); `planar-slab`'s gas file and settings | `IGLOO: bc.txt has 5 header integers (a pure 2D ATLAS file): give IGLOO the boundary file of the one-cell-layer mesh` — refused from the first header in `read_cdp_bc_file` |
| `refuse-bc-copies` | own `INPUT/`: `phase.txt` `A 3` (three families), `two-fam-bc`'s two-copy `bc.txt`, `common/properties.dat` — 2500 records where one copy (1250) or one per family (3750) is readable | `IGLOO: bc.txt record count is neither one copy nor one copy per family` — refused in `read_cdp_bc_file`, which prints the three counts |
| `refuse-bc-copy-order` | own `INPUT/`: `two-fam-bc`'s file with two face-3 headers of copy 2 exchanged (`python3 ../../two-fam-bc/make_fixture.py --swap`, run from the case directory), `two-mat`'s `phase.txt`/`properties.dat` | `IGLOO: bc.txt family copies do not repeat the faces of copy 1` — the codes still line up, only the header integers differ (`family 2 of block 1 ... at record 51`) |
| `refuse-chimera-counts` | own `INPUT/` (`python3 make_chimera.py`, run from this directory): a box of two overlapping blocks joined by chimera (102) faces (`tools/make_box_case.py --nx 8 --ny 2 --nz 2 --blocks 2 --overlap 1 --families 2`), the counts line of block 1's first chimera record holding three integers | `IGLOO: chimera donor counts line is not two non-negative integers` — refused at the first scan of the file |
| `refuse-chimera-donor-line` | the same box, the first donor line of that record holding three integers and a weight | `IGLOO: chimera donor line is not four integers and a weight` |
| `refuse-chimera-donor-cell` | the same box, that donor naming cell i = 99 of block 2 (5 cells along x) | `IGLOO: chimera donor list names a cell outside the donor block` |
| `refuse-chimera-weights` | the same box, that donor's weight 0.9 (its layer sums to 0.9) | `IGLOO: chimera donor list is not normalised` |
| `refuse-chimera-copies` | the same box, copy 2 of that record writing the counts `2 0` against copy 1's `1 1` (the same number of payload lines) | `IGLOO: bc.txt family copies do not repeat the chimera donor counts of copy 1` |
| `refuse-periodic-not-mutual` | own `INPUT/` (`python3 make_periodic.py`, run from this directory): a box with faces 3/4 a periodic pair (`tools/make_box_case.py --nx 8 --ny 2 --nz 2 --periodic-y`) whose face-4 records are all header-only 400, so no face-3 facelet's partner names it back | `IGLOO: a 101/201 connection is not mutual: the partner cell does not name this cell` — refused when the field mollifier's link table is built |
| `refuse-ord2-orient` | own `INPUT/` (`python3 make_ord2.py`, run from this directory): a box of two blocks joined by 101 faces (`tools/make_box_case.py --nx 8 --ny 2 --nz 2 --blocks 2`) at `gas-order = 2`, block 1's first 101 record carrying the orientation `-1 0 0 1`: the face corners it names lie one cell from the partner's | `IGLOO: 101 orientation disagrees with the partner cells` — refused when the dual cells of the interface vertices are joined (`allocation.f90`) |
| `refuse-ord2-slab` | the same generator on a one-layer box (`--nz 1`), that record carrying `1 1 0 1` (d12 = 1) | `IGLOO: 101 orientation mixes the slab direction on a one-layer mesh` |
| `refuse-drag-token` | `drag = no-such-drag-law` | `IGLOO: unknown drag model` |
| `refuse-heat-token` | `heat = no-such-nusselt-law` | `IGLOO: unknown heat model` |
| `refuse-breakup-token` | `breakup = no-such-breakup-model` | `IGLOO: unknown breakup model` |
| `refuse-evaporation-token` | `evaporation = no-such-evaporation-model` | `IGLOO: unknown evaporation model` |
| `refuse-phase-token` | own `INPUT/phase.txt` (`A 1 no-such-key=1`), drag-token's `bc.txt`/`solfile.tec` symlinked | `IGLOO: unknown per-material key in the phase file` |
| `refuse-phase-token-real` | own `INPUT/phase.txt` (`A 1 alpha-e=one`) | `IGLOO: non-numeric per-material value in the phase file` |
| `refuse-breakup-token-material` | own `INPUT/phase.txt` (`A 1 breakup=Foo`), drag-stokes's `bc.txt`/`solfile.tec` symlinked, no `breakup` key in `[IGLOO-Models]` | `IGLOO: unknown breakup model` |
| `refuse-tab-method` | tab-e2e (`INPUT/` symlinked) with `method = 3` | `IGLOO: TAB method must be 1 or 2` |
| `refuse-gas-order` | `gas-order = 3` | `IGLOO: gas-order must be 1 or 2` |
| `refuse-out-file` | `out-file = nonsense` | `IGLOO: unknown out-file token` (accepted: E, S, E+S, ALL -- `ALL` is what hydra's MI2 cases write; `two-mat` runs it) |
| `refuse-ode-solver` | `ode-solver = H-radau5` | `IGLOO: unknown ode-solver` |
| `refuse-boiling-temperature-both` | tc-box (`INPUT/` symlinked) with `boiling-temperature` and `Tboil` both in `[IGLOO-Properties]` | `give boiling-temperature or its alias Tboil, not both` |
| `refuse-solid-hfus` | own `INPUT/phase.txt` `A 1 solidification=on cp-solid=600` (no `h-fus`), drag-stokes's `bc.txt`/`solfile.tec` | `solidification=on requires h-fus > 0` |
| `refuse-solid-cpsolid` | `A 1 solidification=on h-fus=1.07e6` (no `cp-solid`) | `solidification=on requires cp-solid > 0` |
| `refuse-solid-tnuc` | the valid solidification line plus `T-nuc=2400` (above the 2327 K default `T-melt`) | `solidification=on requires T-nuc < T-melt` |
| `refuse-solid-evap` | the valid line plus `[IGLOO-Models] evaporation = CEM` | `solidification=on is exclusive with evaporation` |
| `refuse-solid-comb` | the valid line plus `combustion=Beckstead K-burn=4.5e-7` | `solidification=on with combustion is not supported` |
| `refuse-solid-varcp` | the valid line plus an own 10-row `properties.dat` whose `Cp` column varies (the enthalpy its trapezoid sum, so the table itself is valid) | `solidification=on requires a constant-cp material` |
| `refuse-rho-solid` | the valid line plus `rho-solid=-1` | `IGLOO: rho-solid must be positive` |
| `refuse-wedge-offcentre` | a generated wedge whose k-planes sit at 0 and +1° (`tools/make_wedge_case.py --theta0-deg 0.5 --nx 4 --nr 4`) | `IGLOO: wedge sector must be centred on the azimuth origin (k-planes at -+delthe/2)` — refused at mesh import (`allocation.f90`), not at the properties choke point |
| `refuse-psd-missing` | own `INPUT/bc.txt` (`tools/make_box_case.py --law 0.0 tables/psd.dat`, shared by the six table cases below), no table | `IGLOO: size-distribution file not found` — after trying `tables/psd.dat` and `INPUT/psd.dat` (both printed) |
| `refuse-psd-line-long` | own `bc.txt`: every inlet payload line padded to 600 characters with blanks between columns 7 and 8 | `IGLOO: bc.txt inlet payload line too long` — a line that fills the 512-character buffer |
| `refuse-psd-token-long` | own `bc.txt`: column 8 a path of 256 characters | `IGLOO: bc.txt inlet payload token too long` — a token that fills the 256-character buffer |
| `refuse-psd-law-typo` | own `bc.txt`: column 8 `lognorml` | `IGLOO: unsupported distribution law` — neither a law name nor the name of a file (no `lognorml`, no `INPUT/lognorml`) |
| `refuse-psd-unsorted` | `INPUT/psd.dat` radii 5, 4, 6 µm | `IGLOO: size-distribution file: the radii must increase` |
| `refuse-psd-negative` | a mass fraction of −0.3 | `IGLOO: size-distribution file: a mass fraction is negative` |
| `refuse-psd-three-col` | a data line `5.0e-6 0.3 1.0` | `IGLOO: size-distribution file: a data line must hold exactly two finite numbers` |
| `refuse-psd-empty` | comment and blank lines only | `IGLOO: size-distribution file holds no data rows` |
| `refuse-psd-radius-zero` | a first radius of 0 | `IGLOO: size-distribution file: a radius is not positive` |
| `refuse-psd-mass-zero` | every mass fraction 0 | `IGLOO: size-distribution file: the mass fractions sum to zero` |
| `refuse-comb-brk` | `burn-box` (`INPUT/` symlinked) with `[IGLOO-Models] breakup = TAB` | `combustion=Beckstead with breakup is not supported` |
| `refuse-comb-kburn` | own `INPUT/phase.txt` `A 1 combustion=Beckstead n-burn=1.8` (no `K-burn`, so 0), `burn-box`'s `bc.txt`/`solfile.tec`, `common/properties.dat`; `burn-box`'s settings | `combustion=Beckstead requires K-burn > 0` |
| `refuse-comb-nburn` | as above with `K-burn=4.5e-07 n-burn=3.0` | `combustion=Beckstead requires 0 < n-burn < 3` |
| `refuse-comb-xeff` | as above with `K-burn=4.5e-07 X-eff=1.5` | `combustion=Beckstead requires 0 < X-eff <= 1` |
| `refuse-cap-product` | own `INPUT/phase.txt`: `burn-box`'s line with `xi-cap=0.2` and no `product`; `burn-box`'s `bc.txt`/`solfile.tec`, `common/properties.dat`; `burn-box`'s settings | `xi-cap > 0 requires product = a 0-group material of the phase file` |
| `refuse-cap-varcp` | `combustion/burn-box-cap`'s `phase.txt` (symlinked) and its own ten-row two-zone table, the `OX` zone's `Cp` rising 0.01 J/(kg·K) per K (the enthalpy its integral, so the table itself is valid; `burn-box-cap/make_table.py varcp`) | `xi-cap > 0 requires a constant-cp, constant-density material and product` |
| `refuse-cap-relative` | as above with a relative ten-row table (an `Enthalpy` column, `make_table.py relative`) | `xi-cap > 0 requires the absolute enthalpy datum (an Enthalpy_abs column)` |
| `refuse-brk-sigma` | `tab-e2e` (`INPUT/` symlinked) without `[IGLOO-Properties] sigma` | `Breakup requires sigma` |
| `refuse-brk-mu` | `tab-e2e` without `mu` | `Breakup requires mu` |
| `refuse-sigma-size` | `tab-e2e` with `Le = 1.0` added and `sigma = 0.072 0.072`, `mu = 1.0e-3 1.0e-3`: the vectors are read in the order psat, Mv, Lv, Tboil, cpv, Le, Yinf, sigma, mu and the first present (`Le`, one value) sets the length | `sigma: size mismatch with other property vectors` |
| `refuse-sigma-nm` | `tab-e2e` with `sigma = 0.072 0.072`, `mu = 1.0e-3 1.0e-3` for its one material | `sigma: size /= number of materials` |
| `refuse-sigma-both` | `pe-sigma-table` (`INPUT/` symlinked: its `Sigma` and `Viscosity` columns) with `[IGLOO-Properties] sigma = 0.072` | `IGLOO: sigma given both in [IGLOO-Properties] and as a properties.dat column` |
| `refuse-sigma-column` | `pe-sigma-table`'s case with its own `properties.dat`, the `Sigma` column at −0.01 N/m on the 320 K row (`bc.txt`, `solfile.tec` symlinked to `pilch-erdman-e2e`'s) | `Sigma column: a value that is not positive` |
| `refuse-bc-legacy` | own `INPUT/bc.txt`: `drag-stokes`'s with every code in the legacy single-digit schema (401 → 4, 400 → 1; `python3 make_bc_refusals.py`, run from this directory) | `IGLOO: bc.txt line is not six integers with a 3-digit code` |
| `refuse-bc-empty-payload` | own `INPUT/bc.txt`: `drag-stokes`'s with the first inlet payload line empty (`make_bc_refusals.py`) | `IGLOO: bc.txt inlet payload line is empty` |
| `refuse-body-accel` | `body-force` (`INPUT/` symlinked) with `body-accel = 0.0 -200.0` | `body-accel: expected 3 values (gx gy gz)` |
| `refuse-plane-x` | `INPUT/` symlinked; `plane-x = abc` | `IGLOO: [IGLOO-General] plane-x must be a list of x positions [m]` |
| `refuse-vapour-cc` | mhb98-water-psat's case (`bc.txt`/`phase.txt`/`solfile.tec` symlinked to mhb98-water's), its ATLAS water table with the `Lv` column divided by 1000 (kJ/kg written as J/kg; every other byte unchanged), no `Lv`/`Mv`/`cpv` keys | `Lv column: not the slope of the Psat column` |
| `refuse-vapour-ini-both` | the same case and table (`properties.dat` symlinked to mhb98-water-psat's), `[IGLOO-Properties] Lv = 2.462478e6` kept beside the `Lv` column | `IGLOO: Lv given both in [IGLOO-Properties] and as a properties.dat column` |

All seventy-three exit 128 at setup with nothing injected or integrated. **The two pure-2D cases are RED** on the solver without
their checks: the one-plane mesh runs on a block of zero cells in k and stops at the record count
(`IGLOO: bc.txt record count is neither one copy nor one copy per family`, 880 records against 800 faces), the
5-integer file at the legacy-schema stop (`bcdef = 1`, the first header's sixth value borrowed from the next
line) — neither is their payload. **The two vapour-column cases are RED**
on a reader that does not know the `Lv`, `Mv`, `Cpv` columns: both stop at `Evaporation requires Mv` (exit 128,
not their payload). **The two `Sigma` cases are RED** on a reader that does not know the `Sigma` and
`Viscosity` columns: `refuse-sigma-both` stops at `Breakup requires mu` (its only key is `sigma`),
`refuse-sigma-column` at `Breakup requires sigma` (exit 128, not their payload). **The ten size-distribution cases are
RED** on a reader that splits the inlet line with a list-directed read: a path ends at its slash (`tables`,
`INPUT`), so nine of them stop at `unsupported distribution law` with exit 1 and no payload, and the
600-character line is cut to 512 characters and run as a Dirac (exit 0, output written). The six table
cases are also RED on a loader that resolves the file but refuses every table (`IGLOO: tabulated size
distribution not implemented`, not their payload). **The five chimera cases are RED** on a
reader without the 102 record, which stops each of them at the first chimera payload (`bad BC cell line`,
exit 1, not their payload), and on the reader with its five chimera checks taken out, which runs each of them
to "Stop condition" with exit 0 and writes its output. **`refuse-periodic-not-mutual` is RED** on the link
table of 101 records only, which never looks at the 201 faces: the run integrates to "Stop condition" with
exit 0 and writes its output. **The two gas-order-2 interface cases are RED** on the solver without the interface
pass, which never reads the orientation integers: both run to "Stop condition" with exit 0 and write their output.
**The three oxide-cap cases are RED** on a reader that does not know the `product` token: `refuse-cap-product`
(whose line names no product) runs the burn law to "Stop condition" with exit 0 and writes its output, the other two
stop at `unknown per-material key` (exit 128, not their payload); with its check taken out of `resolve_products`
each of the three runs to "Stop condition" with exit 0. **The three `properties.dat` cases were
RED first**: the reader that took three columns by position ran each of them to "All particles out of
domain" with exit 0 — the two-enthalpy header as an absolute datum with `hOff = 0`, the typed row with
its enthalpy undefined (ORION keeps only the last row's read status), and the decreasing `Psat` column
unread. **The seven solidification cases were
RED first**: the binary before the solidification model stopped each of them with `solidification=on parsed but
not implemented yet`, not their payload. **The four INI-contract
cases were RED first** (ledger O20/O21/O26, recorded on the pre-fix binary): `method = 3` ran to
non-finite states with exit 0 and three output files; `gas-order = 3` silently ran as 2 (the warning was
dead — the value was reassigned before the test); `out-file = nonsense` ran as both with a warning (and
the registry default `E+S` parsed as `E`); `ode-solver = H-radau5` was a plain `stop`. (A fifteenth,
`refuse-dopri5`, was the interim O25 refusal — dropped, see below.) **The five token cases were RED first** too (ledger O26): every model selector in `Lib_Drag`, `Lib_Heat`, `Lib_Breakup` and
`Lib_Evaporation` printed its menu and executed a plain `stop` — exit 0, empty stderr, nothing a
harness or an embedding parent can tell from success (recorded on the pre-fix binary in the session
scratchpad); all twelve `stop`s became `error stop 'IGLOO: unknown …'`.
**`refuse-wedge-offcentre` was RED first** too (2026-09-17, W-plan Q7): the fold band (`axisymFold`,
|θ| ≤ δ/2), the 2.5D gas sample (`sampleGas2D`) and the meridian-frame deposits (`toMeridian`) all take
`refDir` as the sector *centre* — the MOSE/ATLAS layout, k-planes at ∓δ/2, which the production nozzle
solution centres to exactly 0.0 — and on a [0, δ] sector the binary before the tripwire ran the case to
"All particles out of domain" with exit 0, the gas and every deposit rotated by δ/2.
**The two `bc.txt` cases were RED first** too, on a reader that consumed copy 1 only: `bc-copies` ran
its three families of A on copy 1's `krho` 0.34, i.e. `krhoTot` = 1.02 and ṁ = 0.34/(1 − 1.02) · 1.2e-3 =
−2.04e-2 kg/s per parcel — a negative mass flow, 75 exits, exit 0 — and `bc-copy-order` ran like
`two-fam-bc` fed copy 1 for both families (exit 0, six output files).

**The eleven cases of the phase file's burn-law and breakup checks, the property sizes, the two boundary files
and `body-accel` are RED** on the solver with each of their error stops turned into a plain write: eight run on
to "Stop condition" with exit 0 and write their output, and `sigma-size`, `sigma-nm` and `bc-legacy` end on a
later stop with another payload (`mu: size mismatch`, `mu: size /= number of materials`, the record count).
**`refuse-breakup-token-material` is RED** on a reader that knows no `breakup` token on a material line: it stops
at `IGLOO: unknown per-material key in the phase file`, not its payload. **`refuse-plane-x` is RED** on a solver that
does not read `plane-x` (an unknown `[IGLOO-General]` key is ignored): it runs to "Stop condition" with exit 0, no
payload, and writes three `.dat` files.

**Not here.** A mixed-Nk mesh (`allocation.f90`, `af307dd`) needs a two-zone `solfile.tec`; exercised in a
throwaway only. `refuse-dopri5` (the interim `H-dopri5` refusal of ledger O25) was dropped with the
`Lib_INI.f90` refusal and the `DISABLED` flag on `standard/drag-stokes-dopri5` once the OSlo fork carried
the SOLOUT fix (`fix/dopri5-solout`); that case is the gate now.
