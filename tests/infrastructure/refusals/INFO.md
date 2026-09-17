# refusals — setup-time refusals, gated (coverage audit G20)

**Purpose.** IGLOO refuses several inputs with an `error stop` at a single choke point in
`read_cdp_properties` (`src/lib/IO.f90`) — models that are parsed and threaded through but whose
physics is absent (ledger O3/O4), a misleading model combination (evaporation plan F1 guard), and
a `properties.dat` that does not match `phase.txt` (ledger O24 residual, both checks added
2026-09-16). None of these refusals had a test: a regression that turned a refusal into a silent
default — the failure mode the choke point exists to prevent — would have passed the suite.

**Harness.** `igloo_refusal_case(name dir "expected" labels)` in `tests/CMakeLists.txt` runs the
solver and hands its exit code to `tools/check_refusal.py`, which requires ALL of: exit **exactly
128** (ifx `error stop 'msg'` — exactly why `igloo_e2e_case`, whose rc ≥ 128 test means "died on
signal", cannot host a refusal; 0 = a plain `stop` or no refusal, 129–255 = a signal death after
the message); the expected text in **`run_err.txt`** — the `error stop` payload, i.e. the string of
the stop that ended the run, never an `[ERROR]` line merely printed to stdout before some other
death; and nothing set up past the choke point — no "Placing particles" (injection), "Compute
particles dynamics" or "Stop condition" line, no `OUTPUT/*.dat`. **Non-vacuity proven**: the same
tool on drag-stokes's own output fails on every count (exit 0, no payload, injected, integrated,
output written).

| case | fixture (all else = `drag-stokes`) | expected `error stop` payload |
|---|---|---|
| `refuse-p2t` | `INPUT/` symlinked; `[IGLOO-Models] liquid-conduction = P2T` | `liquid-conduction=P2T parsed but not implemented` |
| `refuse-zgr` | `INPUT/` symlinked; `[IGLOO-Models] boiling = ZGR` | `boiling=ZGR parsed but not implemented` |
| `refuse-lk-d2law` | `INPUT/` symlinked; `evaporation = d2-law` + `interface = LK` (F1 guard fires before the `[IGLOO-Properties]` requirements) | `interface=LK needs evaporation in` |
| `refuse-properties-zones` | own `INPUT/`: drag-stokes's `bc.txt`/`solfile.tec`, two-mat's `phase.txt` (A, B), the ONE-zone `common/properties.dat`; `input.ini` also carries two-mat's `[GPB-Phase2]` | `zone count /= number of materials` |
| `refuse-properties-range` | as above with zone B cut to 3000 rows | `zones do not share one temperature range` |
| `refuse-drag-token` | `drag = no-such-drag-law` | `IGLOO: unknown drag model` |
| `refuse-heat-token` | `heat = no-such-nusselt-law` | `IGLOO: unknown heat model` |
| `refuse-breakup-token` | `breakup = no-such-breakup-model` | `IGLOO: unknown breakup model` |
| `refuse-evaporation-token` | `evaporation = no-such-evaporation-model` | `IGLOO: unknown evaporation model` |
| `refuse-evaporation-leb` | `evaporation = LEB` (parsed, declared not implemented) | `IGLOO: evaporation=LEB is not implemented` |
| `refuse-tab-method` | tab-e2e (`INPUT/` symlinked) with `method = 3` | `IGLOO: TAB method must be 1 or 2` |
| `refuse-gas-order` | `gas-order = 3` | `IGLOO: gas-order must be 1 or 2` |
| `refuse-out-file` | `out-file = nonsense` | `IGLOO: unknown out-file token` (accepted: E, S, E+S, ALL -- `ALL` is what hydra's MI2 cases write; `two-mat` runs it) |
| `refuse-ode-solver` | `ode-solver = H-radau5` | `IGLOO: unknown ode-solver` |
| `refuse-wedge-offcentre` | a generated wedge whose k-planes sit at 0 and +1° (`tools/make_wedge_case.py --theta0-deg 0.5 --nx 4 --nr 4`) | `IGLOO: wedge sector must be centred on the azimuth origin (k-planes at -+delthe/2)` — refused at mesh import (`allocation.f90`), not at the properties choke point |

All fifteen exit 128 at setup with nothing injected or integrated (2026-09-17). **The four INI-contract
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

**Not here.** `solidification = on` (O5) has no global key — only the per-material section, whose
name is mid-migration (`[GPB-Phase<i>]` → `[IGLOO-Material<i>]`); add it once that lands. TAB
A mixed-Nk mesh (`allocation.f90`, `af307dd`) needs a two-zone `solfile.tec`; exercised in a
throwaway only. `refuse-dopri5` (the interim `H-dopri5` refusal of ledger O25) was dropped with the
`Lib_INI.f90` refusal and the `DISABLED` flag on `standard/drag-stokes-dopri5` once the OSlo fork carried
the SOLOUT fix (`fix/dopri5-solout`); that case is the gate now.
