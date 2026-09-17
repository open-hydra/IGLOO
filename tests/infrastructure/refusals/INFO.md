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

All ten exit 128 at setup with nothing injected or integrated (2026-09-17). **The five token cases
were RED first** (ledger O26): every model selector in `Lib_Drag`, `Lib_Heat`, `Lib_Breakup` and
`Lib_Evaporation` printed its menu and executed a plain `stop` — exit 0, empty stderr, nothing a
harness or an embedding parent can tell from success (recorded on the pre-fix binary in the session
scratchpad); the eleven `stop`s became `error stop 'IGLOO: unknown … model'`. The twelfth, the
`ode-solver` token in `Lib_INI.f90`, still exits 0 — deferred with the O20/O21 work on that file.

**Not here.** `solidification = on` (O5) has no global key — only the per-material section, whose
name is mid-migration (`[GPB-Phase<i>]` → `[IGLOO-Material<i>]`); add it once that lands. TAB
`method = 3` is **not** refused today (ledger O20: the parser accepts any integer and the model
hits a comment-only stub) — it becomes a `refuse-tab-method` case with the O20 fix, and the
`ode-solver` plain `stop` a `refuse-ode-solver` case with O26's last site. A mixed-Nk mesh
(`allocation.f90`, `af307dd`) needs a two-zone `solfile.tec`; exercised in a throwaway only.
