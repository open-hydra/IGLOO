# INFO — breakup/khrt-stress (A23/S4-E: heavy-shedding robustness, **GREEN**)

## What this is — and is not

A **robustness stress case**, not a validation case. There is no reference paper, no oracle
and no measured quantity compared against theory; `check.py` asserts only that the solver
survives heavy KH shedding with mass intact. Physics validation of KHRT lives in
[breakup/khrt-e2e](../khrt-e2e/INFO.md) (KH-stripping rate vs Reitz-87) and its `check_rt.py`
(RT shatter). Nothing here should ever be cited as evidence that a model is correct.

## Case construction

Same fixture as `breakup/khrt-e2e` — `INPUT/` is a **symlink** to `../khrt-e2e/INPUT`, so the
1 MB `solfile.tec`/`bc.txt`/`properties.dat` set is not duplicated. ⚠ That symlink is relative:
moving *either* case directory breaks it, and the failure looks like a missing gas file.

The single difference from khrt-e2e is `[IGLOO-Models] mShedLim = 0.01` against the production
default of `0.03` (`Lib_INI.f90`): a parent sheds once it has accumulated a 1 % per-drop mass
loss instead of 3 %.

## Why it exists

The A23/S4 work replaced the single `gr%child(ip)` slot with a growable per-parent shed list,
lifting a one-shed-per-parent cap that was a property of the data structure rather than of
Reitz-87. khrt-e2e exercises that path only shallowly — its deepest parent sheds ~22 times —
so list growth under repeated `move_alloc` would otherwise go essentially untested, and any
error in the per-shed strip/create pairing would be at its smallest exactly where the suite
looks.

Measured sensitivity of the fixture (2026-08-04, sprop2, all rc=0, all mass-exact):

| `mShedLim` | children | deepest parcel | wall | peak RSS |
|---|---|---|---|---|
| 0.03 (default, = khrt-e2e) | 241 | 22 | 4.6 s | 28.8 MB |
| **0.01 (this case)** | **701** | **66** | **12.9 s** | **29.8 MB** |
| 0.003 | 2003 | 159 | 31.7 s | 32.3 MB |
| 0.001 | 3673 | 420 | 51.6 s | 34.8 MB |

0.01 was chosen for cost: it drives ~7 doublings of the shed list at a third of the runtime of
0.003. The RSS column is the useful one — a 15× increase in children costs 6 MB, which is what
"the growable lists do not thrash" means concretely.

## What `check.py` gates

1. **≥ 500 children** — a floor; a collapse of the shed path fails loudly rather than silently
   turning this into a duplicate of khrt-e2e.
2. **deepest parcel ≥ 20 sheds** — proves the per-parent list actually grew. If a
   one-shed-per-parent cap ever came back this is the assertion that catches it.
3. **mass conservation to 1e-4** — the point of the case. Every shed both strips the parent and
   creates the child; an error in that pairing compounds with shed count, so it surfaces here
   long before it would at the default `mShedLim`.
4. **the runaway guard did not fire** — `maxShed=1000` (`Lib_Integration.f90`) is a guard, not a
   working limit. Firing it on a merely-aggressive case means either the guard is too low or
   shedding has genuinely diverged; both are findings, so the `[WARNING]` is gated rather than
   tolerated.

The `rc >= 128` signal gate comes from the shared `igloo_e2e_case` command, so a crash or OOM
fails independently of any assertion above.

## Related

- Bug ledger: `tests/BUGS.md` **A23** (A23g, A23h, the one-shed cap).
- Thread-invariance of this same path: `khrt-e2e-threads` (`../khrt-e2e/check_threads.py`).
