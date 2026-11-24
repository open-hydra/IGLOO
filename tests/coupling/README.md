# coupling — the exchange channels

`test_exchange.f90` drives IGLOO as a parent code does: `setup_static` once, then `reset_state` and `solve` per sweep,
the calls of `IGLOO_Mod_Exchange` in between. Each case directory holds an `input.ini` and an `INPUT` link to the
deck it borrows; the driver's argument picks the channel. Its gates are listed at the top of the source; it exits 1
when one fails.

| Case | Deck | Channel |
| :--- | :--- | :--- |
| `exchange-track/` | `combustion/burn-box` | the state the processed mass leaves at (`track_processed`), `gas_rate` |
| `exchange-absorb/` | `standard/drag-stokes` | mass offered to inert drops and absorbed (`add_absorption_rate`) |
| `exchange-remove/` | `standard/drag-stokes` | mass taken off inert drops (`add_removal_rate`) |
| `exchange-births/` | `infrastructure/zero-group-line` | new parcels of a 0-group line (`add_birth_rate`) |

Run: `./tests/test.sh coupling`. What they show: [Coupling & Outputs](../../docs/vv/coupling.md#exchange-channels).
