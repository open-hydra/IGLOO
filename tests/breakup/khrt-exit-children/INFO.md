# INFO — breakup/khrt-exit-children (KH-RT children shed on their parent's exit step; e2e)

## What it verifies
A breakup event is applied only on a final step, and a parcel's last step, the one that crosses the outlet, is one.
A KH shed on that step makes a child at the parent's final position: past the outlet, in the parent's last cell. The
child must leave with the flow the parent lost, the parent's plane record at the outlet must carry the flow over the
step (before the shed), and the child must not record a crossing of the plane it was born beyond.
[khrt-plane](../khrt-plane/INFO.md)'s P5 states that rule, but no child of its run (nor of any other case) is shed
on its parent's exit step; this case makes every parent shed there.

## The case
The 0.15 m box of `make_box_case.py` (60 × 5 × 5 cells, gas 1.2 kg/m³, 1.8e-5 Pa s, 600 K, uniform U = 200 m/s),
water from [khrt-e2e](../khrt-e2e/INFO.md)'s `phase.txt` and `properties.dat` (symlinked). Five drops of d = 0.6 mm,
1e-4 kg/s each, are placed on the outlet's side of the last cell at x = 0.14999999995 m (5e-11 m before the outlet
face, y = 0.005 … 0.045, z = 0.025) at u_p = 0.01 m/s (`[IGLOO-BC]`, the inlet table unused). `breakup =
Reitz-KHRT` with `B1 = 0.1` and `mShedLim = 1e-4`, the other constants at their defaults; Schiller-Naumann drag,
Ranz-Marshall heat; H-sdirk4 at 1e-9; `print-dtime = 1e-12` (a row at every segment end), `out-time = on`, `plane-x
= 0.15`.

`INPUT/`, from this directory:

    python3 ../../tools/make_box_case.py INPUT --ug 200
    ln -s ../../khrt-e2e/INPUT/phase.txt INPUT/phase.txt
    ln -s ../../khrt-e2e/INPUT/properties.dat INPUT/properties.dat

The arithmetic of the first step (IGLOO's formulas): the first segment's `deltat` = max(|ds|, eps)/|v| (eps = 1e-8 m),
its step bound 0.99 `deltat`, SDIRK4's first trial step 1e-2 of it: 9.90e-9 s, 9.96e-11 m, so the step crosses the
outlet and moves less than eps, which makes it final. Over it the drop at We_r 200 strips toward `dStable` = 1.985e-5 m
(s = `dStable`/d = 0.0331) with τ_KH = 5.30e-6 s: 5.41e-3 of the drop's mass, 54 times `mShedLim` and 149 times s³
(the two conditions of the shed), so the shed fires; no RT (the clock's 9.9e-9 s against τ_RT = 1.21e-4 s). With the
default `B1 = 20`, `mShedLim = 0.03` the step strips 2.7e-5: no shed. `B1 = 0.1` is 100 times below Reitz's 10, on
purpose: it makes one step of 9.9e-9 s strip enough to shed; the case gates the bookkeeping of the exit step, not
the stripping rate.

## Gates (`check.py`)
| id | gate | tolerance | measured |
|---|---|---|---|
| X0 | the premise: five injected parcels (IDs 1–5, parent 0), each the parent of at least one child born on its exit step — the child's trajectory zone names it, the child's exit record has `bc` 400 and the parent's exit `t` (the same ES16.8 string), the child's first row the parent's exit X (the same F12.6 string); every exit record has `bc` 400; `run_out.txt` reports as many children of Loop 1 as there are child exit records | exact (the printed strings) | one child per parent, 5 in all, d 0.198461E-04 m, exit at t 9.90000010E-09 s |
| X1 | the outlet plane: each record's `mdot` = its parcel's exit-record `mdot` plus those of its exit-step children (the parent's record carries the flow over the step, before the shed); no exit-step child has a record there; one record per injected parcel | the E13.6 half-ULPs of the exit records plus the record's ES16.8 one | 5 of 5 records, worst 0.431 of the bound |
| X2 | the sum of `mdot` over every exit record = the injected flow, the `[IGLOO-BC]` positions times their `mdot` (5 × 1e-4 kg/s) | the sum of the exit records' E13.6 half-ULPs (2.5e-10 kg/s) | residue 1.1e-10 kg/s |

X0's X clause discriminates nothing at an outlet exit: every exit row and every child's first row prints x =
0.150000. The exit time and the child's zone carry X0, and the RED below fails it.

The child's exit record is written at its first finalize, a step of zero length: it starts outside its cell and
leaves through the outlet at its parent's time. Each child carries 5.41e-3 of its parent's flow (0.541322E-06 against
the parent's 0.994587E-04 kg/s after the shed). No step was rejected: one shed per parent.

## RED
| change | what fails |
|---|---|
| the child of the exit step dropped: the shed capture in the segment finalize skipped on a crossing step (`.and. .not. IamOut`; the shed still lowers the parent's flow) | X0 0 of 5 parents; X1 each record 5.41e-7 kg/s above its exit flow against a bound of 5.1e-11 (1.07e4 times); X2 the exits 2.71e-6 kg/s short against 2.5e-10 (1.08e4 times). The rest of the suite stays green under this change: no other case has such a child |
| `B1` and `mShedLim` at their defaults (20, 0.03) | X0: no child; X1 and X2 green |

## Not covered
The children here are born past the outlet only. A child shed on the exit step through a wall (a parent leaving by
any other boundary) and a step whose first trial is rejected (the shed then fires inside the domain and again on the
crossing step; X0 asks for at least one exit-step child) are not exercised.
