# INFO — solidification/solid-remelt-khrt (a KH-RT drop that freezes, re-melts and breaks up; e2e)

## What it verifies
Model 7 carries the KH-RT clock `[told, tc]` in its state and advances it in the liquid phases only (0 and 1); the
RT breakup fires once `tc` exceeds `τ_RT` while `λ_RT < d`. At a re-melt (the melting plateau reaching `f = 0`) the
segment finalize restarts the clock at the interpolated event time, `breakupStateAtSwitch` returning `[t_rm, 0]`, so the
first RT breakup of the molten drop follows the re-melt by `τ_RT`. Without the restart the clock keeps the state of the
last liquid call before the nucleation, and the first liquid call after the re-melt adds the whole frozen interval.
SG13 ([test_solidification](../INFO.md)) gates the routine's table; this case gates its application end to end.

## The case
The 0.15 m box of `make_box_case.py` (60 × 5 × 5 cells, gas 1.2 kg/m³, 1.8e-5 Pa s, 10 m/s) with the gas temperature
`T = 1500 + 20000 x` K at the cell centres; [solid-brk-khrt](../solid-brk-khrt/INFO.md)'s material (`tests/common/
properties.dat`: ρ 2950, c_l 1250; `solidification=on`, `h-fus` 1.07e6, `cp-solid` 600, `T-melt` 2327 K and `T-nuc`
1861.6 K by default), models (Stokes, Ranz-Marshall, `Reitz-KHRT`, σ 1e-3, μ 1e-3, H-sdirk4 at 1e-11) and drop (30
µm at 35 m/s), injected at `kt · 1525` = 1862.6 K, 1 K above `T-nuc` (the inlet cell's centre at 1525 K). `B0` at its
default 0.61: no KH shed. `print-dtime = 1e-12` and no `print-dcell`: a row at every segment end, so the re-melt and
each RT event have their own row. 25 parcels, one per inlet cell. σ = 1e-3 N/m is a wiring value, about 600 times below
molten alumina's, and the RT breakup fires at We_r 3.9 (the RT branch has no Weber-number floor): the case gates the
re-melt's restart of the clock, not the physics of a breaking alumina drop.

`INPUT/`, from this directory:

    python3 ../../tools/make_box_case.py INPUT --tg 1500 --gradient 20000 --kv 3.5 --kt 1.221377 --rp 15e-6
    cp ../solid-box/INPUT/phase.txt INPUT/phase.txt
    ln -s ../../../common/properties.dat INPUT/properties.dat

The sequence (every parcel alike): KH stripping in the liquid (We_r 11 at injection, above `WeLimit` 6) takes d to
29.9147 µm before the nucleation at 1.098e-5 s, where the recalescence gives `f0 = c_l (T_m − T_nuc)/h_fus` = 0.5437;
the plateau freezes to `f` = 0.699 at x = 0.0425 m, the first cell whose gas is above `T_m`, then melts and re-melts
at `t_rm` = 4.168642e-3 s (x = 0.123263 m). Its slip then 15 m/s (We_r 4.0, no stripping), the drop keeps its d to the
first RT at `t_RT` = 4.308980e-3 s (x = 0.126751 m), `τ_RT` = 1.403378e-4 s after the re-melt; the drop breaks up eight more
times within 0.09 ms, down to 0.254 µm at the outlet (the log's `Oh > 1` warning: Oh = 1 at d = 0.68 µm).

## Gates (`check.py`)
| id | gate | tolerance | measured |
|---|---|---|---|
| R0 | per parcel, rows in time order: a row in phase 2 or 3 (`datcols.solid_phase`) before 0.1 ms; a later row with `f = 0` in phase 0 (the re-melt, `t_rm`); no row before `t_rm` with d below 0.9 d0; the first row after `t_rm` with d below 0.9 of the previous row's (the RT, `t_RT`), every row between them with the re-melt row's d; `λ_RT < d` at both rows (the clock runs throughout); the RT row's d = `λ_RT` of the re-melt row's d at the RT row's `u_p` (one RT breakup, its products of diameter `λ_RT`); every injected parcel | exact on the printed values; the product diameter within the E13.6 half-ULPs of the two d | 25 of 25 parcels; `λ_RT/d` 0.494 at the re-melt, 0.498 at the RT; the product diameter at 0.461 of its bound |
| R1 | `|(t_RT − t_rm) − τ_RT| ≤ 1e-8 s`, `τ_RT = Ctau/Ω_RT` from the RT row's `u_p` and the re-melt row's d: `a = 18 μ_g (u_g − u_p)/(ρ_l d²)` (the case's Stokes law), its component along `u_p`, `force = |a_t| (ρ_l − ρ_g)`, `Ω_RT = (2 force^{3/2}/(3 (3σ)^{1/2} (ρ_l + ρ_g)))^{1/2}`; σ, `Ctau`, `CRT` from `input.ini`, ρ_l from `properties.dat`, the uniform gas from `solfile.tec`, `T-melt` from `phase.txt` | 1e-8 s | worst 1.7e-10 s (25 of 25 parcels) |

R1's tolerance: the clock telescopes (`tc += t − told` at every call), so at the RT it is the elapsed time since the
re-melt; the event is applied on a final step, one that moves less than 1e-8 m (4.0e-10 s at 24.7 m/s); `τ_RT` from
the printed d (E13.6, 1.7e-6 relative, `τ_RT ∝ d^{3/2}`: 3.5e-10 s) and `u_p` (F12.6: 4e-12 s); the two ES16.8 times
(5e-12 s each). The sum, 7.6e-10 s, is 0.08 of the tolerance.

## RED
| change | what fails |
|---|---|
| the finalize's restart of the clock at a phase switch removed (the `breakupStateAtSwitch` of the KH-RT slots, `stateLocal(ind_sb1:ind_sb2)`, in the segment finalize) | the clock keeps `[told, tc]` = [1.27e-5 s, 1.27e-5 s] from the nucleation step's last liquid call; the first RT fires at the re-melt time itself, on a step of zero length (no row of its own), the second one `τ_RT` of the halved drop later: R0 the RT row's d 7.389e-6 m against `λ_RT` 1.496e-5 m (2.5e5 times its bound), R1 4.90e-5 s against `τ_RT` 1.411e-4 s, off by 9.2e-5 s (9.2e3 times the tolerance), every parcel. The rest of the suite stays green under this change |
| the clock restarted at the re-melt step's start instead of the interpolated event time | R1 off by 9.75e-6 s (975 times), the part of the step before the re-melt; R0 holds |

## Not covered
The segment-start catch-up applies the same restart when a segment starts past the re-melt threshold; this re-melt
passes the finalize, so the catch-up route has no e2e witness. A drop that froze fully and melts back through the
solid's melting event, and a re-melt under TAB or ETAB (their oscillator is reset at the nucleation, not at a
re-melt), are not exercised.
