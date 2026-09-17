# swirl-wedge-deposit — the source and euler fields under swirl (W-plan Q7, ledger O23)

**Purpose.** `swirl-wedge` pins the parcel dynamics on the 2.5D wedge; this twin pins what
those parcels *deposit*: the momentum source (`OUTPUT/source.tec`, the per-segment
momentum exchange `Pin − Pout`) and the eulerian field (`OUTPUT/euler1.tec`, the
path-weighted parcel velocity). Both are meridian-plane fields — a cell (x, r) carries
(axial, radial, azimuthal) components for an axisymmetric parent solver — while a parcel on
the wedge sits at an azimuth θ inside the sector (|θ| ≤ δ/2 plus the overshoot before the
fold), so the Cartesian vectors it carries must be rotated by −θ about the axis before they
are deposited: `Lib_Equations::toMeridian`, applied to the segment's (Pin, Pout) at the
segment's mid azimuth (`Lib_Integration.f90`, before `computeSrcField`) and to the euler
moment states `∫v|v|dt` in every `Lib_RHS` system. Same fixture, parcels and ODE settings
as the parent (`INPUT` is a symlink); `out-file = ALL` writes both fields and
`mollify = off` keeps every deposit in the cell that received it, so the fields split
cleanly by radius: P3 (r 0.30–0.35) is alone in its band, P1+P2 (r 0.50–0.53) in theirs.
Registered `swirl-wedge-deposit` (labels `e2e;standard;drag;gas;axisym;euler`).

**Oracle** (`check.py`, all derived; reuses the parent's RK4 cylindrical ODE). The force a
Stokes parcel exerts on the gas is `−ṁ a = ṁ (v − g)/τ`, a Cartesian vector; along the
exact trajectory its meridian components are `F_r = ṁ ∫ v_r/τ dt` (outward — the parcel
drifts out against drag), `F_θ = ṁ ∫ (w − Ωr)/τ dt` (P2 spinning up holds the gas back, P3
over-spun pushes it forward), `F_x = 0` (u ≡ U0), and `E = ṁ (|v_0|² − |v_T|²)/2` (Tp
constant; |v| is fold-invariant, so the energy telescopes exactly). **Not** `ṁ (v_r,in −
v_r,out)`: the difference of the meridian *components* carries the frame terms `w²/r` and
`−v_r w/r`, curvature of the coordinate frame and not force on the gas — those belong to
the gas solver's own geometric source terms. On this fixture that wrong reference has the
opposite sign in r (−6.1e-6 against the true +3.7e-6), which is how the deposit half of O23
was first misread (below). Per-cell euler oracle: the ord2 deposit lands on the dual cell
and `finalizeEUL` projects it onto the four neighbouring geo cells with sub-octant volume
weights (Favre); on this uniform mesh the weights are equal to 3 % and neighbouring dual
averages differ by ~2 %, so the geo value is the time average of (v_r, w) over the
2dx × 2dr window about the geo cell centre to ~1e-4 relative. Cell volumes for the mass
check come from the tec nodes: `dx · δ · (r_j+1² − r_j²)/2`.

**Gates.**

| id | assertion | tolerance |
|---|---|---|
| D0 | the parent's trajectory gates (G0–G6) hold on this run — the euler states ride in the same ODE system | parent's |
| D1 | per band, `ΣFy` vs `F_r` and `ΣFz` vs `F_θ`; `\|ΣFx\| ≤ 1e-12`; `\|F_θ\| ≥ 1e-6` in both bands (non-vacuity); sign: P3 `+`, P1+P2 `−` | 2 % (measured 0.02–0.22 %) |
| D2 | per band, `ΣE` vs the kinetic drop | 2 % (measured 0.01–0.15 %) |
| D3 | every touched cell (source or euler) has centre r ∈ [0.28, 0.56]: mollify is off and no deposit is mis-attributed by a cell | locality |
| D4 | every P3 cell holding ≥ 1 % of the band's peak ρ_p (382 cells; ≥ 150 required): `\|v_p − ⟨v_r⟩\|`, `\|w_p − ⟨w⟩\|`, `\|u_p − U0\|` | 1e-4 / 1e-4 / 1e-9 (floor 9e-7 / 1.5e-5 / 0) |
| D5 | `Σ ρ_p V` over the P3 band = `ṁ T_P3` (tiling invariant, cf. `axis-200` E1) | 5e-3 (measured 1.000127) |

**Measured on the fix binary (2026-09-17).** D1 P3: F_r +2.2592e-6 (oracle +2.2623e-6),
F_θ +3.0059e-6 (+3.0054e-6); P1+P2: F_r +1.4356e-6 (+1.4388e-6), F_θ −6.4123e-6
(−6.4169e-6). D2: +3.6706e-7 (+3.6709e-7), −1.9267e-7 (−1.9295e-7). D4: max|Δv_p| 9.0e-7,
max|Δw_p| 1.45e-5 over 382 cells. 797 cells touched, r ∈ [0.300, 0.540].

**Proven RED (2026-09-17).**
- **Cartesian deposit** (the pre-fix binary `ee4482a`, same INI): D4 red on **325 of 382**
  cells, max|Δv_p| 9.9e-4 (10× tol), median 2.9e-4, max|Δw_p| 2.7e-4 (2.7×). The signature
  is `v_p ≈ 0` in the first cells while `⟨v_r⟩` grows: the parcel leaves θ = 0 on a straight
  line, so its Cartesian `v_y` stays ~0 while `v_r = w sin θ` is real; downstream the
  residual oscillates with the fold cadence and is biased to the overshoot side
  (mean −1.7e-4). D1 was **green** on that binary (0.09–0.4 %): summed over a sector sweep
  the `sin θ` mixing cancels, so the totals cannot see the frame — D4 is the detector.
- **Meridian-component difference** (a throwaway build rotating Pin at the entry azimuth
  and Pout at the exit azimuth, i.e. depositing `ṁ Δ(v_r, w)`): D1 red on every band and
  component — P3 F_r −3.32e-6 vs +2.26e-6 (247 %), F_θ +4.28e-6 vs +3.01e-6 (42 %); P1+P2
  F_r −2.79e-6 vs +1.44e-6 (294 %), F_θ −5.88e-6 vs −6.42e-6 (8 %). D4 stays green (the
  euler states are pointwise). This is the "obvious" fix and it is wrong physics; the
  gate exists so nobody lands it.

**Floor.** Trajectories identical to the parent's (sorted multiset) on the S and ALL runs
and across the fix; the fields are one-parcel-per-thread deterministic (the `!$OMP ATOMIC`
re-association floor is ≤ 5e-15 relative on `source.tec`).

**Suite A/B of the fix** (`ee4482a` → this binary, full serial suite as sorted multisets):
byte-inert everywhere except the two off-plane fixtures — `swirl-wedge/source.tec` (the
fix itself, ≤ 0.5 % per cell, ungated there) and `wedge-fold` (`out-file = e`: the
rotated euler states change the SDIRK4 step sequence, trajectory rows flip by one F12.6
ULP and ID 1 prints two more rows, 349 vs 347; gate green).

**Assumption the gate rests on.** `toMeridian` rotates into the frame at θ = 0 as defined by
`refDir` (the y axis, not user-settable), so the fixture's k-planes at ∓δ/2 make cell frame and
meridian frame coincide — a sector [0, δ] would deposit in a frame δ/2 off the cell centre, and
this gate could not tell. That layout is now refused at mesh import
(`refuse-wedge-offcentre`, RED-proven: the pre-tripwire binary ran it silently); the MOSE
production wedge centres to exactly 0.0. Two more facts the residuals depend on, both
verified in the code: `toMeridian` is fold-invariant (the fold rotates position and
velocity by the same −nδ, so `R(−θ(p))·v` is continuous across it — that is why D4
reaches the floor rather than sawtoothing), and rotating `Pin` and `Pout` by one `pMid`
makes the source deposit `R(−θ_mid)·ṁ∫(v−g)/τ dt`, a midpoint quadrature of the meridian
integral with error O(Δθ_seg²) — which is why D1 holds here and degrades exactly in the
many-sector case below.

**Not pinned here.** Four of the five `Lib_RHS` euler-moment sites (models 2–5: evaporation,
breakup, their combinations) are exercised only on planar fixtures, where `toMeridian` is
the identity — the index splits are the same at all five and the A/B was inert on them, but
only model 1 runs off-plane. The many-sector interior sub-case of O23 (a parcel whose segment
sweeps several sectors before the fold, so the dual cell located by `(x, r cos θ)` sits
inward of its true radius and the deposit lands there): every wedge fixture folds within
≈ 0.6° of the sector edge, so it needs a near-axis parcel with `w ≈ |v|` and a deposit
oracle of its own. The `ord1` deposit path (`gas-order = 1`) is exercised by no gate.
