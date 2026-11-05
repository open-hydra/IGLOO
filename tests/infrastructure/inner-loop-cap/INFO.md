# inner-loop-cap — the inner loop's cap, reached by every parcel

**Purpose.** A cell segment is integrated by a first ODE solver call and, when that call ends on a crossing whose
step moved more than `eps` = 1e-8 m, by refinement calls that integrate the cut step again with a shorter `deltat`
(each about ten times shorter). `[IGLOO-ODE] max-inner-iter` (default 10) caps the calls on one segment: a parcel
that needs more is marked gone with an `Inner loop >` line. This case reaches that cap on every parcel, after
finalized segments, with the Eulerian output on, without touching the code.

**Fixture.** A box from `tools/make_box_case.py --kv 1.0 --rp 1.0e-6 --stretch 1000` (60×5×5 cells over 0.15 m):
the left 30 cells Δx₁ = 0.15/(30 × 1001) = 4.995005e-6 m long, the right 30 cells 1000 Δx₁ = 4.995005e-3 m; the
gas uniform at 10 m/s and 600 K; 25 parcels of 2 µm (ρ_p = 2950 kg/m³, Stokes drag, model 1) injected at the gas
velocity and temperature: no slip and no heat exchange, so every accepted state has x = u_g t. Measured (a print at
each finalize, one thread): every left crossing finalizes on its 4th call, the first right crossing needs 7 calls
(at S = 100 the left needs 5 and the right 7; at S = 10⁴ 3 and 7). With `max-inner-iter = 5` every left crossing
finalizes with one call to spare and every parcel is capped on its first right cell two calls short of the
crossing: entry at x_entry = 30 Δx₁ = 1.4985015e-4 m. `out-file = ALL`, `out-time = on`, `mollify = off`, gas order 1.

**Gates** (`check.py`).
- C0 the cap fires where it should: one `Inner loop >` line per parcel (the IDs of the 25 exit records), every exit
  record code 0 inside the box at X ≥ x_entry − 5e-7 m, every trajectory starts with its injection row and one row
  per left-half crossing, row k at k Δx₁ within 5e-7 m (half the F12.6 quantum) + eps, the 31st at x_entry; the run
  ends normally. Measured: exits at X = 0.005145 m, t = 5.1446014e-4 s.
- C3 the Eulerian mass, Σ ρ_p V_cell over the cells (V_cell from `euler1.tec`'s nodes), equals the exit records'
  Σ ṁ t_entry within 1e-5, t_entry the time of each parcel's 31st row (its entry into the capped cell): the capped
  segment deposits nothing and every finalized one deposits ṁ over its residence (model 1, constant ṁ, mollify off).
  Tolerance: the E13.6 `mdot` column bounds the exit side at 8e-7. Measured −2.955e-7.

**RED.** The key at 10 (the default): no cap, 25 outlet exits (code 400), C0 fails. The reset of
`part%deltaL`, `part%Tstay` and `part%intE` at the start of every segment (in `integrate`'s outer loop, after
`initializeCell`) taken out: the capped segment deposits the previous segment's moments a second time, ṁ Δx₁/u_g per
parcel, C3 +3.333e-2 (= 1/30, measured), C0 holds.

**After the cap.** The cap leaves the parcel with the end of its cut step (past the face) and the time of its last
accepted step, and the cut step's crossing flags set: the cell logic relocates it on that stale crossing, and every
trajectory carries a 32nd row there (X = 0.005145 m); its exit record holds the cut step's position with the latched
time.
