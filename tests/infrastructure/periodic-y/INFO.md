# periodic-y — translational periodic BC (bcdef 201)

**Purpose.** Validates the 201 periodic transport path
(`obj_particles.updateCell → IGLOO_bcBox.periodicTransport`): the reader's
connection/connectionFace parsing, the `T = partner_face_center −
exit_face_center` translation, the velocity-unchanged contract, and the
post-transport relocation.

**Setup.** `standard/body-force` decoupling (kV=kT=1, Stokes exact-linear drag)
with `g_y = −8000` so every particle drifts ≈2.8 box heights across the domain;
faces 3/4 (y-min/y-max) are a 201 pair (`tools/make_box_case.py --periodic-y`,
which writes the 3D connection lines `block i j k face + 4 weights`).

**Gates (check.py).**
1. `v_y(x)` follows the anchored closed form across ALL wraps (velocity mangling
   at any transport breaks the exponential) — worst resid/tol ≈ 0.42;
2. `y(x)` equals the unwrapped closed form **modulo Ly** (measured residual
   ≈ 5·10⁻⁷ m — the transport is an exact ±Ly translation);
3. every particle wraps ≥ 2 times (measured 2–3; non-vacuous) and all 25 exit
   at the outlet;
4. `u≡u_g`, `w≡0`, `T≡T_g` decoupling guards hold through the transports;
5. **S-P, the source smoothed across the periodic pair**: `check.py` runs the case
   again with `mollify = off` (`ref-raw/`), smooths that `source.tec` itself
   (8 passes of the [¼, ½, ¼] filter on the box's uniform volumes, zero flux in x
   and z, the first and last y rows neighbours) and holds `OUTPUT/source.tec` to
   it, field by field, to 1e-12 of the field's scale (measured 1.3e-14). The
   case's 25 parcels start one cell apart and their deposits are uniform in y to
   about 1e-6, so the pair is run again with only face 1's first row injecting
   (`ref-row/`, `ref-row-raw/`): there periodic and zero-flux smoothing differ by
   0.55 of the scale, and the same comparison holds to 1.1e-14.
6. **G3, one trajectory zone per parcel image**: `trajectories-A.dat` holds one
   zone per parcel and image (`ZONE T="Mat A Group  1 ID <id> parent 0[ image
   <k>]"`, `I` its rows); a parcel's images are numbered 0, 1, 2, …, and a wrap
   (a jump in y of more than Ly/2 between consecutive rows) sits at each image
   boundary (70) and inside none. `check.py` also runs the case with
   `print-ds = 2e-3` and `out-time = on` in place of `print-dcell` (`ref-ds/`):
   there each image ends on its exit row on a y face and the next begins on the
   partner face, Ly away, at the same time, every other column alike (70
   boundaries).

**RED (S-P).** With zero flux at the 201 faces the case's own source differs
from the periodic smoothing by 2.5e-8 of the scale (2.5e4 tolerances) and the
first-row run's by 0.36.

**RED (G3).** One trajectory zone per group: no parcel title. Without the image
split: 25 zones with 70 wraps inside them, in the run and in `ref-ds/`.

**Notes.** The `[BCB-*]`/`[per]` ini patches are ATLAS-side only — the solver
reads the 201 codes from `bc.txt`. MOSE's translational-periodic convention
(translate to partner face, velocity unchanged) is what 201 mirrors; the
rotational-periodic variant would need a velocity rotation and is NOT modeled.
