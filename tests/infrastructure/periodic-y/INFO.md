# periodic-y — translational periodic BC (bcdef 201)

**Purpose.** First validation of the 201 periodic transport path
(`obj_particles.updateCell → IGLOO_bcBox.periodicTransport`): the reader's
connection/connectionFace parsing, the `T = partner_face_center −
exit_face_center` translation, the velocity-unchanged contract, and the
post-transport relocation. The path had shipped unexercised (no case had a
periodic pair).

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
4. `u≡u_g`, `w≡0`, `T≡T_g` decoupling guards hold through the transports.

**Notes.** The `[BCB-*]`/`[per]` ini patches are ATLAS-side only — the solver
reads the 201 codes from `bc.txt`. MOSE's translational-periodic convention
(translate to partner face, velocity unchanged) is what 201 mirrors; the
rotational-periodic variant would need a velocity rotation and is NOT modeled.
